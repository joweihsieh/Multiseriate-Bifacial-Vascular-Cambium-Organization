#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})
setwd(tempdir())

# =======================
# USER SETTINGS
# =======================
XENIUM_BASE <- "/home/woodydrylab/FileShare/20260121_Xenium"
TARGET_OUTPUT <- NULL
#TARGET_OUTPUT <- "output-XETG00360__0079031__BIO1_TISSUE_1_and_2__20260115__224443"
#TARGET_OUTPUT <- "output-XETG00360__0079031__BIO2_TISSUE_1__20260115__224443"
#TARGET_OUTPUT <- "output-XETG00360__0079031__BIO2_TISSUE_2__20260115__224443"
#TARGET_OUTPUT <- "output-XETG00360__0079031__BIO2_TISSUE_3__20260115__224443"
#TARGET_OUTPUT <- "output-XETG00360__0079049__BIO1_TISSUE_1__20260115__224443"
#TARGET_OUTPUT <- "output-XETG00360__0079049__BIO1_TISSUE_2_and_3__20260115__224443"
#TARGET_OUTPUT <- "output-XETG00360__0079049__BIO2_TISSUE_1_and_2__20260115__224443"

PLOT_LOWER <- 0.10
PLOT_MAXI  <- 0.75

LCM_THRESHOLDS <- list(
  cambium = list(lower = 0.25, upper = 0.60),
  xylem   = list(lower = 0.25, upper = 0.60),
  phloem  = list(lower = 0.25, upper = 0.45)
)

COL_PATTERN <- "^corr_LCM_"
WRITE_THRESHOLDS_USED_TSV <- TRUE
SLICE_N_ROWS <- 1
BOTTOM_COORD_FILE <- "Selection_Bottom_coordinates.csv"
UPPER_COORD_FILE  <- "Selection_Upper_coordinates.csv"
SHOW_SLICE_SAMPLE_MARKS <- FALSE
SMOOTH_SPAN <- 0.025


# =======================
# HELPERS
# =======================
clean_filename <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

get_lcm_class <- function(cor_colname) {
  x <- sub("^corr_LCM_", "", cor_colname)
  tolower(sub("_.*$", "", x))
}

get_thresholds_for_col <- function(cor_colname, default_lower, default_upper, thr_list) {
  cls <- get_lcm_class(cor_colname)
  if (!is.null(thr_list[[cls]])) {
    return(list(lower = thr_list[[cls]][["lower"]], upper = thr_list[[cls]][["upper"]], cls = cls))
  }
  list(lower = default_lower, upper = default_upper, cls = cls)
}

guess_xy_cols <- function(dt) {
  nms <- names(dt)
  low <- tolower(nms)
  x_idx <- which(low %in% c("x", "x_center", "xcoord", "x_coord", "center_x", "xm"))
  y_idx <- which(low %in% c("y", "y_center", "ycoord", "y_coord", "center_y", "ym"))
  if (length(x_idx) >= 1 && length(y_idx) >= 1) return(list(x = nms[x_idx[1]], y = nms[y_idx[1]]))
  stop("Cannot determine x/y columns from coordinate file.")
}

read_polygon_csv <- function(path) {
  dt <- fread(path)
  xy <- guess_xy_cols(dt)
  poly <- dt[, .(x = as.numeric(get(xy$x)), y = as.numeric(get(xy$y)))]
  poly <- poly[!is.na(x) & !is.na(y)]
  if (nrow(poly) < 3) stop("Polygon file has fewer than 3 valid points: ", path)
  if (!(poly$x[1] == poly$x[nrow(poly)] && poly$y[1] == poly$y[nrow(poly)])) poly <- rbind(poly, poly[1])
  poly
}

point_in_polygon_vec <- function(px, py, vx, vy) {
  n <- length(vx)
  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    xi <- vx[i]; yi <- vy[i]; xj <- vx[j]; yj <- vy[j]
    intersect <- ((yi > py) != (yj > py)) & (px < (xj - xi) * (py - yi) / ((yj - yi) + 1e-30) + xi)
    inside <- xor(inside, intersect)
    j <- i
  }
  inside
}

split_bin_meta_by_selection <- function(bin_meta, sample_dir) {
  bottom_path <- file.path(sample_dir, BOTTOM_COORD_FILE)
  upper_path  <- file.path(sample_dir, UPPER_COORD_FILE)
  if (!file.exists(bottom_path) || !file.exists(upper_path)) return(list(All = copy(bin_meta)))

  poly_bottom <- read_polygon_csv(bottom_path)
  poly_upper  <- read_polygon_csv(upper_path)
  in_bottom <- point_in_polygon_vec(bin_meta$x_center, bin_meta$y_center, poly_bottom$x, poly_bottom$y)
  in_upper  <- point_in_polygon_vec(bin_meta$x_center, bin_meta$y_center, poly_upper$x, poly_upper$y)

  out <- list()
  if (sum(in_bottom) > 0) out[["Bottom"]] <- copy(bin_meta[in_bottom])
  if (sum(in_upper)  > 0) out[["Upper"]]  <- copy(bin_meta[in_upper])
  if (length(out) == 0) out[["All"]] <- copy(bin_meta)
  out
}

get_mid_xslice_profile <- function(bin_meta, cor_colname, y_sel, smooth_span = NULL) {
  dt <- copy(bin_meta)[!is.na(x_center) & (y_center %in% y_sel) & !is.na(get(cor_colname))]
  if (nrow(dt) == 0) return(NULL)

  x_unique_all <- sort(unique(bin_meta$x_center))
  prof_obs <- dt[, .(corr_raw = mean(get(cor_colname), na.rm = TRUE)), by = .(x_center)][order(x_center)]
  prof <- data.table(x_center = x_unique_all)
  prof <- merge(prof, prof_obs, by = "x_center", all.x = TRUE, sort = TRUE)
  prof[, corr_smooth := corr_raw] 

  if (!is.null(smooth_span) && smooth_span > 0) {
    keep_rows <- !is.na(prof$corr_raw)
    if (sum(keep_rows) > 5) {
      loess_fit <- loess(corr_raw ~ x_center, data = prof[keep_rows], span = smooth_span, control = loess.control(surface = "direct"))
      prof[keep_rows, corr_smooth := loess_fit$fitted] 
    }
  }
  list(profile = prof, y_selected_center = mean(y_sel))
}

output_png_figure <- function(plotting_function, output_figure = FALSE, output_path = "temp.png", ...) {
  if (output_figure) {
    png(output_path, pointsize = 10, res = 300, width = 20, height = 18, units = "cm")
    on.exit(dev.off(), add = TRUE)
  }
  plotting_function(...)
}


plot_xenium_with_correlation_and_slice <- function(bin_meta, cor_colname, y_sel, bin_um = 5, lower = 0.10, maxivalue = 0.75,
                                                   reverse_y = TRUE, panel_label = NULL, show_slice_sample_marks = FALSE,
                                                   smooth_span = NULL) {

  x <- bin_meta$x_center; y <- bin_meta$y_center; cor_vector <- bin_meta[[cor_colname]]
  keep <- !is.na(x) & !is.na(y) & !is.na(cor_vector)
  x <- x[keep]; y <- y[keep]; cor_vector <- cor_vector[keep]
  if (length(x) == 0) return(invisible(NULL))

  slice_info <- get_mid_xslice_profile(bin_meta, cor_colname, y_sel = y_sel, smooth_span = smooth_span)
  if (is.null(slice_info) || nrow(slice_info$profile) == 0) return(invisible(NULL))

  y_plot <- if (reverse_y) -y else y
  y_slice_plot <- if (reverse_y) -slice_info$y_selected_center else slice_info$y_selected_center

  cor_for_color <- cor_vector
  cor_for_color[cor_for_color < lower] <- 0
  denom <- (maxivalue - lower); if (denom <= 0) denom <- 1e-6
  color_index <- round(pmax(0, pmin(1, (cor_for_color - lower) / denom)) * 500) + 1
  color_pool <- colorRampPalette(c("#EEF2F9", "#FF0000"))(501)
  plot_order <- order(cor_for_color)

  half <- bin_um / 2
  prof <- slice_info$profile
  xlim_data <- range(x) + c(-half, half)
  ylim_data <- range(y_plot) + c(-half, half)
  y_min <- min(c(prof$corr_smooth, lower), na.rm = TRUE)
  y_max <- max(c(prof$corr_smooth, maxivalue), na.rm = TRUE)

  # =========================================================
  # If this is Cambium, find the highest point on the left and right halves
  # =========================================================
  valid_prof <- prof[!is.na(corr_smooth)]
  peaks_info <- list()
  is_cambium <- grepl("_cambium", cor_colname, ignore.case = TRUE)
  
  if (is_cambium && nrow(valid_prof) > 0) {
    mid_x <- mean(range(valid_prof$x_center))
    
    # Highest point on the left half
    prof_left <- valid_prof[x_center < mid_x]
    if (nrow(prof_left) > 0) {
      idx <- which.max(prof_left$corr_smooth)
      peaks_info[[length(peaks_info) + 1]] <- list(x = prof_left$x_center[idx], y = prof_left$corr_smooth[idx])
    }
    
    # Highest point on the right half
    prof_right <- valid_prof[x_center >= mid_x]
    if (nrow(prof_right) > 0) {
      idx <- which.max(prof_right$corr_smooth)
      peaks_info[[length(peaks_info) + 1]] <- list(x = prof_right$x_center[idx], y = prof_right$corr_smooth[idx])
    }
  }

  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), bg = "black", new = FALSE)
  plot.new()

  dev_w <- par("din")[1]; dev_h <- par("din")[2] 
  fig_left <- 0.085; fig_right <- 0.985
  #top_bottom <- 0.30; top_top <- 0.88 
  #bot_bottom <- 0.075; bot_top <- 0.27
  top_bottom <- 0.4; top_top <- 0.92
  bot_bottom <- 0.075; bot_top <- 0.35

  target_ratio <- ((fig_right - fig_left) * dev_w) / ((top_top - top_bottom) * dev_h)
  data_ratio <- diff(xlim_data) / diff(ylim_data)
  if (data_ratio < target_ratio) {
    x_span_new <- diff(ylim_data) * target_ratio
    xlim_new <- mean(xlim_data) + c(-x_span_new / 2, x_span_new / 2)
    ylim_new <- ylim_data
  } else {
    y_span_new <- diff(xlim_data) / target_ratio
    ylim_new <- mean(ylim_data) + c(-y_span_new / 2, y_span_new / 2)
    xlim_new <- xlim_data
  }
  xticks <- pretty(xlim_data, n = 5); xticks <- xticks[xticks >= xlim_new[1] & xticks <= xlim_new[2]]

  # === TOP PANEL ===
  par(fig = c(fig_left, fig_right, top_bottom, top_top), mar = c(0, 0, 0, 0), bg = "black", new = TRUE)
  plot(NA, NA, xlim = xlim_new, ylim = ylim_new, axes = FALSE, xaxs = "i", yaxs = "i")
  xo <- x[plot_order]; yo <- y_plot[plot_order]; ci <- color_index[plot_order]
  rect(xleft = xo - half, ybottom = yo - half, xright = xo + half, ytop = yo + half, col = color_pool[ci], border = NA)
  abline(h = y_slice_plot, col = "lightblue", lwd = 4)

  if (show_slice_sample_marks) {
    prof_good_x <- prof[!is.na(corr_smooth), x_center]
    if (length(prof_good_x) > 0) points(prof_good_x, rep(y_slice_plot, length(prof_good_x)), pch = 15, cex = 0.25, col = "lightblue")
  }

  # === TOP ANNOTATIONS ===
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), bg = "black", new = TRUE)
  plot.new()
  if (!is.null(panel_label)) text(x = fig_left, y = 0.96, labels = panel_label, col = "white", adj = c(0, 0.5), cex = 2.7)
  cb_left  <- fig_right - 0.30; cb_right <- fig_right - 0.02; cb_y <- 0.98 
  for (i in seq_len(length(color_pool))) segments(cb_left + (i - 1) / length(color_pool) * (cb_right - cb_left), cb_y, cb_left + i / length(color_pool) * (cb_right - cb_left), cb_y, col = color_pool[i], lwd = 12, lend = "butt")
  text(cb_left, cb_y - 0.015, sprintf("%.2f", lower), col = "white", adj = c(0, 1), cex = 2.7)
  text(cb_right, cb_y - 0.015, sprintf("%.2f", maxivalue), col = "white", adj = c(1, 1), cex = 2.7)
  text((cb_left + cb_right) / 2, cb_y - 0.015, "", col = "white", adj = c(0.5, 1), cex = 2.7)

  # === BOTTOM PANEL ===
  par(fig = c(fig_left, fig_right, bot_bottom, bot_top), mar = c(0, 0, 0, 0), bg = "black", new = TRUE)
  plot(NA, NA, xlim = xlim_new, ylim = c(y_min, y_max), axes = FALSE, xaxs = "i", yaxs = "i")
  good <- !is.na(prof$corr_smooth)
  if (any(good)) {
    r <- rle(good); ends <- cumsum(r$lengths); starts <- c(1, head(ends, -1) + 1)
    for (k in seq_along(r$values)) {
      if (r$values[k]) {
        idx <- starts[k]:ends[k]
        lines(prof$x_center[idx], prof$corr_smooth[idx], lwd = 4, col = "lightblue")
        points(prof$x_center[idx], prof$corr_smooth[idx], pch = 16, cex = 0.55, col = "lightblue")
      }
    }
  }
  #abline(h = lower, lty = 2, col = "#6ca6cd") 
  box(col = "lightblue", lwd = 4)

  # === GLOBAL AXES & PRECISION VERTICAL LINES ===
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), bg = "black", new = TRUE)
  plot(0:1, 0:1, type="n", xlim=c(0,1), ylim=c(0,1), axes=FALSE, xaxs="i", yaxs="i", xlab="", ylab="")
  
  # Define global coordinate mapping helpers
  lmx <- function(val) fig_left + (val - xlim_new[1]) / diff(xlim_new) * (fig_right - fig_left)
  lmy <- function(val) bot_bottom + (val - y_min) / (y_max - y_min) * (bot_top - bot_bottom)
  lmy_top <- function(val) top_bottom + (val - ylim_new[1]) / diff(ylim_new) * (top_top - top_bottom)
  
  xt_fig <- lmx(xticks)
  for (i in seq_along(xticks)) {
    segments(x0 = xt_fig[i], y0 = bot_bottom, x1 = xt_fig[i], y1 = bot_bottom - 0.008, col = "lightblue", lwd = 1.5)
    text(xt_fig[i], bot_bottom - 0.015, labels = xticks[i], col = "lightblue", cex = 2.7, adj = c(0.5, 1), font = 2)
  }
  text((fig_left + fig_right) / 2, bot_bottom - 0.055, "", col = "lightblue", cex = 1.5, font = 2)
  yticks <- pretty(c(y_min, y_max), n = 4); yticks <- yticks[yticks >= y_min & yticks <= y_max]
  yt_fig <- lmy(yticks)
  for (i in seq_along(yticks)) {
    segments(x0 = fig_left, y0 = yt_fig[i], x1 = fig_left - 0.008, y1 = yt_fig[i], col = "lightblue", lwd = 1.5)
    text(fig_left - 0.012, yt_fig[i], labels = yticks[i], col = "lightblue", cex = 2.7, adj = c(1, 0.5), font = 2)
  }
  text(fig_left - 0.055, (bot_bottom + bot_top) / 2, "", col = "lightblue", cex = 2.7, srt = 90, font = 2)

  # Draw dashed vertical lines only for Cambium, extending from the peak of the
  # bottom curve up to the horizontal slice line
  #if (length(peaks_info) > 0) {
  #  for (p in peaks_info) {
  #    x_global <- lmx(p$x)
  #    y_start_global <- lmy(p$y)               # Start: the peak of the lower curve
  #    y_end_global <- lmy_top(y_slice_plot)    # End: the horizontal slice line in the top panel
  #    
  #    segments(x0 = x_global, y0 = y_start_global, x1 = x_global, y1 = y_end_global, 
  #             col = "lightblue", lwd = 3, lty = 2)
  #  }
  #}

  invisible(slice_info)
}

# =======================
# MAIN
# =======================
sample_dirs <- list.dirs(XENIUM_BASE, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^.*/output-", sample_dirs)]
if (!is.null(TARGET_OUTPUT)) sample_dirs <- sample_dirs[basename(sample_dirs) == TARGET_OUTPUT]

for (sdir in sample_dirs) {
  grid_dirs <- list.dirs(sdir, full.names = TRUE, recursive = FALSE)
  grid_dirs <- grid_dirs[grepl("/grid05um_out$", grid_dirs)]
  if (length(grid_dirs) == 0) next

  for (gdir in grid_dirs) {
    meta_path <- file.path(gdir, "bin_metadata_LCM.tsv")
    if (!file.exists(meta_path)) next

    bin_um <- 5; bin_meta_all <- fread(meta_path)
    cols_to_plot <- grep(COL_PATTERN, names(bin_meta_all), value = TRUE)
    if (length(cols_to_plot) == 0) next

    region_list <- split_bin_meta_by_selection(bin_meta_all, sdir)
    base_plot_dir <- file.path(gdir, "corr_maps_byLCMthreshold_with_mid_xslice_split")
    dir.create(base_plot_dir, showWarnings = FALSE, recursive = TRUE)

    for (region_name in names(region_list)) {
      bin_meta <- region_list[[region_name]]
      if (nrow(bin_meta) == 0) next
      
      # Ensure all panels use the same absolute horizontal slice line
      y_unique <- sort(unique(bin_meta$y_center))
      y_mid <- mean(range(y_unique))
      y_ord <- order(abs(y_unique - y_mid))
      shared_y_sel <- sort(y_unique[y_ord][seq_len(min(SLICE_N_ROWS, length(y_unique)))])

      for (cc in cols_to_plot) {
        thr <- get_thresholds_for_col(cc, PLOT_LOWER, PLOT_MAXI, LCM_THRESHOLDS)
        out_plot_dir <- file.path(base_plot_dir, region_name, sprintf("%s_L%.2f_U%.2f", clean_filename(thr$cls), thr$lower, thr$upper))
        dir.create(out_plot_dir, showWarnings = FALSE, recursive = TRUE)
        outpng <- file.path(out_plot_dir, paste0(clean_filename(cc), "_", region_name, "_with_mid_xslice_20260331.png"))
        
        output_png_figure(
          plotting_function = plot_xenium_with_correlation_and_slice,
          output_figure = TRUE, output_path = outpng,
          bin_meta = bin_meta, cor_colname = cc, y_sel = shared_y_sel,
          bin_um = bin_um, lower = thr$lower, maxivalue = thr$upper, reverse_y = TRUE,
          #panel_label = region_name, 
          panel_label = NULL,
          show_slice_sample_marks = SHOW_SLICE_SAMPLE_MARKS,
          smooth_span = SMOOTH_SPAN
        )
      }
    }
  }
}
