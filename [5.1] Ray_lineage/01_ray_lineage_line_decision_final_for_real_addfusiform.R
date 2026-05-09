#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


dir_path <- "/home/woodydrylab/FileShare/20260121_Xenium/output-XETG00360__0079049__BIO1_TISSUE_1__20260115__224443/grid05um_out/kmeans_k10_raw_out/ray_match_final_addfu"

if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

setwd(dir_path)


############################################################
# user settings
############################################################
INPUT_TSV  <- "../bin_metadata_with_cluster_raw.tsv"
OUT_PREFIX <- "ray_scan_cluster3_to_5_inner_band"

CLUSTER_LINE   <- 3
CLUSTER_TARGET <- 5
CLUSTER_FUSIFORM <- 6

# scan angle
ANGLE_STEP <- 5
COARSE_ANGLE_SEQ <- seq(-90, 90, by = ANGLE_STEP)
HALF_BIN_WIDTH_DEG <- ANGLE_STEP / 2

# only scan inner side of comp3
INNER_BUFFER <- 0

# distance / band
MAX_SCAN_DIST   <- 500
STEP_SIZE       <- 5
BAND_HALF_WIDTH <- 12

# component split
NEIGHBOR_FACTOR <- 1.5
MIN_COMP_SIZE   <- 3
N_PICK          <- Inf

# overlay line lengths
PCA_AXIS_LEN  <- 0
BEST_AXIS_LEN <- 120

# local density window on step bins
LOCAL_WINDOW_STEPS <- 6

# diagnostic plots
DIAG_BINS <- 40
HIGHLIGHT_COMPS <- c(8, 23, 37, 43, 100, 134)

############################################################
# helper
############################################################
norm_vec <- function(v) {
  s <- sqrt(sum(v^2))
  if (!is.finite(s) || s == 0) return(c(0, 0))
  v / s
}

rotate_vec <- function(v, deg) {
  th <- deg * pi / 180
  c(
    v[1] * cos(th) - v[2] * sin(th),
    v[1] * sin(th) + v[2] * cos(th)
  )
}

find_components <- function(dt_xy, neighbor_dist) {
  n <- nrow(dt_xy)
  if (n == 0) return(integer(0))

  coords <- as.matrix(dt_xy[, .(x_center, y_center)])
  visited <- rep(FALSE, n)
  comp_id <- integer(n)
  cid <- 0L

  neighbor_list <- vector("list", n)

  for (i in seq_len(n)) {
    dx <- coords[, 1] - coords[i, 1]
    dy <- coords[, 2] - coords[i, 2]
    d2 <- dx^2 + dy^2
    neighbor_list[[i]] <- which(d2 <= neighbor_dist^2)
  }

  for (i in seq_len(n)) {
    if (visited[i]) next
    cid <- cid + 1L
    queue <- i
    visited[i] <- TRUE
    comp_id[i] <- cid

    while (length(queue) > 0) {
      cur <- queue[1]
      queue <- queue[-1]
      nbs <- neighbor_list[[cur]]
      nbs <- nbs[!visited[nbs]]
      if (length(nbs) > 0) {
        visited[nbs] <- TRUE
        comp_id[nbs] <- cid
        queue <- c(queue, nbs)
      }
    }
  }

  comp_id
}

get_major_axis <- function(x, y) {
  pts <- cbind(x, y)
  if (nrow(pts) < 2) return(c(1, 0))
  pc <- prcomp(pts, center = TRUE, scale. = FALSE)
  v <- pc$rotation[, 1]
  as.numeric(norm_vec(v))
}

calc_longest_run <- function(x) {
  if (length(x) == 0 || all(x == 0)) return(0L)
  rr <- rle(x)
  max(rr$lengths[rr$values == 1], na.rm = TRUE)
}

calc_max_gap <- function(x) {
  hit_idx <- which(x == 1L)
  if (length(hit_idx) <= 1) return(0L)
  max(diff(hit_idx) - 1L)
}

calc_local_density_binary <- function(x, k = 6L) {
  if (length(x) == 0) return(0)
  if (length(x) < k) return(mean(x))
  vals <- vapply(seq_len(length(x) - k + 1), function(i) {
    mean(x[i:(i + k - 1)])
  }, numeric(1))
  max(vals)
}

calc_local_density_count <- function(x, k = 6L) {
  if (length(x) == 0) return(0)
  if (length(x) < k) return(mean(x))
  vals <- vapply(seq_len(length(x) - k + 1), function(i) {
    mean(x[i:(i + k - 1)])
  }, numeric(1))
  max(vals)
}

save_hist_plot <- function(dt, var_name, out_prefix, bins = 40) {
  p <- ggplot(dt, aes_string(x = var_name)) +
    geom_histogram(bins = bins, fill = "steelblue", color = "white") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Distribution of ", var_name, " across components"),
      x = var_name,
      y = "Count"
    )

  ggsave(paste0(out_prefix, "_hist_", var_name, ".pdf"), p, width = 6, height = 5)
  ggsave(paste0(out_prefix, "_hist_", var_name, ".png"), p, width = 6, height = 5, dpi = 300)
}

save_scatter_plot <- function(dt, xvar, yvar, out_prefix) {
  p <- ggplot(dt, aes_string(x = xvar, y = yvar)) +
    geom_point(
      data = dt[highlight_group == "other"],
      size = 2, alpha = 0.65, color = "grey50"
    ) +
    geom_point(
      data = dt[highlight_group != "other"],
      aes_string(color = "highlight_group"),
      size = 3, alpha = 0.95
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0(yvar, " vs ", xvar),
      x = xvar,
      y = yvar,
      color = "Highlighted comps"
    )

  ggsave(paste0(out_prefix, "_scatter_", xvar, "_vs_", yvar, ".pdf"), p, width = 6, height = 5)
  ggsave(paste0(out_prefix, "_scatter_", xvar, "_vs_", yvar, ".png"), p, width = 6, height = 5, dpi = 300)
}

############################################################
# read data
############################################################
dt <- fread(INPUT_TSV)

required_cols <- c("x_center", "y_center", "cluster_raw")
miss <- setdiff(required_cols, names(dt))
if (length(miss) > 0) {
  stop("Missing required columns: ", paste(miss, collapse = ", "))
}

dt[, cluster_raw := as.character(cluster_raw)]

dt_line   <- dt[cluster_raw == as.character(CLUSTER_LINE)]
dt_target <- dt[cluster_raw == as.character(CLUSTER_TARGET)]
dt_fusiform <- dt[cluster_raw == as.character(CLUSTER_FUSIFORM)]

if (nrow(dt_line) == 0) stop("No cluster 3 bins found.")
if (nrow(dt_target) == 0) stop("No cluster 5 bins found.")

############################################################
# estimate bin size
############################################################
ux <- sort(unique(dt$x_center))
uy <- sort(unique(dt$y_center))
dx <- diff(ux); dx <- dx[dx > 0]
dy <- diff(uy); dy <- dy[dy > 0]

bin_size <- median(c(dx, dy), na.rm = TRUE)
if (!is.finite(bin_size)) bin_size <- 5

neighbor_dist <- bin_size * NEIGHBOR_FACTOR

message("Estimated bin_size = ", round(bin_size, 3))
message("neighbor_dist = ", round(neighbor_dist, 3))

############################################################
# global center
############################################################
center_xy <- c(mean(dt$x_center), mean(dt$y_center))

############################################################
# cluster3 components
############################################################
dt_line[, comp_id := find_components(.SD, neighbor_dist),
        .SDcols = c("x_center", "y_center")]

comp_stat <- dt_line[, .(
  n_bins = .N,
  centroid_x = mean(x_center),
  centroid_y = mean(y_center)
), by = comp_id]

comp_stat <- comp_stat[n_bins >= MIN_COMP_SIZE]

if (nrow(comp_stat) == 0) {
  stop("No cluster 3 components passed MIN_COMP_SIZE.")
}

############################################################
# define each comp3
# start = closest bin to center
# dir0  = PCA major axis direction pointing inward (toward center)
############################################################
comp_info_list <- list()

for (ii in seq_len(nrow(comp_stat))) {
  cid <- comp_stat$comp_id[ii]
  sub <- dt_line[comp_id == cid]

  sub[, r_to_center := sqrt((x_center - center_xy[1])^2 + (y_center - center_xy[2])^2)]

  start_bin <- sub[which.min(r_to_center)][1]
  start_xy <- c(start_bin$x_center, start_bin$y_center)

  major <- get_major_axis(sub$x_center, sub$y_center)
  to_center <- norm_vec(center_xy - start_xy)

  # make dir0 point inward
  if (sum(major * to_center) < 0) {
    major <- -major
  }

  pts <- cbind(sub$x_center, sub$y_center)
  pc <- prcomp(pts, center = TRUE, scale. = FALSE)
  sdev2 <- pc$sdev^2
  elong <- ifelse(length(sdev2) < 2, Inf, sdev2[1] / max(sdev2[2], 1e-8))

  comp_info_list[[length(comp_info_list) + 1L]] <- data.table(
    comp_id = cid,
    n_bins = nrow(sub),
    elongation = elong,
    centroid_x = mean(sub$x_center),
    centroid_y = mean(sub$y_center),
    start_x = start_xy[1],
    start_y = start_xy[2],
    dir0_x = major[1],
    dir0_y = major[2],
    start_r = min(sub$r_to_center)
  )
}

comp_info <- rbindlist(comp_info_list)
setorder(comp_info, -elongation, -n_bins)

if (is.finite(N_PICK)) {
  pick_info <- comp_info[1:min(N_PICK, .N)]
} else {
  pick_info <- comp_info
}

message("Selected comp_id: ", paste(pick_info$comp_id, collapse = ", "))

############################################################
# inner-side band scan
############################################################
projection_point_list <- list()
coarse_trace_list <- list()
ray_line_list <- list()
ray_boundary_list <- list()

n_steps <- floor(MAX_SCAN_DIST / STEP_SIZE)

for (ii in seq_len(nrow(pick_info))) {
  one <- pick_info[ii]
  comp_id <- one$comp_id
  start_xy <- c(one$start_x, one$start_y)
  dir0 <- c(one$dir0_x, one$dir0_y)
  start_r <- one$start_r

  # only target points more inner than this comp
  target_inner <- copy(dt_target)
  target_inner[, r_to_center := sqrt((x_center - center_xy[1])^2 + (y_center - center_xy[2])^2)]
  target_inner <- target_inner[r_to_center < (start_r - INNER_BUFFER)]

  if (nrow(target_inner) == 0) next

  for (coarse_ang in COARSE_ANGLE_SEQ) {
    u <- norm_vec(rotate_vec(dir0, coarse_ang))
    nvec <- c(-u[2], u[1])

    # for plotting center line
    line_x <- start_xy[1] + u[1] * seq(0, MAX_SCAN_DIST, by = 1)
    line_y <- start_xy[2] + u[2] * seq(0, MAX_SCAN_DIST, by = 1)

    ray_line_list[[length(ray_line_list) + 1L]] <- data.table(
      comp_id = comp_id,
      coarse_angle = coarse_ang,
      x = line_x,
      y = line_y,
      point_order = seq_along(line_x)
    )

    # band boundaries
    for (sgn in c(-1, 1)) {
      shift <- nvec * BAND_HALF_WIDTH * sgn
      bd_x <- start_xy[1] + shift[1] + u[1] * seq(0, MAX_SCAN_DIST, by = 1)
      bd_y <- start_xy[2] + shift[2] + u[2] * seq(0, MAX_SCAN_DIST, by = 1)

      ray_boundary_list[[length(ray_boundary_list) + 1L]] <- data.table(
        comp_id = comp_id,
        coarse_angle = coarse_ang,
        boundary_side = ifelse(sgn < 0, "left", "right"),
        x = bd_x,
        y = bd_y,
        point_order = seq_along(bd_x)
      )
    }

    # project target points into current band
    v_x <- target_inner$x_center - start_xy[1]
    v_y <- target_inner$y_center - start_xy[2]

    t_proj <- v_x * u[1] + v_y * u[2]
    perp   <- abs(v_x * nvec[1] + v_y * nvec[2])

    keep <- t_proj >= 0 &
      t_proj <= MAX_SCAN_DIST &
      perp <= BAND_HALF_WIDTH

    coarse_dt <- data.table(
      comp_id = comp_id,
      coarse_angle = coarse_ang,
      step = seq_len(n_steps),
      hit_any = 0L,
      n_points = 0L
    )

    if (any(keep)) {
      kept_dt <- data.table(
        comp_id = comp_id,
        coarse_angle = coarse_ang,
        target_x = target_inner$x_center[keep],
        target_y = target_inner$y_center[keep],
        proj_t = t_proj[keep],
        perp_dist = perp[keep]
      )

      kept_dt[, step_proj := floor(proj_t / STEP_SIZE) + 1L]
      kept_dt <- kept_dt[step_proj >= 1 & step_proj <= n_steps]

      if (nrow(kept_dt) > 0) {
        projection_point_list[[length(projection_point_list) + 1L]] <- kept_dt

        hit_count_dt <- kept_dt[, .N, by = step_proj]
        setnames(hit_count_dt, c("step_proj", "N"), c("step", "n_points"))

        coarse_dt[hit_count_dt, on = "step", `:=`(
          hit_any = 1L,
          n_points = i.n_points
        )]
      }
    }

    coarse_trace_list[[length(coarse_trace_list) + 1L]] <- coarse_dt
  }
}

coarse_trace_dt <- if (length(coarse_trace_list) > 0) rbindlist(coarse_trace_list) else data.table()
projection_point_dt <- if (length(projection_point_list) > 0) rbindlist(projection_point_list) else data.table()
ray_line_dt <- if (length(ray_line_list) > 0) rbindlist(ray_line_list) else data.table()
ray_boundary_dt <- if (length(ray_boundary_list) > 0) rbindlist(ray_boundary_list) else data.table()

if (nrow(coarse_trace_dt) == 0) {
  stop("No inner-side band scan results generated.")
}

############################################################
# angle metrics / score
############################################################
angle_score_dt <- coarse_trace_dt[, .(
  n_hit_steps = sum(hit_any),
  n_hit_points = sum(n_points),
  longest_run = calc_longest_run(hit_any),
  first_hit_step = ifelse(any(hit_any == 1), min(step[hit_any == 1]), NA_integer_),
  last_hit_step  = ifelse(any(hit_any == 1), max(step[hit_any == 1]), NA_integer_),
  max_gap = calc_max_gap(hit_any),
  local_density_binary = calc_local_density_binary(hit_any, LOCAL_WINDOW_STEPS),
  local_density_count = calc_local_density_count(n_points, LOCAL_WINDOW_STEPS)
), by = .(comp_id, coarse_angle)]

angle_score_dt[, span_steps := fifelse(
  is.na(first_hit_step) | is.na(last_hit_step),
  0L,
  last_hit_step - first_hit_step + 1L
)]

angle_score_dt[, run_fill_ratio := fifelse(
  span_steps > 0,
  longest_run / span_steps,
  0
)]

angle_score_dt[, points_per_hit_step := fifelse(
  n_hit_steps > 0,
  n_hit_points / n_hit_steps,
  0
)]

angle_score_dt[, score :=
  0.8 * n_hit_steps +
  0.6 * n_hit_points +
  1.5 * longest_run +
  2.0 * local_density_binary +
  1.2 * points_per_hit_step +
  0.8 * run_fill_ratio -
  0.03 * fifelse(is.na(first_hit_step), n_steps, first_hit_step) -
  0.25 * max_gap
]

best_angle_dt <- angle_score_dt[
  order(comp_id,
        -score,
        -n_hit_points,
        -longest_run,
        -local_density_binary,
        -points_per_hit_step,
        max_gap,
        first_hit_step)
][, .SD[1], by = comp_id]

############################################################
# detailed score table
############################################################
angle_score_detail_dt <- copy(angle_score_dt)
angle_score_detail_dt[, angle_rank := frank(-score, ties.method = "first"), by = comp_id]

best_score_ref_dt <- angle_score_detail_dt[angle_rank == 1, .(
  comp_id,
  best_angle_ref = coarse_angle,
  best_score_ref = score
)]

angle_score_detail_dt <- merge(
  angle_score_detail_dt,
  best_score_ref_dt,
  by = "comp_id",
  all.x = TRUE
)

angle_score_detail_dt[, score_diff_from_best := best_score_ref - score]
angle_score_detail_dt[, score_ratio_to_best := fifelse(
  best_score_ref > 0,
  score / best_score_ref,
  NA_real_
)]
angle_score_detail_dt[, angle_diff_from_best := coarse_angle - best_angle_ref]
angle_score_detail_dt[, is_best_angle := coarse_angle == best_angle_ref]

setorder(angle_score_detail_dt, comp_id, angle_rank, coarse_angle)

fwrite(
  angle_score_detail_dt,
  paste0(OUT_PREFIX, "_angle_score_detailed.tsv"),
  sep = "\t"
)

############################################################
# best / second best per comp
############################################################
second_best_dt <- angle_score_dt[
  order(comp_id,
        -score,
        -n_hit_points,
        -longest_run,
        -local_density_binary,
        -points_per_hit_step,
        max_gap,
        first_hit_step)
][, .SD[2], by = comp_id]

best_angle_dt <- merge(
  best_angle_dt,
  second_best_dt[, .(comp_id, second_best_score = score, second_best_angle = coarse_angle)],
  by = "comp_id",
  all.x = TRUE
)

best_angle_dt[, score_margin := score - second_best_score]
best_angle_dt[, score_ratio_vs_second := fifelse(
  !is.na(second_best_score) & second_best_score > 0,
  score / second_best_score,
  NA_real_
)]

############################################################
# merge back to comp info
############################################################
pick_info <- merge(
  pick_info,
  best_angle_dt[, .(
    comp_id,
    best_angle = coarse_angle,
    best_score = score,
    best_n_hit_steps = n_hit_steps,
    best_n_hit_points = n_hit_points,
    best_longest_run = longest_run,
    best_first_hit_step = first_hit_step,
    best_last_hit_step = last_hit_step,
    best_span_steps = span_steps,
    best_run_fill_ratio = run_fill_ratio,
    best_points_per_hit_step = points_per_hit_step,
    best_local_density_binary = local_density_binary,
    best_local_density_count = local_density_count,
    best_max_gap = max_gap,
    second_best_score,
    second_best_angle,
    score_margin,
    score_ratio_vs_second
  )],
  by = "comp_id",
  all.x = TRUE
)

pick_info[, main_dir_x := mapply(function(dx, dy, a) rotate_vec(c(dx, dy), a)[1], dir0_x, dir0_y, best_angle)]
pick_info[, main_dir_y := mapply(function(dx, dy, a) rotate_vec(c(dx, dy), a)[2], dir0_x, dir0_y, best_angle)]

############################################################
# best metrics table for distributions
############################################################
best_metrics_dt <- copy(pick_info)[, .(
  comp_id,
  n_bins,
  elongation,
  best_angle,
  best_score,
  best_n_hit_steps,
  best_n_hit_points,
  best_longest_run,
  best_first_hit_step,
  best_last_hit_step,
  best_span_steps,
  best_run_fill_ratio,
  best_points_per_hit_step,
  best_local_density_binary,
  best_local_density_count,
  best_max_gap,
  second_best_score,
  second_best_angle,
  score_margin,
  score_ratio_vs_second
)]

best_metrics_dt[, highlight_group := ifelse(
  comp_id %in% HIGHLIGHT_COMPS,
  paste0("comp_", comp_id),
  "other"
)]

############################################################
# summary table
############################################################
angle_summary_dt <- copy(angle_score_dt)

############################################################
# save tables
############################################################
fwrite(pick_info,         paste0(OUT_PREFIX, "_selected_segments.tsv"), sep = "\t")
fwrite(coarse_trace_dt,   paste0(OUT_PREFIX, "_coarse_trace.tsv"), sep = "\t")
fwrite(angle_summary_dt,  paste0(OUT_PREFIX, "_angle_summary.tsv"), sep = "\t")
fwrite(angle_score_dt,    paste0(OUT_PREFIX, "_angle_score.tsv"), sep = "\t")
fwrite(best_angle_dt,     paste0(OUT_PREFIX, "_best_main_angle.tsv"), sep = "\t")
fwrite(best_metrics_dt,   paste0(OUT_PREFIX, "_best_metrics_distribution_table.tsv"), sep = "\t")

if (nrow(projection_point_dt) > 0) {
  fwrite(projection_point_dt, paste0(OUT_PREFIX, "_projected_points.tsv"), sep = "\t")
}

highlight_dt <- best_metrics_dt[highlight_group != "other"]
if (nrow(highlight_dt) > 0) {
  fwrite(highlight_dt, paste0(OUT_PREFIX, "_best_metrics_highlighted_comps.tsv"), sep = "\t")
}

############################################################
# overlay per component
# red    = all hit points across all scanned angles
# red = best-band hit points only
# yellow = original comp3 bins + PCA axis
# green  = best band center line
# grey   = best band boundaries
############################################################
plot_dt <- dt[
  cluster_raw %in% c(
    as.character(CLUSTER_LINE),
    as.character(CLUSTER_TARGET),
    as.character(CLUSTER_FUSIFORM)
  )
]

for (cid in unique(pick_info$comp_id)) {
  sub_ray_all <- ray_line_dt[comp_id == cid]
  sub_pick <- pick_info[comp_id == cid]
  best_ang <- sub_pick$best_angle[1]

  sub_best_ray <- ray_line_dt[comp_id == cid & coarse_angle == best_ang]
  sub_best_bd  <- ray_boundary_dt[comp_id == cid & coarse_angle == best_ang]
  sub_proj_all <- projection_point_dt[comp_id == cid]
  sub_proj_best <- projection_point_dt[comp_id == cid & coarse_angle == best_ang]
  sub_c3_pts <- dt_line[comp_id == cid]

  axis_dt <- sub_pick[, .(
    x1 = start_x,
    y1 = start_y,
    x2 = start_x + dir0_x * PCA_AXIS_LEN,
    y2 = start_y + dir0_y * PCA_AXIS_LEN
  )]

  best_axis_dt <- sub_pick[, .(
    x1 = start_x,
    y1 = start_y,
    x2 = start_x + main_dir_x * BEST_AXIS_LEN,
    y2 = start_y + main_dir_y * BEST_AXIS_LEN
  )]

  p <- ggplot() +
    geom_point(
      data = plot_dt[cluster_raw == as.character(CLUSTER_TARGET)],
      aes(x = x_center, y = y_center),
      size = 0.5, alpha = 1, color = "#E8E8E8"
    ) +
    geom_point(
      data = plot_dt[cluster_raw == as.character(CLUSTER_FUSIFORM)],
      aes(x = x_center, y = y_center),
      size = 0.5, alpha = 1, color = "#00D5B5"
    ) +
    geom_point(
      data = plot_dt[cluster_raw == as.character(CLUSTER_LINE)],
      aes(x = x_center, y = y_center),
      size = 0.5, alpha = 1, color = "#B8B000" 
    )

  if (nrow(sub_ray_all) > 0) {
    p <- p +
      geom_path(
        data = sub_ray_all,
        aes(x = x, y = y, group = coarse_angle),
        linewidth = 0.18, alpha = 0.12, color = "grey80"
      )
  }

  if (nrow(sub_best_bd) > 0) {
    p <- p +
      geom_path(
        data = sub_best_bd,
        aes(x = x, y = y, group = boundary_side),
        linewidth = 0.35, alpha = 0.50, color = "grey65"
      )
  }

  if (nrow(sub_proj_all) > 0) {
    p <- p +
      geom_point(
        data = sub_proj_all,
        aes(x = target_x, y = target_y),
        size = 0.5, alpha = 1, color = "red"
      )
  }

  if (nrow(sub_proj_best) > 0) {
    p <- p +
      geom_point(
        data = sub_proj_best,
        aes(x = target_x, y = target_y),
        size = 0.5, alpha = 1, color = "red"
      )
  }

  if (nrow(sub_best_ray) > 0) {
    p <- p +
      geom_path(
        data = sub_best_ray,
        aes(x = x, y = y, group = comp_id),
        linewidth = 0.95, alpha = 0.95, color = "#2319F5"
      )
  }

  p <- p +
    geom_point(
      data = sub_c3_pts,
      aes(x = x_center, y = y_center),
      size = 0.5, alpha = 1, color = "yellow"
    ) +
    geom_segment(
      data = axis_dt,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      linewidth = 0.75, color = "yellow"
    ) +
    geom_segment(
      data = best_axis_dt,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      linewidth = 0.95, color = "#2319F5"
    ) +
    geom_point(
      data = sub_pick,
      aes(x = start_x, y = start_y),
      size = 2.3, color = "magenta"
    ) +
    coord_equal() +
    scale_y_reverse() +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "black", color = "black"),
      panel.background = element_rect(fill = "black", color = "black")
    ) +
    ggtitle(paste0(
      "comp_", cid,
      " | best angle = ", best_ang, "°",
      " | score=", round(sub_pick$best_score[1], 2),
      " | hit_pts=", sub_pick$best_n_hit_points[1],
      " | hit_steps=", sub_pick$best_n_hit_steps[1],
      " | longest=", sub_pick$best_longest_run[1],
      " | margin=", round(sub_pick$score_margin[1], 2)
    )) +
    theme(
      plot.title = element_text(color = "white", hjust = 0.5, size = 12)
    )

  ggsave(paste0(OUT_PREFIX, "_overlay_comp_", cid, ".pdf"), p, width = 8, height = 8)
  ggsave(paste0(OUT_PREFIX, "_overlay_comp_", cid, ".png"), p, width = 8, height = 8, dpi = 300)
}

############################################################
# heatmap per component
############################################################
for (cid in unique(coarse_trace_dt$comp_id)) {
  sub <- coarse_trace_dt[comp_id == cid]
  best_ang <- pick_info[comp_id == cid]$best_angle[1]

  p_heat <- ggplot(sub, aes(x = coarse_angle, y = step, fill = n_points)) +
    geom_tile() +
    geom_vline(xintercept = best_ang, linewidth = 0.8, color = "#2319F5") +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("comp_", cid, " : inner-side band scan heatmap"),
      x = "Angle relative to inward direction (deg)",
      y = "Step along band",
      fill = "n cluster5 points"
    )

  ggsave(paste0(OUT_PREFIX, "_heatmap_comp_", cid, ".pdf"), p_heat, width = 6, height = 8)
  ggsave(paste0(OUT_PREFIX, "_heatmap_comp_", cid, ".png"), p_heat, width = 6, height = 8, dpi = 300)
}

############################################################
# distribution plots
############################################################
hist_vars <- c(
  "best_score",
  "best_n_hit_steps",
  "best_n_hit_points",
  "best_longest_run",
  "best_run_fill_ratio",
  "best_points_per_hit_step",
  "best_local_density_binary",
  "best_local_density_count",
  "best_max_gap",
  "score_margin",
  "score_ratio_vs_second"
)

for (vv in hist_vars) {
  save_hist_plot(best_metrics_dt, vv, OUT_PREFIX, bins = DIAG_BINS)
}

scatter_pairs <- list(
  c("best_longest_run", "best_n_hit_points"),
  c("best_n_hit_steps", "best_n_hit_points"),
  c("best_score", "best_n_hit_points"),
  c("best_run_fill_ratio", "best_points_per_hit_step"),
  c("best_local_density_binary", "best_n_hit_points"),
  c("score_margin", "best_score")
)

for (pp in scatter_pairs) {
  save_scatter_plot(best_metrics_dt, pp[1], pp[2], OUT_PREFIX)
}

############################################################
# final messages
############################################################
message("Done.")
message("Output files:")
message("  ", paste0(OUT_PREFIX, "_selected_segments.tsv"))
message("  ", paste0(OUT_PREFIX, "_coarse_trace.tsv"))
message("  ", paste0(OUT_PREFIX, "_angle_summary.tsv"))
message("  ", paste0(OUT_PREFIX, "_angle_score.tsv"))
message("  ", paste0(OUT_PREFIX, "_angle_score_detailed.tsv"))
message("  ", paste0(OUT_PREFIX, "_best_main_angle.tsv"))
message("  ", paste0(OUT_PREFIX, "_best_metrics_distribution_table.tsv"))
message("  ", paste0(OUT_PREFIX, "_best_metrics_highlighted_comps.tsv"))

if (nrow(projection_point_dt) > 0) {
  message("  ", paste0(OUT_PREFIX, "_projected_points.tsv"))
}

for (vv in hist_vars) {
  message("  ", paste0(OUT_PREFIX, "_hist_", vv, ".pdf"))
  message("  ", paste0(OUT_PREFIX, "_hist_", vv, ".png"))
}

for (pp in scatter_pairs) {
  message("  ", paste0(OUT_PREFIX, "_scatter_", pp[1], "_vs_", pp[2], ".pdf"))
  message("  ", paste0(OUT_PREFIX, "_scatter_", pp[1], "_vs_", pp[2], ".png"))
}

for (cid in unique(pick_info$comp_id)) {
  message("  ", paste0(OUT_PREFIX, "_overlay_comp_", cid, ".pdf"))
  message("  ", paste0(OUT_PREFIX, "_overlay_comp_", cid, ".png"))
  message("  ", paste0(OUT_PREFIX, "_heatmap_comp_", cid, ".pdf"))
  message("  ", paste0(OUT_PREFIX, "_heatmap_comp_", cid, ".png"))
}


top_comp <- pick_info[which.max(best_score), comp_id]
top_comps <- pick_info[order(-best_score), comp_id][1:20]

data.table(
  comp_id = top_comps,
  overlay_png = paste0(OUT_PREFIX, "_overlay_comp_", top_comps, ".png"),
  heatmap_png = paste0(OUT_PREFIX, "_heatmap_comp_", top_comps, ".png")
)
