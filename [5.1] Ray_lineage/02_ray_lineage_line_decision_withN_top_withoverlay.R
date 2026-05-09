#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setwd("/home/woodydrylab/FileShare/20260121_Xenium/output-XETG00360__0079049__BIO1_TISSUE_1__20260115__224443/grid05um_out/kmeans_k10_raw_out/ray_match3")

############################################################
# user settings
############################################################
INPUT_TSV  <- "../bin_metadata_with_cluster_raw.tsv"
OUT_PREFIX <- "ray_scan_cluster3_to_5_inner_band_permB_topK"

CLUSTER_LINE   <- 3
CLUSTER_TARGET <- 5

# scan angle
ANGLE_STEP <- 5
COARSE_ANGLE_SEQ <- seq(-90, 90, by = ANGLE_STEP)

# only scan inner side of comp3
INNER_BUFFER <- 0

# distance / band
MAX_SCAN_DIST   <- 500
STEP_SIZE       <- 5
BAND_HALF_WIDTH <- 8

# component split
NEIGHBOR_FACTOR <- 1.5
MIN_COMP_SIZE   <- 3
N_PICK          <- Inf

# overlay line lengths
PCA_AXIS_LEN  <- 0
BEST_AXIS_LEN <- 120

# local density window on step bins
LOCAL_WINDOW_STEPS <- 6

# target local direction estimation
TARGET_DIR_RADIUS_FACTOR <- 3

# angle agreement between scanning direction and target local direction
TARGET_ANGLE_TOL <- 30

# score weight for direction agreement
DIR_WEIGHT_MULT <- 1.5

# permutation settings
N_PERM <- 20
SEED <- 123

# top K comparison
TOP_K <- 30

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

angle_diff_axis_deg <- function(dx1, dy1, dx2, dy2) {
  dotv <- dx1 * dx2 + dy1 * dy2
  dotv <- pmax(-1, pmin(1, dotv))
  acos(abs(dotv)) * 180 / pi
}

estimate_local_direction <- function(dt_target, bin_size, radius_factor = 3) {
  dt_target <- copy(dt_target)
  dt_target[, c("local_dir_x", "local_dir_y") := .(NA_real_, NA_real_)]

  if (nrow(dt_target) == 0) return(dt_target)

  dir_radius <- bin_size * radius_factor
  coords <- as.matrix(dt_target[, .(x_center, y_center)])

  for (i in seq_len(nrow(dt_target))) {
    dx <- coords[, 1] - coords[i, 1]
    dy <- coords[, 2] - coords[i, 2]
    d2 <- dx^2 + dy^2

    nb_idx <- which(d2 <= dir_radius^2)

    if (length(nb_idx) >= 3) {
      v <- get_major_axis(dt_target$x_center[nb_idx], dt_target$y_center[nb_idx])
      dt_target[i, `:=`(local_dir_x = v[1], local_dir_y = v[2])]
    }
  }

  dt_target
}

scan_one_comp <- function(one, target_inner, coarse_angles, max_scan_dist, step_size,
                          band_half_width, local_window_steps, target_angle_tol,
                          dir_weight_mult) {

  comp_id <- one$comp_id
  start_xy <- c(one$start_x, one$start_y)
  dir0 <- c(one$dir0_x, one$dir0_y)

  n_steps <- floor(max_scan_dist / step_size)

  angle_rows <- vector("list", length(coarse_angles))

  for (kk in seq_along(coarse_angles)) {
    coarse_ang <- coarse_angles[kk]
    u <- norm_vec(rotate_vec(dir0, coarse_ang))
    nvec <- c(-u[2], u[1])

    v_x <- target_inner$x_center - start_xy[1]
    v_y <- target_inner$y_center - start_xy[2]

    t_proj <- v_x * u[1] + v_y * u[2]
    perp   <- abs(v_x * nvec[1] + v_y * nvec[2])

    keep <- t_proj >= 0 &
      t_proj <= max_scan_dist &
      perp <= band_half_width &
      !is.na(target_inner$local_dir_x) &
      !is.na(target_inner$local_dir_y)

    coarse_dt <- data.table(
      step = seq_len(n_steps),
      hit_any = 0L,
      n_points = 0L,
      mean_dir_weight = 0
    )

    if (any(keep)) {
      kept_dt <- data.table(
        proj_t = t_proj[keep],
        perp_dist = perp[keep],
        dir_x = target_inner$local_dir_x[keep],
        dir_y = target_inner$local_dir_y[keep]
      )

      kept_dt[, angle_diff := angle_diff_axis_deg(dir_x, dir_y, u[1], u[2])]
      kept_dt <- kept_dt[angle_diff <= target_angle_tol]

      if (nrow(kept_dt) > 0) {
        kept_dt[, dir_weight := cos(angle_diff * pi / 180)]
        kept_dt[, step_proj := floor(proj_t / step_size) + 1L]
        kept_dt <- kept_dt[step_proj >= 1 & step_proj <= n_steps]

        if (nrow(kept_dt) > 0) {
          hit_count_dt <- kept_dt[, .(
            n_points = .N,
            mean_dir_weight = mean(dir_weight, na.rm = TRUE)
          ), by = step_proj]
          setnames(hit_count_dt, "step_proj", "step")

          coarse_dt[hit_count_dt, on = "step", `:=`(
            hit_any = 1L,
            n_points = i.n_points,
            mean_dir_weight = i.mean_dir_weight
          )]
        }
      }
    }

    n_hit_steps <- sum(coarse_dt$hit_any)
    n_hit_points <- sum(coarse_dt$n_points)
    longest_run <- calc_longest_run(coarse_dt$hit_any)
    first_hit_step <- if (any(coarse_dt$hit_any == 1)) min(coarse_dt$step[coarse_dt$hit_any == 1]) else NA_integer_
    last_hit_step  <- if (any(coarse_dt$hit_any == 1)) max(coarse_dt$step[coarse_dt$hit_any == 1]) else NA_integer_
    max_gap <- calc_max_gap(coarse_dt$hit_any)
    local_density_binary <- calc_local_density_binary(coarse_dt$hit_any, local_window_steps)
    local_density_count  <- calc_local_density_count(coarse_dt$n_points, local_window_steps)
    mean_dir_weight <- mean(coarse_dt$mean_dir_weight[coarse_dt$hit_any == 1], na.rm = TRUE)
    if (!is.finite(mean_dir_weight)) mean_dir_weight <- 0

    span_steps <- if (is.na(first_hit_step) || is.na(last_hit_step)) 0L else (last_hit_step - first_hit_step + 1L)
    run_fill_ratio <- if (span_steps > 0) longest_run / span_steps else 0
    points_per_hit_step <- if (n_hit_steps > 0) n_hit_points / n_hit_steps else 0

    score <- 
      0.8 * n_hit_steps +
      0.6 * n_hit_points +
      1.5 * longest_run +
      2.0 * local_density_binary +
      1.2 * points_per_hit_step +
      0.8 * run_fill_ratio -
      0.03 * ifelse(is.na(first_hit_step), n_steps, first_hit_step) -
      0.25 * max_gap +
      dir_weight_mult * mean_dir_weight

    angle_rows[[kk]] <- data.table(
      comp_id = comp_id,
      coarse_angle = coarse_ang,
      n_hit_steps = n_hit_steps,
      n_hit_points = n_hit_points,
      longest_run = longest_run,
      first_hit_step = first_hit_step,
      last_hit_step = last_hit_step,
      max_gap = max_gap,
      local_density_binary = local_density_binary,
      local_density_count = local_density_count,
      mean_dir_weight = mean_dir_weight,
      span_steps = span_steps,
      run_fill_ratio = run_fill_ratio,
      points_per_hit_step = points_per_hit_step,
      score = score
    )
  }

  angle_dt <- rbindlist(angle_rows)

  setorder(angle_dt,
           -score,
           -n_hit_points,
           -longest_run,
           -mean_dir_weight,
           -local_density_binary,
           -points_per_hit_step,
           max_gap,
           first_hit_step)

  best_dt <- angle_dt[1]
  second_dt <- if (nrow(angle_dt) >= 2) angle_dt[2] else data.table(
    second_best_score = NA_real_,
    second_best_angle = NA_real_
  )

  best_dt[, second_best_score := if (nrow(angle_dt) >= 2) second_dt$score else NA_real_]
  best_dt[, second_best_angle := if (nrow(angle_dt) >= 2) second_dt$coarse_angle else NA_real_]
  best_dt[, score_margin := score - second_best_score]
  best_dt[, score_ratio_vs_second := fifelse(!is.na(second_best_score) & second_best_score > 0, score / second_best_score, NA_real_)]

  best_dt
}

topk_by_perm <- function(dt, k, dataset_name) {
  dt2 <- copy(dt)
  setorder(dt2, perm_id, -score, comp_id)
  out <- dt2[, head(.SD, min(k, .N)), by = perm_id]
  out[, dataset_type := dataset_name]
  out[, topk_rank := seq_len(.N), by = perm_id]
  out
}

collect_overlay_for_best <- function(one, target_inner, best_angle,
                                     max_scan_dist, step_size,
                                     band_half_width, target_angle_tol) {
  start_xy <- c(one$start_x, one$start_y)
  dir0 <- c(one$dir0_x, one$dir0_y)

  u <- norm_vec(rotate_vec(dir0, best_angle))
  nvec <- c(-u[2], u[1])

  v_x <- target_inner$x_center - start_xy[1]
  v_y <- target_inner$y_center - start_xy[2]

  t_proj <- v_x * u[1] + v_y * u[2]
  perp   <- abs(v_x * nvec[1] + v_y * nvec[2])

  keep <- t_proj >= 0 &
    t_proj <= max_scan_dist &
    perp <= band_half_width &
    !is.na(target_inner$local_dir_x) &
    !is.na(target_inner$local_dir_y)

  kept_dt <- data.table()
  if (any(keep)) {
    kept_dt <- data.table(
      x_center = target_inner$x_center[keep],
      y_center = target_inner$y_center[keep],
      proj_t = t_proj[keep],
      perp_dist = perp[keep],
      dir_x = target_inner$local_dir_x[keep],
      dir_y = target_inner$local_dir_y[keep]
    )
    kept_dt[, angle_diff := angle_diff_axis_deg(dir_x, dir_y, u[1], u[2])]
    kept_dt <- kept_dt[angle_diff <= target_angle_tol]
  }

  main_line_dt <- data.table(
    x = start_xy[1] + u[1] * (0:max_scan_dist),
    y = start_xy[2] + u[2] * (0:max_scan_dist),
    point_order = 0:max_scan_dist
  )

  boundary_list <- list()
  for (bd in c(best_angle - ANGLE_STEP / 2, best_angle + ANGLE_STEP / 2)) {
    u_bd <- norm_vec(rotate_vec(dir0, bd))
    boundary_list[[length(boundary_list) + 1L]] <- data.table(
      boundary_angle = bd,
      x = start_xy[1] + u_bd[1] * (0:max_scan_dist),
      y = start_xy[2] + u_bd[2] * (0:max_scan_dist),
      point_order = 0:max_scan_dist
    )
  }
  boundary_dt <- rbindlist(boundary_list)

  axis_dt <- data.table(
    x1 = one$start_x,
    y1 = one$start_y,
    x2 = one$start_x + one$dir0_x * PCA_AXIS_LEN,
    y2 = one$start_y + one$dir0_y * PCA_AXIS_LEN
  )

  best_axis_dt <- data.table(
    x1 = one$start_x,
    y1 = one$start_y,
    x2 = one$start_x + u[1] * BEST_AXIS_LEN,
    y2 = one$start_y + u[2] * BEST_AXIS_LEN
  )

  list(
    kept_dt = kept_dt,
    main_line_dt = main_line_dt,
    boundary_dt = boundary_dt,
    axis_dt = axis_dt,
    best_axis_dt = best_axis_dt
  )
}

plot_overlay_one <- function(plot_dt_bg, one, best_angle, overlay_obj,
                             out_prefix, title_text,
                             point_color = "red") {

  p <- ggplot() +
    geom_point(
      data = plot_dt_bg[cluster_raw == as.character(CLUSTER_TARGET)],
      aes(x = x_center, y = y_center),
      size = 0.15, alpha = 0.20, color = "#00D5B5"
    ) +
    geom_point(
      data = plot_dt_bg[cluster_raw == as.character(CLUSTER_LINE)],
      aes(x = x_center, y = y_center),
      size = 0.15, alpha = 0.75, color = "#B8B000"
    ) +
    geom_path(
      data = overlay_obj$boundary_dt,
      aes(x = x, y = y, group = boundary_angle),
      linewidth = 0.20, alpha = 0.20, color = "grey70"
    ) +
    geom_path(
      data = overlay_obj$main_line_dt,
      aes(x = x, y = y),
      linewidth = 0.90, alpha = 0.95, color = "limegreen"
    )

  if (nrow(overlay_obj$kept_dt) > 0) {
    p <- p +
      geom_point(
        data = overlay_obj$kept_dt,
        aes(x = x_center, y = y_center),
        size = 0.35, alpha = 0.90, color = point_color
      )
  }

  p <- p +
    geom_segment(
      data = overlay_obj$axis_dt,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      linewidth = 0.70, color = "yellow"
    ) +
    geom_segment(
      data = overlay_obj$best_axis_dt,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      linewidth = 0.95, color = "limegreen"
    ) +
    geom_point(
      data = one,
      aes(x = start_x, y = start_y),
      size = 2.4, color = "magenta"
    ) +
    coord_equal() +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "black", color = "black"),
      panel.background = element_rect(fill = "black", color = "black")
    ) +
    ggtitle(title_text) +
    theme(
      plot.title = element_text(color = "white", hjust = 0.5, size = 14)
    )

  ggsave(paste0(out_prefix, ".pdf"), p, width = 8, height = 8)
  ggsave(paste0(out_prefix, ".png"), p, width = 8, height = 8, dpi = 300)
}
############################################################
# read data
############################################################
set.seed(SEED)

dt <- fread(INPUT_TSV)
dt[, cluster_raw := as.character(cluster_raw)]

required_cols <- c("x_center", "y_center", "cluster_raw")
miss <- setdiff(required_cols, names(dt))
if (length(miss) > 0) {
  stop("Missing required columns: ", paste(miss, collapse = ", "))
}

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

############################################################
# global center and radial distance
############################################################
center_xy <- c(mean(dt$x_center), mean(dt$y_center))
dt[, r_to_center := sqrt((x_center - center_xy[1])^2 + (y_center - center_xy[2])^2)]

# split AFTER r_to_center exists
dt_line   <- copy(dt[cluster_raw == as.character(CLUSTER_LINE)])
dt_target <- copy(dt[cluster_raw == as.character(CLUSTER_TARGET)])

if (nrow(dt_line) == 0) stop("No cluster 3 bins found.")
if (nrow(dt_target) == 0) stop("No cluster 5 bins found.")

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
if (nrow(comp_stat) == 0) stop("No cluster 3 components passed MIN_COMP_SIZE.")

############################################################
# define each comp3
############################################################
comp_info_list <- list()

for (ii in seq_len(nrow(comp_stat))) {
  cid <- comp_stat$comp_id[ii]
  sub <- dt_line[comp_id == cid]

  start_bin <- sub[which.min(r_to_center)][1]
  start_xy <- c(start_bin$x_center, start_bin$y_center)

  major <- get_major_axis(sub$x_center, sub$y_center)
  to_center <- norm_vec(center_xy - start_xy)

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

pick_info <- rbindlist(comp_info_list)
setorder(pick_info, -elongation, -n_bins)

if (is.finite(N_PICK)) {
  pick_info <- pick_info[1:min(N_PICK, .N)]
}

message("Selected comp_id count: ", nrow(pick_info))

############################################################
# run real data
############################################################
message("Running real-data scan ...")

real_best_list <- vector("list", nrow(pick_info))

for (ii in seq_len(nrow(pick_info))) {
  one <- pick_info[ii]

  target_inner <- copy(dt_target[r_to_center < (one$start_r - INNER_BUFFER)])
  if (nrow(target_inner) == 0) next

  target_inner <- estimate_local_direction(
    target_inner,
    bin_size = bin_size,
    radius_factor = TARGET_DIR_RADIUS_FACTOR
  )

  real_best_list[[ii]] <- scan_one_comp(
    one = one,
    target_inner = target_inner,
    coarse_angles = COARSE_ANGLE_SEQ,
    max_scan_dist = MAX_SCAN_DIST,
    step_size = STEP_SIZE,
    band_half_width = BAND_HALF_WIDTH,
    local_window_steps = LOCAL_WINDOW_STEPS,
    target_angle_tol = TARGET_ANGLE_TOL,
    dir_weight_mult = DIR_WEIGHT_MULT
  )
}

real_best_dt <- rbindlist(real_best_list, fill = TRUE)
real_best_dt[, perm_id := 0L]
real_best_dt[, dataset_type := "real"]

############################################################
# permutation B
############################################################
message("Running permutation B negative controls ...")

perm_best_all <- list()

for (pp in seq_len(N_PERM)) {
  message("  permutation ", pp, "/", N_PERM)

  perm_best_list <- vector("list", nrow(pick_info))

  for (ii in seq_len(nrow(pick_info))) {
    one <- pick_info[ii]

    real_target_inner <- copy(dt_target[r_to_center < (one$start_r - INNER_BUFFER)])
    n_target <- nrow(real_target_inner)

    if (n_target == 0) next

    inner_pool_all <- copy(dt[r_to_center < (one$start_r - INNER_BUFFER), .(x_center, y_center)])

    if (nrow(inner_pool_all) < n_target || nrow(inner_pool_all) < 3) next

    samp_idx <- sample.int(nrow(inner_pool_all), size = n_target, replace = FALSE)
    target_inner_perm <- copy(inner_pool_all[samp_idx])

    target_inner_perm <- estimate_local_direction(
      target_inner_perm,
      bin_size = bin_size,
      radius_factor = TARGET_DIR_RADIUS_FACTOR
    )

    perm_best_list[[ii]] <- scan_one_comp(
      one = one,
      target_inner = target_inner_perm,
      coarse_angles = COARSE_ANGLE_SEQ,
      max_scan_dist = MAX_SCAN_DIST,
      step_size = STEP_SIZE,
      band_half_width = BAND_HALF_WIDTH,
      local_window_steps = LOCAL_WINDOW_STEPS,
      target_angle_tol = TARGET_ANGLE_TOL,
      dir_weight_mult = DIR_WEIGHT_MULT
    )
  }

  perm_best_dt <- rbindlist(perm_best_list, fill = TRUE)
  perm_best_dt[, perm_id := pp]
  perm_best_dt[, dataset_type := "permB"]

  perm_best_all[[pp]] <- perm_best_dt
}

perm_best_all_dt <- rbindlist(perm_best_all, fill = TRUE)

############################################################
# TOP K comparison
############################################################
message("Building top-K comparison ...")

real_topk_dt <- topk_by_perm(real_best_dt, TOP_K, "real_topK")
perm_topk_dt <- topk_by_perm(perm_best_all_dt, TOP_K, "permB_topK")

real_topk_dt[, perm_group := "real"]
perm_topk_dt[, perm_group := "permB"]

topk_all_dt <- rbindlist(list(real_topk_dt, perm_topk_dt), fill = TRUE)

############################################################
# save tables
############################################################
fwrite(real_best_dt,     paste0(OUT_PREFIX, "_real_best_metrics.tsv"), sep = "\t")
fwrite(perm_best_all_dt, paste0(OUT_PREFIX, "_permB_best_metrics.tsv"), sep = "\t")
fwrite(real_topk_dt,     paste0(OUT_PREFIX, "_real_topK_metrics.tsv"), sep = "\t")
fwrite(perm_topk_dt,     paste0(OUT_PREFIX, "_permB_topK_metrics.tsv"), sep = "\t")
fwrite(topk_all_dt,      paste0(OUT_PREFIX, "_real_vs_permB_topK_metrics.tsv"), sep = "\t")

############################################################
# overlay output dirs
############################################################
dir.create(paste0(OUT_PREFIX, "_overlay_real_topK"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(OUT_PREFIX, "_overlay_permB_topK"), showWarnings = FALSE, recursive = TRUE)

############################################################
# real topK overlay
############################################################
message("Drawing real-data topK overlays ...")

plot_dt_bg <- dt[cluster_raw %in% c(as.character(CLUSTER_LINE), as.character(CLUSTER_TARGET))]

real_topk_draw_dt <- merge(
  real_topk_dt,
  pick_info,
  by = "comp_id",
  all.x = TRUE
)

for (ii in seq_len(nrow(real_topk_draw_dt))) {
  one <- real_topk_draw_dt[ii]

  target_inner <- copy(dt_target[r_to_center < (one$start_r - INNER_BUFFER)])
  if (nrow(target_inner) == 0) next

  target_inner <- estimate_local_direction(
    target_inner,
    bin_size = bin_size,
    radius_factor = TARGET_DIR_RADIUS_FACTOR
  )

  overlay_obj <- collect_overlay_for_best(
    one = one,
    target_inner = target_inner,
    best_angle = one$coarse_angle,
    max_scan_dist = MAX_SCAN_DIST,
    step_size = STEP_SIZE,
    band_half_width = BAND_HALF_WIDTH,
    target_angle_tol = TARGET_ANGLE_TOL
  )

  out_prefix_one <- file.path(
    paste0(OUT_PREFIX, "_overlay_real_topK"),
    sprintf("real_topK_rank%03d_comp%03d_angle%+03d", one$topk_rank, one$comp_id, as.integer(one$coarse_angle))
  )

  ttl <- paste0(
    "REAL top", TOP_K,
    " rank=", one$topk_rank,
    " | comp_", one$comp_id,
    " | angle=", one$coarse_angle,
    " | score=", round(one$score, 2)
  )

  plot_overlay_one(
    plot_dt_bg = plot_dt_bg,
    one = one,
    best_angle = one$coarse_angle,
    overlay_obj = overlay_obj,
    out_prefix = out_prefix_one,
    title_text = ttl,
    point_color = "red"
  )
}

############################################################
# permutation-B topK overlay
# re-run permutations with same seed/order to reproduce sampled points
############################################# ###############
message("Drawing permutation-B topK overlays ...")

perm_topk_need_dt <- merge(
  perm_topk_dt[, .(perm_id, comp_id, coarse_angle, topk_rank, score)],
  pick_info,
  by = "comp_id",
  all.x = TRUE
)

set.seed(SEED)

for (pp in seq_len(N_PERM)) {
  message("  overlay permutation ", pp, "/", N_PERM)

  need_this_perm <- perm_topk_need_dt[perm_id == pp]
  if (nrow(need_this_perm) == 0) {
    for (ii in seq_len(nrow(pick_info))) {
      one <- pick_info[ii]
      real_target_inner <- copy(dt_target[r_to_center < (one$start_r - INNER_BUFFER)])
      n_target <- nrow(real_target_inner)
      if (n_target == 0) next
      inner_pool_all <- copy(dt[r_to_center < (one$start_r - INNER_BUFFER), .(x_center, y_center)])
      if (nrow(inner_pool_all) < n_target || nrow(inner_pool_all) < 3) next
      samp_idx <- sample.int(nrow(inner_pool_all), size = n_target, replace = FALSE)
    }
    next
  }

  for (ii in seq_len(nrow(pick_info))) {
    one <- pick_info[ii]

    real_target_inner <- copy(dt_target[r_to_center < (one$start_r - INNER_BUFFER)])
    n_target <- nrow(real_target_inner)
    if (n_target == 0) next

    inner_pool_all <- copy(dt[r_to_center < (one$start_r - INNER_BUFFER), .(x_center, y_center)])
    if (nrow(inner_pool_all) < n_target || nrow(inner_pool_all) < 3) next

    samp_idx <- sample.int(nrow(inner_pool_all), size = n_target, replace = FALSE)
    target_inner_perm <- copy(inner_pool_all[samp_idx])

    target_inner_perm <- estimate_local_direction(
      target_inner_perm,
      bin_size = bin_size,
      radius_factor = TARGET_DIR_RADIUS_FACTOR
    )

    need_row <- need_this_perm[comp_id == one$comp_id]
    if (nrow(need_row) == 0) next

    best_angle <- need_row$coarse_angle[1]

    overlay_obj <- collect_overlay_for_best(
      one = one,
      target_inner = target_inner_perm,
      best_angle = best_angle,
      max_scan_dist = MAX_SCAN_DIST,
      step_size = STEP_SIZE,
      band_half_width = BAND_HALF_WIDTH,
      target_angle_tol = TARGET_ANGLE_TOL
    )

    out_prefix_one <- file.path(
      paste0(OUT_PREFIX, "_overlay_permB_topK"),
      sprintf("perm%03d_topK_rank%03d_comp%03d_angle%+03d",
              pp, need_row$topk_rank[1], one$comp_id, as.integer(best_angle))
    )

    ttl <- paste0(
      "permB top", TOP_K,
      " | perm=", pp,
      " | rank=", need_row$topk_rank[1],
      " | comp_", one$comp_id,
      " | angle=", best_angle,
      " | score=", round(need_row$score[1], 2)
    )

    plot_overlay_one(
      plot_dt_bg = plot_dt_bg,
      one = one,
      best_angle = best_angle,
      overlay_obj = overlay_obj,
      out_prefix = out_prefix_one,
      title_text = ttl,
      point_color = "orange"
    )
  }
}

############################################################
# topK score distribution
############################################################
p_topk_score <- ggplot(topk_all_dt, aes(x = score, color = perm_group, fill = perm_group)) +
  geom_density(alpha = 0.20) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Real vs permutation-B top ", TOP_K, ": best_score distribution"),
    x = "best_score",
    y = "Density"
  )

ggsave(paste0(OUT_PREFIX, "_real_vs_permB_topK_score_density.pdf"), p_topk_score, width = 7, height = 5)
ggsave(paste0(OUT_PREFIX, "_real_vs_permB_topK_score_density.png"), p_topk_score, width = 7, height = 5, dpi = 300)

############################################################
# topK margin distribution
############################################################
p_topk_margin <- ggplot(topk_all_dt, aes(x = score_margin, color = perm_group, fill = perm_group)) +
  geom_density(alpha = 0.20) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Real vs permutation-B top ", TOP_K, ": score_margin distribution"),
    x = "score_margin",
    y = "Density"
  )

ggsave(paste0(OUT_PREFIX, "_real_vs_permB_topK_margin_density.pdf"), p_topk_margin, width = 7, height = 5)
ggsave(paste0(OUT_PREFIX, "_real_vs_permB_topK_margin_density.png"), p_topk_margin, width = 7, height = 5, dpi = 300)

############################################################
# topK score vs margin
############################################################
p_topk_scatter <- ggplot(topk_all_dt, aes(x = score, y = score_margin, color = perm_group)) +
  geom_point(alpha = 0.55, size = 2) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Real vs permutation-B top ", TOP_K, ": best_score vs score_margin"),
    x = "best_score",
    y = "score_margin"
  )

ggsave(paste0(OUT_PREFIX, "_real_vs_permB_topK_score_vs_margin.pdf"), p_topk_scatter, width = 7, height = 5)
ggsave(paste0(OUT_PREFIX, "_real_vs_permB_topK_score_vs_margin.png"), p_topk_scatter, width = 7, height = 5, dpi = 300)

############################################################
# permutation topK summary statistics
############################################################
perm_topk_summary_dt <- perm_topk_dt[, .(
  mean_topK_score = mean(score, na.rm = TRUE),
  median_topK_score = median(score, na.rm = TRUE),
  mean_topK_margin = mean(score_margin, na.rm = TRUE),
  median_topK_margin = median(score_margin, na.rm = TRUE),
  kth_score = min(score, na.rm = TRUE)
), by = perm_id]

real_topk_summary_dt <- real_topk_dt[, .(
  mean_topK_score = mean(score, na.rm = TRUE),
  median_topK_score = median(score, na.rm = TRUE),
  mean_topK_margin = mean(score_margin, na.rm = TRUE),
  median_topK_margin = median(score_margin, na.rm = TRUE),
  kth_score = min(score, na.rm = TRUE)
), by = perm_id]

fwrite(perm_topk_summary_dt, paste0(OUT_PREFIX, "_permB_topK_summary.tsv"), sep = "\t")
fwrite(real_topk_summary_dt, paste0(OUT_PREFIX, "_real_topK_summary.tsv"), sep = "\t")

############################################################
# empirical p-values for topK summaries
############################################################
real_mean_topK_score <- real_topk_summary_dt$mean_topK_score[1]
real_median_topK_score <- real_topk_summary_dt$median_topK_score[1]
real_mean_topK_margin <- real_topk_summary_dt$mean_topK_margin[1]
real_kth_score <- real_topk_summary_dt$kth_score[1]

p_mean_score <- (sum(perm_topk_summary_dt$mean_topK_score >= real_mean_topK_score, na.rm = TRUE) + 1) / (nrow(perm_topk_summary_dt) + 1)
p_median_score <- (sum(perm_topk_summary_dt$median_topK_score >= real_median_topK_score, na.rm = TRUE) + 1) / (nrow(perm_topk_summary_dt) + 1)
p_mean_margin <- (sum(perm_topk_summary_dt$mean_topK_margin >= real_mean_topK_margin, na.rm = TRUE) + 1) / (nrow(perm_topk_summary_dt) + 1)
p_kth_score <- (sum(perm_topk_summary_dt$kth_score >= real_kth_score, na.rm = TRUE) + 1) / (nrow(perm_topk_summary_dt) + 1)

topk_empirical_p_dt <- data.table(
  TOP_K = TOP_K,
  N_PERM = N_PERM,
  real_mean_topK_score = real_mean_topK_score,
  perm_mean_topK_score = mean(perm_topk_summary_dt$mean_topK_score, na.rm = TRUE),
  p_mean_topK_score = p_mean_score,
  real_median_topK_score = real_median_topK_score,
  perm_median_topK_score = mean(perm_topk_summary_dt$median_topK_score, na.rm = TRUE),
  p_median_topK_score = p_median_score,
  real_mean_topK_margin = real_mean_topK_margin,
  perm_mean_topK_margin = mean(perm_topk_summary_dt$mean_topK_margin, na.rm = TRUE),
  p_mean_topK_margin = p_mean_margin,
  real_kth_score = real_kth_score,
  perm_mean_kth_score = mean(perm_topk_summary_dt$kth_score, na.rm = TRUE),
  p_kth_score = p_kth_score
)

fwrite(topk_empirical_p_dt, paste0(OUT_PREFIX, "_topK_empirical_p.tsv"), sep = "\t")

############################################################
# histogram of permutation topK summaries with real line
############################################################
plot_perm_hist <- function(dt, varname, real_value, filename_prefix, title_text) {
  p <- ggplot(dt, aes_string(x = varname)) +
    geom_histogram(bins = 30, fill = "grey70", color = "white") +
    geom_vline(xintercept = real_value, color = "red", linewidth = 1) +
    theme_bw(base_size = 12) +
    labs(
      title = title_text,
      x = varname,
      y = "Count"
    )

  ggsave(paste0(filename_prefix, ".pdf"), p, width = 7, height = 5)
  ggsave(paste0(filename_prefix, ".png"), p, width = 7, height = 5, dpi = 300)
}

plot_perm_hist(
  perm_topk_summary_dt,
  "mean_topK_score",
  real_mean_topK_score,
  paste0(OUT_PREFIX, "_permB_topK_mean_score_hist"),
  paste0("Permutation-B top ", TOP_K, ": mean topK score")
)

plot_perm_hist(
  perm_topk_summary_dt,
  "mean_topK_margin",
  real_mean_topK_margin,
  paste0(OUT_PREFIX, "_permB_topK_mean_margin_hist"),
  paste0("Permutation-B top ", TOP_K, ": mean topK margin")
)

plot_perm_hist(
  perm_topk_summary_dt,
  "kth_score",
  real_kth_score,
  paste0(OUT_PREFIX, "_permB_topK_kth_score_hist"),
  paste0("Permutation-B top ", TOP_K, ": K-th score")
)

############################################################
# final messages
############################################################
message("Done.")
message("Top-K comparison finished with TOP_K = ", TOP_K)
message("Output files:")
message("  ", paste0(OUT_PREFIX, "_real_best_metrics.tsv"))
message("  ", paste0(OUT_PREFIX, "_permB_best_metrics.tsv"))
message("  ", paste0(OUT_PREFIX, "_real_topK_metrics.tsv"))
message("  ", paste0(OUT_PREFIX, "_permB_topK_metrics.tsv"))
message("  ", paste0(OUT_PREFIX, "_real_vs_permB_topK_metrics.tsv"))
message("  ", paste0(OUT_PREFIX, "_real_vs_permB_topK_score_density.png"))
message("  ", paste0(OUT_PREFIX, "_real_vs_permB_topK_margin_density.png"))
message("  ", paste0(OUT_PREFIX, "_real_vs_permB_topK_score_vs_margin.png"))
message("  ", paste0(OUT_PREFIX, "_permB_topK_summary.tsv"))
message("  ", paste0(OUT_PREFIX, "_real_topK_summary.tsv"))
message("  ", paste0(OUT_PREFIX, "_topK_empirical_p.tsv"))
message("  ", paste0(OUT_PREFIX, "_permB_topK_mean_score_hist.png"))
message("  ", paste0(OUT_PREFIX, "_permB_topK_mean_margin_hist.png"))
message("  ", paste0(OUT_PREFIX, "_permB_topK_kth_score_hist.png"))
message("  ", paste0(OUT_PREFIX, "_overlay_real_topK/"))
message("  ", paste0(OUT_PREFIX, "_overlay_permB_topK/"))

message("Top-K empirical summary:")
print(topk_empirical_p_dt)





