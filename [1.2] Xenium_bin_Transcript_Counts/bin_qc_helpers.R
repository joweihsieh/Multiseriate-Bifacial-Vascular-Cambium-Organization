guess_xy_cols <- function(dt) {
  nms <- names(dt)
  low <- tolower(nms)

  x_idx <- which(low %in% c("x", "x_center"))
  y_idx <- which(low %in% c("y", "y_center"))

  if (length(x_idx) == 0 || length(y_idx) == 0) {
    stop("Cannot find x/y columns in polygon csv.")
  }

  list(x = nms[x_idx[1]], y = nms[y_idx[1]])
}

read_polygon_csv <- function(path) {
  dt <- data.table::fread(path)
  xy <- guess_xy_cols(dt)

  poly <- dt[, .(
    x = as.numeric(get(xy$x)),
    y = as.numeric(get(xy$y))
  )]

  if (nrow(poly) < 3) {
    stop("Polygon file has fewer than 3 points: ", path)
  }

  if (!(poly$x[1] == poly$x[nrow(poly)] && poly$y[1] == poly$y[nrow(poly)])) {
    poly <- rbind(poly, poly[1])
  }

  poly
}

point_in_polygon_vec <- function(px, py, vx, vy) {
  n <- length(vx)
  inside <- rep(FALSE, length(px))

  j <- n
  for (i in seq_len(n)) {
    xi <- vx[i]
    yi <- vy[i]
    xj <- vx[j]
    yj <- vy[j]

    intersect <- ((yi > py) != (yj > py)) &
      (px < (xj - xi) * (py - yi) / ((yj - yi) + 1e-30) + xi)

    inside <- xor(inside, intersect)
    j <- i
  }

  inside
}

split_bin_meta_by_selection <- function(
  bin_meta,
  sample_dir,
  bottom_coord_file = "Selection_Bottom_coordinates.csv",
  upper_coord_file = "Selection_Upper_coordinates.csv"
) {
  bottom_path <- file.path(sample_dir, bottom_coord_file)
  upper_path <- file.path(sample_dir, upper_coord_file)

  if (!file.exists(bottom_path) || !file.exists(upper_path)) {
    return(list(ALL = data.table::copy(bin_meta)))
  }

  poly_bottom <- read_polygon_csv(bottom_path)
  poly_upper <- read_polygon_csv(upper_path)

  in_bottom <- point_in_polygon_vec(
    bin_meta$x_center, bin_meta$y_center,
    poly_bottom$x, poly_bottom$y
  )

  in_upper <- point_in_polygon_vec(
    bin_meta$x_center, bin_meta$y_center,
    poly_upper$x, poly_upper$y
  )

  list(
    Bottom = data.table::copy(bin_meta[in_bottom]),
    Upper = data.table::copy(bin_meta[in_upper])
  )
}

calc_distribution <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(list(
      min = NA_real_,
      Q1 = NA_real_,
      median = NA_real_,
      mean = NA_real_,
      Q3 = NA_real_,
      max = NA_real_
    ))
  }

  q <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

  list(
    min = unname(q[1]),
    Q1 = unname(q[2]),
    median = unname(q[3]),
    mean = mean(x, na.rm = TRUE),
    Q3 = unname(q[4]),
    max = unname(q[5])
  )
}

plot_hist <- function(x, outfile, title, xlab_txt = "Transcripts per bin") {
  x <- x[is.finite(x)]
  x <- x[x > 0]

  if (length(x) < 10) {
    return(invisible(NULL))
  }

  xmax <- as.numeric(quantile(x, 0.99, na.rm = TRUE))

  png(outfile, width = 1200, height = 900)

  par(
    bg = "black",
    col = "white",
    col.axis = "white",
    col.lab = "white",
    col.main = "white"
  )

  hist(
    x[x <= xmax],
    breaks = 100,
    col = "white",
    border = "white",
    main = title,
    xlab = xlab_txt,
    ylab = "Bin count"
  )

  dev.off()
}

aggregate_to_superbins <- function(dt_region, target_bin_um = 50) {
  if (nrow(dt_region) == 0) {
    return(data.table::data.table(
      super_x = integer(),
      super_y = integer(),
      total_counts = numeric(),
      n_subbins = integer(),
      x_center = numeric(),
      y_center = numeric(),
      row_ids = list()
    ))
  }

  dt2 <- data.table::copy(dt_region)

  x0 <- min(dt2$x_center, na.rm = TRUE)
  y0 <- min(dt2$y_center, na.rm = TRUE)

  dt2[, super_x := floor((x_center - x0) / target_bin_um)]
  dt2[, super_y := floor((y_center - y0) / target_bin_um)]

  agg <- dt2[, .(
    total_counts = sum(total_counts, na.rm = TRUE),
    n_subbins = .N,
    x_center = mean(x_center, na.rm = TRUE),
    y_center = mean(y_center, na.rm = TRUE),
    row_ids = list(row_id)
  ), by = .(super_x, super_y)]

  agg[]
}

calc_superbin_gene_counts <- function(agg_dt, mat) {
  n <- nrow(agg_dt)

  if (n == 0) {
    return(numeric())
  }

  g_counts <- numeric(n)

  for (i in seq_len(n)) {
    rows <- agg_dt$row_ids[[i]]

    if (length(rows) == 0) {
      g_counts[i] <- 0
      next
    }

    gene_sum <- Matrix::colSums(mat[rows, , drop = FALSE])
    g_counts[i] <- sum(gene_sum > 0)
  }

  g_counts
}

get_sample_dirs <- function(xenium_base) {
  sample_dirs <- list.dirs(xenium_base, full.names = TRUE, recursive = FALSE)
  sample_dirs[grepl("^.*/output-", sample_dirs)]
}
