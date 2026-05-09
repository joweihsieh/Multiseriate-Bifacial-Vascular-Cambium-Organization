#!/usr/bin/env Rscript
# use this!!
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(42)

############################################################
# 1. basic settings
############################################################
OUTDIR <- "theoretical_patterns_final_cambium_xylem"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

NX <- 160
NY <- 60
N_GENES <- 300
K_VALUES <- 2:6

make_grid <- function() {
  dt <- CJ(ix = 1:NX, iy = 1:NY)
  dt[, x := seq(-1, 1, length.out = NX)[ix]]
  dt[, y := seq(-0.5, 0.5, length.out = NY)[iy]]
  dt[]
}

############################################################
# 2. Define unified gene-expression anchors
############################################################
base_expr <- rnorm(N_GENES, 0, 0.5)

# Ensure that K=2 cleanly separates the phloem side from the xylem side
anchor_P <- base_expr + c(rep(20, 60), rep(0, 240)) 
anchor_C <- base_expr + c(rep(0, 60), rep(8, 60), rep(0, 180))
anchor_X <- base_expr + c(rep(0, 120), rep(10, 60), rep(0, 120))

# Salt-and-pepper noise axis, used specifically for the central region of the multiseriate model
vec_SP   <- c(rep(0, 180), rep(12, 60), rep(0, 60))

############################################################
# 3. Data-generation function
############################################################
generate_data <- function(model_name) {
  dt <- make_grid()
  dt[, wavy_x := pmin(pmax(x + 0.05 * sin(2 * pi * y), -1), 1)]

  if (model_name == "Uniseriate") {
    # Ideal gradient: cambium initial cells form a single central layer
    dt[, w_P := pmax(0, -wavy_x)]
    dt[, w_X := pmax(0,  wavy_x)]
    dt[, w_C := pmax(0, 1 - abs(wavy_x))]
    
  } else if (model_name == "Multiseriate") {
    # Multiseriate model: cambium initial cells occupy a broad central region
    dt[, w_P := pmax(0, pmin(1, -wavy_x * 2.5 - 0.5))]
    dt[, w_X := pmax(0, pmin(1,  wavy_x * 2.5 - 0.5))]
    dt[, w_C := pmax(0, 1 - w_P - w_X)]
    
  } else if (model_name == "Segregate") {
    # One-directional trajectory: pure cambium to pure xylem, ignoring phloem
    progression <- (dt$wavy_x + 1) / 2
    dt[, w_P := 0]
    dt[, w_C := pmax(0, 1 - progression)]
    dt[, w_X := pmax(0, progression)]
  }

  # Ensure that the weights sum to 1
  dt[, sum_w := w_P + w_C + w_X]
  dt[, `:=`(w_P = w_P/sum_w, w_C = w_C/sum_w, w_X = w_X/sum_w)]

  # Set plotting colors: blue for phloem, orange for cambium, and yellow for xylem
  cP <- c(46, 120, 199)
  cC <- c(244, 177, 131)  # Orange, used for the uniseriate and multiseriate models
  cX <- c(242, 183, 0)    # Yellow
  
  # Use light yellow as the left-side starting color specifically for the segregated model
  if (model_name == "Segregate") {
    c_Start <- c(255, 242, 178) # Light yellow
  } else {
    c_Start <- cC               # Keep orange for the other models
  }

  dt[, plot_col := grDevices::rgb(
    (w_P*cP[1] + w_C*c_Start[1] + w_X*cX[1])/255,
    (w_P*cP[2] + w_C*c_Start[2] + w_X*cX[2])/255,
    (w_P*cP[3] + w_C*c_Start[3] + w_X*cX[3])/255
  )]

  # Build the expression matrix
  expr <- matrix(0, nrow = nrow(dt), ncol = N_GENES)
  dt[, sp_state := sample(c(-1, 1), .N, replace = TRUE)]

  for(i in 1:nrow(dt)) {
    base_sig <- dt$w_P[i] * anchor_P + dt$w_C[i] * anchor_C + dt$w_X[i] * anchor_X
    
    # === Multiseriate-specific salt-and-pepper noise ===
    if (model_name == "Multiseriate") {
      # Scale the strong binary noise by w_C so that it occurs only in the central region
      sp_strength <- dt$w_C[i] * 1.5 
      base_sig <- base_sig + sp_strength * dt$sp_state[i] * vec_SP
    }
    
    # Add a small global cell-level baseline noise term
    expr[i, ] <- base_sig + rnorm(N_GENES, 0, 1.0)
  }

  dt[, model := model_name]
  list(meta = dt, expr_scaled = scale(expr))
}

############################################################
# 4. Run simulations and clustering
############################################################
cat("Generating models...\n")
models <- c("Uniseriate", "Multiseriate", "Segregate")
all_sims <- lapply(models, generate_data)
names(all_sims) <- models

# 1. Spatial ground-truth gradient PNG
all_meta <- rbindlist(lapply(all_sims, function(x) x$meta))
p_base <- ggplot(all_meta, aes(x, y, fill = plot_col)) +
  geom_tile() + scale_fill_identity() + coord_fixed() +
  facet_wrap(~ model, ncol = 1) + theme_bw(base_size = 12) +
  labs(title = "Ground Truth Spatial Identity") +
  theme(panel.grid = element_blank())
ggsave(file.path(OUTDIR, "00_base_patterns.png"), p_base, width = 8, height = 12, dpi = 300)

# 2. Run K-means from K=2 to K=6
cat("Running K-means for K=2 to 6...\n")
all_cluster_dt <- rbindlist(lapply(names(all_sims), function(nm) {
  rbindlist(lapply(K_VALUES, function(k) {
    km <- kmeans(all_sims[[nm]]$expr_scaled, centers = k, nstart = 50)
    out <- copy(all_sims[[nm]]$meta)
    out[, cluster := factor(km$cluster)]
    out[, k := k]
    out[, model := nm]
    out
  }))
}), fill = TRUE)

# 3. Plot comparisons across all K values
cat("Plotting...\n")
p_all <- ggplot(all_cluster_dt, aes(x, y, fill = cluster)) +
  geom_tile() + coord_fixed() +
  facet_grid(model ~ k) + theme_bw(base_size = 12) +
  labs(title = "K-means Clustering Patterns across Models (K = 2 to 6)") +
  theme(
    panel.grid = element_blank(), 
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )
ggsave(file.path(OUTDIR, "01_k2_to_k6_compare.png"), p_all, width = 18, height = 10, dpi = 300)

cat("Done! Check PNG outputs in:", OUTDIR, "\n")