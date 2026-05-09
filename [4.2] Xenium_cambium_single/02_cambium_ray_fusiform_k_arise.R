library(data.table)
library(ggplot2)
library(readxl)

setwd("/Users/jo-wei/Library/CloudStorage/Dropbox/YCL/Spatial_Transcriptomes/Xenium/20260121_Data/20260208_grid_based/5um_k10_cambium_one_two_arise")

dt <- as.data.table(
  read_excel("/Users/jo-wei/Library/CloudStorage/Dropbox/YCL/Spatial_Transcriptomes/Xenium/20260121_Data/code/Cambium_Ray_pairs_start.xlsx")
)

dt[, KStart_cambium_one := as.numeric(KStart_cambium_one)]
dt[, Kstart_cambium_two := as.numeric(Kstart_cambium_two)]
dt[, k := as.numeric(k)]

dt0 <- copy(dt)

dt[, row_id := .I]

dt_long <- dt[, {
  sample_vec <- trimws(unlist(strsplit(sample_id, ",")))
  label_vec  <- trimws(unlist(strsplit(label, ",")))
  
  n_max <- max(length(sample_vec), length(label_vec))
  if (length(sample_vec) == 1 && n_max > 1) sample_vec <- rep(sample_vec, n_max)
  if (length(label_vec)  == 1 && n_max > 1) label_vec  <- rep(label_vec, n_max)
  
  data.table(
    sample_id_split = sample_vec,
    label_split = label_vec
  )
}, by = .(
  row_id, sample_id, label, um, k, pair_cambium, pair_ray,
  KStart_cambium_one, cambium, Kstart_cambium_two,
  cambium_ray, cambium_fusiform
)]

dt_long[, sample_show := label_split]

plot_dt <- dt_long[, {
  ks <- 2:k[1]
  
  state <- sapply(ks, function(one_k) {
    if (one_k < KStart_cambium_one[1]) {
      "No"
    } else if (one_k < Kstart_cambium_two[1]) {
      "Cambium"
    } else {
      "Cambium: ray and fusiform"
    }
  })
  
  data.table(k_plot = ks, state = state)
}, by = .(
  sample_id_split, sample_show,
  KStart_cambium_one, Kstart_cambium_two
)]

#sample_order <- unique(
#  dt_long[order(KStart_cambium_one, Kstart_cambium_two), sample_show]
#)

sample_order <- sort(unique(plot_dt$sample_show))

plot_dt[, sample_show := factor(sample_show, levels = rev(sample_order))]

p <- ggplot(plot_dt, aes(x = k_plot, y = sample_show, fill = state)) +
  geom_tile(color = "white", linewidth = 0.8, height = 0.85) +
  scale_fill_manual(values = c(
    "No" = "grey90",
    "Cambium" = "#4DBBD5",
    "Cambium: ray and fusiform" = "#E64B35"
  )) +
  scale_x_continuous(
    breaks = 2:max(dt_long$k),
    expand = c(0, 0)
  ) +
  labs(
    x = "K",
    y = NULL,
    fill = NULL,
    title = "Emergence and splitting of cambium across K"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("cambium_emergence_split_timeline.pdf", p, width = 8, height = 5)
ggsave("cambium_emergence_split_timeline.png", p, width = 8, height = 5, dpi = 300)
