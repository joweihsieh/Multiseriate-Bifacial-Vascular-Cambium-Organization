
setwd("/home/woodydrylab/FileShare/20260117_MultiLayer_Unifacial")
options(jupyter.rich_display = FALSE)

library(dplyr)

conde_color <- read.csv("/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/Conde_annotation.csv")
conde <- read.csv("/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/UMAP_csv/Conde_projection_UMAP.csv")

conde_color$cluster <- as.integer(conde_color$cluster)

# Colors
conde_color <- conde_color %>%
  mutate(
    Color = case_when(
      cluster == 12 ~ "#26BDF8",      
      cluster == 14 ~ "#079EDB",      
      cluster == 18 ~ "#235889",      
      cluster %in% c(4, 16) ~ "#9DC3E6",
      cluster %in% c(1, 2, 8, 20) ~ "#0070C0",      
      cluster == 5 ~ "#DEED17",
      cluster %in% c(6, 19, 22) ~ "#FAB40A",
      cluster %in% c(0, 3, 7, 13) ~ "#FFE699",
      cluster %in% c(9, 10, 11, 15, 17, 21, 23) ~ "#A6A6A6",
      TRUE ~ Color
    )
  )

# merge
conde_cluster <- conde %>%
  inner_join(
    conde_color %>% select(Barcode, cluster, Color),
    by = "Barcode"
  )


plot_order <- seq_len(nrow(conde_cluster))

png(
  "integrated_Conde_WeiLi_annotation_Conde_only_recolor_20260411.png",
  pointsize = 10,
  res = 300,
  width = 15,
  height = 15,
  units = "cm"
)

par(mai = c(0.7, 0.7, 0.9, 0.5))

plot(
  conde_cluster$umap_1[plot_order],
  conde_cluster$umap_2[plot_order],
  col = conde_cluster$Color[plot_order],
  pch = 20,
  cex = 0.75,
  main = "Conde",
  xlab = "UMAP.1",
  ylab = "UMAP.2",
  las = 1,
  asp = 1
)

dev.off()


###########

setwd("/home/woodydrylab/FileShare/20260117_MultiLayer_Unifacial")
options(jupyter.rich_display = FALSE)

library(dplyr)

WeiLi <- read.csv("/home/chingweilu/Multi_Bilayer_Unifacial/Adjusting_UMAP_neighbor/UMAP_csv/WeiLi_projection_UMAP.csv")




png(
  "integrated_WeiLi_WeiLi_annotation_WeiLi_only_recolor_20260411.png",
  pointsize = 10,
  res = 300,
  width = 15,
  height = 15,
  units = "cm"
)

par(mai = c(0.7, 0.7, 0.9, 0.5))

plot(
  WeiLi$umap_1,
  WeiLi$umap_2,
  col = "black",   
  pch = 20,
  cex = 0.75,
  main = "WeiLi",
  xlab = "UMAP.1",
  ylab = "UMAP.2",
  las = 1,
  asp = 1
)

dev.off()