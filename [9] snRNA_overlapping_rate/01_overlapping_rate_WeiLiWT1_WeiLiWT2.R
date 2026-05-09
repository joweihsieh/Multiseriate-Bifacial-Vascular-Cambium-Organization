###### change filenames for each pair of your interest

#install.packages("Seurat")
library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)

library(igraph)
library(magrittr)
library(Matrix)
library(MASS)
library(Seurat)
library(ggplot2)


oriPar = par(no.readonly=T)
setwd("/home/woodydrylab/FileShare/20260117_MultiLayer_Unifacial/overlapping")


#Setup the Seurat Object
WeiLiWT1.data = Read10X(data.dir = '/home/woodydrylab/DiskArray/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/cellranger_reanalysis_TenX_PtrWT1forWOX1/outs/filtered_feature_bc_matrix')
#/home/woodydrylab/DiskArray/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/cellranger_reanalysis_TenX_PtrWT2forWOX2/outs/filtered_feature_bc_matrix
WeiLiWT2.data = Read10X(data.dir = '/home/woodydrylab/DiskArray/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/cellranger_reanalysis_TenX_PtrWT2forWOX2/outs/filtered_feature_bc_matrix')
#/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/celescope_out/NuWeiLiWT2WTBio2/outs/filtered

WeiLiWT1 = CreateSeuratObject(counts = WeiLiWT1.data, project = 'WeiLiWT1', min.cells = 3, min.features = 200)
WeiLiWT2 = CreateSeuratObject(counts = WeiLiWT2.data, project = 'WeiLiWT2', min.cells = 3, min.features = 200)


WeiLiWT2$orig.ident <- "WeiLiWT2"
WeiLiWT1$orig.ident <- "WeiLiWT1"

WeiLiWT1 <- RenameCells(WeiLiWT1, add.cell.id = "WeiLiWT1")
WeiLiWT2 <- RenameCells(WeiLiWT2, add.cell.id = "WeiLiWT2")

WeiLiWT1$orig.ident <- "WeiLiWT1"
WeiLiWT2$orig.ident <- "WeiLiWT2"

#Normalize the data
WeiLiWT1 = NormalizeData(WeiLiWT1, normalization.method = 'LogNormalize', scale.factor = 10000)
WeiLiWT2 = NormalizeData(WeiLiWT2, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
WeiLiWT1 = FindVariableFeatures(WeiLiWT1, selection.method = 'vst', nfeatures = 2000)
WeiLiWT2 = FindVariableFeatures(WeiLiWT2, selection.method = 'vst', nfeatures = 2000)

#Find integration anchors and integrate data
integration_anchors = FindIntegrationAnchors(object.list = list(WeiLiWT1,WeiLiWT2),
                                             anchor.features = 2000,
                                             scale = TRUE,
                                             reduction = 'cca',
                                             l2.norm = TRUE,
                                             k.anchor = 5)

saveRDS(integration_anchors,'RDS_integration_anchors_WeiLiWT1_WeiLiWT2.rds')

integration_anchors = readRDS('RDS_integration_anchors_WeiLiWT1_WeiLiWT2.rds')
Combined_object = IntegrateData(anchorset = integration_anchors)

#Run the standard workflow for visualization and clustering
# stopifnot(DefaultAssay(Combined_object) == 'integrated')
Combined_object = ScaleData(Combined_object)
Combined_object = RunPCA(Combined_object, npcs = 30)
Combined_object = RunUMAP(Combined_object, reduction = 'pca', dims = 1:30)
Combined_object = FindNeighbors(Combined_object, reduction = 'pca', dims = 1:30, k.param = 3)
Combined_object = FindClusters(Combined_object, resolution = 0.5)


saveRDS(Combined_object,'RDS_Combined_object_WeiLiWT1_WeiLiWT2.rds')



#Combined_object <- readRDS('RDS_Combined_object_WeiLiWT1_WeiLiWT2.rds')
png("UMAP_WeiLiWT1_WeiLiWT2.png", width = 2400, height = 1800, res = 300)
print(
  DimPlot(
    Combined_object,
    reduction = "umap",
    group.by = "orig.ident",
    pt.size = 0.3
  ) +
    ggtitle("Integrated UMAP by dataset") +
    theme_classic()
)
dev.off()



Combined_object <- readRDS("RDS_Combined_object_WeiLiWT1_WeiLiWT2.rds")

umap_coords <- Embeddings(Combined_object, "umap")
idents <- Combined_object$orig.ident

group_colors <- c("WeiLiWT1" = "#000000", "WeiLiWT2" = "#C59738")
cell_colors <- group_colors[as.character(idents)]

png("UMAP_WeiLiWT1_WeiLiWT2_20260411.png", pointsize = 10, width = 20, height = 15, units = "cm", res = 300)

plot(NA,
     xlim = range(umap_coords[, 1]),
     ylim = range(umap_coords[, 2]),
     xlab = '', ylab = '', axes = FALSE, main = "WeiLiWT1_WeiLiWT2")

points(umap_coords[, 1], umap_coords[, 2], col = cell_colors, pch = 16, cex = 0.6)

dev.off()


##########################


library(igraph)
library(magrittr)
library(Matrix)
library(MASS)
library(Seurat)
library(ggplot2)

setwd("/home/woodydrylab/FileShare/20260117_MultiLayer_Unifacial/overlapping")


oriPar = par(no.readonly=T)

############################################################################## Step 1. Generating RDS 


getMSTsubtreeCenter = function(projection, sp1, sp2){
    print('Calculate the distance between each pair of cells')
    distMatrix = as.matrix(dist(projection)) %>% Matrix(sparse=T)
    stopifnot(sum(distMatrix==0) == nrow(projection))

    print('Create the graph from adjacent matrix')
    graphFull = graph_from_adjacency_matrix(distMatrix,
                                            mode='undirected',weighted=T)

    print('Construct the MST')
    graphMST = mst(graphFull)

    print('Remove inter-species edges')
    edgeVname = attr(E(graphMST),'vnames')
    delEdge = edgeVname %>% grep(sp1,.,value=T) %>% grep(sp2,.,value=T)
    graphCutMST = delete_edges(graphMST,delEdge)

    print('Extract the subgraph centers')
    subgraphCenter = c()
    candidateVertices = attr(V(graphCutMST),'name')
    while(length(candidateVertices)>0){
        pickedVertex = candidateVertices[1]
        pickedVertices = attr(subcomponent(graphCutMST,pickedVertex),'name')
        pickedGraph = induced_subgraph(graphCutMST,pickedVertices)
        pickedCloseness = closeness(pickedGraph)
        if(length(pickedCloseness)==1){
            subgraphCenter %<>% c(names(pickedCloseness))
        }else{
            subgraphCenter %<>% c(names(which.max(pickedCloseness)))
        }
        candidateVertices %<>% setdiff(pickedVertices)
    }
    return(subgraphCenter)
}



runUMAPandSaveSubtreeCenter <- function(rdsFilePath, prefix, sp1, sp2) {
    # Read the RDS file
    ClaPAIR_combinedObject = readRDS(rdsFilePath)

    # Run UMAP
    ClaPAIR_combinedObject <- RunUMAP(
        object = ClaPAIR_combinedObject,
        reduction = "pca",
        dims = 1:30,
        seed.use = 42,
        min.dist = 0.3, # 0.3 # c(0.001, 0.5)
        n.neighbors = 30, # 30L # c(5, 50)
        umap.method = "uwot", metric = "cosine"
    )

    # Get UMAP projections
    ClaPAIR_projectionUMAP = ClaPAIR_combinedObject@reductions$umap@cell.embeddings

    print('Calculating MST and saving subtree center')
    ClaPAIR_subtreeCenter = getMSTsubtreeCenter(ClaPAIR_projectionUMAP, sp1, sp2)
    saveRDS(ClaPAIR_subtreeCenter, paste0('RDS_', prefix, '_subtreeCenter.rds'))
}

runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_WeiLiWT1_WeiLiWT2.rds',
    prefix = 'WeiLiWT1_WeiLiWT2',
    sp1 = "WeiLiWT1", 
    sp2 = "WeiLiWT2" 
)

############################################################################## Step 2. maps


runUMAPandSaveSubtreeCenter <- function(rdsFilePath) {
    # Read the RDS file
    ClaPAIR_combinedObject = readRDS(rdsFilePath)

    # Run UMAP
    ClaPAIR_combinedObject <- RunUMAP(
        object = ClaPAIR_combinedObject,
        reduction = "pca",
        dims = 1:30,
        seed.use = 42,
        min.dist = 0.3, # 0.3 # c(0.001, 0.5)
        n.neighbors = 30, # 30L # c(5, 50)
        umap.method = "uwot", metric = "cosine"
    )

    # Get UMAP projections
    ClaPAIR_projectionUMAP = ClaPAIR_combinedObject@reductions$umap@cell.embeddings

    print('Calculating MST and saving subtree center')
    return(ClaPAIR_projectionUMAP)
}



centerContourPlot = function(subgraphCenter = ClaPAIR_subtreeCenter,
                             projectionUMAP = ClaPAIR_projectionUMAP,
                             lowerPercentile = 5,
                             higherPercentile = 40,
                             main = 'PtrTens12_rep',
                             plot_name = 'PtrTens12_rep',
                             sp1 = "TW_",
                             df_path = NULL){
    
    Ptr_projectionUMAP = projectionUMAP %>% extract(grepl(sp1,rownames(.)),)
    Ptr_density_map = kde2d(Ptr_projectionUMAP[,1],
                            Ptr_projectionUMAP[,2],n=500,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
    get_territory_density = function(Coor){
        out = Ptr_density_map$z[max(which(Ptr_density_map$x < Coor[1])),
                                max(which(Ptr_density_map$y < Coor[2]))]
        return(out)
    }
    Ptr_territory_density = Ptr_projectionUMAP %>% apply(1,get_territory_density)
    norm_factor = 1/sum(Ptr_territory_density)
    # print(norm_factor)
    
    Ptr_SubgraphCenter = subgraphCenter %>% extract(grepl(sp1,.))
    Sp2_SubgraphCenter = subgraphCenter %>% extract(!grepl(sp1,.))
    
    num_Ptr =  sum(grepl(sp1,rownames(projectionUMAP)))
    num_Sp2 =  sum(!grepl(sp1,rownames(projectionUMAP)))
    num_total = nrow(projectionUMAP)
    num_center_Ptr = length(Ptr_SubgraphCenter)
    num_center_Sp2 = length(Sp2_SubgraphCenter)

    get_subgraphcenter_density = function(partSubgraphCenter){

        allCenterCoordinate = projectionUMAP[subgraphCenter,]
        partCenterCoordinate = projectionUMAP[partSubgraphCenter,]
        density_map = kde2d(partCenterCoordinate[,1],
                            partCenterCoordinate[,2],n=500,
                            h=apply(allCenterCoordinate,2,bandwidth.nrd)/2,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
        return(density_map)
    }
    Ptr_subgraphcenter_density = get_subgraphcenter_density(Ptr_SubgraphCenter)
    Sp2_subgraphcenter_density = get_subgraphcenter_density(Sp2_SubgraphCenter)
    stopifnot(Ptr_subgraphcenter_density$x==Sp2_subgraphcenter_density$x)
    stopifnot(Ptr_subgraphcenter_density$y==Sp2_subgraphcenter_density$y)
    
    merge_subgraphcenter_density = Ptr_subgraphcenter_density
    merge_subgraphcenter_density$z =
        (num_Ptr/num_total)*num_center_Ptr*Ptr_subgraphcenter_density$z +
        (num_Sp2/num_total)*num_center_Sp2*Sp2_subgraphcenter_density$z
    
    merge_subgraphcenter_density$z %<>% multiply_by(norm_factor)
    
    message(
      "Plot density max:",
      round(max(merge_subgraphcenter_density$z), 5)
    )
    message(
      "Plot density Q75:",
      round(quantile(merge_subgraphcenter_density$z, 0.75), 5)
    )
    message(
      "Plot density Q50:",
      round(quantile(merge_subgraphcenter_density$z, 0.50), 5)
    )
    message(
      "Plot density Q25:",
      round(quantile(merge_subgraphcenter_density$z, 0.25), 5)
    )
    message(
      "Plot density min:",
      round(min(merge_subgraphcenter_density$z), 5)
    )
    
     territory_map = kde2d(projectionUMAP[,1],
                          projectionUMAP[,2],n=500,h=0.02,
                          lims=c(min(projectionUMAP[,1])-0.5,
                                 max(projectionUMAP[,1])+0.5,
                                 min(projectionUMAP[,2])-0.5,
                                 max(projectionUMAP[,2])+0.5))
    
    png(paste0(plot_name,'.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    {
        plot(NA,
             xlim=range(projectionUMAP[,'umap_1']),
             ylim=range(projectionUMAP[,'umap_2']),
             xlab='',ylab='',axes=F,main=main)
        Levels = seq(0,1,length.out=500)
        lowerRank = 500 * lowerPercentile/100
        higherRank = 500 * higherPercentile/100
        .filled.contour(merge_subgraphcenter_density$x,
                        merge_subgraphcenter_density$y,
                        merge_subgraphcenter_density$z,
                        levels=Levels,
                        col=c(colorRampPalette(c('#AED2F5','#FECB71'))(lowerRank),
                              colorRampPalette(c('#FECB71','#F8696B'))(higherRank-lowerRank),
                              colorRampPalette(c('#F8696B','#713020'))(500-higherRank)))
         .filled.contour(territory_map$x,
                        territory_map$y,
                        ifelse(territory_map$z>0,1,0),
                        levels=c(0,0.5),col=c('white',NA))
        contour(territory_map$x,territory_map$y,ifelse(territory_map$z>0,1,0),
                levels=0.5,lwd=0.8,drawlabels=F,add=T)
    }
    dev.off()
}



##################### 

ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_WeiLiWT1_WeiLiWT2.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_WeiLiWT1_WeiLiWT2_subtreeCenter.rds')


centerContourPlot(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'WeiLiWT1_WeiLiWT2',
    plot_name = 'WeiLiWT1_WeiLiWT2',
    sp1 = "WeiLiWT1" # WeiLiWT2
    )



############################################################################## Step 3. pie



centerContourPlot_pie = function(subgraphCenter = ClaPAIR_subtreeCenter,
                             projectionUMAP = ClaPAIR_projectionUMAP,
                             lowerPercentile = 5,
                             higherPercentile = 40,
                             main = 'PtrVert2',
                             plot_name = 'PtrVert2',
                             sp1 = "Ptr_",
                             df_path = paste0('Overlap_heatpie_', "PtrVert2", '_on_Ptr.csv')){
    
    Ptr_projectionUMAP = projectionUMAP %>% extract(grepl(sp1,rownames(.)),)
    Ptr_density_map = kde2d(Ptr_projectionUMAP[,1],
                            Ptr_projectionUMAP[,2],n=500,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
    get_territory_density = function(Coor){
        out = Ptr_density_map$z[max(which(Ptr_density_map$x < Coor[1])),
                                max(which(Ptr_density_map$y < Coor[2]))]
        return(out)
    }
    Ptr_territory_density = Ptr_projectionUMAP %>% apply(1,get_territory_density)
    norm_factor = 1/sum(Ptr_territory_density)
    # print(norm_factor)
    
    Ptr_SubgraphCenter = subgraphCenter %>% extract(grepl(sp1,.))
    Sp2_SubgraphCenter = subgraphCenter %>% extract(!grepl(sp1,.))
    
    num_Ptr =  sum(grepl(sp1,rownames(projectionUMAP)))
    num_Sp2 =  sum(!grepl(sp1,rownames(projectionUMAP)))
    num_total = nrow(projectionUMAP)
    num_center_Ptr = length(Ptr_SubgraphCenter)
    num_center_Sp2 = length(Sp2_SubgraphCenter)

    get_subgraphcenter_density = function(partSubgraphCenter){
        allCenterCoordinate = projectionUMAP[subgraphCenter,]
        partCenterCoordinate = projectionUMAP[partSubgraphCenter,]
        density_map = kde2d(partCenterCoordinate[,1],
                            partCenterCoordinate[,2],n=500,
                            h=apply(allCenterCoordinate,2,bandwidth.nrd)/2,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
        return(density_map)
    }
    Ptr_subgraphcenter_density = get_subgraphcenter_density(Ptr_SubgraphCenter)
    Sp2_subgraphcenter_density = get_subgraphcenter_density(Sp2_SubgraphCenter)
    stopifnot(Ptr_subgraphcenter_density$x==Sp2_subgraphcenter_density$x)
    stopifnot(Ptr_subgraphcenter_density$y==Sp2_subgraphcenter_density$y)
    
    merge_subgraphcenter_density = Ptr_subgraphcenter_density
    merge_subgraphcenter_density$z =
        (num_Ptr/num_total)*num_center_Ptr*Ptr_subgraphcenter_density$z +
        (num_Sp2/num_total)*num_center_Sp2*Sp2_subgraphcenter_density$z
    
    merge_subgraphcenter_density$z %<>% multiply_by(norm_factor)
    
    
    territory_map = kde2d(projectionUMAP[,1],
                          projectionUMAP[,2],n=500,h=0.02,
                          lims=c(min(projectionUMAP[,1])-0.5,
                                 max(projectionUMAP[,1])+0.5,
                                 min(projectionUMAP[,2])-0.5,
                                 max(projectionUMAP[,2])+0.5))
    
    message(
      "Cell region density max:",
      round(max(merge_subgraphcenter_density$z[territory_map$z > 0]), 5)
    )
    message(
      "Cell region density Q75:",
      round(quantile(
        merge_subgraphcenter_density$z[territory_map$z > 0], 0.75), 5)
    )
    message(
      "Cell region density Q50:",
      round(quantile(
        merge_subgraphcenter_density$z[territory_map$z > 0], 0.50), 5)
    )
    message(
      "Cell region density Q25:",
      round(quantile(
        merge_subgraphcenter_density$z[territory_map$z > 0], 0.25), 5)
    )
    message(
      "Cell region density min:",
      round(min(merge_subgraphcenter_density$z[territory_map$z > 0]), 5)
    )

    {
        n_col_levels <- 500
        lower_rank <- n_col_levels * lowerPercentile / 100
        higher_rank <- n_col_levels * higherPercentile / 100
        filled_col <-
            c(colorRampPalette(c('#AED2F5','#FECB71'))(lower_rank),
              colorRampPalette(c('#FECB71','#F8696B'))(higher_rank-lower_rank),
              colorRampPalette(c('#F8696B','#713020'))(500-higher_rank))
        col_counts <-
            cut(
                merge_subgraphcenter_density$z[territory_map$z > 0],
                breaks = seq(0, 1, length.out = n_col_levels + 1)
            ) %>%
            table() %>%
            as.vector()
        col_counts_df <- data.frame(filled_col, col_counts)
#Densities from 0 to 1 are divided into 500 bins with different color shading, with proportions of different densities shown in a pie chart in each panel.  
#remove #B4D1EA  
#total is 60   
#sum(col_counts_df$col_counts)
#75898
#(75898-60)/75898

        agg_filled_col <- sapply(
            seq(n_col_levels / 5),
            function(i) {
                out <- col_counts_df$filled_col[(i - 1) * 5 + 3]
                return(out)
            }
        )
        agg_col_counts <- sapply(
            seq(n_col_levels / 5),
            function(i) {
                out <- sum(col_counts_df$col_counts[1:5 + (i - 1) * 5])
                return(out)
            }
        )
        col_counts_df <- data.frame(
            filled_col = agg_filled_col,
            col_counts = agg_col_counts
        )
        
        col_counts_df$col_props <-
            col_counts_df$col_counts / sum(col_counts_df$col_counts)
        
        if (!is.null(df_path)) write.csv(col_counts_df, df_path, row.names = FALSE)
        
        col_counts_df$factor_filled_col <-
            factor(col_counts_df$filled_col, levels = col_counts_df$filled_col)
        ggplot(
            col_counts_df,
            aes(x = "", y = col_props, fill = factor_filled_col)
        ) +
            geom_bar(stat = "identity", colour = "white", size = 0.05) +
            scale_fill_manual("legend", values = setNames(filled_col, filled_col)) +
            coord_polar("y", start = 0) +
            theme_void() +
            theme(legend.position="none")
        ggsave(
            paste0("Heathist_", plot_name, ".png"),
            width = 20, height = 15, units = "cm", dpi = 300
        )

    }
}



##################### 


ClaPAIR_projectionUMAP = runUMAPandSaveSubtreeCenter(
    rdsFilePath = 'RDS_Combined_object_WeiLiWT1_WeiLiWT2.rds')
ClaPAIR_subtreeCenter = readRDS('RDS_WeiLiWT1_WeiLiWT2_subtreeCenter.rds')


centerContourPlot_pie(
    subgraphCenter = ClaPAIR_subtreeCenter,
    projectionUMAP = ClaPAIR_projectionUMAP,
    lowerPercentile = 5,
    higherPercentile = 40,
    main = 'WeiLiWT1_WeiLiWT2',
    plot_name = 'WeiLiWT1_WeiLiWT2',
    sp1 = "WeiLiWT1",
    df_path = paste0('Overlap_heatpie_', "WeiLiWT1", '_on_WeiLiWT2.csv')
    )


###### color bar

plot_standalone_colorbar <- function(lowerPercentile = 5,
                                     higherPercentile = 40,
                                     plot_name = "Standalone_ColorBar") {
  
  Levels <- seq(0, 1, length.out = 500)
  lowerRank <- 500 * lowerPercentile / 100
  higherRank <- 500 * higherPercentile / 100
  
  filled_col <- c(colorRampPalette(c('#AED2F5','#FECB71'))(lowerRank),
                  colorRampPalette(c('#FECB71','#F8696B'))(higherRank-lowerRank),
                  colorRampPalette(c('#F8696B','#713020'))(500-higherRank))
  
  png(paste0(plot_name, ".png"), width = 4, height = 12, units = "cm", res = 300, bg = "transparent")
  
  par(mar = c(2, 0.5, 2, 4))
  
  image(x = 1, y = Levels, z = t(as.matrix(Levels)), 
        col = filled_col, axes = FALSE, xlab = "", ylab = "")
  
  box(lwd = 1.5)             
  axis(4, las = 1, cex.axis = 1.2) 
  
  dev.off()
  
  message("Standalone color bar saved as: ", paste0(plot_name, ".png"))
}

plot_standalone_colorbar()

######

ovelaps <- read.csv("Overlap_heatpie_WeiLiWT1_on_WeiLiWT2.csv")
percent <- (sum(ovelaps$col_counts)-(ovelaps[ovelaps$filled_col=="#B4D1EA","col_counts"]))/sum(ovelaps$col_counts)*100
print(percent)
#99.66033




