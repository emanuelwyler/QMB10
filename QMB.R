
library(VoltRon)


library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(igraph)
library(reshape2)
library(DESeq2)
library(apeglm)
library(dendextend)
library(stringr)
library(gridExtra)
library(cowplot)
#library(sf)
#library(concaveman)
library(data.table)
library(pals)
library(presto)
library(purrr)
library(ggpubr)
library(RColorBrewer)
library(lme4)

#Colors for cell types
celltype_cols <- c(`Adipocytes/Sebaceous Glands` = "#F8766D", `ciliated epithelial cells` = "#F27C56", `Arrector pili muscle` = "#EC823A",
                   `Endothelial cells` = "#DC8D00", `Melanocytes` = "#D89000", `MerkelCells` = "#CE9500", `Neuronal cells_1` = "#C39A00", `Neuronal cells_2` = "#B79F00",
                   `SchwannCells` = "#AAA300", `Chondrocytes` = "#93AA00",
                   `Fibroblasts: Dermal Papilla` = "#007700", `Fibroblasts: Myofibroblasts_1` = "#00BA38", `Fibroblasts: Mesenchymal_1` = "#55FF8A",
                   `Fibroblasts: Mesenchymal_2` = "#00807E", `Fibroblasts: Pro-inflammatory` = "#00C0BD", `Fibroblasts: Dermal Sheath (cyc)` = "#00FFE4", 
                   `Fibroblasts: Papillary_1` = "#0078AF", `Fibroblasts: Papillary_2` = "#00B5EE", `Fibroblasts: Papillary_3` = "#00F4FF",
                   `Fibroblasts: Dermal Sheath` = "#0957CE", `Fibroblasts: undefined` = "#BFDDFF",
                   `Keratinocytes: Stratum Basale` = "#A600BB", `Keratinocytes: Stratum corneum` = "#E16FF7", `Keratinocytes: IRS` = "#FFC7FF",
                   `Keratinocytes: Matrix` = "#BB0092", `Keratinocytes: Bulge region` = "#FD61D1", `Keratinocytes: ORS/companion layer` = "#FFB5CF",
                   `Keratinocytes: Stratum spinosum/Granulosum` = "#B42049", 
                   Undefined1 ="gray90", 
                   Undefined2 ="gray88", Undefined3 ="gray86", Undefined4 ="gray84", Undefined5 ="gray82", Undefined6 ="gray80", 
                   Undefined7 ="gray78", Undefined8 ="gray76", Undefined9 ="gray74", Undefined10 ="gray72", Undefined11 ="gray70", UndefinedCorePap="gray30")

time_levels <- c("0dpi", "2dpi", "3dpi", "4dpi")

setwd("~/temp/QMB10/")

#download the fully annotated object from GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9378020
#HSV1all_ann <- readRDS("HSV1all_anncore_annot_3.rds")
#Load pseudobulk data
pseudobulk_all <- readRDS("pseudobulk_all.rds")

#######################
#Clustering infected papillary fibroblasts
#This was used in figure 6 as well as in supplementary figures S20 and S21

#Variable genes: These are genes that in papillary fibroblasts at some points vary during infection, either up or down, and includes viral genes
#we then cluster using these genes only 
variable_genes <- pseudobulk_all %>% filter(grepl("Papillary", celltype)) %>% filter(padj < 0.01) %>% pull (var="gene_name") %>% unique()
variable_genes_seu <- gsub("_", "-", variable_genes)

#seu <- VoltRon::as.Seurat(HSV1all_ann, cell.assay = "Xenium", type = "image")
#seu$HSV1load = HSV1all_ann@metadata@cell$HSV1Count/HSV1all_ann@metadata@cell$Count*100
#t <- 2
#select_cells <- seu@meta.data %>% filter(grepl("Papillary", celltype)) %>% filter(HSV1Count >= t) %>% rownames()
#seu <- subset(seu, cells=select_cells)
#saveRDS(seu, "seu.rds")

seu <- readRDS("seu.rds")

all.genes <- rownames(seu)
seu <- NormalizeData(seu)
seu <- ScaleData(seu, features = all.genes)
seu <- RunPCA(seu, features = variable_genes_seu)
ElbowPlot(seu)
seu <- RunUMAP(seu, dims = 1:10)
seu <- FindNeighbors(seu, dims = 1:10)
#Select a resolution that somehow does not make not too many and not too few clusters
seu <- FindClusters(seu, resolution = 0.24)
#Rename clusters by increasing viral load (new column infcluster)
loadsortedvector <- seu@meta.data %>% group_by(seurat_clusters) %>% dplyr::summarize(loadmean=mean(HSV1load)) %>% arrange(loadmean) %>% pull(var="seurat_clusters")
mutateframe <- data.frame(loadsortedvector, paste0("infcluster", seq(1,length(loadsortedvector),1)))
colnames(mutateframe) <- c("seurat_clusters", "infcluster")
df <- left_join(seu@meta.data, mutateframe, by="seurat_clusters")
seu@meta.data$infcluster <- df$infcluster

#Make parula colors for color infection clusters
cluster_colors=parula(length(loadsortedvector))
names(cluster_colors)=paste0("infcluster", seq(1,length(loadsortedvector),1))

#Plot a UMAP colored by cluster, a boxplot with viral load, and (not used) and expression dot plot  
p1 <- DimPlot(seu, group.by="infcluster", reduction = "umap", cols=cluster_colors)
p2 <- ggplot()+geom_boxplot(data=seu@meta.data, aes(x=factor(infcluster), y=HSV1load, fill=infcluster))+
  theme_bw()+
  scale_fill_manual(values = cluster_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.ticks=element_blank(), panel.grid.minor=element_blank())
#PDGFRA shows host cell shut off (downregulated in later clusters), BBC3as shows antisense transcription, CD79A is a proxy for read-through (read-in from the RPS19 gene)
p3 <- DotPlot(seu, group.by="infcluster", features = c("HSV1-UL54", "HSV1-LAT", "BBC3as", "PDGFRA", "CD79A")) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
p1+p2+p3
ggsave("inf_clusters_pap.pdf", useDingbats=FALSE, width=18, height = 6)

#Feature plots with some genes for supplementary figure S20
FeaturePlot(seu, c("HSV1-US1", "IER3", "HSV1-LAT", "DUX4", "BBC3as", "CD79A"), pt.size = 0.05, ncol = 3)
ggsave("inf_clusters_featureplots_pap.pdf", useDingbats=FALSE, width=15, height = 10)

#Transfer information to the 3dpi subset (for Fig. 6B)
the_name="3dpi"
figure_subset <- readRDS(paste0("figure_subset_", the_name, ".rds"))

df <- figure_subset@metadata@cell
df$cell <- rownames(df)
df2 <- HSV1all_ann@metadata@cell
df2$cell <- rownames(df2)
df <- left_join(df, df2 %>% select("cell", "inf_clusters_pap"), by="cell")
rownames(df) <- df$cell
identical(rownames(df), rownames(figure_subset@metadata@cell))
figure_subset@metadata@cell <- df
vrSpatialPlot(figure_subset, group.by = 'inf_clusters_pap', colors = cluster_colors_2, pt.size = 0.5, alpha = 1, plot.segments = TRUE, background.color = "black")+coord_fixed()
ggsave("figure_subset_3dpi_inf_clusters_pap.pdf", width=6, height=2)

#Next step: determine the "infection front". These are cells that are in infcluster1 (i.e. the cluster with the lowest amount viral reads, I just got infected)
#But we term these cells only "infection front" if there are in a group of at least 5 members, and every member of the group is less than 150 µm away from
#at least one other group member. By that, individual cells with some few random counts are not "infection" front, but their number is small.

#Make dataframe df_metadata_coords with cell information and coordinates
df_metadata_coords <- subset(HSV1all_ann@metadata@cell, annotation != "undefined")
df_metadata_coords$cell <- rownames(df_metadata_coords)
coords <- as.data.frame(vrCoordinates(HSV1all_ann[['Assay2']]))
coords$cell <- rownames(coords)
df_metadata_coords <- left_join(df_metadata_coords, coords, by="cell")
rownames(df_metadata_coords) <- df_metadata_coords$cell

#These values are of course adjustable, they sounded reasonable to us given how far cells are away from each other and how the infection front looks like
the_proxdist=150
the_groupsize=5
df_distances <- list()
#Cycle through all individual organoids
for (the_annotation in unique(df_metadata_coords$annotation)) {
  print(the_annotation)
  
  #Select only papillary fibroblasts in the annotation we are handling
  df_metadata_coords_temp <- df_metadata_coords %>% filter(annotation == the_annotation) %>% filter(grepl("Papillary", celltype))
  all_cells <- df_metadata_coords_temp %>% dplyr::select(c("x", "y"))
  #Make distance matrix of papillary fibroblasts
  d <- dist(df_metadata_coords_temp %>% dplyr::select(c("x", "y")))
  distance_table <- as.data.frame(as.matrix(d))
  
  #define the cells that make up the base cells, from which then distance regions are expanded (i.e. infcluster 1)
  points_coords <- df_metadata_coords_temp %>%
    filter(inf_clusters_pap=="infcluster1") %>%
    dplyr::select(c("x", "y", "cell"))
  dist_mat <- as.numeric(as.matrix(dist(points_coords %>% dplyr::select(c("x", "y")))))
  
  #Calculate adjacency matrix, everything that is closer than the_proxdist µm to each other is put into a group
  adj2 <- matrix(dist_mat < the_proxdist, nrow = nrow(points_coords))
  g <- graph_from_adjacency_matrix(adj2)
  points_coords$infareagroup <- igraph::components(g)$membership
  df_metadata_coords_temp <- left_join(df_metadata_coords_temp, points_coords %>% dplyr::select(c("cell", "infareagroup")), by="cell") %>% replace(., is.na(.), 0)
  rownames(df_metadata_coords_temp) <- df_metadata_coords_temp$cell
  
  #Groups with less than the_groupsize cells are filtered out
  small_clusts <- df_metadata_coords_temp %>% group_by(infareagroup) %>% tally() %>% filter(n<the_groupsize) %>% pull(infareagroup)
  df_metadata_coords_temp <- df_metadata_coords_temp %>% mutate(infareagroup=if_else(infareagroup %in% small_clusts, 10000, infareagroup))
  df_metadata_coords_temp$infareagroup <- paste0(df_metadata_coords_temp$infareagroup, "_", the_annotation)
  
  #metadata_basecells contains the infcluster1 papillary fibroblasts in the groups defined before, i.e. the infection fron
  metadata_basecells <- df_metadata_coords_temp %>% filter(!grepl("^0_|^10000_", infareagroup, perl=TRUE))
  
  #For every cell, calculate distance to the closest basecell
  expr2 <- list()
  i <- 1
  print(length(rownames(all_cells)))
  for (the_cell in rownames(all_cells)) {
    if (i %% 500 == 0) {
      print(i)
    }
    df <- distance_table %>% dplyr::select(all_of(the_cell)) %>% rownames_to_column(var="cell")
    colnames(df) <- c("cell", "distance")
    expr2[[the_cell]] <- merge(df, metadata_basecells, by="cell") %>% slice_min(order_by=distance, n=1)
    i = i+1
  }
  #Transfer distance information to df_metadata_coords_temp, which is the dataframe with the metadata of the annotation we are working on 
  distances <- do.call(rbind,expr2) %>% replace(., is.na(.), 0)
  distances$cell <- rownames(distances)
  df_metadata_coords_temp <- left_join(df_metadata_coords_temp, distances %>% dplyr::select(c("cell", "distance")), by="cell")
  df_metadata_coords_temp <- df_metadata_coords_temp %>% mutate(infectionfront = case_when(!grepl("^0_|^10000_", infareagroup, perl=TRUE) & inf_clusters_pap == "infcluster1" ~ "infectionfront",
                                                                                           distance <50 & inf_clusters_pap == "pap no virus" ~ "uninfected_next_to_infectionfront",
                                                                                           distance <50 & inf_clusters_pap == "infcluster2" ~ "infected_next_to_infectionfront",
                                                                                           inf_clusters_pap %in% c("infcluster3", "infcluster4", "infcluster5", "infcluster6", "infcluster7") ~ "later_infected",
                                                                                           .default = inf_clusters_pap))
  df_distances[[the_annotation]] <- df_metadata_coords_temp
}

df_distances_temp <- do.call(rbind, df_distances)

df_metadata <- HSV1all_ann@metadata@cell
df_metadata$cell <- rownames(df_metadata)
df_metadata <- left_join(df_metadata, df_distances_temp %>% dplyr::select(c("cell", "inf_clusters_pap", "infectionfront")), by="cell") %>% replace(., is.na(.), "not pap")
rownames(df_metadata) <- df_metadata$cell
saveRDS(df_metadata$infectionfront, "df_metadata_infectionfront.rds")
HSV1all_ann@metadata@cell$infectionfront <- df_metadata$infectionfront

cluster_colors_3 <- c("not pap"="#EEEEEE", "pap no virus" = "#FFBBBB", uninfected_next_to_infectionfront="#00FF00", infectionfront="turquoise1", infected_next_to_infectionfront = "blue", infcluster1 = "#352A87", infcluster2 = "#056EDE", infcluster3 = "#089BCE", 
                      infcluster4 = "#33B7A0", infcluster5 = "#A3BD6A", infcluster6 = "#F9BD3F", 
                      infcluster7 = "#F9FB0E", later_infected="#FF5555")

for (the_timepoint in time_levels) {
  vrSpatialPlot(subset(HSV1all_ann, timepoint==the_timepoint), group.by = 'infectionfront', colors = cluster_colors_3[unique(subset(HSV1all_ann, timepoint==the_timepoint)@metadata@cell$infectionfront)], pt.size = 0.5, alpha = 1, plot.segments = TRUE, background.color = "black")
  ggsave(paste0("HSV1all_infectionfront_", the_timepoint, "_new.pdf"), width=10, height=8)
}
vrSpatialPlot(HSV1all_ann, group.by = 'infectionfront', colors = cluster_colors_3[unique(HSV1all_ann@metadata@cell$infectionfront)], pt.size = 0.5, alpha = 1, plot.segments = TRUE, background.color = "black")+coord_fixed()
ggsave("HSV1all_infectionfront.pdf", width=11, height=8)

#Transfer information to the 3dpi subset (for Figure S20B)
the_name="3dpi"
figure_subset <- readRDS(paste0("figure_subset_", the_name, ".rds"))
df <- figure_subset@metadata@cell
df$cell <- rownames(df)
df2 <- HSV1all_ann@metadata@cell
df2$cell <- rownames(df2)
df <- left_join(df, df2 %>% select("cell", "infectionfront"), by="cell")
rownames(df) <- df$cell
figure_subset@metadata@cell <- df
vrSpatialPlot(figure_subset, group.by = 'infectionfront', colors = cluster_colors_3[unique(df$infectionfront)], pt.size = 0.5, alpha = 1, plot.segments = TRUE, background.color = "black")+coord_fixed()
ggsave("figure_subset_3dpi_infection_front_pap.pdf", width=6, height=2)




######################
#Differential expression at the infection front <-  Fig. S22C

#Assign infection front category 
table(HSV1all_ann@metadata@cell$infectionfront)

#Make count table with cell types and distgroups
expr <- list()
expr_stats <- list()
count_data <- vrData(HSV1all_ann)
for (the_timepoint in c("0dpi", "3dpi", "4dpi")) {
  for (the_params in c("pap no virus", "uninfected_next_to_infectionfront", "infectionfront", "infected_next_to_infectionfront")) {
    for (the_annotation in unique(subset(HSV1all_ann@metadata@cell, timepoint==the_timepoint)$annotation)) {
      namestring <- paste0(the_timepoint,'|', the_annotation,'|',the_params)
      cells <- HSV1all_ann@metadata@cell %>% filter(infectionfront == the_params) %>% filter(annotation==the_annotation) %>% rownames()
      expr_stats[[namestring]] <- length(cells)
      if (length(cells) > 10) {
        expr[[namestring]] <- rowSums(count_data[,cells])
      }
    }
  }
}

counts_infectionfront <- do.call(cbind,expr)

colData_infectionfront <- data.frame(
  timepoint=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\1',names(expr), perl=TRUE)),
  sample=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\2',names(expr), perl=TRUE)),
  infectionfront=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\3',names(expr), perl=TRUE)),
  group=names(expr),
  row.names=names(expr)
)

res <- list()

#Calculate differential expression
#First, compare all pap no virus and all infectionfront groups in 3 dpi/4 dpi to 0 dpi
for (the_timepoint in c("3dpi", "4dpi")) {
  for (the_params in c("pap no virus", "uninfected_next_to_infectionfront", "infectionfront", "infected_next_to_infectionfront")) {
    take_row <- rowSums(counts_infectionfront) > 0
    take_col <- (colSums(counts_infectionfront) > 0) &
      (colData_infectionfront[,'timepoint']=="0dpi" | (colData_infectionfront[,'timepoint']==the_timepoint & colData_infectionfront[,'infectionfront']==the_params))
    dds <- DESeqDataSetFromMatrix(countData=counts_infectionfront[take_row,take_col],
                                  colData=colData_infectionfront[take_col,,drop=FALSE],
                                  design=~timepoint)
    dds$timepoint <- relevel(dds$timepoint, ref = "0dpi")
    dds <- DESeq(dds, quiet=TRUE)
    for (the_comparison in grep("Intercept", resultsNames(dds), value=TRUE, invert = TRUE)) {
      res[[paste0(the_params, "_", the_timepoint, "_vs_0dpi")]] <- lfcShrink(dds,
                                                                             coef=the_comparison,
                                                                             type='apeglm',
                                                                             format="DataFrame",
                                                                             quiet = TRUE) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene_name') %>%
        dplyr::mutate(timepoint=the_timepoint, params=the_params, contrast=paste0(the_params, "_", the_timepoint, "_vs_0dpi"))
      
    }
  }
}
#Second, compare all infection front groups to pap no virus groups in 3 dpi and 4dpi
for (the_timepoint in c("3dpi", "4dpi")) {
  for (the_params in c("uninfected_next_to_infectionfront", "infectionfront", "infected_next_to_infectionfront")) {
    take_row <- rowSums(counts_infectionfront) > 0
    take_col <- (colSums(counts_infectionfront) > 0) &
      ((colData_infectionfront[,'timepoint']==the_timepoint & colData_infectionfront[,'infectionfront']=="pap no virus") | (colData_infectionfront[,'timepoint']==the_timepoint & colData_infectionfront[,'infectionfront']==the_params))
    #compare all distgroups against each other by cycling through reference levels
    dds <- DESeqDataSetFromMatrix(countData=counts_infectionfront[take_row,take_col],
                                  colData=colData_infectionfront[take_col,,drop=FALSE],
                                  design=~infectionfront)
    dds$infectionfront <- relevel(dds$infectionfront, ref = "pap no virus")
    dds <- DESeq(dds, quiet=TRUE)
    for (the_comparison in grep("Intercept", resultsNames(dds), value=TRUE, invert = TRUE)) {
      res[[paste0(the_params, "_", the_timepoint, "_vs_pap no virus_", the_timepoint)]] <- lfcShrink(dds,
                                                                                                     coef=the_comparison,
                                                                                                     type='apeglm',
                                                                                                     format="DataFrame",
                                                                                                     quiet = TRUE) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene_name') %>%
        dplyr::mutate(timepoint=the_timepoint, params=the_params, contrast=paste0(the_params, "_", the_timepoint, "_vs_pap no virus_", the_timepoint))
      
    }
  }
}

pseudobulk_infectionfront <- do.call(rbind,res)
saveRDS(pseudobulk_infectionfront, "pseudobulk_infectionfront.rds")
pseudobulk_infectionfront <- readRDS("new_pap_clustering/pseudobulk_infectionfront.rds")

contrasts_to_use <- c("infectionfront_3dpi_vs_0dpi", 
                      "uninfected_next_to_infectionfront_3dpi_vs_0dpi",
                      "pap no virus_3dpi_vs_0dpi")

df <- pseudobulk_infectionfront %>% filter(contrast %in% contrasts_to_use)
write.table(df, "pvalues_FigS20C.txt", sep="\t", row.names = FALSE, quote=FALSE)

genes <- pseudobulk_infectionfront %>%
  dplyr::filter(!grepl("HSV", gene_name)) %>% 
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  #only induced genes
  dplyr::filter(log2FoldChange>0) %>%
  dplyr::filter(!(is.na(padj)) & (padj < .05)) %>%
  dplyr::group_by(contrast, timepoint) %>%
  dplyr::slice_min(order_by=padj,n=50,with_ties=FALSE) %>%
  dplyr::pull(gene_name) %>%
  unique()
#Clustering to get order of genes and groups
df <- pseudobulk_infectionfront %>%
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  dplyr::filter(gene_name %in% genes) %>%
  dplyr::mutate(group=paste(timepoint,contrast,sep='_'))  %>%
  dplyr::select(group,gene_name,log2FoldChange) %>%
  spread(group,log2FoldChange) %>%
  tibble::column_to_rownames('gene_name')
df[is.na(df)] <- 0
hc <- hclust(dist(df))
gene.order <- row.names(df)[order.hclust(hc)]
#plot
minmax <- as.integer(max(max(df), abs(min(df))))+1
#Make circles later in Adobe Illustrator
ggplot(pseudobulk_infectionfront %>%
         dplyr::filter(contrast %in% contrasts_to_use) %>%
         dplyr::filter(gene_name %in% genes) %>%
         dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
         dplyr::mutate(gene=factor(gene_name,levels=gene.order)),
       aes(y=gene,x=factor(contrast, levels=contrasts_to_use),size=-log10(padj),color=log2FoldChange)) +
  geom_point(shape=16) +
  #scale_size_continuous(name='adj. p-value',
  #                      breaks=c(2,4,6),
  #                      labels=c("1e-2","1e-4", "<1e-6"),
  #                      limits=c(0,6)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-minmax,minmax),oob=scales::squish) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')
h=length(gene.order)/5+3.5
w=length(colnames(df))*0.7+2
ggsave("DE_infection_front.pdf", width=w, height = h, units="in", useDingbats=FALSE)
