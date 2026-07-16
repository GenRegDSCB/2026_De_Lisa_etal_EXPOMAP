library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)


##### Supplementary Fig. 2A - Heatmap miRNA - tissue/cell levels - chromatine ----

df <- read.delim("c", row.names = 1, check.names = F)
df <- filter(df, to_exclude == "no")

#### HEATMAP chromatine state 
df_states <- select(df, c("Esophagus":"Rectal_smooth_muscle"))

states.c <- c("#C82226", "#E34A33", "#9CD569",  "#128A48", "#006837","#E3F399", "#FFFFBF", "#48A0B2",  "#5E4FA2","#9E0142", "#C52C4B",  "#FDCA78", "#656565", "#CFCFCF", "#FFFFFF" )


ht_states <- Heatmap(as.matrix(df_states),
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 7),
                     cluster_columns = F,
                     cluster_rows = T,
                     na_col = "white",
                     col= states.c,
                     name = "Chromatin states",
                     border = "black",
                     column_names_rot = 45,
                     column_split = c(rep("mucosa", 11), rep("muscle", 4)),
                     heatmap_legend_param = list(border = "black"),
                     clustering_method_rows = "ward.D2",
                     clustering_method_columns = "ward.D2",
                     clustering_distance_columns = "euclidean",
                     clustering_distance_rows = "euclidean")

ht_states

ro <- rev(row_order(ht_states))

ht_states <- Heatmap(as.matrix(df_states[ro,]),
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 7),
                     cluster_columns = F,
                     cluster_rows = F,
                     na_col = "white",
                     col= states.c,
                     name = "Chromatin states",
                     border = "black",
                     column_names_rot = 80,
                     column_split = c(rep("mucosa", 11), rep("muscle", 4)),
                     heatmap_legend_param = list(border = "black"),
                     clustering_method_rows = "ward.D2",
                     clustering_method_columns = "ward.D2",
                     clustering_distance_columns = "euclidean",
                     clustering_distance_rows = "euclidean")

ht_states

##### Barplot log2FC & median
df <- df[ro,]

df_fc <- c(df$log2FC)
df_median_VSL <- c(df$Median_VSL)
df_median_LF <- c(df$Median_LF)

row_ha  = rowAnnotation(log2FC = anno_barplot(df_fc),
                        annotation_name_rot = 80,
                        annotation_name_gp = gpar(fontsize = 8))

### Heatmap tissue count

df_tissue <- select(df, c("small_intestine":"colon"))

count.c <- colorRampPalette(brewer.pal(9, "Greens"))(9)

df2 <- (log10(df_tissue+1))
df2[is.na(df2)] <- NA


ht_mirna_tissue_count <- Heatmap(as.matrix(df2),
                                 row_names_gp = gpar(fontsize = 6),
                                 column_names_gp = gpar(fontsize = 7),
                                 #right_annotation = row_ha,
                                 cluster_columns = F,
                                 cluster_rows = F,
                                 name = "Tissue counts",
                                 na_col = "white",
                                 col= count.c,
                                 border = "black",
                                 column_names_rot = 80,
                                 heatmap_legend_param = list(border = "black"),
                                 clustering_method_rows = "ward.D2",
                                 clustering_method_columns = "ward.D2",
                                 clustering_distance_columns = "euclidean",
                                 clustering_distance_rows = "euclidean")

### Heatmap cells count 
df_cells <- select(df, c("Esophagus_epithelial_cell":"Smooth_muscle_cell_colon"))

cell.c <- colorRampPalette(brewer.pal(9, "Blues"))(9)

df3 <- (log10(df_cells+1))
df3[is.na(df3)] <- NA

ht_cells <- Heatmap(as.matrix(df3),
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize = 7),
                    right_annotation = row_ha,
                    cluster_columns = F,
                    cluster_rows = F,
                    name = "Cell counts",
                    na_col = "white",
                    col= cell.c,
                    border = "black",
                    column_names_rot = 80,
                    heatmap_legend_param = list(border = "black"),
                    clustering_method_rows = "ward.D2",
                    clustering_method_columns = "ward.D2",
                    clustering_distance_columns = "euclidean",
                    clustering_distance_rows = "euclidean")


ht_states + ht_mirna_tissue_count + ht_cells


##### Supplementary Fig. 2B-C - Alpha diversities (Shannon and Inverse Simpson) ----

toplot <- "a"
sdat <- read.delim(toplot, check.names = F, row.names = 1)
sdat <- filter(sdat, to_exclude_metag == "no")

sdat$Group <- factor(sdat$Group, levels=c("VSL", "LF"))

ch <- colorRampPalette(brewer.pal(11,"RdBu"))(11)
colorh2 <- c(ch[10],  ch[2])

ggplot(sdat, aes(x=Group, y=diversity_shannon, col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+ 
  theme_bw() +
  labs(y="alpha diversity - Shannon", x="Group") +
  scale_color_manual(values=colorh2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p")

ggplot(sdat, aes(x=Group, y=diversity_inverse_simpson, col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+ 
  theme_bw() +
  labs(y="alpha diversity - Inverse Simpson", x="Group") +
  scale_color_manual(values=colorh2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p")

