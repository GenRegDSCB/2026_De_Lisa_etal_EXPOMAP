library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(ggfortify)
library(circlize)
library(ComplexHeatmap)

##### Figure 2A - Volcano plot DE ----

DE_results <- read.delim("b", check.names = F)

#mean of mean to assign the size
DE_results <- DE_results %>% rowwise() %>% mutate(Mean=mean(c_across(c("Median.LF", "Median.VSL"))))

#differential expression
DE_results$diffexpressed <- "NO"
DE_results$diffexpressed[DE_results$log2FC > 0 & DE_results$pvalue < 0.05] <- "Higher"
DE_results$diffexpressed[DE_results$log2FC < 0 & DE_results$pvalue < 0.05] <- "Lower"

#write top 5 gene names
DE_results$delabel <- NA
DE_results$delabel[DE_results$diffexpressed != "NO"] <- DE_results$ID[DE_results$diffexpressed != "NO"]


#assign colors to categories
class.c <- (colorRampPalette(brewer.pal(11,"Spectral"))(11))
de_col <- c("Higher"=class.c[2],
            "Lower"=class.c[10])


ggplot(data=DE_results, aes(x=log2FC, y=-log10(pvalue), size=log10(Mean), col=diffexpressed)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed") +
  scale_color_manual(values =  de_col)  +
  geom_text_repel(aes(label = ifelse(pvalue < 0.01, delabel, "")))


#### Figure 2B - Boxplot miRNA levels ----
top10_nc <- read.delim("b", check.names = F, row.names = 1)
top10_nc <- as.data.frame(t(top10_nc))

top10_nc <- dplyr:::select(top10_nc, c("miR-3185", "miR-3929", "miR-4518-5p", "miR-425-3p", "miR-4436b-3p", "miR-874-3p"))

top10_nc <- as.data.frame(t(top10_nc))
top10_nc <- rownames_to_column(top10_nc, "ID")

top10_nc <- gather(top10_nc, Variable, Level, -ID)
top10_nc$Group <- ifelse(grepl("LF",top10_nc$Variable),"LF","VSL")
top10_nc$Group <- factor(top10_nc$Group, levels=c("VSL", "LF"))

top10_nc$ID <- factor(top10_nc$ID, levels=rev(unique(top10_nc$ID)))

class.c <- (colorRampPalette(brewer.pal(11,"Spectral"))(11))
class_col <- c("VSL"=class.c[10],
               "LF"=class.c[2])

ggplot(top10_nc, aes(x=Group, y=log10(Level+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  theme_bw() +
  facet_wrap(~ID, nrow  = 2, scales = "free")+
  labs(y="log10(Normalized Counts)", x="miRNA") +
  scale_color_manual(values =  class_col) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p") 


#### Figure 2C - Heatmap DE miRNAs ----
df <- read.delim("a", check.names = F, row.names = 1)

#to exclude: VSL-10M
df <- filter(df, to_exclude_miRNA == "no")

mirna_fc <- read.delim("b", row.names = 1, check.names = F)

mirna <- mirna_fc[, rownames(df)]

colorh <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(256))
group.c <- colorRampPalette(brewer.pal(9, "Set1"))(9)
sex.c <- colorRampPalette(brewer.pal(8, "Dark2"))(8)
bmi.c <- colorRamp2(c(1, 5, 30),  hcl_palette = "Reds", reverse = TRUE)
ipaq.c <- colorRampPalette(brewer.pal(9, "Oranges"))(9)
predimed.c <- colorRamp2(c(0, 5, 10), hcl_palette = "Greens", reverse = TRUE)
smoking.c <- colorRampPalette(brewer.pal(11, "Spectral"))(11)
floor.c <- colorRampPalette(brewer.pal(11, "Spectral"))(11)
passive.c <- colorRampPalette(brewer.pal(11, "Spectral"))(11)
toluene.c <- colorRamp2(c(0, 5, 20), hcl_palette = "YlOrRd", reverse = TRUE)
benzene.c <- colorRamp2(c(0, 5, 150), hcl_palette = "YlGnBu", reverse = TRUE)


group.cv <- c("VSL"=group.c[2], "LF"=group.c[1])
sex.cv <- c("M"=sex.c[1], "F"=sex.c[2])
smoking.cv <- c("no"=smoking.c[8], "former"=smoking.c[10], "yes"=smoking.c[11])
floor.cv <- c("<=2"=floor.c[1], ">2"=floor.c[3])
ipaq.cv <- c("Low"=ipaq.c[3], "Moderate"=ipaq.c[5], "High"=ipaq.c[7])
passive.cv <-c("no"=passive.c[7], "yes"=passive.c[2])


df$Group <- factor(df$Group, levels = c("VSL", "LF"))
df$IPAQ.score <- factor(df$IPAQ.score, levels = c("Low", "Moderate", "High"))
df$Smoking <- factor(df$Smoking, levels = c("no", "former", "yes"))

# Definition of the annotations
column_ha = HeatmapAnnotation(Class=df$Group,
                              Sex=df$Sex,
                              BMI=df$BMI,
                              Smoking=df$Smoking,
                              Passive_Smoke=df$living.smokers,
                              Floor=df$floor1,
                              PREDIMED=df$PREDIMED.score,
                              IPAQ=df$IPAQ.score,
                              "Toluene (personal)"=df$Toluene_personal,
                              "Benzene (personal)"=df$Benzene_personal,
                              col = list(Class = group.cv,
                                         Sex=sex.cv,
                                         BMI=bmi.c,
                                         IPAQ = ipaq.cv,
                                         PREDIMED = predimed.c,
                                         Smoking = smoking.cv,
                                         Floor = floor.cv,
                                         Passive_Smoke = passive.cv,
                                         Toluene=toluene.c,
                                         Benzene=benzene.c
                              ),
                              border = TRUE,
                              simple_anno_size = unit(0.3, "cm"),
                              annotation_name_gp = gpar(fontsize = 8),
                              annotation_legend_param = list(border="black"))


mirna1 <- t(mirna)
mirna1 <- scale(log10(mirna1+1))
mirna_fc <- ifelse(mirna_fc$log2FC < 0, "decreasing", "increasing")

ht_mirna <- Heatmap(t(mirna1),
                    row_names_gp = gpar(fontsize = 5),
                    column_names_gp = gpar(fontsize = 5),
                    top_annotation = column_ha,
                    cluster_columns = T,
                    cluster_rows = T,
                    name = "Z-score",
                    na_col = "white",
                    col=colorRamp2(breaks=c(-4,-2,0,2,4), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
                    border = "black",
                    column_names_rot = 45,
                    row_split = mirna_fc,
                    heatmap_legend_param = list(border = "black"),
                    clustering_method_rows = "ward.D2",
                    clustering_method_columns = "ward.D2",
                    clustering_distance_columns = "spearman",
                    clustering_distance_rows = "spearman")
ht_mirna


#### Figure 2D - Scatterplot Sp. correlations ----

pdata <- "a"

cdata <- "b"

metadata <- read.delim(pdata, check.names = F, row.names = 1)
miRNAs <-read.delim(cdata, row.names = 1, check.names = F)

metadata <- filter(metadata, to_exclude_miRNA == "no")

miRNAs <- miRNAs[, row.names(metadata)]

toplot <- cbind(metadata, t(miRNAs))

group.c <- colorRampPalette(brewer.pal(9, "Set1"))(9)
group.cv <- c("VSL"=group.c[2], "LF"=group.c[1])

ggplot(toplot, aes(x = log10(toplot$Xylene_personal+1),
                   y = log10(toplot$`miR-641`+1),
                  col=toplot$Group)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm")  +
  scale_color_manual(values = group.cv)+
  stat_cor(method="spearman") +
  labs(x = "log10(Xylene personal)", y="log10(miR-641)")

#### Figure 2E - Dotplot enrichment----

enriched_GOBP <- read.delim("c", check.names = F)

ggplot(enriched_GOBP, aes(x=-log10(adj.p.val), y=GS, size=log10(gene.tested), col=coef), size = 0.2) + 
  geom_dotplot(binaxis = 'x',
               stackdir = 'center',
               stackratio = 0.2,
               binwidth = 0.02) +
  theme_bw() +
  geom_point() +
  facet_wrap(~`DE miRNA class`, nrow=2, scale="free_y", space="free_y")+
  scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(11,"RdBu"))(11)), limits = c(-1,1)) +
  labs(y="GO Term", x="-log10(Adj. p-value)", size="log10(N° target genes)", col="Coefficient") 










