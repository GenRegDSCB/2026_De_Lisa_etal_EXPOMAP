library(RColorBrewer)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(caret)
library(ggthemes)
library(ggpubr)
library(ggrepel)
library(ggplot2)

## Figure 4A - Richness plot ----

toplot <- "a"
sdat <- read.delim(toplot, check.names = F, row.names = 1)
sdat <- filter(sdat, to_exclude_metag == "no")

sdat$Group <- factor(sdat$Group, levels=c("VSL", "LF"))

ch <- colorRampPalette(brewer.pal(11,"RdBu"))(11)
colorh2 <- c(ch[10],  ch[2])

ggplot(sdat, aes(x=Group, y=richness, col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0)+
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+ 
  theme_bw() +
  labs(y="Richness", x="Group") +
  scale_color_manual(values=colorh2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p")

## Figure 4B - Beta Diversity ----

cdat <- read.delim("d", check.names = F, row.names = 1)

sdat <- read.delim("a", check.names = F, row.names = 1)
sdat <- filter(sdat, sdat$to_exclude_metag == "no")

cdat <- cdat[, row.names(sdat)]

sdat$Group <- factor(sdat$Group, levels = c("VSL", "LF"))

#Computing beta diversity 
beta_dist <- vegdist(t(cdat), index = "bray")

mds <- metaMDS(beta_dist)
mds_data <- as.data.frame(mds$points)
mds_data$SampleID <- rownames(mds_data)
mds_data <- cbind(mds_data, sdat)

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = factor(Group))) +
  geom_point(cex=2)+
  labs(color="Class")+
  scale_color_manual(values=colorh2) +
  theme_bw() +
  stat_ellipse()

## Figure 4C - Corrplot correlations -----

library(Hmisc)
library(corrplot)
library(RColorBrewer)

metadata <- read.delim("a", check.names = F, row.names = 1)
metadata <- filter(metadata, to_exclude_metag == "no")

cov <- select(metadata, "BMI","MET_tot", "sedentariety.wd.minutes", "redmeat.freq", "wine.freq", "PREDIMED.score", "food.diversity", "Toluene_outdoor", "Benzene_outdoor", "Ethylbenzene_outdoor", "Xylene_outdoor", "Benzene_personal", "Toluene_personal", "Ethylbenzene_personal", "Xylene_personal", "microb.load", "richness", "diversity_inverse_simpson", "diversity_shannon","evenness_simpson")

#pairwise correlation
res_p <- rcorr(as.matrix(cov), type="spearman")

colcor <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(256))

corrplot(res_p$r,
         type="lower",
         p.mat = res_p$P,
         tl.cex=0.5,
         tl.col="black",
         tl.srt = 65,
         col=colcor,
         order = "hclust",
         hclust.method="ward.D2",
         #addCoef.col=“black”,
         number.cex=0.7,
         diag = F,
         insig = "label_sig",
         pch.cex = 0.7,
         sig.level = c(0.001, 0.01, 0.05))

## Figure 4D - Boxplot selected species ----

toplot <- read.delim("d", check.names = F)

#top3 up and down
selected_species <- c("GGB34228_SGB72916|t__SGB72916", "GGB9759_SGB15370|t__SGB15370", "GGB9719_SGB15273|t__SGB15273", "Roseburia_hominis|t__SGB4659", "Bifidobacterium_bifidum|t__SGB17256", "Roseburia_faecis|t__SGB4925")

toplot <- filter(toplot, ID %in% selected_species)

toplot$ID <- factor(toplot$ID, levels=c(selected_species))

toplot <- gather(toplot, Subject, Level, -ID)

toplot$Group <- ifelse(grepl("LF",toplot$Subject),"LF","VSL")

toplot$Group <- factor(toplot$Group, levels=c("VSL", "LF"))

class.c <- (colorRampPalette(brewer.pal(11,"Spectral"))(11))

class_col <- c("VSL"=class.c[10],
               "LF"=class.c[2])

ggplot(toplot, aes(x=Group, y=log2(Level+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  scale_color_manual(values = class_col) +
  labs(x="Class",
       y="Relative abundance") +
  theme_bw()+
  facet_wrap(~ID, nrow = 2, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=8))+
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p")

## Figure 4E - Heatmap differentially abundant species ----

##### HEATMAP DA species ----

#table containing SIAMCAT outputs, zoe ranking and correlations
ra_zoe_corr <- read.delim("d", check.names = F, row.names = 1)
df <- read.delim("a", check.names = F, row.names = 1)

df <- filter(df, to_exclude_metag == "no")

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
richness.c <- colorRamp2(c(150, 250, 400), hcl_palette = "Purples", reverse = TRUE)


group.cv <- c("VSL"=group.c[2], "LF"=group.c[1])
sex.cv <- c("M"=sex.c[1], "F"=sex.c[2])
smoking.cv <- c("no"=smoking.c[8], "former"=smoking.c[10], "yes"=smoking.c[11])
floor.cv <- c("<=2"=floor.c[1], ">2"=floor.c[3])
ipaq.cv <- c("Low"=ipaq.c[3], "Moderate"=ipaq.c[5], "High"=ipaq.c[7])
passive.cv <-c("no" = passive.c[7], "yes" = passive.c[2])


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
                              Toluene=df$Toluene_personal,
                              Benzene=df$Benzene_personal,
                              Richness=df$richness,
                              col = list(Class = group.cv,
                                         Sex=sex.cv,
                                         BMI=bmi.c,
                                         IPAQ = ipaq.cv,
                                         PREDIMED = predimed.c,
                                         Smoking = smoking.cv,
                                         Floor = floor.cv,
                                         Passive_Smoke = passive.cv,
                                         Toluene=toluene.c,
                                         Benzene=benzene.c,
                                         Richness = richness.c
                              ),
                              border = TRUE,
                              simple_anno_size = unit(0.3, "cm"),
                              annotation_name_gp = gpar(fontsize = 8),
                              annotation_legend_param = list(border="black"))

daspecies <- ra_zoe_corr[, row.names(df)]
daspecies <- scale(t(log10(daspecies+1)))

ht <- Heatmap(t(daspecies),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 4),
              top_annotation = column_ha,
              #right_annotation = row_ha,
              cluster_columns = T,
              cluster_rows = T,
              name = "Z-score",
              col=colorRamp2(breaks=c(-4,-2,0,2,4), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
              border = "black",
              column_names_rot = 55,
              heatmap_legend_param = list(border = "black"),
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              clustering_distance_columns = "spearman",
              clustering_distance_rows = "spearman")


#### Heatmap log2FC ----
df_ra1 <- select(ra_zoe_corr, "fc")

ht_ra_fc <- Heatmap(as.matrix(df_ra1),
                        row_names_gp = gpar(fontsize = 6),
                        column_names_gp = gpar(fontsize = 8),
                        cluster_columns = F,
                        cluster_rows = T,
                        name = "log2FC",
                        col=colorRamp2(breaks=c(-2,0,2),colors=c(colorh[64], "white", colorh[192])),
                        width = 3,
                        border = "black",
                        column_names_rot = 45,
                        heatmap_legend_param = list(border = "black"),
                        clustering_method_rows = "ward.D2",
                        clustering_method_columns = "ward.D2",
                        cell_fun = function(j,i,x,y,w,h,fill) {
                          if(ra_zoe_corr$`p.val`[i] < 0.001) {
                            grid.text("***", x, y, vjust = 0.75)
                          } else if(ra_zoe_corr$`p.val`[i] < 0.01) {
                            grid.text("**", x, y, vjust = 0.75)
                          }else if(ra_zoe_corr$`p.val`[i] < 0.05) {
                            grid.text("*", x, y, vjust = 0.75)
                          }})

#### Barplot correlation with Toluene levels ----

#p=1 equals to p<0.05 and coherence
pval.cv <- ra_zoe_corr$`p<0.05`
pval.cv <- c(ifelse(ra_zoe_corr$`p<0.05` == 0, "white", "grey"))

row_ha  = rowAnnotation(Rho = anno_barplot(ra_zoe_corr$Rho,
                                           gp = gpar(fill = pval.cv)),
                        annotation_name_rot = 50,
                        annotation_name_gp = gpar(fontsize = 6))

#### Heatmap ZOE ranking ----
df_ra2 <- select(ra_zoe_corr, "ZOE MB Health-rank")
df_ra2[is.na(df_ra2)] <- -1
ht_ra_zoe <- Heatmap(as.matrix(df_ra2),
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 8),
                     right_annotation = row_ha, 
                     cluster_columns = F,
                     cluster_rows = T,
                     name = "ZOE Ranking",
                     col= colorRamp2(c(0,1), hcl_palette = "RdPu", reverse = TRUE),
                     width = 2,
                     border = "black",
                     column_names_rot = 45,
                     heatmap_legend_param = list(border = "black"),
                     clustering_method_rows = "ward.D2",
                     clustering_method_columns = "ward.D2")

ht + ht_ra_fc + ht_ra_zoe

## Figure 4F - Volcano plot ----

data <- read.delim("a", check.names = F, row.names = 1)

rownames(data) <- str_replace(rownames(data), ":.*$", "")

data$diffexpressed <- "NO"
data$diffexpressed[data$log2FC > 0] <- "Higher"
data$diffexpressed[data$log2FC < 0] <- "Lower"

class.c <- (colorRampPalette(brewer.pal(11,"Spectral"))(11))
de_col <- c("Higher"=class.c[2],
            "Lower"=class.c[10])

ggplot(data, aes(x=log2FC, y=-log10(pvalue), col=diffexpressed)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed") +
  scale_color_manual(values =  de_col) +
  labs(y="-log10(p-value)") +
  geom_text_repel(aes(label = ifelse(pvalue < 0.05, rownames(data), "")))


## Figure 4G - Boxplot selected pathways ----

toplot <- read.delim("a", check.names = F, row.names = 1)

toplot <- gather(pdat_filtered, Pathway, Level, -Group)

toplot$Group <- factor(toplot$Group, levels=c("VSL", "LF"))

selected_path <- c("PWY-6435: 4-hydroxybenzoate biosynthesis III (plants)", "PWY1ZNC-1: assimilatory suLF-ate reduction IV")

class.c <- (colorRampPalette(brewer.pal(11,"Spectral"))(11))
class_col <- c("VSL"=class.c[10],
               "LF"=class.c[2])

toplot <- filter(toplot, Pathway %in% selected_path)

ggplot(toplot, aes(x=Group, y=log10(Level+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0 ) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  scale_color_manual(values = class_col) +
  facet_wrap(~Pathway, scales = "free")+
  labs(x="Class",
       y="Relative abundance") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=8))+
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p", y.position = 2.8)

## Fig. 4H - TSEA boxplot ----

toplot <- read.delim("a", check.names = F)
toplot$Group <- factor(out$Group, levels=c("VSL", "LF"))

toplot <- filter(toplot, toplot$Signature %in% c("Butyrate producers"))

ggplot(toplot, aes(x=Group, y=log10(Score1+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0 ) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  theme_bw() +
  facet_wrap(~Signature, scales = "free", nrow = 1)+
  scale_color_manual(values =  class_col) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p")

