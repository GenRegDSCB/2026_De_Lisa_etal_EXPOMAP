path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

library("tidyverse")
library("ggpubr")
library("RColorBrewer")
library("gtsummary")

#### Table 1 ----
sdata <- read.delim("a", check.names=F)

sdata$Group  <- factor(sdata$Group, levels=c("VSL", "LF"))

table1 <- 
  sdata %>%
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})"),
              include = c(Group, Age, Sex, BMI, Smoking, living.smokers, Year.building, floor1),
              type = all_dichotomous() ~ "categorical", 
              by= Group) |> add_p()


table1 %>%
  as_gt() %>%
  gt::gtsave(filename = "1.0_Table1_Metadata.docx")


### Figure 1 ----

pdata <- "a"

data <- read.delim(pdata, check.names = F)

#### Figure 1B - BARPLOT PREDIMED ----

pred <- dplyr::select(data, c("ID", "Group", "PREDIMED_adherence"))
pred <- dplyr::filter(pred, !is.na(PREDIMED_adherence))

pred$Group <- factor(pred$Group, levels = c("VSL", "LF"))
pred$PREDIMED_adherence <- factor(pred$PREDIMED_adherence, levels = c("Low", "Medium-High", "High"))

ggplot(pred, aes(x = Group, fill = PREDIMED_adherence)) +
  geom_bar(position = "fill") +
  theme_bw() +
  scale_fill_brewer(palette = "YlGn") +
  labs(x="Group", y="Fraction", fill = "PREDIMED score", title = "PREDIMED score (p > 0.05)")
  
#### Figure 1C - BARPLOT IPAQ ----

ipaq <- dplyr::select(data, c("ID", "Group", "IPAQ.score"))
ipaq <- dplyr::filter(ipaq, !is.na(IPAQ.score))

ipaq$Group <- factor(ipaq$Group, levels = c("VSL", "LF"))
ipaq$IPAQ.score <- factor(ipaq$IPAQ.score, levels = c("Low", "Moderate", "High"))

ggplot(ipaq, aes(x = Group, fill = IPAQ.score)) +
  geom_bar(position = "fill") +
  theme_bw() +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(x="Residential Area", y="Fraction", fill = "IPAQ score", title = "IPAQ (p-value = 0.1332)")

#### Figure 1D - Boxplot SEDENTARIETY ----

sed <- dplyr::select(data, c("ID", "Group", "sedentariety.wd.minutes"))

sed <- gather(sed, Variable, Level, -ID, -Group)

sed$Group <- factor(sed$Group, levels=c("VSL", "LF"))

class.c<-colorRampPalette(brewer.pal(11,"Spectral"))(11)
class_col <- c("VSL"=class.c[10],
               "LF"=class.c[2])

sed$Variable <- factor(sed$Variable, levels=c("sedentariety.wd.minutes"))

ggplot(sed, aes(x=Group, y=(Level/60), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0 ) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  theme_bw() +
  facet_wrap(~Variable, labeller = as_labeller(c(
    "sedentariety.wd.minutes" = "Sedentariety")))+
  labs(y="Hours", x="Residential Area", title = "Sedentariety Levels (during working days)")+
  scale_color_manual(values =  class_col) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p")


#### Figure 1E - Boxplot BTEX (with RfC) ----

#outdoor / personal
a <- dplyr::select(data, c("ID","Group","Benzene_outdoor","Toluene_outdoor", "Ethylbenzene_outdoor", "Xylene_outdoor"))

a <- gather(a, Variable, Level, -ID, -Group)
a$Group <- factor(a$Group, levels=c("VSL", "LF"))

class.c<-colorRampPalette(brewer.pal(11,"Spectral"))(11)
class_col <- c("VSL"=class.c[10],
               "LF"=class.c[2])

a$Variable <- factor(a$Variable, levels=c("Benzene_outdoor","Toluene_outdoor", "Ethylbenzene_outdoor", "Xylene_outdoor"))

treshold <- c(0.03, 5, 1, 0.1)
treshold <- log10(treshold+1)

thr_df <- data.frame(Variable = levels(a$Variable), thr = treshold)
thr_df$Variable <- factor(thr_df$Variable, levels=c("Benzene_outdoor","Toluene_outdoor", "Ethylbenzene_outdoor", "Xylene_outdoor"))

ggplot(a, aes(x=Group, y=log10(Level+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  theme_bw() +
  facet_wrap(~Variable, nrow = 1)+
  scale_color_manual(values =  class_col) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p") +
  geom_hline(data = thr_df, mapping = aes(yintercept = thr), color = "red", linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE)


#### Figure 1F - Heatmap Exposure + Food intake ----

pdata <- "a"

data <- read.delim(pdata, check.names = F, row.names = 1)

library("circlize")
library("ComplexHeatmap")
library(RColorBrewer)

colorh <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(256))
group.c <- colorRampPalette(brewer.pal(9, "Set1"))(9)
sex.c <- colorRampPalette(brewer.pal(8, "Dark2"))(8)
bmi.c <- colorRamp2(c(1, 5, 30),  hcl_palette = "Reds", reverse = TRUE)
ipaq.c <- colorRampPalette(brewer.pal(9, "Oranges"))(9)
predimed.c <- colorRamp2(c(0, 5, 10), hcl_palette = "Greens", reverse = TRUE)
smoking.c <- colorRampPalette(brewer.pal(11, "Spectral"))(11)
floor.c <- colorRampPalette(brewer.pal(11, "Spectral"))(11)
passive.c <- colorRampPalette(brewer.pal(11, "Spectral"))(11)

group.cv <- c("VSL"=group.c[2], "LF"=group.c[1])
sex.cv <- c("M"=sex.c[1], "F"=sex.c[2])
smoking.cv <- c("no"=smoking.c[8], "former"=smoking.c[10], "yes"=smoking.c[11])
floor.cv <- c("<=2"=floor.c[1], ">2"=floor.c[3])
ipaq.cv <- c("Low"=ipaq.c[3], "Moderate"=ipaq.c[5], "High"=ipaq.c[7])
passive.cv <-c("no"=passive.c[7], "yes"=passive.c[2])


data$Group <- factor(data$Group, levels = c("VSL", "LF"))
data$IPAQ.score <- factor(data$IPAQ.score, levels = c("Low", "Moderate", "High"))
data$Smoking <- factor(data$Smoking, levels = c("no", "former", "yes"))

# Definition of the annotations
column_ha = HeatmapAnnotation(Class=data$Group,
                              Sex=data$Sex,
                              BMI=data$BMI,
                              Smoking=data$Smoking,
                              Passive_Smoke=data$living.smokers,
                              Floor=data$floor1,
                              PREDIMED=data$PREDIMED.score,
                              IPAQ=data$IPAQ.score,
                              col = list(Class = group.cv,
                                         Sex=sex.cv,
                                         BMI=bmi.c,
                                         IPAQ = ipaq.cv,
                                         PREDIMED = predimed.c,
                                         Smoking = smoking.cv,
                                         Floor = floor.cv,
                                         Passive_Smoke = passive.cv),
                              border = TRUE,
                              simple_anno_size = unit(0.3, "cm"),
                              annotation_name_gp = gpar(fontsize = 8),
                              annotation_legend_param = list(border="black"))

df1 <- dplyr::select(data, c("Benzene_personal","Toluene_personal","T-to-B_personal", "Ethylbenzene_personal", "Xylene_personal", "Benzene_outdoor","Toluene_outdoor", "T-to-B_outdoor", "Ethylbenzene_outdoor", "Xylene_outdoor"))


df1[is.na(df1)] <- NA
df1 <- t(scale(log10(df1+1)))

split_vec <- data$Group[match(colnames(df1), rownames(data))]
split_vec <- factor(split_vec, levels = c("VSL", "LF"))

ht_sensors <- Heatmap(df1,
                      row_names_gp = gpar(fontsize = 8),
                      column_names_gp = gpar(fontsize = 6),
                      top_annotation = column_ha,
                      cluster_columns = F,
                      cluster_rows = T,
                      name = "Z-score",
                      na_col = "white",
                      col=colorRamp2(breaks=c(-4,-2,0,2,4), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
                      border = "black",
                      column_names_rot = 45,
                      column_split = split_vec,
                      heatmap_legend_param = list(border = "black"),
                      clustering_method_rows = "ward.D2",
                      clustering_method_columns = "ward.D2",
                      clustering_distance_columns = "spearman",
                      clustering_distance_rows = "spearman")
## Heatmap Food Intake 

df2 <- dplyr::select(data, c("icecream.freq", "maize.freq", "pizza.freq", "freshcheese.freq", "syrupfruit.freq", "wine.freq", "redmeat.freq", "sausages.freq", "salad.freq", "tea.freq", "hotdrink.freq"))

df2[is.na(df2)] <- NA
df2 <- t(scale(log10(df2+1)))

ht_diet <- Heatmap(df2,
                   row_names_gp = gpar(fontsize = 8),
                   #row_labels = rownames,
                   column_names_gp = gpar(fontsize = 6),
                   cluster_columns = F,
                   cluster_rows = T,
                   name = "Z-score",
                   na_col = "white",
                   col=colorRamp2(breaks=c(-4,-2,0,2,4), colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
                   border = "black",
                   column_names_rot = 45,
                   heatmap_legend_param = list(border = "black"),
                   clustering_method_rows = "ward.D2",
                   clustering_method_columns = "ward.D2",
                   clustering_distance_columns = "spearman",
                   clustering_distance_rows = "spearman")

ht_sensors %v% ht_diet


