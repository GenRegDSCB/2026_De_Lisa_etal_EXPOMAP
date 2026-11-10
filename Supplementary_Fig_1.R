##### Supplementary Table 1 ----
library("gtsummary")
library(tidyverse)

sdata <- read.delim("0_Metadata_metrics.tsv", check.names=F)

sdata$Group  <- factor(sdata$Group, levels=c("VSL", "LF"))

table1 <- 
  sdata %>%
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})"),
              include = c(Group, pa.intense.minutes, MET_intense, pa.moderate.days, pa.moderate.minutes, MET_moderate, walking.days, walking.minutes, MET_walking, MET_tot, IPAQ.score, sedentariety.wd.minutes, sedentariety.we.minutes, oil.consumption, vegetables.pd, fruit.pd, redmeat.pd, butter.pd, legumes.pw, fish.pw, sweets.pw, dried.fruit.pw, whitemeat.preference, saute.pw, PREDIMED.score, pasta.freq, rice.freq, whole.grains.freq, pizza.freq, maize.freq, soup.freq, redmeat.freq, whitemeat.freq, legumes.freq, agedcheese.freq, freshcheese.freq, fish.freq, shellfish.freq, ham.freq, sausages.freq, cookedveg.freq, tomatoes.freq, salad.freq, orangeveg.freq, potatoes.freq, eggs.freq, fruit.freq, syrupfruit.freq, driedfruit.freq, seeds.freq, wine.freq, beer.freq, liquors.freq, softdrinks.freq, fruitjuice.freq, coffee.freq, milk.freq, sugarmilk.freq, vegandrink.freq, tea.freq, hotdrink.freq, biscuits.freq, sweetsnack.freq, sweets.freq, cereals.freq, chocolate.freq, sugars.freq, icecream.freq, food.diversity, Benzene_personal, Toluene_personal, TtoB_personal, Ethylbenzene_personal, Xylene_personal, Total_PERSONAL, Benzene_outdoor, Toluene_outdoor, TtoB_outdoor, Ethylbenzene_outdoor, Xylene_outdoor, Total_OUTDOOR),
              type = list(redmeat.freq ~ "continuous",
                          pa.moderate.days ~ "continuous",
                          walking.days ~ "continuous",
                          redmeat.freq ~ "continuous",
                          pizza.freq ~ "continuous",
                          legumes.freq ~ "continuous",
                          agedcheese.freq ~ "continuous",
                          freshcheese.freq ~ "continuous",
                          eggs.freq ~ "continuous",
                          seeds.freq ~ "continuous",
                          all_continuous() ~ "continuous"), 
              by= Group) |> add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2))


table1 %>%
  as_gt() %>%
  gt::gtsave(filename = "Supplementary_Table1_Metadata.docx")

#### Supp. Fig. 1A - Boxplot Food intake ----

data <- read.delim("0_Metadata_metrics.tsv", check.names=F)

diet <- dplyr::select(data, c("ID", "Group", "icecream.freq", "maize.freq", "pizza.freq", "syrupfruit.freq", "wine.freq", "redmeat.freq", "sausages.freq", "salad.freq", "tea.freq", "hotdrink.freq", "freshcheese.freq"))

diet <- gather(diet, Variable, Level, -ID, -Group)

diet$Group <- factor(diet$Group, levels=c("VSL", "LF"))

class.c<-colorRampPalette(brewer.pal(11,"Spectral"))(11)
class_col <- c("VSL"=class.c[10],
               "LF"=class.c[2])

diet$Variable <- factor(diet$Variable, levels=c("maize.freq", "freshcheese.freq", "redmeat.freq", "sausages.freq", "salad.freq", "icecream.freq", "pizza.freq", "syrupfruit.freq", "tea.freq", "hotdrink.freq", "wine.freq"))

ggplot(diet, aes(x=Group, y=log10(Level+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0 ) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  theme_bw() +
  facet_wrap(~Variable, nrow = 2, labeller = as_labeller(c(
    "redmeat.freq" = "Red Meat", 
    "maize.freq" = "Maize",
    "freshcheese.freq" = "Fresh Cheese",
    "pizza.freq" = "Pizza",
    "icecream.freq" = "Icecream",
    "syrupfruit.freq" = "Fruit Syrup",
    "wine.freq" = "Wine",
    "sausages.freq" = "Cured meat",
    "salad.freq" = "Salad",
    "tea.freq" = "Hot Tea",
    "hotdrink.freq" = "Tisane"
  )), scales = "free")+
  labs(y="log10(Food Consumption Frequency)",
       x="Residential Area", 
       title = "Food Consumption Frequency per week")+
  scale_color_manual(values =  class_col) +
  geom_pwc(method = "wilcox_test", 
           label = "p.signif", hide.ns = "p", y.position = 0.8)


#### Supp. Fig. 1B - Corrplt BTEX ----

library(Hmisc)
library(corrplot)
library(RColorBrewer)

cov <- select(data, "Toluene_outdoor", "Benzene_outdoor", "Ethylbenzene_outdoor", "Xylene_outdoor", "Benzene_personal", "Toluene_personal", "Ethylbenzene_personal", "Xylene_personal")

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


#### Supp. Fig. 1C - Forest plot GLM miRNAs ----

out <- read.delim("1_GLM_DE_miRNA_Group_Sex_Covs.txt", check.names = F)

#Retrieve the log2FC values for each miRNA tested 
fc <- read.delim(input_miRNA_log2FC, check.names = F)

out <- left_join(out, fc, by = "miRNA")

out$Significance <- "ns"
out$Significance <- ifelse(out$p.value < 0.05, "p<0.05", out$Significance)
out$Significance <- ifelse(out$q < 0.05, "q<0.05", out$Significance)

sig.cv <- c("q<0.05"="firebrick4","p<0.05"="purple2","ns"="grey30")

class.c <- (colorRampPalette(brewer.pal(11,"Spectral"))(11))
out$color <- ifelse(out$log2FC < 0, class.c[10], class.c[2])

out$label <- paste0(
  "<span style='color:", out$color, "'>",
  out$miRNA,
  "</span>"
)                    

ggplot(out, aes(y=estimate, x= reorder(label, estimate, FUN = mean), shape=covariate, col=Significance)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(width = 0.8), 
                  size = 0.2,
                  linewidth = 0.2) +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual("Significance", values = sig.cv) +
  scale_shape_manual(values = c(
    "BMI" = 15,
    "food.diversity" = 17,
    "PREDIMED.score" = 18,
    "sedentariety.wd.minutes" = 3,
    "smoke_exp" = 4
  )) +
  theme_bw() +
  theme(axis.text.x = element_markdown(size=6, angle = 65, hjust = 1))

#### Suppl. Fig. 1D ----

#S.F. 1D is the Elbow plot computed during clustering analysis (see Fig. 2C)

