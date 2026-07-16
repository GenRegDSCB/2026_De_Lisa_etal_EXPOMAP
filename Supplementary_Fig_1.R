path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

pdata <- "a"

data <- read.delim(pdata, check.names = F)

#### Supp. Fig. 1A - Boxplot Food intake ----

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




