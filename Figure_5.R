library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

##### Fig. 5A - SCFA ----

a <- read.delim("a", check.names = F)

a <- gather(a, Variable, Level, -ID, -Group)

a$Group <- factor(a$Group, levels=c("VSL", "LF"))

class.c<-colorRampPalette(brewer.pal(11,"Spectral"))(11)

class_col <- c("VSL"=class.c[10],
               "LF"=class.c[3])

a$Variable <- factor(a$Variable, levels=c("Acetic acid", "Propionic acid", "Butyric acid"))

ggplot(a, aes(x=Group, y=log10(Level+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  theme_bw() +
  facet_wrap(~Variable, scale = "free", nrow =1) +
  labs(y="log10(Levels - ug/mL)", x="Residential Area")+
  scale_color_manual(values =  class_col) +
  geom_pwc(method="wilcox.test", label = "p.signif", hide.ns=T)

##### Fig. 5B - Scatter plot Propionate vs B. Bifidum ----

scfa <- read.delim("b", row.names = 1, check.names = F)

metadata <- read.delim("a", check.names = F, row.names = 1)
metadata <- filter(metadata, to_exclude_metag == "no")

scfa <- scfa[row.names(scfa)%in%row.names(metadata), ]
metadata <- metadata[row.names(scfa), ]

bac <- read.delim("c", check.names = F, row.names = 1)
bac <- bac[,row.names(scfa)]

toplot <- cbind(scfa, t(bac))
toplot <- select(toplot, c("Group", "Propionate", "Bifidobacterium_bifidum|t__SGB17256", "GGB9646_SGB15123|t__SGB15123", "Campylobacter_hominis|t__SGB19429", "Bifidobacterium_catenulatum|t__SGB17241"))

toplot <- gather(toplot, Species, Level, -Propionate, -Group)

toplot$Group <- factor(toplot$Group, levels = c("VSL", "LF"))

group.c <- colorRampPalette(brewer.pal(9, "Set1"))(9)
group.cv <- c("VSL"=group.c[2], "LF"=group.c[1])

ggplot(toplot, aes(x = log10(toplot$Level+1),
                   y =log10(toplot$Propionate+1))) +
  geom_point(aes(col=toplot$Group)) +
  theme_bw() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = group.cv)+
  facet_wrap(~Species, nrow = 2, scales = "free") +
  stat_cor(method="spearman")
