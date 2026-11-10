library(tidyverse)
library(ggplot2)
library(RColorBrewer)

## Supplementary Fig. 3A - Forest plot GLM SGBs ----

out <- read.delim("1_GLM_DE_SGB_Group_Sex_Covs.txt", check.names = F)

#Retrieve the log2FC values for each SGB tested 
fc <- read.delim(input_SGB_log2FC, check.names = F)

out <- left_join(out, fc, by = "Species")

out$Significance <- "ns"
out$Significance <- ifelse(out$p.value < 0.05, "p<0.05", out$Significance)
out$Significance <- ifelse(out$q < 0.05, "q<0.05", out$Significance)

sig.cv <- c("q<0.05"="firebrick4","p<0.05"="purple2","ns"="grey30")

class.c <- (colorRampPalette(brewer.pal(11,"Spectral"))(11))
out$color <- ifelse(out$log2FC < 0, class.c[10], class.c[2])

out$label <- paste0(
  "<span style='color:", out$color, "'>",
  out$Species,
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

## Supplementary Fig. 3B ----

#S.F. 3B is the Elbow plot computed during clustering analysis (see Fig. 4E)

## Supplementary Figure 3 - TSEA boxplots ----

toplot <- read.delim("0_TSEA_scoring.tsv", check.names = F)
toplot$Group <- factor(out$Group, levels=c("VSL", "LF"))

#select the enrichment type 
toplot <- filter(toplot, toplot$Type == "xx")

#select the signature to plot
toplot <- filter(toplot, toplot$Signature %in% c("xx"))

ggplot(toplot, aes(x=Group, y=log10(Score1+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0 ) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  theme_bw() +
  facet_wrap(~Signature, scales = "free", nrow = 1)+
  scale_color_manual(values =  class_col) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p")