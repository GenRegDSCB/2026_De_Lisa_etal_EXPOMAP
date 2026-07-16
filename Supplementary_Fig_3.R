## Supplementary Figure 3 - TSEA boxplots ----

toplot <- read.delim("a", check.names = F)
toplot$Group <- factor(out$Group, levels=c("VSL", "LF"))

#select the enrichment to plot
toplot <- filter(toplot, toplot$Signature %in% c("xx"))

ggplot(toplot, aes(x=Group, y=log10(Score1+1), col=Group))+
  geom_boxplot(width = 0.2, fill="white", outlier.alpha = 0 ) +
  geom_jitter(width = 0.1, shape=21, alpha=0.75, size=1)+
  theme_bw() +
  facet_wrap(~Signature, scales = "free", nrow = 1)+
  scale_color_manual(values =  class_col) +
  geom_pwc(method = "wilcox_test", label = "p.signif", hide.ns = "p")