prodigal -p meta -i <viral_input> -d <viral_output>

diamond blastx -d <diamond_nr_database> \
-q <prodigal_genecalls> \
-o <diamond_output> \
--outfmt 6 qseqid evalue bitscore stitle -k 1 \
--more-sensitive

```{r}
library(ggfittext)
library(gggenes)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

virus <- read.csv("<VIRSorter_diamond_output>", sep = "\t")

colourCount = length(unique(virus$COG_cat))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot(virus, aes(xmin = Start, xmax = End, y = Genome, fill = COG_cat, forward = Strand, label = role)) +
  geom_gene_arrow(arrowhead_height = unit(10, "mm"), arrowhead_width = unit(3, "mm"), arrow_body_height = unit(10, "mm") ) +
  geom_gene_label(align = "middle") +
  facet_wrap(~ Genome, scales = "free", ncol = 1) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(colourCount)) +
  theme_genes() +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=5))

ggsave("<gggenes_output>", width = 420, height = 297, units = "mm", limitsize = FALSE, dpi = 300)
```
