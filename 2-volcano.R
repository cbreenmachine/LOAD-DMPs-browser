library(data.table)
library(tidyverse)
library(plotly)

#TODO
#Want tooltip to popup for each point, where the tooltip gives chrom, and any possible genes
#Restrict to significant points only for tractability

global_size = 24

# show_col(pal_jama("default", alpha=0.95)(7))
colors <- pal_jama("default", alpha=1)(4)
# show_col(colors)

c_gene <- colors[4]
c_sig <- colors[2]
c_not_sig <- colors[3]


#--> Load in data (DO NOT CHANGE NAMES)
df <- fread("DMPs.bed") %>% dplyr::mutate(chrom = factor(chrom)) %>% dplyr::rename(start = chromStart, end = chromEnd)
genes.df <- fread("app.genes.csv") %>% dplyr::mutate(gene = factor(gene), chrom = factor(chrom))

# df$gene <- NA

gg <- 'BRCA1'
tmp <- dplyr::filter(genes.df, gene==gg)
chrom <- tmp$chrom
a <-tmp$start
b <- tmp$end
# ix <- which(df$chr == chrom & df$start > a & df$end < b)
df[chr==chrom & start>a & end<b, gene:=gg]



ix <- which(sample(df$lfdr) > 0.2)
ix.to.ignore <- ix[1:floor(length(ix) * 0.997)]

x_begin <- 1.2
x_end <- 2.2
yy <- 2.5

p <- df[!ix.to.ignore, ] %>%
  mutate(Significance = ifelse(lfdr < 0.05, "LFDR < 0.05", "LFDR > 0.05")) %>%
  ggplot(aes(
    x = LOAD.minus.control, y = -log10(lfdr), 
    color = Significance)) +
  geom_point()+
  theme_minimal() +
  scale_color_manual(values = c("#20854EFF", "grey")) +
  xlab("Effect size") +
  ylab(expression('-log'[10]*'(LFDR)')) +
  xlim(c(-2.5,2.5)) +
  ylim(c(0, 7)) +
  ggtitle("Volcano plot of DMPs")

text_size = 4.5

p.annotated <- p + annotate("segment", x = x_begin, y = yy, xend = x_end, yend = yy, 
                            color="black",
                            arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = x_begin, size=text_size, hjust=0, y = yy+0.3, label = "Hypermethylated") +
  annotate("segment", x = -x_begin, y = yy, xend = -x_end, yend = yy, 
           color="black",
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = -x_end, y = yy+0.3, size=text_size, hjust=0, label = "Hypomethylated") +
  theme(legend.position=c(.2,.75),
        plot.background = element_rect(fill = "white", colour = "white"),
        text = element_text(size=18),
        plot.title = element_text(hjust = 0),
        legend.background = element_rect(color="grey", fill=alpha("white", 1)))
