library(data.table)
library(tidyverse)
library(EnsDb.Hsapiens.v86)


df <- fread("DMPs.bed", select = c("chrom", "chromStart", "chromEnd", "lfdr", "LOAD.minus.control")) # DO NOT CHANGE

DMPs.df <- df %>%
  dplyr::transmute(chrom = paste0("chr",str_pad(str_remove(chrom, "chr"), width=2, pad = "0")), 
                   position = chromStart, signed.logged.lfdr = -log10(lfdr) * sign(LOAD.minus.control)) %>%
  dplyr::arrange(chrom, position)




ensdb <- EnsDb.Hsapiens.v86 # DO NOT CHANGE

#--> Subset the genes and get in [chrom, gene] format
genes <- genes(ensdb)

# Just protein coding
genes <- genes[genes$gene_biotype %in% 'protein_coding']

# https://support.bioconductor.org/p/84355/#84356
genes.collapsed <- unlist(range(split(genes, ~symbol)))

# Pad: str_pad(seqnames, width=2, pad="0")
genes.df <- cbind(data.frame(genes.collapsed), gene = names(genes.collapsed)) %>%
  dplyr::filter(seqnames %in% 1:22) %>%
  dplyr::transmute(chrom = paste0("chr", str_pad(seqnames, width=2, pad="0")),
                   start, end, gene, strand) %>%
  dplyr::mutate(chrom = factor(chrom), gene = factor(gene)) %>%
  dplyr::arrange(chrom, start)


#--> Check our work
head(DMPs.df)
head(genes.df)

#--> Write out only what we need
fwrite(x = DMPs.df, file = "app.DMPs.csv")
fwrite(x = genes.df, file = "app.genes.csv")
