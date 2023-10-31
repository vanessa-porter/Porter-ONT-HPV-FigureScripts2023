.libPaths("/path.to.libraries")

library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

## --------------------------------------------------------------------------------------
## Load in the expression matrices
## --------------------------------------------------------------------------------------

htmcp_mat <- read.delim("/path/to/gdc_gene_expression/HTMCP_GDC_TPM-unstranded_matrix.txt", header = T)
tcga_mat <- read.delim("/path/to/gdc_gene_expression/TCGA_GDC_TPM-unstranded_matrix.txt", header = T)
allgenes <- read.delim("/path/to/gene_expression/event_locations_1Mb_window_genes.bed", header = F)
colnames(allgenes) <- c("chr", "start", "end", "sample", "event", "hpv.sites", "HPV.type", "gene.chr", "gene.start", "gene.end", "gene.id", "strand", "locus", "gene.type", "chr.locus")

# select the samples in the cohort
samples <- read.delim("/path/to/sourceTables/fig1a-2.txt", header = T)
tcga_samples <- samples$library[grep("TCGA", samples$library)]
htmcp_samples <- samples$library[grep("HTMCP", samples$library)]
htmcp_samples <- gsub("-",".",htmcp_samples)
tcga_samples <- gsub("-",".",tcga_samples)

# subset the matrix
htmcp_mat_sub <- htmcp_mat[,colnames(htmcp_mat) %in% c("Gene",htmcp_samples)]
tcga_mat_sub <- tcga_mat[,colnames(tcga_mat) %in% c("Gene",tcga_samples)]
all_mat_sub <- cbind(htmcp_mat_sub,tcga_mat_sub[,-1])
colnames(all_mat_sub)

# save genes
write.table(allgenes[,1:11],"tables/event_locations_1Mb_genes.bed", quote = F, sep = "\t", col.names = T, row.names = F)

## --------------------------------------------------------------------------------------
## Test specific genes
## --------------------------------------------------------------------------------------

# put the gene name in to test
gene <- "MYC"
test_sample <- gsub("-", ".",allgenes$sample[allgenes$gene.id == gene])
expr_mat <- all_mat_sub

p1 <- expr_mat %>%
  filter(Gene == gene) %>%
  gather(sample, tpm,-Gene) %>%
  mutate(colour = ifelse(sample %in% test_sample, "test", "others")) %>%
  ggplot(aes(x = colour, y = log10(tpm), colour = colour)) +
  geom_boxplot(outlier.shape = NA, position = "dodge") +
    stat_compare_means(method = "wilcox.test") +
  geom_jitter(aes(colour = colour, x = colour), height = 0, width = 0.2, size =3, alpha=0.5) +
    labs(x = NULL, y = "log10(TPM)") +
  scale_colour_manual(values = c("grey", "dark red")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 13),
          axis.title = element_text(size=14), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none")

ggsave(plot = p1, filename = paste0("/path/to/figurePanels/figure1/",gene, "_integrated_vs_not_integrated.pdf"), width = 2.7, height = 2.7)

genes <- c("NOTUM","TUBB2B","MAML3","NR4A1","NR4A3",
           "RAB40B","AP1M1","ARHGEF19","ATG9A","DNTTIP2")
test_samples <- c("HTMCP.03.06.02179","HTMCP.03.06.02179","HTMCP.03.06.02176","HTMCP.03.06.02109","HTMCP.03.06.02428",
                  "HTMCP.03.06.02179","HTMCP.03.06.02261","HTMCP.03.06.02175","HTMCP.03.06.02175","HTMCP.03.06.02054")

plist <- NULL

for (i in 1:length(genes)){
  p2 <- expr_mat %>%
    filter(Gene == genes[i]) %>%
    gather(sample, tpm,-Gene) %>%
    mutate(colour = ifelse(sample == test_samples[i], "test", "others")) %>%
    ggplot(aes(x = Gene, y = log10(tpm))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = colour), height = 0, width = 0.2, size =3, alpha=0.5) +
    labs(x = NULL, y = "log10(TPM)") +
    scale_colour_manual(values = c("grey", "dark red")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 13),
          axis.title = element_text(size=14), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none") 
  plist[[i]] <- p2
}

pgrid <- plot_grid(plotlist=plist, nrow = 2)

ggsave(plot = p2, filename = paste0("/path/to/figurePanels/figure6/", gene, "_single_integrated.pdf"), width = 2.7, height = 2.7)
ggsave(plot = pgrid, filename = "/path/to/figurePanels/figure6/outlierGeneExpressionBoxplot.pdf", width = 6.4, height = 3.8)







