.libPaths("/path/to/libraries")

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(stats))
suppressMessages(library(tidyverse))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))

### -------------------------------------------------------------------------------
### READ IN FILES 
### -------------------------------------------------------------------------------

sum <- read.delim("/path/to/tables/intSitesSummaryAll.txt", header = T)
genes <- read.delim("/path/to/biomart_ensembl100_GRCh38.sorted.bed.gz", header = F)
hpvgene <- read.delim("/path/to/all_hpv_integration_events_distance_genes.bed", header = F)
colnames(hpvgene) <- c("chr", "start", "end", "sample", "event", "sites", "integration.type", "n.sites", "hpv.type",
                       "gene.chr", "gene.start", "gene.end","gene.id","gene.strand","locus","distance", "chrlocus")
htmcp_mat <- read.delim("/path/to/gdc_gene_expression/HTMCP_GDC_TPM-unstranded_matrix.txt", header = T)
tcga_mat <- read.delim("/path/to/dc_gene_expression/TCGA_GDC_TPM-unstranded_matrix.txt", header = T)

### -------------------------------------------------------------------------------
### MODIFY DATAFRAME
### -------------------------------------------------------------------------------
expr_sample <- htmcp_mat[,c(1,3)]
expr_sample <- expr_sample[expr_sample$HTMCP.03.06.02002 > 0,]

colnames(htmcp_mat)
colnames(tcga_mat)

### -------------------------------------------------------------------------------
### MODIFY DATAFRAME
### -------------------------------------------------------------------------------

# add the distance between the middle of the integration event and the middle of the gene
hpvgene$event.middle <- hpvgene$start+(hpvgene$end-hpvgene$start)
hpvgene$gene.middle <- hpvgene$gene.start+(hpvgene$gene.end-hpvgene$gene.start)
hpvgene$distance <- ifelse(hpvgene$gene.start < hpvgene$start & hpvgene$gene.end < hpvgene$start, -(hpvgene$start - hpvgene$gene.start), ifelse(
    hpvgene$gene.start > hpvgene$end, hpvgene$gene.start - hpvgene$end, 0))

# genes type - filter for protein coding
hpvgene$gene.type <- genes$V7[match(hpvgene$gene.id, genes$V4)]
hpvgene <- hpvgene %>% filter(gene.type == "protein_coding")
genes <- genes %>% filter(V7 == "protein_coding")

# take out genes with no expression 
max <- apply(htmcp_mat[,-1], 1, max) 
htmcpExpr <- htmcp_mat[max > 1 & htmcp_mat$Gene %in% genes$V4,] # Only genes expressed in at least 1 sample (RPKM > 5 in at least 1 sample)

max <- apply(tcga_mat[,-1], 1, max) 
tcgaExpr <- tcga_mat[max > 1 & tcga_mat$Gene %in% genes$V4,] # Only genes expressed in at least 1 sample (RPKM > 5 in at least 1 sample)

# melt the gene expression
htmcpEDF <- melt(htmcpExpr)
htmcpEDF$variable <- gsub("\\.", "-", htmcpEDF$variable)
htmcpEDF$value <- log10(htmcpEDF$value+0.001)
tcgaEDF <- melt(tcgaExpr)
tcgaEDF$variable <- gsub("\\.", "-", tcgaEDF$variable)
tcgaEDF$value <- log10(tcgaEDF$value+0.001)

### -------------------------------------------------------------------------------
### CREATE FUNCTION
### -------------------------------------------------------------------------------

genehtmcp <- hpvgene[grep("HTMCP", hpvgene$sample),]
genehtmcp <- genehtmcp[genehtmcp$gene.id %in% htmcpExpr$Gene,]
genetcga <- hpvgene[grep("TCGA", hpvgene$sample),]
genetcga <- genetcga[genetcga$gene.id %in% tcgaExpr$Gene,]

find_outlier <- function(x) {
    return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

expr_significance <- function(gene,sample,expr_mat){
    g <- gene
    s <- sample
    expr <- expr_mat
    
    # pick 1000 random genes
    #random <- sample(expr$gene, 1000)
    #rFC <- numeric()
    #for (i in 1:1000) {
    #    rg <- random[i]
    #    exprR <- expr[expr$gene == rg,]
    #    med <- median(exprR$value)
    #    med <- ifelse(med == 0, 0.00001,med)
    #    sample <- exprR$value[exprR$variable == s]
    #    sample <-  ifelse(sample == 0, 0.00001,sample)
    #    fc <- log2(sample/med)
    #    rFC <- c(rFC,fc)
    #}
    
    # get the expression fold change from the median of the gene
    exprG <- expr[expr$Gene == g,]
    med <- median(exprG$value)
    sam <- exprG$value[exprG$variable == s]
    fc <- log2(sam/med)
    exprG$outlier <- ifelse(find_outlier(exprG$value), "outlier", "not")
    outlier <- exprG$outlier[exprG$variable == s]
    exprG$zscores <- scale(exprG$value)
    z <- exprG$zscores[exprG$variable == s]
    
    # put the values together
    #n <- ifelse(fc > 0, sum(rFC > fc), sum(rFC < fc))
    #p <- n/1000
    #p <- ifelse(p == 0, 0.001,p)
    df <- data.frame(sample.expr=sam, median.expr=med, log2FC=fc, is.outlier=outlier, zscore = z)
    return(df)
}

# run the function
pvalH2 <- list()
for (i in 1:nrow(genehtmcp)) {
    df <- expr_significance(gene = genehtmcp$gene.id[i], sample = genehtmcp$sample[i], expr_mat = htmcpEDF)
    pvalH2[[i]] <- df
    print(i)
}

pvalHdf <- dplyr::bind_rows(pvalH2)
genehtmcpP <- cbind(genehtmcp,pvalHdf)
genehtmcpPFilt <- genehtmcpP[genehtmcpP$sample.expr > 0,]


# run the function
pvalT <- list()
for (i in 1:nrow(genetcga)) {
  df <- expr_significance(gene = genetcga$gene.id[i], sample = genetcga$sample[i], expr_mat = tcgaEDF)
  pvalT[[i]] <- df
  print(i)
}


pvalTdf <- dplyr::bind_rows(pvalT)
genetcgaP <- cbind(genetcga,pvalTdf)
genetcgaPFilt <- genetcgaP[genetcgaP$sample.expr > 0 ,]

# put together
allgenes <- rbind(genehtmcpPFilt,genetcgaPFilt)
allgenes <- allgenes[complete.cases(allgenes),]

### -------------------------------------------------------------------------------
### ADD ASE DATA
### -------------------------------------------------------------------------------

ase_files1 <- read.delim("/path/to/htmcp/ase/output/results_files", header = F)
ase_files2 <- read.delim("/path/to/tcga/ase/output/results_files", header = F)
ase_files <- c(ase_files1$V1, ase_files2$V1)
samples <- gsub("/path/to/htmcp/ase/|/path/to/tcga/ase/output/results_files|/mBASED/MBASED_expr_gene_results.txt", "", ase_files)
ase_dfs <- lapply(ase_files, read.delim, header = TRUE, sep = "\t")
names(ase_dfs) <- samples
ase <- bind_rows(ase_dfs, .id = "sample_name")
ase$aseResults <- ifelse(ase$padj < 0.05 & ase$majorAlleleFrequency > 0.65, "ASE", "BAE")

allgenes$ase_result <- ase$aseResults[match(paste0(allgenes$sample, allgenes$gene.id), paste0(ase$sample_name, ase$gene))]
allgenes$ase_MAF <- ase$majorAlleleFrequency[match(paste0(allgenes$sample, allgenes$gene.id), paste0(ase$sample_name, ase$gene))]
allgenes$event.loci <- paste0(allgenes$sample, ":", allgenes$sites)

upgenes <- allgenes[allgenes$is.outlier == "outlier" & allgenes$log2FC > 0,]
downgenes <- allgenes[allgenes$is.outlier == "outlier" & allgenes$log2FC < 0,]
outliergenes <- allgenes[allgenes$is.outlier == "outlier",]

testgene <- "NR4A3"
test_sample <- allgenes$sample[allgenes$gene.id == testgene]

p1_ase <- ase %>%
    filter(gene == testgene) %>%
    mutate(colour = ifelse(sample_name %in% test_sample, "test", "others")) %>%
    ggplot(aes(x = gene, y = majorAlleleFrequency)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(fill = colour, colour = aseResults), height = 0, width = 0.2, size =4, shape = 21, stroke = 1, alpha = 0.8) +
    scale_colour_manual(values = c("black", "grey"), na.value = "white")+
    labs(x = NULL, y = "RNA major allele frequency") +
    scale_fill_manual(values = c("grey", "dark red")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 13, colour = "black"),
          axis.title = element_text(size=14, colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none") 
p1_ase

ggsave(plot = p1_ase, filename = paste0("figurePanels/figure6/",testgene, "_MAF_ASE.pdf"), width = 2.7, height = 2.7)

### -------------------------------------------------------------------------------
### FIND OUTLIERS
### -------------------------------------------------------------------------------

# find the recurring genes
t <- table(upgenes$gene.id)
upgenes_reocurring <-  t[t > 2]
upgenes_reocurring <- upgenes_reocurring[order(upgenes_reocurring, decreasing = T)]


### -------------------------------------------------------------------------------
### LOOK AT NR4A3 REGULATED GENES
### -------------------------------------------------------------------------------

NR4A3_upgenes <- read.delim("/path/to/gene_expression/NR4A3_regulated_genes_up.txt", header = F)

test_sample <- "HTMCP.03.06.02428"
gene <- NR4A3_upgenes$V1[1:50]

p2 <- expr_mat %>%
    filter(Gene %in% gene) %>%
    gather(sample, tpm,-Gene) %>%
    mutate(colour = ifelse(sample %in% test_sample, "test", "others")) %>%
    arrange(colour) %>%
    ggplot(aes(x = Gene, y = log10(tpm))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = colour), height = 0, width = 0.2, size =2, alpha=0.5) +
    labs(x = NULL, y = "log10(TPM)") +
    scale_colour_manual(values = c("grey", "dark red")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 10),
          axis.title = element_text(size=14), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none") 
p2

### -------------------------------------------------------------------------------
### PLOT THE EXPRESSION CHANGES
### -------------------------------------------------------------------------------

genehtmcpPFilt <- genehtmcpPFilt2
genehtmcpPFilt$log2FC.L <- "> 3"
genehtmcpPFilt$log2FC.L[genehtmcpPFilt$log2FC <= 3] <- NA

labs <- outliergenes[outliergenes$log2FC > 3 | abs(outliergenes$log2FC) > 0.5 & outliergenes$ase_result == "ASE",]
labs <- labs[complete.cases(labs[,1:3]),]

p_outliers <- ggplot(outliergenes,aes(x=distance,y=event.loci,fill=log2FC, colour = ase_result)) +
    geom_point(size=4, shape = 21, stroke = 1) +
    scale_colour_manual(values = c("black", "grey"), na.value = "white")+
    geom_text_repel(data=labs, 
                    aes(x=distance,y=event.loci, label = gene.id), max.overlaps = 20, colour = "black", size = 3) +
    xlim(-1000000,1000000)+
    facet_grid(integration.type ~ ., scales = "free_y", space = "free_y") +
    scale_fill_distiller(palette = "RdBu", limits = c(-6,6)) + 
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(size = 13, colour = "black"),
          axis.title = element_text(size=14, colour = "black", face = "bold"))
p_outliers
ggsave(plot = p_outliers, filename = "figurePanels/figure6/outlier_genes_ase_expression.pdf", width = 5, height = 6)
          
### -------------------------------------------------------------------------------
### SAVE TABLES
### -------------------------------------------------------------------------------

write.table(outliergenes, file = "tables/outlierGeneExpression.txt", quote = F, col.names = T, row.names = F, sep = "\t")

