# ecDNA ChIP P Values at GeneHancers
# Written by Vanessa Porter

## --------------------------------------------------------------------------------
# Load libraries  
## --------------------------------------------------------------------------------
.libPaths("/path/to/libraries")
library(plyr)
library(dplyr)
library(reshape2)
library(data.table)

## --------------------------------------------------------------------------------
# User functions
## --------------------------------------------------------------------------------

list_n_item <- function(string, split,  n){
  list <- strsplit(string, split)
  sapply(list, `[`, n)
}

## --------------------------------------------------------------------------------
# Pick GeneHancers
## --------------------------------------------------------------------------------

# read in the genehancers
gh <- read.delim("/path/to/htmcp/rproj/ecDNA_epigenome/GeneHancer_version_4-4.bed", header = F)
gh$id <- gsub("genehancer_id=", "", list_n_item(gh$V7, ";", 1))

# remove the genehancers that are within the ecDNA
enh_files <- grep("genehancer_overlap.bed", list.files("/path/to/htmcp/rproj/ecDNA_epigenome/samples", recursive = T, full.names = T), value = T)
enh_files_sub <- enh_files[-c(8,10,14)]
smGH <- NULL
for (i in 1:length(enh_files_sub)){
  smGH[[i]] <- read.delim(enh_files_sub[i], header = F)
}
smGH <- bind_rows(smGH)
smGH$id <- gsub("genehancer_id=", "", list_n_item(smGH$V7, ";", 1))
ghRM <- gh[!gh$id %in% smGH$id,]

# remove the Y chromosome (female)
ghRM <- ghRM[ghRM$V1 != "chrY",]

# pick 1000 random GeneHancers
ghRM <- ghRM[!duplicated(ghRM$id),]
random <- ghRM[sample(nrow(ghRM),size = 1000,replace=FALSE),]
length(unique(random$id))

# save the bed file for coverage
write.table(random, file = "/path/to/gh1000sample.bed", sep = "\t", col.names = F, row.names = F, quote = F)

### Calculate normalized coverage (RPM counts) on the enhancers using deeptools (in terminal):
# source deeptools
# multiBigwigSummary BED-file -b /path/to/htmcp/hg38_chip/h3k27ac/*.bw -o geneHancer_sample_cov_h3k27ac.npz --BED gh1000sample.bed --outRawCounts geneHancer_sample_counts_h3k27ac.txt --smartLabels -p 40

## --------------------------------------------------------------------------------
# Prepare the function
## --------------------------------------------------------------------------------

chip_pval <- function(cov, smplCov, sample){
  ###### Prepare the coverage matrices
  colnames(cov) <- gsub("X.", "", colnames(cov))
  colnames(cov) <- gsub("_h3k27ac.", "", colnames(cov))
  colnames(smplCov) <- gsub("X.", "", colnames(smplCov))
  colnames(smplCov) <- gsub("_h3k27ac.", "", colnames(smplCov))
  
  # add the ID to the dataframe and make the rownames
  cov$id <- gh$id[match(paste0(cov$.chr., ":", cov$start.), paste0(gh$V1, ":", gh$V2))]
  smplCov$id <- gh$id[match(paste0(smplCov$.chr., ":", smplCov$start.), paste0(gh$V1, ":", gh$V2))]
  
  # add the length of the enhancer and convert to RPKM
  cov$size <- cov$end. - cov$start.
  cov_long <- melt(cov[,4:ncol(cov)], id.vars = c("id", "size"))
  cov_long$RPKM <- cov_long$value / (cov_long$size/1000)
  smplCov$size <- smplCov$end. - smplCov$start.
  smplCov_long <- melt(smplCov[,4:ncol(smplCov)], id.vars = c("id", "size"))
  smplCov_long$RPKM <- smplCov_long$value / (smplCov_long$size/1000)
  
  # mark the event sample
  cov_long$event_sample <- ifelse(cov_long$variable == sample, T, F)
  smplCov_long$event_sample <- ifelse(smplCov_long$variable == sample, T, F)
  
  # add a small number to the RPKM as to not divide by zero
  cov_long$RPKM <- cov_long$RPKM + 0.00001
  smplCov_long$RPKM <- smplCov_long$RPKM + 0.00001
  
  ###### Calculate the foldchange and pval for each event
  
  # Make dfs for the control and test enhancers' log2FC
  fc_mat <- data.frame(id = unique(cov_long$id)) 
  fc_mat$log2FC <- NA
  fc_mat_test <- data.frame(id = unique(smplCov_long$id)) 
  fc_mat_test$log2FC <- NA
  
  # test samples
  for (id in fc_mat_test$id) {
    samp_cov <- smplCov_long[smplCov_long$event_sample == T & smplCov_long$id == id,"RPKM"]
    mean <- mean(smplCov_long[smplCov_long$event_sample == F & smplCov_long$id == id,"RPKM"])
    fc <- log2(samp_cov / mean)
    fc_mat_test[fc_mat_test$id == id, "log2FC"] <- fc
  }
  
  # control samples
  for (id in fc_mat$id) {
    samp_cov <- cov_long[cov_long$event_sample == T & cov_long$id == id,"RPKM"]
    mean <- mean(cov_long[cov_long$event_sample == F & cov_long$id == id,"RPKM"])
    fc <- log2(samp_cov / mean)
    fc_mat[fc_mat$id == id, "log2FC"] <- fc
  }
  
  # calculate the pvalues
  fc_mat_test$pvalue <- NA
  
  for (t in 1:nrow(fc_mat_test)) {
    fc <- fc_mat_test$log2FC[t]
    num_higher <- nrow(fc_mat[fc_mat$log2FC > fc,])
    p <- ifelse(num_higher == 0, 0.001, num_higher / 1000)
    fc_mat_test$pvalue[t] <- p
  }
  
  # BH correction
  fc_mat_test$padj <- p.adjust(fc_mat_test$pvalue, method = "BH")
  
  ## plot
  
  # add colour variable
  smplCov_long$colID <- ifelse(smplCov_long$variable == sample, "sample", "others")
  
  # plot
  plot <- ggplot(smplCov_long %>% arrange(colID), aes(x = id, y = RPKM))+
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(height = 0, width = 0.2, size = 2.5, alpha=0.7, aes(colour = colID)) +
    theme_bw() + 
    scale_color_manual(values = c("grey", "red")) +
    labs(x = NULL, y = "Normalized H3K27ac Coverage", fill = NULL, colour = NULL) +
    theme(axis.text.y = element_text(colour = "black", size = 10),
          axis.text.x = element_text(colour = "black", size = 10, angle = 40, hjust = 1, vjust = 1),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 12, colour = "black"),
          panel.grid = element_blank(), 
          axis.ticks.y = element_line(),
          axis.line = element_line(),
          legend.position = "none")
  
  # return the relevant values
  return(list(fc_mat_test,plot))
}

## --------------------------------------------------------------------------------
# Run the function on the two ecDNA samples
## --------------------------------------------------------------------------------

cov <- read.delim("/path/to/geneHancer_sample_counts_h3k27ac.txt", header = T)

# HTMCP-02179 
smplCov1 <- read.delim("/path/to/ecDNA_epigenome/samples/HTMCP-03-06-02179/gh_h3k27ac_matrix.txt", header = T)
sample1 <- "HTMCP.03.06.02179"
pvaldf1 <- chip_pval(cov = cov, smplCov = smplCov1, sample = sample1)
pvaldf1[[1]]["sample"]  <- sample1
#write.table(pvaldf1, "/path/to/genehancer_chip_cov/F46074/gh_chip_pvalue_results.txt", sep = "\t", col.names = F, row.names = F, quote = F)

# HTMCP-02182
smplCov2 <- read.delim("/path/to/ecDNA_epigenome/samples/HTMCP-03-06-02182/gh_h3k27ac_matrix.txt", header = T)
sample2 <- "HTMCP.03.06.02182"
pvaldf2 <- chip_pval(cov = cov, smplCov = smplCov2, sample = sample2)
pvaldf2[[1]]["sample"] <- sample2
#write.table(pvaldf1, "/path/to/genehancer_chip_cov/F46074/gh_chip_pvalue_results.txt", sep = "\t", col.names = F, row.names = F, quote = F)

# HTMCP-02320
smplCov3 <- read.delim("/path/to/ecDNA_epigenome/samples/HTMCP-03-06-02320/gh_h3k27ac_matrix.txt", header = T)
sample3 <- "HTMCP.03.06.02320"
pvaldf3 <- chip_pval(cov = cov, smplCov = smplCov3, sample = sample3)
pvaldf3[[1]]["sample"] <- sample3
#write.table(pvaldf1, "/path/to/genehancer_chip_cov/F46074/gh_chip_pvalue_results.txt", sep = "\t", col.names = F, row.names = F, quote = F)

# HTMCP-02170
smplCov4 <- read.delim("/path/to/ecDNA_epigenome/samples/HTMCP-03-06-02170/gh_h3k27ac_matrix.txt", header = T)
sample4 <- "HTMCP.03.06.02170"
pvaldf4 <- chip_pval(cov = cov, smplCov = smplCov4, sample = sample4)
pvaldf4[[1]]["sample"] <- sample4
#write.table(pvaldf1, "/path/to/genehancer_chip_cov/F46074/gh_chip_pvalue_results.txt", sep = "\t", col.names = F, row.names = F, quote = F)

# HTMCP-02320 (HPV52 ecDNA)
smplCov5 <- read.delim("/path/to/ecDNA_epigenome/samples/HTMCP-03-06-02040/gh_h3k27ac_matrix.txt", header = T)
sample5 <- "HTMCP.03.06.02040"
pvaldf5 <- chip_pval(cov = cov, smplCov = smplCov5, sample = sample5)
pvaldf5[[1]]["sample"] <- sample5
#write.table(pvaldf2, "/path/to/genehancer_chip_cov/F47779/gh_chip_pvalue_results.txt", sep = "\t", col.names = F, row.names = F, quote = F)

### save
pval_df <- rbind(pvaldf1[[1]], pvaldf2[[1]], pvaldf3[[1]], pvaldf4[[1]], pvaldf5[[1]])
pval_df$sample <- gsub('\\.', '-', pval_df$sample)
all_plot <- plot_grid(pvaldf2[[2]], pvaldf4[[2]], pvaldf5[[2]], nrow = 1, rel_widths = c(3.5,7,3))

write.table(pval_df,"/path/to/tables/eccDNA_chip_pvalues.txt", sep = "\t", col.names = T, row.names = F, quote = F)
ggsave(all_plot, filename = "/path/to/figurePanels/figure3/chip_ecDNA_cov.pdf", width = 7, height = 3)


