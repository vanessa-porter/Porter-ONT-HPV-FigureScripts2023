### ----------------------------------------------------------
### MULTI-INTEGRATION DNA METHYLATION IN THE HPV GENOME
### ----------------------------------------------------------

## Set Library Path
.libPaths("/path/to/libraries")

## Load Libraries
library(reshape2)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggsci)
library(cowplot)
library(valr)
library(pheatmap)
library(cutr)
library(pheatmap)
library(DSS)

### ----------------------------------------------------------
### LOAD DATA
### ----------------------------------------------------------

# Set file paths
htmcp_path <- "/path/to/htmcp/call_integration/output"
tcga_path <- "/path/to/tcga/call_integration/output"

### HPV METHYLATION

# Find files with HPV methylation by read
hpvMethFiles1 <- grep("hpv_reads_methylation.tsv", list.files(htmcp_path, recursive = TRUE, full.names = TRUE), value = TRUE)
hpvMethFiles2 <- grep("hpv_reads_methylation.tsv", list.files(tcga_path, recursive = TRUE, full.names = TRUE), value = TRUE)
hpvMethFiles <- c(hpvMethFiles1, hpvMethFiles2)

# Extract library names from file paths
lib <- gsub(paste0(htmcp_path, "/|", tcga_path, "/|/methylation/hpv_reads_methylation.tsv"), "", hpvMethFiles)

# Import data
hpvMeth <- lapply(hpvMethFiles, read.delim, header = FALSE, sep = "\t")
names(hpvMeth) <- lib
hpvMeth <- dplyr::bind_rows(hpvMeth, .id = "sample")

### BREAKPOINT PAIRS

# Find files with breakpoint pair information
hpvReadFiles1 <- grep("hpvSizeReads.txt", list.files(htmcp_path, recursive = TRUE, full.names = TRUE), value = TRUE)
hpvReadFiles2 <- grep("hpvSizeReads.txt", list.files(tcga_path, recursive = TRUE, full.names = TRUE), value = TRUE)
hpvReadFiles <- c(hpvReadFiles1, hpvReadFiles2)

# Extract library names from file paths
lib2 <- gsub(paste0(htmcp_path, "/|", tcga_path, "/|/hpv_size/hpvSizeReads.txt"), "", hpvReadFiles)

# Import data
hpvRead <- lapply(hpvReadFiles, read.delim, header = TRUE, sep = "\t")
names(hpvRead) <- lib2
hpvRead <- dplyr::bind_rows(hpvRead, .id = "sample") 

### EVENT LOCATIONS
events <- read.delim("tables/intSitesSummaryAll.txt", header = TRUE)
events$id <- paste0(events$Patient, "/", events$event)

### EVENT CATEGORIES
summary <- read.delim("tables/integrationTypes.txt", header = T)
summary$id <- paste0(summary$sample, "/", summary$event)

### HETEROLOGOUS EVENTS
integrants <- read.delim("tables/integrantBreakpairs.txt", header = T)
integrants$id <- paste0(integrants$sample, ":", integrants$event, ":", integrants$bp_pair_name)

### ----------------------------------------------------------
### ADD INFO TO THE DATAFRAME
### ----------------------------------------------------------

hpvMeth$bp_pair_name <- hpvRead$bp_pair_name[match(hpvMeth$V5, hpvRead$read)]
hpvMeth$bp_pair <- hpvRead$bp_pair[match(hpvMeth$V5, hpvRead$read)]

# remove reads that aren't in a bp pair
#hpvMeth <- hpvMeth[!is.na(hpvMeth$bp_pair),]

# add the event and event type
hpvMeth$event <- hpvRead$event[match(hpvMeth$V5, hpvRead$read)]
hpvMeth$integration.type <- summary$type[match(paste0(hpvMeth$sample, "/", hpvMeth$event), summary$id)]
hpvMeth$id <- paste0(hpvMeth$sample, ":", hpvMeth$event, ":", hpvMeth$bp_pair_name)
hpvMeth$integrant.category <- integrants$category[match(hpvMeth$id, integrants$id)]
hpvMeth$integrant.size.category <- integrants$size_category[match(hpvMeth$id, integrants$id)]

# subset to multibreakpoint types
mb <- hpvMeth
#mb <- hpvMeth[hpvMeth$integration.type == "multi-breakpoint_integration",]

### ----------------------------------------------------------
### CALCULATE THE METHYLATION FREQUENCY VALUES
### ----------------------------------------------------------

# name the missing columns 
colnames(mb)[2:12] <- c("chromosome","strand","start","end","read_name","log_lik_ratio",
                  "log_lik_methylated","log_lik_unmethylated","num_calling_strands","num_motifs","sequence")

# determine if the call is methylated
mb$methylation.call <- ifelse(mb$log_lik_ratio > 0, "methylated", "unmethylated")
#mb <- mb[complete.cases(mb),]

# calculate the methylation frequency per breakpoint pair across the HPV genome
freq <- mb %>%
  dplyr::group_by(id, sample,event,bp_pair_name, integration.type,integrant.category,integrant.size.category, chromosome, start, methylation.call) %>%
  dplyr::summarise(n.calls = n()) %>%
  dplyr::ungroup() %>%
  spread(methylation.call,  n.calls) %>%
  dplyr::mutate(methylated = replace_na(methylated, 0),
         unmethylated = replace_na(unmethylated, 0))

freq$total.calls <- freq$methylated + freq$unmethylated
freq$perc.methylated <- freq$methylated / freq$total.calls

# filter for calls with at least 4 reads
freq <- freq[freq$total.calls > 4,]

# only hpv chromosomes
freq <- freq[grep("HPV",freq$chromosome),]
freqSub <- freq[complete.cases(freq),]

### ----------------------------------------------------------
### PLOT ACROSS THE HPV GENOME
### ----------------------------------------------------------

# HPV Gene Info
hpvAll <- read.delim("/path/to/ref/hpv_gff/HPV_all.gff", header = F)
hpv <- hpvAll[hpvAll$V3 == "gene",]
hpvLengths <- hpvAll[hpvAll$V4 == 1 & hpvAll$V3 == "region",]
hpv_info <- strsplit(hpv$V9, ";")
hpv_info <- lapply(hpv_info, function(x){
  x <- x[grep("Name=", x)]
  x <- gsub("Name=", "", x)
  return(x)
})
hpv$gene <- unlist(hpv_info)
hpv <- hpv[,-c(2,3,6,7,8,9)]
colnames(hpv) <- c("genome", "start", "end","gene")

# Add Reading Frame
hpv$reading.frame <- ifelse(hpv$gene %in% c("E7","E1","E5","L2"), 1,
                            ifelse(hpv$gene %in% c("E6","E2","L1"), 2, 3))
hpv$position <- hpv$start + (hpv$end - hpv$start)/2

# Find and remap HPV genes
freqSub$hpv.gene <- NA 
freqSub$E6.start <- NA
for (i in 1:nrow(freqSub)) {
  gene <- hpv$gene[hpv$start < freqSub$start[i] & hpv$end > freqSub$start[i] & hpv$genome == freqSub$chromosome[i]]
  
  freqSub$hpv.gene[i] <- ifelse(length(gene) > 0, paste(gene, collapse = ","),
                                     ifelse(freqSub$start[i] < hpv$start[hpv$genome == freqSub$chromosome[i] & hpv$gene == "E6"] | freqSub$start[i] > hpv$end[hpv$genome == freqSub$chromosome[i] & hpv$gene == "L1"],
                                            "LCR", NA))
  
  freqSub$E6.start[i] <- ifelse(freqSub$start[i] < hpv$end[hpv$genome == freqSub$chromosome[i] & hpv$gene == "L1"], 
                                     freqSub$start[i] - hpv$start[hpv$genome == freqSub$chromosome[i] & hpv$gene == "E6"],
                                     -(hpv$start[hpv$genome == freqSub$chromosome[i] & hpv$gene == "E6"] + 
                                         (hpvLengths$V5[hpvLengths$V1 == freqSub$chromosome[i]] - freqSub$start[i])))
}

### ----------------------------------------------------------
### FIND THE MOST VARIABLY METHYLATED POSITIONS ACROSS THE HPV GENOME
### ----------------------------------------------------------

# make a matrix of the percent methylated
sub <- freqSub[,c("id","chromosome","start","perc.methylated")]
sub$position <- paste0(sub$chromosome, ":", sub$start)
wide <- pivot_wider(sub, names_from = id, values_from = perc.methylated)

# write a function to find variable positions in an HPV genome
variable_methyl <- function(HPV.type, n.positions){
  wideSub <- wide[wide$chromosome == HPV.type,]
  toRemove <- apply(wideSub,2, function(x) all(is.na(x)))
  wideSub <- wideSub[,!toRemove]
  mat <- as.matrix(wideSub[,-c(1:3)])
  rownames(mat) <- wideSub$position
  v <- names(sort(apply(mat,1,var, na.rm=T), decreasing=T)[1:n.positions])
  var <- mat[row.names(mat) %in% v, ] 
  # remove cols with zero variance
  colVar <- apply(var,2,var, na.rm=T)
  varRm <- c(names(colVar[colVar==0])[!is.na(names(colVar[colVar==0]))], names(colVar[is.na(colVar)]))
  var <- var[,!colnames(var) %in% varRm]
  return(var)
}

# variable HPV16 methylation
var16 <- variable_methyl("HPV16", 10)
var18 <- variable_methyl("HPV18", 10)
var45 <- variable_methyl("HPV45", 10)

#get the info for the variable positions across the genome
varProbes16 <- rownames(var16)
varPos16 <- wide[wide$position %in% varProbes16, 1:3]
varPos16$gene <- freqSub$hpv.gene[match(varPos16$position, paste0(freqSub$chromosome, ":", freqSub$start))]


varProbes18 <- rownames(var18)
varPos18 <- wide[wide$position %in% varProbes18, 1:3]
varPos18$gene <- freqSub$hpv.gene[match(varPos18$position, paste0(freqSub$chromosome, ":", freqSub$start))]

### ----------------------------------------------------------
### PLOT WHERE THE VARIABLE CPGS ARE
### ----------------------------------------------------------

## Create HPV Position Plots
hpvPlots <- lapply(c("HPV16", "HPV18"), function(i) {
  plot1 <- ggplot(hpv %>% filter(genome == i), aes(xmin = start, xmax = end, ymin = reading.frame - 0.4, ymax = reading.frame + 0.4)) + 
    geom_rect(aes(fill = gene), colour = "black") +
    scale_fill_lancet() +
    geom_text(mapping = aes(x = position, y = reading.frame, label = gene), size = 3, fontface = "bold")+
    scale_x_continuous(breaks = c(seq(0,8000,by=1000)), limits = c(0,8000)) +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none")
  
  #get the info for the variable positions across the genome
  var <- variable_methyl(i, 10)
  varProbes <- rownames(var)
  varPos <- wide[wide$position %in% varProbes, 1:3]
  
  # plot
  plot2 <- ggplot(varPos, 
                  aes(x = start, y = 1)) +
    geom_point(size = 5, shape= 4) +
    theme_bw() + 
    xlim(0,8000) +
    labs(y = NULL, x = "position")+
    theme(axis.text.x = element_text(colour = "black", size = 11),
          axis.text.y = element_blank(), 
          axis.title = element_text(colour = "black", size = 12, face = "bold"),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 12, colour = "black"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_line())
  
  return(list(plot1, plot2))
})

# hpv16
grid1 <- plot_grid(plotlist = hpvPlots[[1]], rel_heights = c(1,0.8), align = "v", ncol = 1, labels = c("HPV16"), label_size = 14)

# hpv18
grid2 <- plot_grid(plotlist = hpvPlots[[2]], rel_heights = c(1,0.8), align = "v", ncol = 1, labels = c("HPV18"), label_size = 14)

# put together
gridplot <- plot_grid(grid1, grid2, rel_heights = c(1,1), align = "v", ncol = 1)

ggsave(gridplot, filename = "figurePanels/figure5/hpv_variable_cpgs_position.pdf", width = 5.5, height = 4.5, units = "in")

### ----------------------------------------------------------
### PLOT THE METHYLATED VALUES AT THE VARIABLE POSITIONS IN A HEATMAP
### ----------------------------------------------------------

### make the annotation matrix for rows and columns 

# HPV16
anno_row16 <- as.data.frame(varPos16[,4])
rownames(anno_row16) <- varPos16$position

anno_row18 <- as.data.frame(varPos18[,4])
rownames(anno_row18) <- varPos18$position

freqSub$integration.type <- gsub("multi-breakpoint_integration", "multibreakpoint_integration", freqSub$integration.type)

### make the anno colours for the heatmap
ann_colors <- list(integrant.size.category = c(less1="#d3e7d4", over1="#709775", over2="#415d43", over3="#111d13"),
                   integrant.category = c(full="#003d5b", heterologous="#d1495b", incomplete="#fbfefb", partial="#30638e"),
                   integration.type = c(deletion_integration="#ded0e5", repeat_integration="#3292b4", duplication_integration="#1c6868", 
                                        ecDNA_integration="#512469", multibreakpoint_integration="#a6104c", unmatched_integration="#D3D3D3",
                                        translocation_integration = "#f3d86a"))


### write a function to make the heatmap

methylHeat <- function(var_mat, anno){
  mat <- t(var_mat)
  
  # identify columns in which there is not enough variance to calculate distance matrices
  giveNAs = which(is.na(as.matrix(dist(mat))),arr.ind=TRUE)
  tab = sort(table(c(giveNAs)),decreasing=TRUE)
  checkNA = sapply(1:length(tab),function(i){
    sum(is.na(as.matrix(dist(mat[-as.numeric(names(tab[1:i])),]))))
  })
  rmv = names(tab)[1:min(which(checkNA==0))]
  
  # remove those columns
  var_mat_rm = var_mat[,-as.numeric(rmv)]
  var_mat_rm <- var_mat_rm[order(rownames(var_mat_rm)),order(colnames(var_mat_rm))]
  
  # fix the annotation data
  anno_row <- anno
  anno_row$col <- "-"
  anno_row <- anno_row[rownames(var_mat_rm),]
  anno_col <- freqSub[freqSub$id %in% colnames(var_mat_rm), c("id", "integration.type", "integrant.category", "integrant.size.category")]
  anno_col <- as.data.frame(anno_col[!duplicated(anno_col),])
  rownames(anno_col) <- anno_col$id
  anno_col <- anno_col[order(colnames(var_mat_rm)),]
  colnames(var_mat_rm) == rownames(anno_col)
  rownames(var_mat_rm) == rownames(anno_row)
  anno_col <- anno_col %>% dplyr::select(integration.type, integrant.category, integrant.size.category)
  
  # plot the heatmap
  pheatmap(mat = var_mat_rm, cluster_cols = T, cluster_rows = T, show_rownames = FALSE, show_colnames = FALSE,
           annotation_col = anno_col, annotation_colors = ann_colors,
           color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(101), 
           #annotation_row = anno_row,
           method = "ward.D2", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean") 
}

### run function
heat16 <- methylHeat(var_mat = var16, anno=anno_row16)
heat18 <- methylHeat(var18, anno=anno_row18)

# get the cluster info

# rows
var16order <- rownames(var16)[order(rownames(var16))]
clust16row <- var16order[heat16$tree_row$order]
var18order <- rownames(var18)[order(rownames(var18))]
clust18row <- var18order[heat18$tree_row$order]

# cols
clust16col <- heat16$tree_col$labels[heat16$tree_col$order]
clust18col <- heat18$tree_col$labels[heat18$tree_col$order]
clust16 <- cutree(heat16$tree_col, k=3)
clust18 <- cutree(heat18$tree_col, k=3)

anno_col16 <- freqSub[freqSub$id %in% clust16col, c("id", "integration.type", "integrant.category", "integrant.size.category")]
anno_col16 <- anno_col16[!duplicated(anno_col16),]
anno_col16$cluster <- clust16[match(anno_col16$id, names(clust16))]

anno_col18 <- freqSub[freqSub$id %in% clust18col, c("id", "integration.type", "integrant.category", "integrant.size.category")]
anno_col18 <- anno_col18[!duplicated(anno_col18),]
anno_col18$cluster <- clust18[match(anno_col18$id, names(clust18))]

# add the methylation
var16t <- t(var16)
var16t <- var16t[rownames(var16t) %in% anno_col16$id,]
var16t <- as.data.frame(var16t[order(rownames(var16t)),])
anno_col16 <- anno_col16[order(anno_col16$id),]
anno_col16 <- cbind(anno_col16, var16t)

var18t <- t(var18)
var18t <- var18t[rownames(var18t) %in% anno_col18$id,]
var18t <- as.data.frame(var18t[order(rownames(var18t)),])
anno_col18 <- anno_col18[order(anno_col18$id),]
anno_col18 <- cbind(anno_col18, var18t)

# save
write.table(anno_col16, file = "tables/hpv16MethylationClusters.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(anno_col18, file = "tables/hpv18MethylationClusters.txt", quote = F, col.names = T, row.names = F, sep = "\t")

### ----------------------------------------------------------
### PLOT ACROSS THE HPV GENOME
### ----------------------------------------------------------

freqSub$event.id <- paste0(freqSub$sample, "/", freqSub$event)

plot1 <- ggplot(freqSub %>% filter(integration.type == "multi-breakpoint_integration"), 
       aes(x = E6.start, y = id, colour = perc.methylated)) +
  geom_point(size = 2) +
  theme_bw() + 
  #facet_grid(event.id ~ ., scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 0, linetype = 2, size = 1) +
  scale_color_distiller(palette = "RdBu") + 
  labs(x = "adjusted position", y = "sample", colour = "% methylated") +
  theme(axis.text.x = element_text(colour = "black", size = 11),
        axis.text.y = element_blank(), 
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line = element_line())

ggsave(plot1, filename = "multibreak_HPVmethyl.pdf", width = 8, height = 35, units = "in")

### ----------------------------------------------------------
### DIFFERENTIAL METHYLATION OF THE HPV GENOME
### ----------------------------------------------------------

# Create a list of the methylation calls in all HPV16 integrants
freqDSS <- freqSub[,c("id", "chromosome", "start", "total.calls", "methylated")]
colnames(freqDSS)[2:5] <- c('chr', 'pos', 'N', 'X')
freqDSS_integrated_HPV16 <- split(freqDSS[freqDSS$chr == "HPV16",], f = freqDSS$id[freqDSS$chr == "HPV16"])
freqDSS_integrated_HPV16 <- lapply(freqDSS_integrated_HPV16, subset, select = c('chr', 'pos', 'N', 'X'))

# Create a list of the methylation calls in the unintegrated HPV16 genomes
freqDSSU <- freq[freq$sample %in% c("TCGA-C5-A2M2", "TCGA-C5-A905", "HTMCP-03-06-02036"),c("id", "chromosome", "start", "total.calls", "methylated")]
colnames(freqDSSU)[2:5] <- c('chr', 'pos', 'N', 'X')
freqDSS_unintegrated_HPV16 <- split(freqDSSU, f = freqDSSU$id)
freqDSS_unintegrated_HPV16 <- lapply(freqDSS_unintegrated_HPV16, subset, select = c('chr', 'pos', 'N', 'X'))

# final list and specifying the groups
integrated <- names(freqDSS_integrated_HPV16)
unintegrated <- names(freqDSS_unintegrated_HPV16)
freqDSS_list <- c(freqDSS_integrated_HPV16,freqDSS_unintegrated_HPV16)

# run DSS 
DSObject <- makeBSseqData(freqDSS_list, names(freqDSS_list))
dmlTest.sm = DMLtest(DSObject, group1 = integrated, group2 = unintegrated, 
                     smoothing=TRUE)
dmr <- callDMR(dmlTest.sm, p.threshold = 0.01)

### ----------------------------------------------------------
### SAVE TABLES
### ----------------------------------------------------------

# variable positions
varPosSave <- rbind(varPos16, varPos18)


write.table(varPosSave, file = "tables/variableMethylPositions.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(dmr, file = "tables/dmrHPV16EpisomalIntegrated.txt", quote = F, col.names = T, row.names = F, sep = "\t")


