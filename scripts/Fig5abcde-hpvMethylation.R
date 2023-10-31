### ----------------------------------------------------------
### LOAD LIBRARIES
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

### ----------------------------------------------------------
### READ IN DATAFRAMES
### ----------------------------------------------------------

# Event locations
events <- read.delim("/path/to/tables/intSitesSummaryAll.txt", header = TRUE)
events$id <- paste0(events$library, "/", events$event)
events <- events[events$event != "none",]

# Event categories 
summary <- read.delim("/path/to/tables/integrationTypes.txt", header = T)
summary$id <- paste0(summary$sample, "/", summary$event)

### ----------------------------------------------------------
### HPV GENE POSITION INFO
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

## Create HPV Position Plots
plot1 <- lapply(c("HPV16", "HPV18", "HPV45"), function(i) {
  ggplot(hpv %>% filter(genome == i), aes(xmin = start, xmax = end, ymin = reading.frame - 0.4, ymax = reading.frame + 0.4)) + 
    geom_rect(aes(fill = gene), colour = "black") +
    scale_fill_lancet() +
    geom_text(mapping = aes(x = position, y = reading.frame, label = gene), size = 3, fontface = "bold")+
    scale_x_continuous(breaks = c(seq(0,8000,by=1000)), limits = c(0,8000)) +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          axis.line.x = element_line(size = 0.5),
          axis.text.x = element_text(size = 12, colour = "black"),
          legend.position = "none")
})

### ----------------------------------------------------------
### INTEGRATION EVENT METHYLATION DATA
### ----------------------------------------------------------

# Set file paths
htmcp_path <- "/path/to/htmcp/call_integration/output"
tcga_path <- "/path/to/tcga/call_integration/output"

# Find event methylation files
eventMethFiles1 <- grep("event_methyl_freq.tsv", list.files(htmcp_path, recursive = TRUE, full.names = TRUE), value = TRUE)
eventMethFiles2 <- grep("event_methyl_freq.tsv", list.files(tcga_path, recursive = TRUE, full.names = TRUE), value = TRUE)
eventMethFiles <- c(eventMethFiles1, eventMethFiles2)

# Extract library names from file paths
ontLib <- gsub(paste0(htmcp_path, "/|", tcga_path, "/|/methylation/event_methyl_freq.tsv"), "", eventMethFiles)

# Import ont sites
eventMeth <- lapply(eventMethFiles, read.delim, header = TRUE, sep = "\t")
names(eventMeth) <- ontLib
eventMeth <- dplyr::bind_rows(eventMeth, .id = "library")

# Subset HPV data
eventMethHPV <- eventMeth[grep("HPV", eventMeth$chromosome),]

# Find and remap HPV genes
eventMethHPV$hpv.gene <- NA 
eventMethHPV$E6.start <- NA
for (i in 1:nrow(eventMethHPV)) {
  gene <- hpv$gene[hpv$start < eventMethHPV$start[i] & hpv$end > eventMethHPV$start[i] & hpv$genome == eventMethHPV$chromosome[i]]
  
  eventMethHPV$hpv.gene[i] <- ifelse(length(gene) > 0, paste(gene, collapse = ","),
                                     ifelse(eventMethHPV$start[i] < hpv$start[hpv$genome == eventMethHPV$chromosome[i] & hpv$gene == "E6"] | eventMethHPV$start[i] > hpv$end[hpv$genome == eventMethHPV$chromosome[i] & hpv$gene == "L1"],
                                            "LCR", NA))
  
  eventMethHPV$E6.start[i] <- ifelse(eventMethHPV$start[i] < hpv$end[hpv$genome == eventMethHPV$chromosome[i] & hpv$gene == "L1"], 
                                     eventMethHPV$start[i] - hpv$start[hpv$genome == eventMethHPV$chromosome[i] & hpv$gene == "E6"],
                                     -(hpv$start[hpv$genome == eventMethHPV$chromosome[i] & hpv$gene == "E6"] + 
                                         (hpvLengths$V5[hpvLengths$V1 == eventMethHPV$chromosome[i]] - eventMethHPV$start[i])))
}

# Reorder factor levels
eventMethHPV$chromosome <- factor(eventMethHPV$chromosome, levels = c("HPV16", "HPV18", "HPV45", "HPV82", "HPV52", "HPV58", "HPV33", "HPV31","HPV68", "HPV73","HPV35"))

### ----------------------------------------------------------
### UNINTEGRATED HPV METHYLATION
### ----------------------------------------------------------

# unintegrated samples
unint_samples <- summary$sample[summary$type == "no_detected_integration"]

# Find event methylation files
hpvMethFiles1 <- grep("hpv_methylation_frequency.tsv", list.files(htmcp_path, recursive = TRUE, full.names = TRUE), value = TRUE)
hpvMethFiles2 <- grep("hpv_methylation_frequency.tsv", list.files(tcga_path, recursive = TRUE, full.names = TRUE), value = TRUE)
hpvFiles <- c(hpvMethFiles1, hpvMethFiles2)

# Extract library names from file paths
hpvMethNames <- gsub(paste0(htmcp_path, "/|", tcga_path, "/|/methylation/hpv_methylation_frequency.tsv"), "", hpvFiles)

# Import
hpvMeth <- lapply(hpvFiles, read.delim, header = FALSE, sep = "\t")
names(hpvMeth) <- hpvMethNames
hpvMeth <- dplyr::bind_rows(hpvMeth, .id = "library")
hpvMeth <- hpvMeth[hpvMeth$library %in% unint_samples,c(1:3,5:8)]
colnames(hpvMeth) <- c("library","chromosome","start","methylated","unmethylated","total.calls","perc.methylated")
hpvMeth$integration.type <- "no_integration_detected"
hpvMeth$id <- hpvMeth$library
hpvMeth$event <- "event0"

# Find and remap HPV genes
hpvMeth$hpv.gene <- NA 
hpvMeth$E6.start <- NA
for (i in 1:nrow(hpvMeth)) {
  gene <- hpv$gene[hpv$start < hpvMeth$start[i] & hpv$end > hpvMeth$start[i] & hpv$genome == hpvMeth$chromosome[i]]
  
  hpvMeth$hpv.gene[i] <- ifelse(length(gene) > 0, paste(gene, collapse = ","),
                                     ifelse(hpvMeth$start[i] < hpv$start[hpv$genome == hpvMeth$chromosome[i] & hpv$gene == "E6"] | hpvMeth$start[i] > hpv$end[hpv$genome == hpvMeth$chromosome[i] & hpv$gene == "L1"],
                                            "LCR", NA))
  
  hpvMeth$E6.start[i] <- ifelse(hpvMeth$start[i] < hpv$end[hpv$genome == hpvMeth$chromosome[i] & hpv$gene == "L1"], 
                                hpvMeth$start[i] - hpv$start[hpv$genome == hpvMeth$chromosome[i] & hpv$gene == "E6"],
                                     -(hpv$start[hpv$genome == hpvMeth$chromosome[i] & hpv$gene == "E6"] + 
                                         (hpvLengths$V5[hpvLengths$V1 == hpvMeth$chromosome[i]] - hpvMeth$start[i])))
}


### ----------------------------------------------------------
### GET THE METHYLATION DATA ADJACENT TO HPV INTEGRATION
### ----------------------------------------------------------

# Get a list of files containing "methylregions.bed"
files <- c(
  grep("methylregions.bed", list.files("/path/to/htmcp/call_integration", recursive = TRUE, full.names = TRUE), value = TRUE),
  grep("methylregions.bed", list.files("/path/to/tcga/call_integration", recursive = TRUE, full.names = TRUE), value = TRUE)
)

# Extract sample names from file paths
names1 <- gsub("/path/to/|htmcp|tcga|/call_integration/output/|methylation/regions/|/methylregions.bed", "", files)

# Read data from files into a list
mrList <- lapply(files, read.delim, header = FALSE)
names(mrList) <- names1

# Combine all data into one data frame and add an "id" column
mr <- dplyr::bind_rows(mrList, .id = "id")

# Add event categories and other information to the methyl regions data frame
mr$integration.type <- summary$type[match(mr$id, summary$id)]
mr$hpv.type <- summary$HPV.type[match(mr$id, summary$id)]

# Subset event locations data frame to only the events in the methyl regions data frame
eventsSub <- events[events$id %in% mr$id,]

# Add 5' and 3' positions to the event locations data frame
pos <- events %>%
  dplyr::filter(event != "none") %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(up_pos = min(pos),
            down_pos = max(pos))

for (i in 1:nrow(pos)){
  pos$HPVdir[i] <- ifelse(events$HPVpos[events$pos == pos$up_pos[i]] > events$HPVpos[events$pos == pos$down_pos[i]], "forward", "reverse") 
}
pos$HPVdir[pos$id == "HTMCP-03-06-02185/event1"] <- "reverse"
pos$HPVdir[pos$id == "HTMCP-03-06-02179/event3"] <- "reverse"

# Add 5' and 3' positions and HPV direction to the methyl regions data frame
mr$up_hpv_pos <- pos$up_pos[match(mr$id, pos$id)]
mr$down_hpv_pos <- pos$down_pos[match(mr$id, pos$id)]
mr$HPVdir <- pos$HPVdir[match(mr$id, pos$id)]

# Add HPV distance to the methyl regions data frame
mr$hpv.dist <- ifelse(grepl("5", mr$V12), mr$V2 - mr$up_hpv_pos, mr$V2 - mr$down_hpv_pos)

#################

# switch the direction to match the translation direction
mrTln <- mr
mrTln$hpv.dist <- case_when(
  mr$HPVdir == "forward" & mr$V12 ==  "5_upstream" ~ -(abs(mr$hpv.dist)),
  mr$HPVdir == "reverse" & mr$V12 ==  "5_upstream" ~ abs(mr$hpv.dist),
  mr$HPVdir == "forward" & mr$V12 ==  "3_downstream" ~ abs(mr$hpv.dist),
  mr$HPVdir == "reverse" & mr$V12 ==  "3_downstream" ~ -(abs(mr$hpv.dist)),
  mr$HPVdir == "forward" & mr$V12 ==  "5_downstream" ~ -(abs(mr$hpv.dist)),
  mr$HPVdir == "reverse" & mr$V12 ==  "5_downstream" ~ abs(mr$hpv.dist),
  mr$HPVdir == "forward" & mr$V12 ==  "3_upstream" ~ abs(mr$hpv.dist),
  mr$HPVdir == "reverse" & mr$V12 ==  "3_upstream" ~ -(abs(mr$hpv.dist)),
  TRUE ~ mr$hpv.dist
)

## bin the methlation into 500 bp bins
mrTln$bins <- smart_cut(mrTln$hpv.dist,
                        seq(from = -5000, to = 5000, by = 500), 
                        labels = ~paste0(.y[1],':',.y[2]-1), 
                        simplify = FALSE)

allbins <- smart_cut(seq(from = -5000, to = 5000, by = 500),
                     seq(from = -5000, to = 5000, by = 500), 
                     labels = ~paste0(.y[1],':',.y[2]-1), 
                     simplify = FALSE)

# make an average value for each bin
mrAve <- mrTln %>%
  dplyr::group_by(id,bins,V12, integration.type, hpv.type) %>%
  dplyr::summarise(average_methyl = mean(V8))

mrAve$bins <- factor(mrAve$bins, levels = levels(allbins))

### ----------------------------------------------------------
### PLOT THE METHYLATION DATA ADJACENT TO HPV INTEGRATION
### ----------------------------------------------------------

mrPlots <- lapply(c("ecDNA_integration","deletion_integration"), function(t){
  # Filter mrAve data for ecDNA and deletion integration types
  mrAveSub <- mrAve %>% filter(integration.type == t)
  
  # Order and factor the id variable for ecDNA data
  order <- mrAveSub %>% 
    dplyr::filter(bins %in% allbins[11:13]) %>% 
    dplyr::group_by(id, hpv.type) %>%
    dplyr::summarise(ave = mean(average_methyl)) %>%
    dplyr::arrange(hpv.type, ave) %>% 
    dplyr::distinct(id)
  mrAveSub$id <- factor(mrAveSub$id, levels = order$id)
  
  # Remove NA values for id variable
  mrAveSub <- mrAveSub[!is.na(mrAveSub$id) & !is.na(mrAveSub$bins),]
  
  # Create a plot for data
  plot <- ggplot(mrAveSub, aes(x = bins, y = id, fill = average_methyl)) +
    geom_tile() +
    scale_fill_distiller(palette = "RdBu") +
    theme_bw() +
    facet_grid(hpv.type ~ ., scales = "free", space = "free") +
    labs(fill = "average Methylation", y = "events", x = "distance from HPV integration") +
    theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5,size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          panel.grid.major.y = element_blank())
  return(list(order,plot))
}
)

### ----------------------------------------------------------
### PLOT METHYLATION ACROSS THE HPV GENOME
### ----------------------------------------------------------

# get the order for the ecDNA and deletion plots
ecOrder <- as.data.frame(mrPlots[[1]][1])[,1]
delOrder <- as.data.frame(mrPlots[[2]][1])[,1]

# Add the event categories and other info to the methyl regions
eventMethHPV$id <- paste0(eventMethHPV$library, "/", eventMethHPV$event)
eventMethHPV$integration.type <- summary$type[match(eventMethHPV$id, summary$id)]

# match the unintegrated samples
hpvMeth <- hpvMeth[,colnames(eventMethHPV)]

# add to the table 
eventMethHPV <- rbind(eventMethHPV,hpvMeth)

# define the integration types
types <- c("ecDNA_integration","deletion_integration", "duplication_integration","translocation_integration","no_integration_detected")

# plot the HPV methylation for each integration type
hpvPlot <- lapply(types, function(type) {
  if(type == "ecDNA_integration"){
    eventMethHPV$id = factor(eventMethHPV$id,levels = ecOrder)
  } else if (type == "deletion_integration") {
    eventMethHPV$id = factor(eventMethHPV$id,levels = delOrder)
  }
  
  ggplot(eventMethHPV %>% dplyr::arrange(chromosome) %>% filter(integration.type == type & !is.na(id)), 
         aes(x = E6.start, y = id, colour = perc.methylated)) +
    geom_point(size = 2) +
    facet_grid(chromosome ~ integration.type, scales = "free", space = "free") +
    theme_bw() + 
    geom_vline(xintercept = 0, linetype = 2, size = 1) +
    scale_color_distiller(palette = "RdBu") + 
    labs(x = "adjusted position", y = "sample", colour = "% methylated") +
    theme(axis.text.x = element_text(colour = "black", size = 11),
          #axis.text.y = element_blank(), 
          axis.title = element_text(colour = "black", size = 12, face = "bold"),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 12, colour = "black"),
          panel.grid.minor = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.line = element_line())
})

### ----------------------------------------------------------
### SAVE THE PLOTS 
### ----------------------------------------------------------

# eccDNA integration
ggsave(mrPlots[[1]][[2]], filename = "figurePanels/figure5/ecDNA_methylRegions.pdf", width = 6, height = 6, units = "in")
ggsave(hpvPlot[[1]], filename = "figurePanels/figure5/ecDNA_HPVmethyl.pdf", width = 6, height = 6, units = "in")

# deletion integration
ggsave(mrPlots[[2]][[2]], filename = "figurePanels/figure5/deletion_methylRegions.pdf", width = 6, height = 6, units = "in")
ggsave(hpvPlot[[2]], filename = "figurePanels/figure5/deletion_HPVmethyl.pdf", width = 6, height = 6, units = "in")

# duplication integration
ggsave(hpvPlot[[3]], filename = "figurePanels/figure5/duplication_HPVmethyl.pdf", width = 6, height = 6, units = "in")

# translocation integration
ggsave(hpvPlot[[4]], filename = "figurePanels/figure5/translocation_HPVmethyl.pdf", width = 6, height = 6, units = "in")

# unintegrated 
ggsave(hpvPlot[[5]], filename = "figurePanels/figure5/unintegrated_HPVmethyl.pdf", width = 6, height = 3, units = "in")

### ----------------------------------------------------------
### SAVE THE TABLES 
### ----------------------------------------------------------

# methyl regions
mrAveSave1 <- mrAve[mrAve$integration.type %in% c("ecDNA_integration"),]
mrAveSave2 <- mrAve[mrAve$integration.type %in% c("deletion_integration"),]

ecOrderDf <- data.frame(name=ecOrder, order=1:length(ecOrder))
delOrderDf <- data.frame(name=delOrder, order=1:length(delOrder))

mrAveSave1$id.number <- ecOrderDf$order[match(mrAveSave1$id, ecOrderDf$name)]
mrAveSave2$id.number <- delOrderDf$order[match(mrAveSave2$id, delOrderDf$name)]

mrAveSave <- rbind(mrAveSave1,mrAveSave2)
mrAveSave <- mrAveSave[complete.cases(mrAveSave),]

write.table(mrAveSave, file = "tables/methylationAdjacentRegions.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")

# methyl HPV
eventMethHPVSave1 <- eventMethHPV[eventMethHPV$integration.type %in% c("ecDNA_integration"),]
eventMethHPVSave2 <- eventMethHPV[eventMethHPV$integration.type %in% c("deletion_integration"),]
eventMethHPVSave3 <- eventMethHPV %>% dplyr::arrange(chromosome) %>% filter(integration.type == "duplication_integration" & !is.na(id))
eventMethHPVSave4 <- eventMethHPV %>% dplyr::arrange(chromosome) %>% filter(integration.type == "translocation_integration" & !is.na(id))
eventMethHPVSave5 <- eventMethHPV %>% dplyr::arrange(chromosome) %>% filter(integration.type == "no_integration_detected" & !is.na(id))

eventMethHPVSave <- rbind(eventMethHPVSave1, eventMethHPVSave2, eventMethHPVSave3, eventMethHPVSave4,eventMethHPVSave5)

write.table(eventMethHPVSave, 
            file = "tables/methylationHPVRegions.txt",
            quote = F, col.names = T, row.names = F, sep = "\t")
