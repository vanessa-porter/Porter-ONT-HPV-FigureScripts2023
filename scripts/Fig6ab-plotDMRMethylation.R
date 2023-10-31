### ----------------------------------------------------------
### PLOTTING THE DMR DIRECTION OF METHYLATION
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
### READ IN THE FILES
### ----------------------------------------------------------
files1 <- c(grep("diff_meth.csv", list.files("/path/to/htmcp/call_integration/output", 
                                                  recursive = T, full.names = T),value = T),
            grep("diff_meth.csv", list.files("/path/to/tcga/call_integration/output", 
                                                  recursive = T, full.names = T),value = T))

files2 <- c(grep("hpv_phaseblock_hp.txt", list.files("/path/to/htmcp/call_integration/output", 
                                             recursive = T, full.names = T),value = T),
            grep("hpv_phaseblock_hp.txt", list.files("/path/to/tcga/call_integration/output", 
                                             recursive = T, full.names = T),value = T))

# get the sample names
names1 <- gsub("/path/to/|htmcp|tcga|/call_integration/output/|/methylation/diff_meth.csv", "", files1)
names2 <- gsub("/path/to/|htmcp|tcga|/call_integration/output/|/event_phase/|/hpv_phaseblock_hp.txt", "", files2)

events <- strsplit(names2, "event")

### ----------------------------------------------------------
### FUNCTION
### ----------------------------------------------------------

### MAKE THE FUNCTION
generatePlotDF <- function(sample, event){
  dmrPath <- grep(sample, files1, value = T)
  pbPath <- grep(event, files2, value = T)
  
  print(sample)
  print(event)
  
  # dmrs
  dmr <- read.delim(dmrPath, header = T)
  
  # HPV phaseblock
  hpv <- read.delim(pbPath, header = F)
  
  if(ncol(hpv) == 10){
    colnames(hpv) <- c("hpvChr", "hpvStart", "hpvEnd", "hpvSites", "id", "pbChr", "pbStart", "pbEnd", "pbSize", "HP")
    hpv$midpoint <- (hpv$hpvEnd+hpv$pbStart)/2
    hpv$hpvDistStart <- hpv$pbStart - hpv$midpoint
    hpv$hpvDistEnd <- hpv$pbEnd - hpv$midpoint
    
    # subset dmrs to the phase block with HPV
    dmrHPV <- dmr %>% filter(chr == hpv$pbChr[1], start > hpv$pbStart[1], end < hpv$pbEnd[1])
    
    # convert to distance from HPV
    dmrHPV$hpvDistStart <- dmrHPV$start-hpv$midpoint[1]
    dmrHPV$hpvDistEnd <- dmrHPV$end-hpv$midpoint[1]
    
    if (hpv$HP[1] == "HP1"){
      colnames(dmrHPV)[6:7] <- c("hpv", "non.HPV")
    } else{
      colnames(dmrHPV)[6:7] <- c("non.HPV", "hpv")
    }
    
    # plot the HPV haplotype vs the other haplotype
    dmrHPV$diff.Methy.Dir <- ifelse(dmrHPV$hpv > dmrHPV$non.HPV, abs(dmrHPV$diff.Methy), -(abs(dmrHPV$diff.Methy)))
    
    # subset to 200kb on either side
    dmrHPV <- dmrHPV[dmrHPV$hpvDistStart > -500000 & dmrHPV$hpvDistEnd < 500000,]
    
    # determine the missing regions of the phase blocks
    if(hpv$hpvDistStart[1] > -500000){
      up <- data.frame(start=-500000,end=hpv$hpvDistStart[1])
    } else{
      up <- NULL
    }
    
    if(hpv$hpvDistEnd[1] < 500000){
      down <- data.frame(start=hpv$hpvDistEnd[1],end=500000)
    } else{
      down <- NULL
    }
    
    ## average methylation in DMRs
    
    h1 <- mean(dmrHPV$non.HPV)
    h2 <- mean(dmrHPV$hpv)
    meth <- ifelse(h1 < h2, abs(h2-h1), -abs((h1-h2)))
    
    return(list(dmrHPV, hpv, up, down, meth)) 
  }
  
}


### RUN THE FUNCTION
allEvents <- lapply(events, function(x){
  generatePlotDF(sample = x[1], event = x[2])
})

### NAME THE EVENTS
lab <- lapply(events, function(x){
  event <- paste0("event", x[2])
  l <- paste0(x[1], ":", event)
  return(l)
})

lab <- as.data.frame(t(rbind(lab)))

names(allEvents) <- lab$lab

### ----------------------------------------------------------
### GET THE FILES FOR THE DMR DENSITY
### ----------------------------------------------------------

files3 <- c(grep("densityHPVBins.txt", list.files("/path/to/htmcp/call_integration/output", 
                                                  recursive = T, full.names = T),value = T),
            grep("densityHPVBins.txt", list.files("/path/to/cga/call_integration/output", 
                                                  recursive = T, full.names = T),value = T))

files4 <- c(grep("dmrHotspotIntersectHPV.bed", list.files("/path/to/htmcp/call_integration/output", 
                                                          recursive = T, full.names = T),value = T),
            grep("dmrHotspotIntersectHPV.bed", list.files("/path/to/tcga/call_integration/output", 
                                                          recursive = T, full.names = T),value = T))

files5 <- c(grep("densityDMRHotspotsRegions.sorted.bed", list.files("/path/to/htmcp/call_integration/output", 
                                                                    recursive = T, full.names = T),value = T),
            grep("densityDMRHotspotsRegions.sorted.bed", list.files("/path/to/tcga/call_integration/output", 
                                                                    recursive = T, full.names = T),value = T))

# get the sample names
names1 <- gsub("/path/to/|htmcp|tcga|/call_integration/output/|/methylation/densityHPVBins.txt", "", files3)
names2 <- gsub("/path/to/|htmcp|tcga|/call_integration/output/|/methylation/dmrHotspotIntersectHPV.bed", "", files4)

dsHPV <- list()
for (f in 1:length(files3)) {
  file <- files3[f]
  df <- read.delim(file, header = T)
  dsHPV[[f]] <- df
}
names(dsHPV) <- names1
dsHPV <- bind_rows(dsHPV, .id = "sample")
dsHPV$region <- paste0(dsHPV$sample, ":", dsHPV$sites)

intHPV <- list()
for (f in 1:length(files4)) {
  file <- files4[f]
  df <- read.delim(file, header = F)
  intHPV[[f]] <- df
}
names(intHPV) <- names2
intHPV <- bind_rows(intHPV, .id = "sample")
intHPV$region <- paste0(intHPV$sample, ":", intHPV$V4)
intHPV <- intHPV[intHPV$V5 != ".",]
intHPV$DMR.size <- intHPV$V7 - intHPV$V6
summary(intHPV$DMR.size)
intHPV$test <- "DMR_HPV"

allDMR <- list()
for (f in 1:length(files5)) {
  file <- files5[f]
  df <- read.delim(file, header = F)
  allDMR[[f]] <- df
}
names(allDMR) <- names2
allDMR <- bind_rows(allDMR, .id = "sample")
allDMR$DMR.size <- allDMR$V3 - allDMR$V2 
summary(allDMR$DMR.size)
allDMR$test <- "all"
dmrSizeDF <- rbind(allDMR, intHPV[,c(1:4,11,12)])

### ----------------------------------------------------------
### PLOT THE DMR HOTSPOT SIZES
### ----------------------------------------------------------

dsHPVsub <- dsHPV[dsHPV$bin_type == "HPV",]
dsHPVsub <- dsHPVsub %>% arrange(dsNorm)

dsHPV$region <- factor(dsHPV$region, levels = dsHPVsub$region)

myPalette <- colorRampPalette(brewer.pal(11, "Greys"))
heat <- ggplot(dsHPV %>% filter(region %in% intHPV$region), aes(x = as.factor(bin_num), y = region, fill = dsNorm)) + 
  geom_tile() +
  scale_fill_gradientn(colours = myPalette(100), limits = c(0, 1)) +
  geom_vline(xintercept = 301, colour = "red", linetype=2) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
heat

pSize <- ggplot(dmrSizeDF, aes(y = DMR.size, x = test)) +
  geom_violin() +
  #geom_jitter(height = 0, width = 0.1, size =2, alpha=0.5) +
  theme_bw() +
  labs(y = "DMR hotspot size", x = NULL) +
  theme(axis.text.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, colour = "black", face = "bold"))

pSize <- ggplot(dmrSizeDF, aes(x = DMR.size/1000000, colour = test)) +
  geom_density(size = 1) +
  theme_bw() +
  scale_colour_manual(values = c("grey", "#ef233c")) +
  labs(x = "DMR hotspot size (Mb)", colour = NULL) +
  xlim(0,30)+
  geom_vline(xintercept = median(dmrSizeDF$DMR.size[dmrSizeDF$test == "DMR_HPV"])/1000000, linetype = 2, colour = "#ef233c")+
  geom_vline(xintercept = median(dmrSizeDF$DMR.size[dmrSizeDF$test == "all"])/1000000, linetype = 2, colour = "grey")+
  theme(axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
        axis.title = element_text(size = 12, colour = "black", face = "bold")) 

pSize

ggsave(plot = pSize, filename = "figurePanels/figure6/DMRHPVHotspotSize.pdf", width = 4, height = 5, units = "in")

### ----------------------------------------------------------
### Read in the files
### ----------------------------------------------------------

types <- read.delim("/path/to/tables/integrationTypes.txt", header = T)
sum <- read.delim("/path/to/tables/intSitesSummaryAll.txt", header = T)

# get the total number of each event
tot <- types %>%
  group_by(type) %>%
  summarise(n = n())

# match the HPV DMRs with events
intHPV$event <- sum$event[match(paste0(intHPV$sample,gsub("\\,.*","", intHPV$V4)), paste0(sum$library,sum$HPV.site))]
intHPV$integration_type <- types$type[match(paste0(intHPV$sample,intHPV$event), paste0(types$sample,types$event))]
intHPV$total_type <- tot$n[match(intHPV$integration_type, tot$type)]
intHPV$integration_type <- gsub("_integration", "", intHPV$integration_type)
intHPV$HPV.type <- sum$HPVchr[match(intHPV$sample, sum$library)]

type_perc <- intHPV %>%
  group_by(integration_type) %>%
  summarise(perc=n()/unique(total_type))

type_perc$integration_type <- factor(type_perc$integration_type, levels = type_perc$integration_type[order(type_perc$perc, decreasing = T)])

# event types that effect methylation
intHPV$integration_type <- ifelse(is.na(intHPV$integration_type), "unmatched", intHPV$integration_type)
intHPV$integration_type <- factor(intHPV$integration_type, 
                                  levels = names(summary(as.factor(intHPV$integration_type))[rev(order(summary(as.factor(intHPV$integration_type))))]))


ptypes1 <- ggplot(intHPV, aes(x = integration_type)) +
  geom_bar() +
  theme_bw() + 
  #ylim(0,13) +
  labs(x = "integration type of HPV-overlapped DMR hotspots") +
  theme(axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text = element_text(size = 12, colour = "black"))

ptypes2 <- ggplot(type_perc %>% filter(integration_type != "unmatched"), aes(x = integration_type, y = perc)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  labs(x = "integration type of HPV-overlapped DMR hotspots", y = "percent of events") +
  theme(axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.text = element_text(size = 12, colour = "black"))

ggsave(plot = ptypes2, filename = "figurePanels/figure6/DMRHPVHotspotTypesPerc.pdf", width = 4, height = 4, units = "in")
ggsave(plot = ptypes1, filename = "figurePanels/figure6/DMRHPVHotspotTypes.pdf", width = 8, height = 6, units = "in")

### ----------------------------------------------------------
### SUBSET FOR THE SIGNIFICANT EVENTS
### ----------------------------------------------------------

intHPV$event.id <- paste0(intHPV$sample, ":", intHPV$event, "_", intHPV$V1, "_", intHPV$V2)
subEvents <- allEvents[-grep(pattern = "chrX",names(allEvents))]

sigEvents <- subEvents[names(subEvents) %in% intHPV$event.id]

# find order
methyl <- lapply(sigEvents, function(x){
  methyl <- x[[5]]
  methylDF <- data.frame(aveMethyl = methyl)
  return(methylDF)
})
m <- bind_rows(methyl, .id = "id")
m <- arrange(m, aveMethyl)

sigEvents <- sigEvents[m$id]

### ----------------------------------------------------------
### ADD OTHER INFO
### ----------------------------------------------------------

m$status <- ifelse(m$aveMethyl > 0.25, "methylated", ifelse(m$aveMethyl < -0.25, "unmethylated", "variable"))
m$int.type <- intHPV$integration_type[match(m$id, intHPV$event.id)]
m$event <- intHPV$event[match(m$id, intHPV$event.id)]
m$sample <- intHPV$sample[match(m$id, intHPV$event.id)]
m$HPV.type <- intHPV$HPV.type[match(m$id, intHPV$event.id)]
m$DMR.size <- intHPV$DMR.size[match(m$id, intHPV$event.id)]

### ----------------------------------------------------------
### Read in the files
### ----------------------------------------------------------

files6 <- c(grep("dmrIntersectHPV.bed", list.files("/path/to/htmcp/call_integration/output", 
                                                   recursive = T, full.names = T),value = T),
            grep("dmrIntersectHPV.bed", list.files("/path/to/tcga/call_integration/output", 
                                                   recursive = T, full.names = T),value = T))

names3 <- gsub("/path/to/|htmcp|tcga|/call_integration/output/|/methylation/dmrIntersectHPV.bed", "", files6)

dmrs <- list()
for (f in 1:length(files6)) {
  file <- files6[f]
  df <- read.delim(file, header = F)
  df <- df[df$V4 != ".",]
  dmrs[[f]] <- df
}
names(dmrs) <- names3
dmrs <- bind_rows(dmrs, .id = "sample")

dmrsFilt <- dmrs[abs(dmrs$V8) < 10000000,]
dmrsFilt$hpv_midpoint <- (dmrsFilt$V5+dmrsFilt$V6)/2
dmrsFilt$region <- paste0(dmrsFilt$sample, ":", dmrsFilt$V7)
dmrsFilt$start <- dmrsFilt$V2 - dmrsFilt$hpv_midpoint
dmrsFilt$end <- dmrsFilt$V3 - dmrsFilt$hpv_midpoint

dmrsFilt$event.id <- intHPV$event.id[match(dmrsFilt$region, intHPV$region)]
dmrsFilt <- dmrsFilt[dmrsFilt$event.id %in% m$id,]
dmrsFilt$event.id <- factor(dmrsFilt$event.id, levels = m$id)

pdmrs <- ggplot(dmrsFilt, aes(xmin = start/1000000, xmax = end/1000000, ymin = as.integer(event.id)-1, ymax = as.integer(event.id)+1)) +
  geom_rect(size=0.01, colour = "black") +
  #geom_rect(aes(xmin = start, xmax = end, ymin = as.integer(as.factor(region))-1, ymax = as.integer(as.factor(region))+1, fill = "#ef233c")) +
  facet_wrap(region ~ ., scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(x = "distance from HPV integration (Mb)")+
  geom_vline(xintercept = 0, linetype = 1, colour = "#ef233c", size = 0.5)+
  xlim(-10,10)+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank())
pdmrs

ggsave(plot = pdmrs, filename = "figurePanels/figure6/DMRsAroundHPV10Mb.pdf", width = 8, height = 10, units = "in")
ggsave(plot = pdmrs, filename = "figurePanels/figure6/DMRsAroundHPV10Mb.png", width = 8, height = 11, units = "in")

### ----------------------------------------------------------
### MAKE THE FIGURE SHOWING DMR METHYLATION AROUND HPV
### ----------------------------------------------------------

p <- NULL
for (i in 1:length(sigEvents)) {
  sub <- sigEvents[[i]]
  
  if(is.null(sub[[3]]) & is.null(sub[[4]])){
    plot <- ggplot() +
      geom_rect(data=sub[[1]], aes(xmin = hpvDistStart/1000, xmax = hpvDistEnd/1000, ymin = 0, ymax = 1, fill = diff.Methy.Dir, colour = diff.Methy.Dir), size = 0.5) +
      theme_bw() + 
      xlim(c(-500,500))+
      scale_fill_distiller(palette = "RdBu", limits = c(-1,1)) +
      scale_colour_distiller(palette = "RdBu", limits = c(-1,1)) +
      geom_vline(xintercept = 0, linetype = 2) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), "cm"))
  } else if(is.null(sub[[3]])){
    plot <- ggplot() +
      geom_rect(data=sub[[1]], aes(xmin = hpvDistStart/1000, xmax = hpvDistEnd/1000, ymin = 0, ymax = 1, fill = diff.Methy.Dir, colour = diff.Methy.Dir), size = 0.5) +
      theme_bw() + 
      geom_rect(data=sub[[4]], aes(xmin=(start/1000), xmax = (end/1000),ymin = 0, ymax = 1), fill = "grey")+
      xlim(c(-500,500))+
      scale_fill_distiller(palette = "RdBu", limits = c(-1,1)) +
      scale_colour_distiller(palette = "RdBu", limits = c(-1,1)) +
      geom_vline(xintercept = 0, linetype = 2) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), "cm"))
  } else if(is.null(sub[[4]])){
    plot <- ggplot() +
      geom_rect(data=sub[[1]], aes(xmin = hpvDistStart/1000, xmax = hpvDistEnd/1000, ymin = 0, ymax = 1, fill = diff.Methy.Dir, colour = diff.Methy.Dir), size = 0.5) +
      theme_bw() + 
      geom_rect(data=sub[[3]], aes(xmin=(start/1000), xmax = (end/1000),ymin = 0, ymax = 1), fill = "grey")+
      xlim(c(-500,500))+
      scale_fill_distiller(palette = "RdBu", limits = c(-1,1)) +
      scale_colour_distiller(palette = "RdBu", limits = c(-1,1)) +
      geom_vline(xintercept = 0, linetype = 2) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), "cm"))
  } else{
    plot <- ggplot() +
      geom_rect(data=sub[[1]], aes(xmin = hpvDistStart/1000, xmax = hpvDistEnd/1000, ymin = 0, ymax = 1, fill = diff.Methy.Dir, colour = diff.Methy.Dir), size = 0.5) +
      theme_bw() + 
      geom_rect(data=sub[[3]], aes(xmin=(start/1000), xmax = (end/1000),ymin = 0, ymax = 1), fill = "grey")+
      geom_rect(data=sub[[4]], aes(xmin=(start/1000), xmax = (end/1000),ymin = 0, ymax = 1), fill = "grey")+
      xlim(c(-500,500))+
      scale_fill_distiller(palette = "RdBu", limits = c(-1,1)) +
      scale_colour_distiller(palette = "RdBu", limits = c(-1,1)) +
      geom_vline(xintercept = 0, linetype = 2) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), "cm"))
  } 
  
  p[[i]] <- plot
}

pFinal <- plot_grid(plotlist = p, align = "v", ncol = 1)

ggsave(plot=pFinal, filename = "figurePanels/figure6/dmrMethylDirSigEvents.pdf", height = 10, width = 8, units = "in")

sub <- sigEvents[[2]]
plot2 <- ggplot() +
  geom_rect(data=sub[[1]], aes(xmin = hpvDistStart/1000, xmax = hpvDistEnd/1000, ymin = 0, ymax = 1, fill = diff.Methy.Dir, colour = diff.Methy.Dir), size = 0.5) +
  theme_bw() + 
  geom_rect(data=sub[[3]], aes(xmin=(start/1000), xmax = (end/1000),ymin = 0, ymax = 1), fill = "grey")+
  geom_rect(data=sub[[4]], aes(xmin=(start/1000), xmax = (end/1000),ymin = 0, ymax = 1), fill = "grey")+
  xlim(c(-500,500))+
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1)) +
  scale_colour_distiller(palette = "RdBu", limits = c(-1,1)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave(plot=plot2, filename = "figurePanels/figure6/dmrMethylDirSigEventsLegendAxis.pdf", height = 2, width = 8, units = "in")

### ----------------------------------------------------------
### COVARIANTS OF METHYLATION
### ----------------------------------------------------------

cov <- gather(m, covariant, value, c(3,4,6,7), factor_key=TRUE)

m$id <- factor(m$id, levels = m$id)

ggplot() +
  geom_tile(data=m, aes(x = id, y = 1, fill = aveMethyl)) +
  #geom_tile(data=cov, aes(x = id, y = covariant, fill = value)) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

### ----------------------------------------------------------
### DMR P VALUES
### ----------------------------------------------------------

files7 <- c(grep("/pvalue.txt", list.files("/path/to/htmcp/call_integration/output", 
                                                   recursive = T, full.names = T),value = T),
            grep("/pvalue.txt", list.files("/path/to/tcga/call_integration/output", 
                                                   recursive = T, full.names = T),value = T))

names <- gsub("/path/to/|htmcp|tcga|/call_integration/output/|/event_phase|/dmr_permute/pvalue.txt", "", files7)

dmrP <- list()
for (f in 1:length(files7)) {
  file <- files7[f]
  df <- read.delim(file, header = T)
  dmrP[[f]] <- df
}
names(dmrP) <- names
dmrP <- bind_rows(dmrP, .id = "id")
# adjusted p value
dmrP$padj <- p.adjust(dmrP$pval)
dmrP$id <- gsub("/", ":",dmrP$id)

# add pvalue to table
m$pval <- dmrP$pval[match(m$id,dmrP$id)]

# add number 
m$number <- as.character(1:nrow(m))

### ----------------------------------------------------------
### DMR DENSITY COMPARISONS
### ----------------------------------------------------------

files8 <- c(grep("/plottingDensityValues.txt", list.files("/path/to/htmcp/call_integration/output", 
                                           recursive = T, full.names = T),value = T),
            grep("/plottingDensityValues.txt", list.files("/path/to/tcga/call_integration/output", 
                                           recursive = T, full.names = T),value = T))
names <- gsub("/path/to/|htmcp|tcga|/call_integration/output/|/event_phase|/dmr_permute/plottingDensityValues.txt", "", files8)

# Import data
methDensity <- lapply(files8, read.delim, header = TRUE, sep = "\t")
names(methDensity) <- names
methDensity <- dplyr::bind_rows(methDensity, .id = "id")

# add info
methDensity$id <- gsub("/", ":", methDensity$id)
methDensity$pvalue <- m$pval[match(methDensity$id, m$id)]
methDensity$number <- m$number[match(methDensity$id, m$id)]

methDensity$number <- factor(methDensity$number, levels = 1:nrow(m))

# plot
plotD <- ggplot(methDensity, aes(x = region, y = dmr.density)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = sample), height = 0, width = 0.2, size =2, alpha=0.5) + 
  facet_wrap(vars(number), ncol = 3, scales = "free_y") +
  theme_minimal() +
  labs(x = NULL, y = "DMR Density") +
  #annotate("text", x=1.3, y=max(de$dmr.density[de$region == "hpv_region"]), 
  #         label= paste0("pvalue = ", p), 
  #         colour = "grey30", size = 5) + 
  scale_colour_manual(values = c("#e63946", "#1d3557")) +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.title = element_text(size=14), 
        axis.ticks.y = element_line(),
        axis.line = element_line(), 
        legend.position = "none")

ggsave(plot=plotD, filename = "figurePanels/figure6/dmrDensityBoxplots.pdf", height = 11, width = 5, units = "in")


### ----------------------------------------------------------
### SAVE TABLES
### ----------------------------------------------------------

# integration event characteristics
write.table(m, file = "tables/integrationEventsOnDMRHotspots.txt", quote = F, col.names = T, row.names = F, sep = "\t")
