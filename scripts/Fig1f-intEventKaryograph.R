.libPaths("/path/to/libraries")

## --------------------------------------------
# Load libraries
## --------------------------------------------

library(dplyr)
library(ggplot2)
library(ggsci)
library(cowplot)
library(ggpubr)

## --------------------------------------------
# Load files
## --------------------------------------------

chromSize <- read.delim("/path/to/refs/hg38_no_alt_TCGA_HTMCP_HPVs_chromSizes.txt", header = F)
centPos <-  read.delim("/path/to/refs/hg38_centromere_positions_merged.bed", header = F)

## ---------------------------------------------------------------------------
## POSITION CHROMOSOME INFO
## ---------------------------------------------------------------------------

# subset to the main chromosomes
chromSize <- chromSize[chromSize$V1 %in% c(paste0("chr", 1:22), "chrX"),]

# Rename columns to chromosome and size
colnames(chromSize) <- c("chr","size")

# Reorder levels for plotting
chromSize$chr <- factor(chromSize$chr,levels=c(paste0("chr", 1:22), "chrX"))
chromSize <- chromSize[chromSize$chr %in% c(paste0("chr", 1:22), "chrX"),]

# Divide by 1Mb to clean up axis
chromSize$size <- chromSize$size/1000000

# centromere mapping
colnames(centPos) <- c("chr", "start", "end")
centPos$chr <- factor(centPos$chr,levels=c(paste0("chr", 1:22), "chrX"))
centPos$centre <- centPos$start + ((centPos$end - centPos$start)/2)
centPos$centre <- centPos$centre/1000000

## --------------------------------------------
# Load files
## --------------------------------------------

sum <- read.delim("/path/to/tables/integrationTypes.txt")
sites <- read.delim("/path/to/tables/intSitesSummaryAll.txt")

htmcp_path <- "/path/to/htmcp/call_integration/output"
tcga_path <- "/path/to/tcga/call_integration/output"

### HPV METHYLATION

# Find files with HPV methylation by read
Files1 <- grep("hpv_integration_events_distance.bed", list.files(htmcp_path, recursive = TRUE, full.names = TRUE), value = TRUE)
Files2 <- grep("hpv_integration_events_distance.bed", list.files(tcga_path, recursive = TRUE, full.names = TRUE), value = TRUE)
Files <- c(Files1, Files2)

# Extract library names from file paths
sample <- gsub(paste0(htmcp_path, "/|", tcga_path, "/|/events/hpv_integration_events_distance.bed"), "", Files)
sample <- gsub("MaM_486","TCGA-C5-A2LY",sample)
sample <- gsub("MaM_485","TCGA-C5-A2LX",sample)

# Import data
dist <- lapply(Files, read.delim, header = FALSE, sep = "\t")
names(dist) <- sample
dist <- dplyr::bind_rows(dist, .id = "sample")
dist$hpv_site <- unlist(lapply(dist$V4,function(x){strsplit(x, ",")[[1]][1]}))
dist$event <- sites$event[match(paste0(dist$sample, dist$hpv_site), paste0(sites$library, sites$HPV.site))]
dist$integration.type <- sum$type[match(paste0(dist$sample, dist$event), paste0(sum$sample, sum$event))]
dist$nsites <- sum$nsites[match(paste0(dist$sample, dist$event), paste0(sum$sample, sum$event))]
dist$HPV.type <- sum$HPV.type[match(paste0(dist$sample, dist$event), paste0(sum$sample, sum$event))]

# save bed files
multiBed <- dist[dist$nsites > 2,c(2,3,4,1,7,5,10)]
multiBed <- multiBed[complete.cases(multiBed),]
ecBed <- dist[dist$integration.type == "ecDNA_integration" | dist$integration.type == "complex_ecDNA_integration" & dist$nsites == 2,c(2,3,4,1,7,5)]
ecBed <- ecBed[complete.cases(ecBed),]
delBed <- dist[dist$integration.type == "deletion_integration",c(2,3,4,1,7,5,10)]
delBed <- delBed[complete.cases(delBed),]
dupBed <- dist[dist$integration.type == "duplication_integration",c(2,3,4,1,7,5,10)]
dupBed <- dupBed[complete.cases(dupBed),]
traBed <- dist[dist$integration.type == "translocation_integration",c(2,3,4,1,7,5,10)]
traBed <- traBed[complete.cases(traBed),]
repBed <- dist[dist$integration.type == "repeat_integration",c(2,3,4,1,7,5,10)]
repBed <- repBed[complete.cases(repBed),]

bedL <- list(multiBed, ecBed, delBed, dupBed, traBed, repBed)

write.table(dist[,c(2,3,4,1,7,5,10)],"/path/to/tables/event_locations.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(multiBed,"/path/to/tables/multi-breakpoint_event_locations.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(ecBed,"/path/to/tables/ecDNA_event_locations.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(delBed,"/path/to/tables/deletion_event_locations.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(dupBed,"/path/to/tables/duplication_event_locations.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(traBed,"/path/to/tables/translocation_event_locations.bed", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(repBed,"/path/to/tables/repeat_event_locations.bed", quote = F, sep = "\t", col.names = F, row.names = F)

## --------------------------------------------
# edit dataframes
## --------------------------------------------

bedL2 <- lapply(bedL, function(x){
    colnames(x)[1:3] <- c("chr", "start", "end")
    x$start <- x$start/1000000
    x$end <- x$end/1000000
    x$midpoint <- round(x$start + ((x$end - x$start)/2))
    x <- x %>% arrange(chr, midpoint)
    x$position <- paste0(x$chr, ":", x$midpoint)
    return(x)
})
names(bedL2) <- c("multibreakpoint", "ecDNA", "deletion", "duplication", "translocation", "repeat")

bedDF <- dplyr::bind_rows(bedL2, .id = "integration.type")

hs <- bedDF %>%
    group_by(integration.type, position) %>%
    summarise(n = n())

# split the chromosome name and bin position
hsPlot <- separate(hs, col = position, into = c("chr", "pos"), sep = ":", remove = T)

# scale to fit the plot - i.e. make the maximum width 0.65
max <- max(hsPlot$n)
hsPlot$percMax <- hsPlot$n/max
hsPlot$percMax <- hsPlot$percMax * 0.65

# Change to factor and reorder levels
hsPlot$chr <- factor(hsPlot$chr, levels=c(paste0("chr", 1:22), "chrX"))
hsPlot$pos <- as.numeric(hsPlot$pos)

## --------------------------------------------
# make figure
## --------------------------------------------

mb <- sum[sum$type == "multi-breakpoint_integration",]

ph <- ggplot(mb, aes(x = nsites)) +
    geom_histogram(colour = "black") +
    theme_minimal() +
    labs(x = "# of breakpoints", y = "count") +
    theme(axis.line.x.bottom = element_line(),
          axis.line.y.left = element_line(),
          axis.ticks = element_line(),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, colour = "black"))

ggsave(ph, filename = "multibreakEventsNsitesHist.pdf", width = 4, height = 4, units = "in")

## --------------------------------------------
# make figure
## --------------------------------------------

adL <- data.frame(xmin = c(9.7, 9.7, 10.3), xmax = c(10.35, 9.75,10.35), ymin = c(240,238,238), ymax = c(243,243,243),
                  fill = c("ase","ase","ase"))

adW <- data.frame(x = c(11.5,10), y = c(240,250),
                  label = c("# of integration events", as.character(c(max))))

p1 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSize %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 4) +
    # centromeres
    geom_point(data = centPos %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, y = centre), 
               size = 4, colour = "black") +
    # ASE genes
    geom_rect(data = hsPlot %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+3, fill = integration.type),
              size = 1) +
    # legend bars
    geom_rect(data = adL, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "black", size = 0.15) +
    # legend text
    geom_text(data = adW, 
              aes(x = x, y = y, label = label))+
    ylim(0, 250) +
    scale_fill_lancet()+
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank())+
    labs(x=NULL,y="Chromosome Size (Mb)")

# very annoying but you have to filter all the dataframes or else the factor levels won't match the integer value
chromSizeFilt <- chromSize %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
chromSizeFilt$chr <- factor(chromSizeFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
centPosFilt <- centPos %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
hsPlotFilt <- hsPlot %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
hsPlotFilt$chr <- factor(hsPlotFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))

p2 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSizeFilt, aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 4) +
    # centromeres
    geom_point(data = centPosFilt, aes(x = chr, y = centre), 
               size = 4, colour = "black") +
    # ASE genes
    geom_rect(data = hsPlotFilt, 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+3, fill = integration.type),
              size = 1) +
    ylim(0, 250) +
    scale_fill_lancet()+
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank())+
    labs(x=NULL,y="Chromosome Size (Mb)")

# put them together
plot <- plot_grid(p1, p2, align = "v", axis = "l", nrow = 2)
plot
# save plot
ggsave(plot, filename = "figurePanels/figure1/position_events_karyograph.pdf", width = 10, height = 7, units = "in")

