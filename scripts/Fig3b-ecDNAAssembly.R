#### ecDNA ASSEMBLIES 
#### VANESSA PORTER

## --------------------------------------------------------------------------------
## Load libraries
## --------------------------------------------------------------------------------
.libPaths("/path/to/libraries")
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(wesanderson)
library(ggsci)
library(scales)
library(DataEditR)
library(cowplot)
library(biomaRt)

## --------------------------------------------------------------------------------
## Load in HPV type gff files for ecDNA (HPV16, HPV45, HPV52)
## --------------------------------------------------------------------------------

# HPV16
hpv16_all <- read.delim("//path/to/ref/hpv_gff/HPV16.gff", header = F, comment.char = "#", stringsAsFactors = F)
hpv16 <- hpv16_all[hpv16_all$V3 == "gene",]
hpv16_size <- hpv16_all[1,]
hpv16_info <- strsplit(hpv16$V9, ";")
hpv16_info <- lapply(hpv16_info, function(x){
  x <- x[grep("Name=", x)]
  x <- gsub("Name=", "", x)
  return(x)})
hpv16$gene <- unlist(hpv16_info)
hpv16 <- hpv16[,-c(2,3,6,8,9)]
colnames(hpv16) <- c("genome", "start", "end", "strand","gene")

# HPV45
hpv45_all <- read.delim("/path/to/ref/hpv_gff/HPV45.gff", header = F, comment.char = "#", stringsAsFactors = F)
hpv45 <- hpv45_all[hpv45_all$V3 == "gene",]
hpv45_size <- hpv45_all[hpv45_all$V3 == "region",]
hpv45_info <- strsplit(hpv45$V9, ";")
hpv45_info <- lapply(hpv45_info, function(x){
  x <- x[grep("Name=", x)]
  x <- gsub("Name=", "", x)
  return(x)})
hpv45$gene <- unlist(hpv45_info)
hpv45 <- hpv45[,-c(2,3,6,8,9)]
colnames(hpv45) <- c("genome", "start", "end", "strand","gene")

# HPV52
hpv52_all <- read.delim("/path/to/ref/hpv_gff/HPV52.gff", header = F, comment.char = "#", stringsAsFactors = F)
hpv52 <- hpv45_all[hpv52_all$V3 == "gene",]
hpv52_size <- hpv52_all[hpv52_all$V3 == "region",]
hpv52_info <- strsplit(hpv52$V9, ";")
hpv52_info <- lapply(hpv52_info, function(x){
  x <- x[grep("Name=", x)]
  x <- gsub("Name=", "", x)
  return(x)})
hpv52$gene <- unlist(hpv52_info)
hpv52 <- hpv52[,-c(2,3,6,8,9)]
colnames(hpv52) <- c("genome", "start", "end", "strand","gene")

## --------------------------------------------------------------------------------
## Read in the aligned paf file of the assembly
## --------------------------------------------------------------------------------

# Fix the contig paf file
paf16 <- read.delim("/path/to/tcga/integration_events/F24953/event1.txt-flye-asm/assembly.hybrid.paf", header = F, stringsAsFactors = F)
paf16 <- paf16[,c(1,3:6,8,9)]
colnames(paf16) <- c("query.name", "query.start", "query.end", "strand","target.name", "target.start", "target.end")
paf16 <- paf16 %>%
  mutate(temp = target.start) %>%
  mutate(target.start=ifelse(strand=='-', target.end, target.start)) %>%
  mutate(target.end=ifelse(strand=='-', temp, target.end)) %>%
  dplyr::select(1:7)

paf45 <- read.delim("/path/to/htmcp/integration_events/F46074/event1.txt-flye-asm/assembly.hybrid.paf", header = F, stringsAsFactors = F)
paf45 <- paf45[,c(1,3:6,8,9)]
colnames(paf45) <- c("query.name", "query.start", "query.end", "strand","target.name", "target.start", "target.end")
paf45 <- paf45 %>%
  mutate(temp = target.start) %>%
  mutate(target.start=ifelse(strand=='-', target.end, target.start)) %>%
  mutate(target.end=ifelse(strand=='-', temp, target.end)) %>%
  dplyr::select(1:7)

paf52 <- read.delim("/path/to/htmcp/integration_events/F47779/event1.txt-flye-asm/assembly.hybrid.paf", header = F, stringsAsFactors = F)
paf52 <- paf52[,c(1,3:6,8,9)]
colnames(paf52) <- c("query.name", "query.start", "query.end", "strand","target.name", "target.start", "target.end")
paf52 <- paf52 %>%
  mutate(temp = target.start) %>%
  mutate(target.start=ifelse(strand=='-', target.end, target.start)) %>%
  mutate(target.end=ifelse(strand=='-', temp, target.end)) %>%
  dplyr::select(1:7)

## --------------------------------------------------------------------------------
## HPV16 ecDNA
## --------------------------------------------------------------------------------

# Values for rearrangement
start <- min(paf16$query.start[paf16$target.name == "HPV16"])
end <- max(paf16$query.end[paf16$target.name == "HPV16"])
hpv.start <- paf16$target.start[paf16$query.start == start]
hpv.end <- paf16$target.end[paf16$query.end == end]
size <- abs(end - start)
hpv.size <- hpv16_size$V5

# genes that are entirely removed
rm <- hpv16$gene[hpv16$start > hpv.end & hpv16$end < hpv.start] # none

# genes that are fine as is
stay <- hpv16$gene[hpv16$start < hpv.end & hpv16$end < hpv.end | hpv16$start > hpv.start & hpv16$end > hpv.start ] # E6, E7, L1, E5, L2, E1

# genes that are partially removed 
part <- hpv16$gene[!hpv16$gene %in% c(rm, stay)] # E4 and L2

# remove the genes that are cut out
hpv16 <- hpv16[!hpv16$gene %in% rm,]

# change the start/end on the genes kept as is
hpv16$new.start <- ifelse(hpv16$start < hpv.start & hpv16$start > hpv.end, 1, 
                          ifelse(hpv16$start < hpv.end, (hpv.size - hpv.start + hpv16$start), 
                                 hpv16$start - hpv.start))

hpv16$new.end <- ifelse(hpv16$end > hpv.end & hpv16$end < hpv.start, size, 
                        ifelse(hpv16$end < hpv.end, (hpv.size - hpv.start + hpv16$end), 
                               hpv16$end - hpv.start))

# there's 2 genes that are split in half, so we need to add the other half and fix the start and end. 
E2 <- hpv16[hpv16$gene == "E2",]
E4 <- hpv16[hpv16$gene == "E4",]
E2$new.end <- size
E4$new.end <- size
hpv16$new.start[hpv16$gene == "E2"] <- 1
hpv16$new.start[hpv16$gene == "E4"] <- 1
hpv16 <- rbind(hpv16, E2, E4)

# arrange the genes on the contig
hpv16$contig.start <- hpv16$new.start + start
hpv16$contig.end <- hpv16$new.end + start

hpv16$contig <- ifelse(hpv16$gene %in% c("E7","E1","E5","L2"), "anno3",
                       ifelse(hpv16$gene %in% c("E6","E2","L1"), "anno2", "anno1"))

hpv16 <- hpv16[,c("contig", "contig.start", "contig.end", "strand", "gene", "new.start", "new.end")]
colnames(hpv16) <- colnames(paf16)
circ16 <- rbind(paf16, hpv16)
circ16 <- arrange(circ16, query.start)

circ16 <- rbind(circ16, E2, E4) %>% arrange(query.start)

# set the plotting values
circ16$ystart <- ifelse(circ16$query.name == "contig_1", 1, 
                             ifelse(circ16$query.name == "anno3", 2, 
                                    ifelse(circ16$query.name == "anno2", 3, 4)))
circ16$yend <- ifelse(circ16$query.name == "contig_1", 1, 
                             ifelse(circ16$query.name == "anno3", 2, 
                                    ifelse(circ16$query.name == "anno2", 3, 4)))

circ16$size <- ifelse(circ16$query.name == "contig_1", 4, 5)

# colours
ann_colors <- list(target.name = c(chr13="#2C2E43", chr3="#2C2E43",
                                   HPV16="#FB9300", HPV52 ="#FB9300", HPV45 ="#FB9300",
                                   E6="#925E9FFF",E7="#FDAF91FF", E1="#0099B4FF",
                                   E2="#ED0000FF",E4="#42B540FF", E5="#ADB6B6FF", 
                                   L1="#00468BFF", L2="#AD002AFF", 
                                   A="#2C2E43", C="#F0E5CF", E="#8599ad"))

ggplot()+
  geom_segment(data = circ16, aes(x=query.start, xend=query.end, y=ystart, yend=yend, colour=target.name, size = size)) +
  scale_colour_manual(values = ann_colors[["target.name"]]) +
  theme_void() + 
  coord_polar() + 
  ylim(-3,5) + 
  scale_size(range = c(3, 5)) +
  geom_vline(xintercept = start, linetype = "longdash", colour = "#FB9300", size = 0.5) +
  geom_vline(xintercept = end, linetype = "longdash", colour = "#FB9300", size = 0.5) +
  theme(legend.position = "none")

## --------------------------------------------------------------------------------
## HPV45 ecDNA
## --------------------------------------------------------------------------------

# Values for rearrangement
start <- min(paf45$query.start[paf45$target.name == "HPV45"])
end <- max(paf45$query.end[paf45$target.name == "HPV45"])
hpv.start <- paf45$target.end[paf45$query.end == end]
hpv.end <- paf45$target.start[paf45$query.start == start]
size <- abs(end - start)
hpv.size <- hpv45_size$V5

# genes that are entirely removed
rm <- hpv45$gene[hpv45$start > hpv.end & hpv45$end < hpv.start] # E2 

# genes that are fine as is
stay <- hpv45$gene[hpv45$start < hpv.end & hpv45$end < hpv.end | hpv45$start > hpv.start & hpv45$end > hpv.start ] # E6, E7, L1

# genes that are partially removed 
part <- hpv45$gene[!hpv45$gene %in% c(rm, stay)] # E1 and L2

# remove the genes that are cut out
hpv45 <- hpv45[!hpv45$gene %in% rm,]

# change the start/end on the genes kept as is
hpv45$new.start <- ifelse(hpv45$start < hpv.start & hpv45$start > hpv.end, 1, 
                          ifelse(hpv45$start < hpv.end, (hpv.size - hpv.start + hpv45$start), 
                                 hpv45$start - hpv.start))

hpv45$new.end <- ifelse(hpv45$end > hpv.end & hpv45$end < hpv.start, size, 
                        ifelse(hpv45$end < hpv.end, (hpv.size - hpv.start + hpv45$end), 
                               hpv45$end - hpv.start))

# Since this is on the negative strand, need to reverse the order
hpv45$new.start <- size - hpv45$new.start + 1
hpv45$new.end <- size - hpv45$new.end + 1

# arrange the genes on the contig
hpv45$contig.start <- hpv45$new.start + start
hpv45$contig.end <- hpv45$new.end + start

hpv45$contig <- ifelse(hpv45$gene %in% c("E7","E1","E5","L2"), "anno3",
                       ifelse(hpv45$gene %in% c("E6","E2","L1"), "anno2", "anno1"))

colnames(hpv45) <- colnames(paf45)

# import gene information to map the TP63 exons
grch38 <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
tp63 <- getBM(attributes = c("ensembl_transcript_id_version","ensembl_exon_id","exon_chrom_start", "exon_chrom_end", "5_utr_start", "5_utr_end"),
              mart = grch38, 
              filters = "ensembl_transcript_id_version", 
              values = c("ENST00000264731.8", "ENST00000354600.10", "ENST00000418709.6",
                         "ENST00000437221.5", "ENST00000320472.9", "ENST00000440651.6",
                         "ENST00000392460.7", "ENST00000449992.5", "ENST00000392463.6",
                         "ENST00000392461.7", "ENST00000456148.1", "ENST00000434928.5"))

# retrieve the exons that are within the ecDNA
tp63 <- tp63[tp63$exon_chrom_start > 189745513 & tp63$exon_chrom_end < 189809230,]
tp63 <- tp63 %>% arrange(exon_chrom_start)
tp63_exons <- tp63[3:nrow(tp63),]
tp63_exons <- tp63_exons[!duplicated(tp63_exons$exon_chrom_start),3:4]
utr <- tp63[1,5:6]
colnames(utr) <- colnames(tp63_exons)
tp63_exons <- rbind(tp63_exons, utr)

# plot the exons within the assembly coordinates
chrstart <- min(paf45$query.start[paf45$target.name == "chr3"])
chr.start <- paf45$target.start[paf45$query.start == chrstart]
tp63_exons$contig.start <- tp63_exons$exon_chrom_start - chr.start + chrstart
tp63_exons$contig.end <- tp63_exons$exon_chrom_end - chr.start + chrstart
tp63_exons$query.name <- "anno2"
tp63_exons$target.name <- "exons"
tp63_exons$strand <- "+"
tp63_exons <- tp63_exons[,c(5,3,4,7,6,1,2)]
colnames(tp63_exons) <- colnames(paf45)

# put the df together
circ45 <- rbind(paf45, hpv45[,1:7], tp63_exons)
circ45 <- arrange(circ45, query.start)

# set the plotting values
circ45$ystart <- ifelse(circ45$query.name == "contig_1", 1, 
                        ifelse(circ45$query.name == "anno3", 2, 
                               ifelse(circ45$query.name == "anno2", 3, 4)))
circ45$yend <- ifelse(circ45$query.name == "contig_1", 1, 
                      ifelse(circ45$query.name == "anno3", 2, 
                             ifelse(circ45$query.name == "anno2", 3, 4)))

circ45$size <- ifelse(circ45$query.name == "contig_1", 4, 5)

# colours
ann_colors <- list(target.name = c(chr13="grey", chr3="grey", exons="#2C2E43",
                                   HPV16="#FB9300", HPV52 ="#FB9300", HPV45 ="#FB9300",
                                   E6="#925E9FFF",E7="#FDAF91FF", E1="#0099B4FF",
                                   E2="#ED0000FF",E4="#42B540FF", E5="#ADB6B6FF", 
                                   L1="#00468BFF", L2="#AD002AFF", 
                                   A="#2C2E43", C="#F0E5CF", E="#8599ad"))

ggplot()+
  geom_segment(data = circ45, aes(x=query.start, xend=query.end, y=ystart, yend=yend, colour=target.name, size = size)) +
  geom_segment(data = data.frame(start=c(1,66499),end=c(61087,69079)), 
               aes(x=start, xend=end, y=2, yend=2), colour="#2C2E43", size = 0.5, linetype = "F1") +
  geom_segment(data = data.frame(start=41568,end=61087), 
               aes(x=start, xend=end, y=3, yend=3), colour="#2C2E43", size = 0.5, linetype = "F1") +
  geom_segment(data = data.frame(start=60181,end=60435), 
               aes(x=start, xend=end, y=2, yend=2), colour="#2C2E43", size = 5) +
  scale_colour_manual(values = ann_colors[["target.name"]]) +
  theme_void() + 
  coord_polar() + 
  ylim(-3,5) + 
  scale_size(range = c(3, 5)) +
  geom_vline(xintercept = start, linetype = "longdash", colour = "#FB9300", size = 0.5) +
  geom_vline(xintercept = end, linetype = "longdash", colour = "#FB9300", size = 0.5) +
  theme(legend.position = "none")

## --------------------------------------------------------------------------------
## HPV52 ecDNA
## --------------------------------------------------------------------------------

# Values for rearrangement
start <- min(paf52$query.start[paf52$target.name == "HPV52"])
end <- max(paf52$query.end[paf52$target.name == "HPV52"])
hpv.start <- paf52$target.start[paf52$query.start == start]
hpv.end <- paf52$target.end[paf52$query.end == end]
size <- abs(end - start)
hpv.size <- hpv52_size$V5

# genes that are entirely removed
rm <- hpv52$gene[hpv52$start > hpv.end & hpv52$end < hpv.start] # E2 and L2

# genes that are fine as is
stay <- hpv52$gene[hpv52$start < hpv.end & hpv52$end < hpv.end | hpv52$start > hpv.start & hpv52$end > hpv.start ] # E6 and E7

# genes that are partially removed 
part <- hpv52$gene[!hpv52$gene %in% c(rm, stay)] # E1 and E2

# remove the genes that are cut out
hpv52 <- hpv52[!hpv52$gene %in% rm,]

# change the start/end on the genes kept as is
hpv52$new.start <- ifelse(hpv52$start < hpv.start & hpv52$start > hpv.end, 1, 
                          ifelse(hpv52$start < hpv.end, (hpv.size - hpv.start + hpv52$start), 
                                 hpv52$start - hpv.start))

hpv52$new.end <- ifelse(hpv52$end > hpv.end & hpv52$end < hpv.start, size, 
                        ifelse(hpv52$end < hpv.end, (hpv.size - hpv.start + hpv52$end), 
                               hpv52$end - hpv.start))

hpv52$contig.start <- hpv52$new.start + start
hpv52$contig.end <- hpv52$new.end + start

hpv52$contig <- ifelse(hpv52$gene %in% c("E7","E1","E5","L2"), "anno3",
                       ifelse(hpv52$gene %in% c("E6","E2","L1"), "anno2", "anno1"))

hpv52 <- hpv52[,c("contig", "contig.start", "contig.end", "strand", "gene", "new.start", "new.end")]

colnames(hpv52) <- colnames(paf52)
circ52 <- rbind(paf52, hpv52)
circ52 <- arrange(circ52, query.start)
circ52$target.name[c(1,10)] <- "C"
circ52$target.name[2] <- "A"
circ52$target.name[3] <- "E"
circ52$ystart <- c(rep(1.5, 4), 4, 1.5, 4, 3, 3, 1.5)
circ52$yend <- c(rep(1.5, 4), 4, 1.5, 4, 3, 3, 1.5)
circ52$size <- c(rep(4, 4), 5, 4, 5, 5, 5, 4)

# colours
mypal = pal_lancet("lanonc")(8)
show_col(mypal)
ann_colors <- list(target.name = c(chr13="black", HPV16="#FB9300", E6="#925E9FFF",E7="#FDAF91FF", E1="#0099B4FF",
                            E2="#ED0000FF",E4="#42B540FF", E5="#ADB6B6FF", L1="#00468BFF", L2="#AD002AFF", 
                            A="#2C2E43", C="#F0E5CF", E="#8599ad", HPV52 = "#FB9300"))


# ecDNA configuration figure 
p.circ52 <- ggplot()+
  geom_segment(data = circ52, aes(x=query.start, xend=query.end, y=ystart, yend=yend, colour=target.name, size = size)) +
  scale_colour_manual(values = ann_colors[["target.name"]]) +
  theme_void() + 
  coord_polar() + 
  ylim(-3,5) + 
  scale_size(range = c(3, 5)) +
  geom_vline(xintercept = 73317, linetype = "longdash", colour = "#FB9300", size = 0.5) +
  geom_vline(xintercept = 76821, linetype = "longdash", colour = "#FB9300", size = 0.5) +
  geom_vline(xintercept = 9431, linetype = "longdash", colour = "grey", size = 0.5) +
  geom_vline(xintercept = 34667, linetype = "longdash", colour = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = 1, size = 0.5) +
  geom_vline(xintercept = 30000, linetype = 1, size = 0.5) +
  geom_vline(xintercept = 60000, linetype = 1, size = 0.5) +
  geom_text(data = data.frame(intercepts=c(0,30000,60000),labels=c("0","30kb", "60kb")), 
            aes(x = intercepts, label = labels), 
            y = 6, angle = 0, size = 5) +
  geom_text(data = data.frame(intercepts=c(4500,20000,52000),labels=c("C", "A", "E")), 
            aes(x = intercepts, label = labels), 
            y = 2.1, angle = 0, size = 5, fontface = "bold") + 
  theme(legend.position = "none")

# linear configuration
lin52 <- circ52[circ52$target.name %in% c("A", "C","E"),5:10] %>% arrange(target.start)
lin52$text.x <- (lin52$target.start + lin52$target.end) / 2

p.lin52 <- ggplot()+
  geom_hline(yintercept = 1.5, linetype = 1, colour = "grey", size = 1) +
  geom_segment(data = lin52, aes(x=target.start, xend=target.end, y=ystart, yend=yend, colour=target.name, size = size)) +
  theme_void() +
  scale_colour_manual(values = ann_colors[["target.name"]]) +
  geom_text(data = data.frame(intercepts=c(74445483,74479405,74505476,74511381, 74530690),labels=c("A", "B", "C", "D", "E")), 
            aes(x = intercepts, label = labels), 
            y = 1.505, angle = 0, size = 5, fontface = "bold") + 
  geom_text(data = data.frame(intercepts=c(74491452),labels=c("117 kb")), 
            aes(x = intercepts, label = labels), 
            y = 1.515, angle = 0, size = 5) + 
  theme(legend.position = "none")

plot_grid(p.lin52, p.circ52, ncol = 2)

