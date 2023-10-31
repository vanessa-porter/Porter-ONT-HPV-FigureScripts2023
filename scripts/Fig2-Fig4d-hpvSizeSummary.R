.libPaths("/path/to/libraries")

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pafr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))

### -------------------------------------------------------------------------------
### READ IN FILES 
### -------------------------------------------------------------------------------

### CATEGORIES

f <- grep("hpvSizeCategories.txt", c(list.files("/path/to/htmcp/call_integration/output", recursive = T, full.names = T),
                                  list.files("/path/to/tcga/call_integration/output", recursive = T, full.names = T)), value = T)
n <- gsub("/path/to/|tcga|htmcp|/call_integration/output/|/hpv_size/hpvSizeCategories.txt", "", f)

list <- list()
for (i in 1:length(f)) {
    df <- read.delim(f[i], header = T)
    list[[i]] <- df
}
names(list) <- n
cat <- bind_rows(list, .id = "sample")

cat$sample <- gsub("MaM_485", "TCGA-C5-A2LX", cat$sample)
cat$sample <- gsub("MaM_486", "TCGA-C5-A2LY", cat$sample)

### -------------------------------------------------------------------------------
### ADD THE INTEGRATION CATEGORIES
### -------------------------------------------------------------------------------

int <- read.delim("/path/to/tables/integrationTypes.txt", header = T)
sum <- read.delim("/path/to/tables/intSitesSummaryAll.txt", header = T)

cat$bp1 <- gsub("\\_.*","",cat$bp_pair)
sum$bp <- paste0(sum$chr, ":", sum$pos)

# add the info to the df
cat$event <- sum$event[match(cat$bp1, sum$bp)]
cat$integration_type <- int$type[match(paste0(cat$sample, cat$event), paste0(int$sample, int$event))]
cat$integration_type[is.na(cat$integration_type)] <- "unmatched_integration"
cat$HPV.type <- sum$HPVchr[match(paste0(cat$sample, cat$event), paste0(sum$library, sum$event))]

cat <- cat[complete.cases(cat),]

write.table(cat, file = "/path/to/tables/integrantBreakpairs.txt", quote = F, col.names = T, row.names = F, sep = "\t")
catC <- cat[cat$status == "complete",]

### -------------------------------------------------------------------------------
### FIGURES
### -------------------------------------------------------------------------------

ann_colors <- list(HPV.clade = c(A7="#B55F8F", A9="#253083", Other="#A9A9A9"),
                   HPV.type = c(HPV16="#3953A4", HPV18="#9768ad", HPV45="#CC138C", HPV82="#369797",
                                HPV52="#0B8DCD", HPV31="#d2ecf9", HPV73="#737474", HPV68="#f3b2d4", HPV97="black", HPV58="#8bafba", HPV59="#ae3030"))

# hpv type
p1 <- ggplot(catC %>% filter(HPV.type %in% c("HPV16", "HPV18")), aes(x = HPV.type, colour = HPV.type, y = max_nHPV))+
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(height = 0, width = 0.2, size =3, alpha=0.5) +
    theme_minimal() +
    xlab("HPV type") + 
    ylab("max # of HPV genomes in integrant") +
    scale_colour_manual(values = c(ann_colors[["HPV.type"]])) +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none") +
    stat_compare_means(method = "wilcox")

ggsave(p1, file="/path/to/figurePanels/figure2/fig2i.pdf", height = 3, width = 3)

catC$category_simple <- ifelse(catC$category %in% c("partial", "full"), "single", "heterologous")


catC$integration_type <- factor(catC$integration_type, levels = c("multi-breakpoint_integration","repeat_integration",
                                                                  "deletion_integration","duplication_integration","ecDNA_integration",
                                                                  "translocation_integration"))

# integration type
p2 <- ggplot(catC %>% filter(category_simple == "heterologous"), aes(x = integration_type))+
    geom_bar() + 
    theme_minimal() +
    xlab("integration type") + 
    ylab("# of heterologous integrants") +
    #scale_colour_manual(values = c(ann_colors[["HPV.type"]])) +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black", angle = 60, hjust = 1, vjust = 1),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line())

ggsave(p2, file="/path/to/figurePanels/figure2/fig2f.pdf", height = 4.5, width = 3.2)

catC$all <- "all"
catC$HPV18_vs_all <- ifelse(catC$HPV.type == "HPV18", "HPV18", "Other_HPV")
# integration type
p3 <- ggplot(catC, aes(x = HPV18_vs_all, fill = category_simple))+
    geom_bar(position="fill", colour = "black") + 
    theme_minimal() +
    xlab("integration type") + 
    ylab("% of integrants") +
    labs(fill = NULL) +
    scale_fill_manual(values = c("#219ebc","#ffb703")) +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black", angle = 60, hjust = 1, vjust = 1),
          legend.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line())

ggsave(p3, file="/path/to/figurePanels/figure2/fig2h.pdf", height = 3.2, width = 4.5)

# hpv type
p4 <- ggplot(catC, aes(x = integration_type, y = max_nHPV))+
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(height = 0, width = 0.2, size =2, alpha=0.5, colour = "dark grey") +
    theme_minimal() +
    xlab("integration type") + 
    ylab("max # of HPV genomes in integrant") +
    #scale_colour_manual(values = c(ann_colors[["HPV.type"]])) +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black", angle = 60, hjust = 1, vjust = 1),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none") 

ggsave(p4, file="/path/to/figurePanels/figure2/fig2g.pdf", height = 5.5, width = 5.5)

# ngroups 

p5 <- ggplot(catC, aes(x = as.factor(ngroups), fill = category_simple))+
    geom_bar(colour = "black") + 
    theme_minimal() +
    xlab("# of integrant structures") + 
    ylab("# of breakpoint pairs") +
    labs(fill = NULL) +
    scale_fill_manual(values = c("#219ebc","#ffb703")) +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black"),
          legend.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line())

ggsave(p5, file="/path/to/figurePanels/figure2/fig2c.pdf", height = 2.5, width = 6)

# incomplete 

catI <- cat[cat$status == "incomplete",]
mean(catI$max_nHPV)
table(catI$HPV.type)

p6 <- ggplot(catI, aes(x = max_length/1000))+
    geom_histogram(colour = "black", size = 0.5) + 
    theme_minimal() +
    xlab("max length of incomplete HPV integrant (bp)") + 
    ylab("count") +
    labs(fill = NULL) +
    theme(panel.grid.minor = element_blank(), 
          axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(size = 12,colour = "black"),
          legend.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          axis.line = element_line())

ggsave(p6, file="/path/to/figurePanels/figure2/fig2j.pdf", height = 4, width = 4)

### -------------------------------------------------------------------------------
### GET STATS ON THE INTEGRANTS
### -------------------------------------------------------------------------------

# Longest integrant
head(cat %>% arrange(desc(max_nHPV)), 1)
head(cat %>% filter(status == "complete") %>% arrange(desc(max_nHPV)), 1)

# complete integrants
table(catC$integration_type)
table(catC$HPV.type)

# percent heterologous
table(catC$category)
table(catC$category_simple)
table(catC$HPV18_vs_all)
table(catC$category)/sum(table(catC$category))

fisher.test(matrix(c(nrow(catC[catC$HPV.type != "HPV18" & catC$category == "heterologous",]),
                     nrow(catC[catC$HPV.type == "HPV18" & catC$category == "heterologous",]),
                     nrow(catC[catC$HPV.type != "HPV18" & catC$category != "heterologous",]),
                     nrow(catC[catC$HPV.type == "HPV18" & catC$category != "heterologous",])),byrow = T, ncol=2))

# heterologous 
catCH <- catC %>% filter(category == "heterologous")
table(catCH$integration_type)
table(catCH$HPV.type)
table(catCH$integration_type)/sum(table(catCH$integration_type))

# heterologous per HPV type
table(catC$category[catC$HPV.type == "HPV16"])/sum(table(catC$category[catC$HPV.type == "HPV16"]))
table(catC$category[catC$HPV.type == "HPV18"])/sum(table(catC$category[catC$HPV.type == "HPV18"]))
table(catC$category[catC$HPV.type != "HPV18"])/sum(table(catC$category[catC$HPV.type != "HPV18"]))

# length per HPV type
catC %>%
  group_by(HPV.type) %>%
  summarise(mean = mean(max_nHPV))

catC %>%
  group_by(HPV.type) %>%
  summarise(mean = mean(max_nHPV))

# number of incomplete breakpoints
table(cat$category)
table(cat$HPV.type[cat$category == "incomplete"])

# largest integrants
max(cat$max_nHPV[cat$category == "incomplete"])
max(cat$max_length[cat$category == "incomplete"])
max(catC$max_nHPV)
max(catC$max_length)

### -------------------------------------------------------------------------------
### MULTI-BREAKPOINT INTEGRATION
### -------------------------------------------------------------------------------

# prepare the dataframe to include teh multi-breakpoint events
catM <- catC[catC$integration_type %in% c("multi-breakpoint_integration", "complex_ecDNA_integration"),]
catM <- separate(catM, col = bp_pair, into = c("bp1", "bp2"), sep = "_", remove = F)
catM <- separate(catM, col = bp1, into = c("bp1_chr", "bp1_pos"), sep = ":", remove = T)
catM <- separate(catM, col = bp2, into = c("bp2_chr", "bp2_pos"), sep = ":", remove = T)
catM$event.id <- paste0(catM$sample, "/", catM$event)

plot_bp_connections <- function(event){
    # For each event, find the unique breakpoints, and find the unique chromosomes, and then plot them out 
    event <- catM[catM$event.id == event,]
    chrs <- unique(c(event$bp1_chr, event$bp2_chr))
    sites <- unique(c(paste0(event$bp1_chr,":",event$bp1_pos), paste0(event$bp2_chr,":",event$bp2_pos)))
    sitesDF <- data.frame(sites = sites, chr = gsub("\\:.*", "", sites))
    
    # dataframes to be filled
    intSites <- data.frame(chr = gsub("\\:.*", "", sites), site = as.numeric(gsub(".*:", "", sites)), position = NA)
    chrSize <- data.frame(chr = chrs, start = 0, size = NA, intStart = NA, intEnd = NA)
    connect <- event[,c("bp1_chr", "bp1_pos", "bp2_chr", "bp2_pos", "category")]
    connect$bp1_position <- NA
    connect$bp2_position <- NA
    
    for (c in chrs) {
        subsite <- sitesDF$sites[sitesDF$chr == c]
        
        # find the minimum position on the chromosome 
        positions <- as.numeric(gsub(paste0(c,":"), "", subsite))
        minsite <- min(positions)
        maxsite <- max(positions)
        
        # add the length to the chrsize
        len <- maxsite - minsite
        chrSize$intStart[chrSize$chr == c] <- minsite
        chrSize$intEnd[chrSize$chr == c] <- maxsite
        chrSize$size[chrSize$chr == c] <- len/1000
        
        # subtract that value from the positions on that chromosome to zero them
        intSites$position[intSites$chr == c] <- (intSites$site[intSites$chr == c] - minsite)/1000
        
        # adjust the positions on the connections
        connect$bp1_position[connect$bp1_chr == c] <- (as.numeric(connect$bp1_pos[connect$bp1_chr == c]) - minsite)/1000
        connect$bp2_position[connect$bp2_chr == c] <- (as.numeric(connect$bp2_pos[connect$bp2_chr == c]) - minsite)/1000
    }
    
    # count how many connections for each int site
    count <- table(c(connect$bp1_pos, connect$bp2_pos))
    intSites$nConnect <- count[match(intSites$site, names(count))]
    intSites$connect <- ifelse(intSites$nConnect == 1, "one", ifelse(intSites$nConnect == 2, "two", "three_plus"))
    
    # seperate heterologous integrants
    connectA <- connect[connect$category != "heterologous",]
    connectB <- connect[connect$category == "heterologous",]
    
    # make plot
    ggplot() +
        # chromosome bars
        geom_segment(data = chrSize, aes(y = chr, yend = chr, x = start, xend = size), 
                     lineend = "round", color = "lightgrey", size = 3) + 
        geom_point(data = intSites, aes(y = chr, x = position, colour = connect), 
                   size = 3) +
        geom_curve(data = connectA, aes(x = bp1_position, xend = bp2_position, y = bp1_chr, yend = bp2_chr), linetype = 2, size = 0.5, curvature = -0.3, colour = "#fdb714") +
        geom_curve(data = connectB, aes(x = bp1_position, xend = bp2_position, y = bp1_chr, yend = bp2_chr), linetype = 2, size = 0.5, curvature = 0.4, colour = "#219EBC") +
        scale_colour_manual(values = c("black","#BE1E2D", "#ec8992")) +
        theme_bw() +
        labs(x = "position in event (kb)") +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.line = element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size = 10, colour = "black"),
              axis.title.x = element_text(size = 12, colour = "black", face = "bold"),
              axis.title.y = element_blank(), legend.position = "none")
}

    
plot_bp_connections("HTMCP-03-06-02238/event2")
plot_bp_connections(event = "HTMCP-03-06-02175/event2")
plot_bp_connections(event = "HTMCP-03-06-02058/event1")
plot_bp_connections(event = "TCGA-C5-A1BL/event4")
plot_bp_connections(event = "HTMCP-03-06-02267/event1")
plot_bp_connections(event = "HTMCP-03-06-02128/event2")    
plot_bp_connections(event = "HTMCP-03-06-02149/event1")   
plot_bp_connections(event = "HTMCP-03-06-02428/event1")  

    

    









    