.libPaths("/path/to/libraries")

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(cowplot))

## -------------------------------------------------------------------------
## Import Data
## -------------------------------------------------------------------------

# read bed files
dupR <- read.delim("/path/to/htmcp/call_integration/output/HTMCP-03-06-02068/bam/hpv_reads.bed", header = F)
ecR <- read.delim("/path/to/tcga/call_integration/output/TCGA-C5-A1MP/bam/hpv_reads.bed", header = F)
delR <- read.delim("/path/to/htmcp/call_integration/output/HTMCP-03-06-02063/events/event1.bed", header = F)

# depth files
dupD <- read.delim("/path/to/htmcp/call_integration/output/HTMCP-03-06-02068/depth/allreads_event_depth.txt", header = T)
ecD <- read.delim("/path/to/tcga/call_integration/output/TCGA-C5-A1MP/depth/allreads_event_depth.txt", header = T)
delD <- read.delim("/path/to/htmcp/call_integration/output/HTMCP-03-06-02063/depth/allreads_event1_depth.txt", header = T)

## -------------------------------------------------------------------------
## Make alignment figure - DUP
## -------------------------------------------------------------------------

## HPV read alignment 

# add a y value to the reads in order of size
dupR <- dupR[dupR$V1 == "chr1",]
dupR$length <- dupR$V3 - dupR$V2
dupR <- dupR %>%
  mutate(portion = 
           case_when(
             V2 < 209381334-5 ~ "up",
             V3 > 209387173+5 ~ "down",
             TRUE  ~  "down"
             )) 
dupR <- dupR %>% arrange(portion,desc(length))
yvals <- dupR %>%
  group_by(portion) %>%
  summarise(y = 1:n())
dupR$y <- yvals$y

# size of the aligned region
max(dupR$V3) - min(dupR$V2)

## figure

a1 <- ggplot(dupR, aes(ymin = y, ymax = y+0.75, xmin = V2, xmax = V3)) + 
  geom_rect(fill = "grey") + 
  scale_y_continuous(trans = "reverse") + 
  theme_void() 
  theme(panel.grid = element_blank(), 
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.y = element_blank())

a2 <- ggplot(dupD, aes(y = bam.all_reads.sorted.bam, x = POS)) + 
  geom_bar(stat = "identity", fill = "#5d5d5d", colour = "#5d5d5d") +
  theme_minimal() +
  xlim(c(min(dupR$V2), max(dupR$V3)))+
  labs(x = NULL, y = "read depth") +
    theme(axis.text = element_text(colour = "black", size = 10),
          axis.title = element_text(colour = "black", size = 10, face = "bold"),
          panel.grid = element_blank(), 
          axis.ticks.y = element_line(),
          axis.line = element_line())

p1 <- plot_grid(a2, a1, rel_heights = c(1,4), ncol = 1, align = 'v', axis = "lr")
ggsave(filename = "dup_align.pdf", plot = p1, height = 5, width = 5, units = "in")

## -------------------------------------------------------------------------
## Make alignment figure - ecDNA
## -------------------------------------------------------------------------

## HPV read alignment 

mid <- (73397788 - 73342042) + 73397788
# add a y value to the reads in order of size
ecR <- ecR[ecR$V1 == "chr13",]
ecR$length <- ecR$V3 - ecR$V2
ecR <- ecR %>%
  mutate(portion = 
           case_when(
             V2 < 73342042+10 ~ "up",
             TRUE  ~  "down"
           )) 
ecR <- ecR %>% arrange(portion,desc(length))
yvals1 <- ecR %>%
  group_by(portion) %>%
  summarise(y = 1:n())
ecR$y <- yvals1$y

# size of the aligned region
max(ecD$POS) - min(ecD$POS)

## figure

b1 <- ggplot(ecR, aes(ymin = y, ymax = y+0.75, xmin = V2, xmax = V3)) + 
  geom_rect(fill = "grey") + 
  scale_y_continuous(trans = "reverse") + 
  theme_void() +
  xlim(c(min(ecD$POS), max(ecD$POS)))
b1

b2 <- ggplot(ecD, aes(y = bam.all_reads.sorted.bam, x = POS)) + 
  geom_bar(stat = "identity", fill = "#5d5d5d", colour = "#5d5d5d") +
  theme_minimal() +
  labs(x = NULL, y = "read depth") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10, face = "bold"),
        panel.grid = element_blank(), 
        axis.ticks.y = element_line(),
        axis.line = element_line())
b2

p2 <- plot_grid(b2, b1, rel_heights = c(1,4), ncol = 1, align = 'v', axis = "lr")
p2
ggsave(filename = "ecDNA_align.pdf", plot = p2, height = 5, width = 5, units = "in")

## -------------------------------------------------------------------------
## Make alignment figure - deletion
## -------------------------------------------------------------------------

## HPV read alignment 

mid <- (73397788 - 73342042) + 73397788
# add a y value to the reads in order of size
delR <- delR[delR$V1 == "chr2",]
delR$length <- delR$V3 - delR$V2
delR <- delR %>%
  mutate(portion = 
           case_when(
             V2 < 141176425+10 ~ "up",
             TRUE  ~  "down"
           )) 

delR <- delR %>% arrange(portion,desc(length))
yvals2 <- delR %>%
  group_by(portion) %>%
  summarise(y = 1:n())
delR$y <- yvals2$y

# size of the aligned region
max(delD$POS) - min(delD$POS)

## figure

c1 <- ggplot(delR, aes(ymin = y, ymax = y+0.75, xmin = V2, xmax = V3)) + 
  geom_rect(fill = "grey") + 
  scale_y_continuous(trans = "reverse") + 
  theme_void() 
c1

c2 <- ggplot(delD, aes(y = bam.all_reads.sorted.bam, x = POS)) + 
  geom_bar(stat = "identity", fill = "#5d5d5d", colour = "#5d5d5d") +
  theme_minimal() +
  labs(x = NULL, y = "read depth") +
  xlim(c(min(delR$V2), max(delR$V3))) +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10, face = "bold"),
        panel.grid = element_blank(), 
        axis.ticks.y = element_line(),
        axis.ticks.x = element_line(),
        axis.line = element_line())
c2

p3 <- plot_grid(c2, c1, rel_heights = c(1,4), ncol = 1, align = 'v', axis = "lr")
p3
ggsave(filename = "del_align.pdf", plot = p3, height = 5, width = 5, units = "in")
