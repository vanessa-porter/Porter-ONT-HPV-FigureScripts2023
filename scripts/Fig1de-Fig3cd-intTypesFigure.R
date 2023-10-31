.libPaths("/path/to/libraries")

#Note: these packages need to be installed.
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggupset))
suppressMessages(library(viridis))
suppressMessages(library(FSA))
suppressMessages(library(ggpubr))

## -------------------------------------------------------------------------
## Import Integration Type Data
## -------------------------------------------------------------------------

# Define the directories for the first and second datasets
dir1 <- "/path/to/htmcp/call_integration/output"
dir2 <- "/path/to/tcga/call_integration/output"

# List all files containing "integration_type.txt" in both directories
all_files <- c(
  list.files(dir1, recursive = TRUE, full.names = T, pattern = "integration_type.txt"),
  list.files(dir2, recursive = TRUE, full.names = T, pattern = "integration_type.txt")
)

# Extract sample names from the file paths
names <- gsub(paste0(dir1, "/|", dir2, "/|intType/|/integration_type.txt"), "", all_files)

# Read data from all files and store in a list
data_list <- lapply(all_files, read.delim, header=F)

# Set names for the list elements based on sample names
names(data_list) <- names

# Combine the dataframes in the list into a single dataframe, adding an "event" column
int <- bind_rows(data_list, .id = "event")

# Extract the "event" column and separate it into "sample" and "event" columns
int$event <- as.character(int$event)
int <- separate(int, event, into = c("sample", "event"), sep = "/")
colnames(int)[3] <- "type"

## -------------------------------------------------------------------------
## Import Integration Summary Data
## -------------------------------------------------------------------------

# load in all the events files
ontSitesFiles <- c(grep("events/summary.txt", list.files(dir1, recursive = T, full.names = T), value = T),
                   grep("events/summary.txt", list.files(dir2, recursive = T, full.names = T), value = T))

# Extract sample names from the file paths
names2 <- gsub(paste0(dir1, "/|", dir2, "/|/events/summary.txt"), "", ontSitesFiles)

# import ont sites
sum_list <- lapply(ontSitesFiles, read.delim, header=T)

# Set names for the list elements based on sample names
names(sum_list) <- names2

# Combine the dataframes in the list into a single dataframe, adding an "library" column
sum <- bind_rows(sum_list, .id = "library")

# remove samples with no integration
#sum <- sum[complete.cases(sum),]

## -------------------------------------------------------------------------
## Fix the integration event calls for select samples
## -------------------------------------------------------------------------

#find the missing events
missing <- data.frame(event = unique(paste0(sum$library, ":", sum$event))[!unique(paste0(sum$library, ":", sum$event)) %in% paste0(int$sample, ":", int$event)])
missing <- missing %>% separate(event, c("sample", "event"), sep = ":")
missing$nsites <- summary1$n[match(paste0(missing$sample, ":", missing$event), paste0(summary1$library, ":", summary1$event))]

missing <- missing %>%
  mutate(type = case_when(
    is.na(nsites) ~ "no_detected_integration"))

# add missing calls to df
int <- rbind(int[,1:3], missing[,c(1,2,4)])
int$HPV.type <- sum$HPVchr[match(int$sample, sum$library)]
int$type[is.na(int$type)] <- "no_detected_integration"

# fix a sample with atypical deletion read mapping
int$type[int$sample == "TCGA-C5-A1MQ" & int$event == "event1"] <- "deletion_integration"

# fix a non-repeat sample
int$type[int$sample == "HTMCP-03-06-02267" & int$event == "event1"] <- "multi-breakpoint_integration" 

# add nsites
# number of integration events
summary1 <- sum[complete.cases(sum),] %>% 
  group_by(library, event) %>%
  summarize(n = n())

int$nsites <- summary1$n[match(paste0(int$sample, ":", int$event), paste0(summary1$library, ":", summary1$event))]
int$nsites[is.na(int$nsites)] <- 0

# add HPV type
int$HPV.type <- sum$HPVchr[match(paste0(int$sample, ":", int$event), paste0(sum$library, ":", sum$event))]

# Recall "complex ecDNA event" 
int$type[int$type == "complex_ecDNA_integration" & int$nsites > 2] <- "multi-breakpoint_integration"
int$type[int$type == "complex_ecDNA_integration" & int$nsites == 2] <- "ecDNA_integration"

# only integration
intOnly <- int[!int$type %in% c("no_detected_integration", "unmatched_integration"),]

# add to summary 
sum$integration_type <- int$type[match(paste0(sum$library, sum$event), paste0(int$sample, int$event))]

## -------------------------------------------------------------------------
## MAKE SUMMARY TABLES
## -------------------------------------------------------------------------

summary <- summary1 %>% 
  group_by(library) %>%
  summarize(events = n(), sites = sum(n))

summary2 <- sum[complete.cases(sum),] %>% 
  group_by(library, event) %>%
  summarize(chromosomes = paste(unique(chr), collapse = ","))

summary3 <- sum[complete.cases(sum),] %>% 
  group_by(library, event, chr) %>%
  summarize(n.chr = n())

summary3 <- summary3 %>% 
  group_by(library, event) %>%
  summarize(n.chr = n())

summary4 <- sum[complete.cases(sum),] %>% 
  group_by(library, event, integration_type, chr) %>%
  summarize(size = max(pos) - min(pos))

summary$HPV.type <- sum$HPVchr[match(summary$library, sum$library)]
summary$chromosomes <- summary2$chromosomes[match(summary$library, summary2$library)]
summary$n.chromosomes <- summary3$n.chr[match(summary$library, summary3$library)]


## -------------------------------------------------------------------------
## Figures - Fig 1ef
## -------------------------------------------------------------------------

intOnly$HPV.type.simp <- ifelse(!intOnly$HPV.type %in% c("HPV16", "HPV45", "HPV18"), "Other_HPV", intOnly$HPV.type)
intOnly <- intOnly %>%
  mutate(HPV.clade = case_when(
    HPV.type %in% c("HPV16", "HPV52", "HPV31", "HPV33", "HPV35", "HPV53", "HPV58")  ~ "A9",
    HPV.type %in% c("HPV18", "HPV45", "HPV59", "HPV68")  ~ "A7",
    HPV.type %in% c("HPV82", "HPV51", "HPV26", "HPV30", "HPV69", "HPV73")  ~ "Other"))

## bubble plot
ntypes <- intOnly %>%
  group_by(HPV.type.simp, type) %>%
  summarise(n = n())
nHPV <- intOnly %>%
  group_by(HPV.type.simp) %>%
  summarise(n = n())
nAll <- intOnly %>%
  group_by(type) %>%
  summarise(n = n())

ntypes$nHPV <- nHPV$n[match(ntypes$HPV.type.simp, nHPV$HPV.type.simp)]
ntypes$percent <- ntypes$n / ntypes$nHPV

ntypes$type <- factor(ntypes$type, levels = rev(c("translocation_integration", "duplication_integration", "deletion_integration",
                                              "ecDNA_integration", "multi-breakpoint_integration", "repeat_integration")))

p <- ggplot(ntypes, aes(x = HPV.type.simp, y = type, size = percent)) +
  geom_point(colour="black",pch=21, fill = "grey") + 
  theme_minimal() +
  scale_size(range = c(2,15)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_text(size = 12, face = "bold", colour = "black"))

p
ggsave(filename = "figurePanels/figure1/intType_HPVs.pdf", plot = p, height = 5, width = 8, units = "in")

## Integration types bar plot
nAll <- nAll %>% arrange(desc(n))
nAll$type <- gsub("_integration|_integrated", "", nAll$type)
nAll$type <- gsub("linear_dup", "duplication", nAll$type)
nAll$type <- factor(nAll$type, levels = nAll$type)

p1 <- ggplot(nAll, aes(x = type, y = n)) +
  geom_bar(stat = "identity") + 
  theme_minimal() +
  labs(y = "# of events", x = NULL, fill = NULL) +
  theme(axis.text.y = element_text(colour = "black", size = 13),
        axis.text.x = element_text(colour = "black", size = 13, angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(colour = "black", size = 14, face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.ticks.y = element_line(),
        axis.line = element_line())
p1
ggsave(filename = "figurePanels/figure1/intType_num.pdf", plot = p1, height = 4, width = 6, units = "in")

## UPSET PLOT DIFF TYPES

p2 <- intOnly %>%
  distinct(sample, type) %>%
  group_by(sample) %>%
  summarise(list = list(type)) %>%
  ggplot(aes(x=list)) +
  geom_bar() +
  scale_x_upset() +
  theme_minimal() +
  theme(axis.text.y = element_text(colour = "black", size = 13),
        axis.title = element_blank(),
        legend.text = element_text(size = 12, colour = "black"),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.ticks.y = element_line(),
        axis.line = element_line(),
        plot.margin = margin(1,1,1.5,6, "cm"))
p2
ggsave(filename = "figurePanels/figure1/intType_upset.pdf", plot = p2, height = 4, width = 6, units = "in")

list <- intOnly %>%
  distinct(sample, type) %>%
  group_by(sample) %>%
  summarise(list = list(type))

## save environment
save.image("intTypesImage.RData")

## -------------------------------------------------------------------------
## Stats - Fig 1ef
## -------------------------------------------------------------------------

# number of samples
length(unique(int$sample))
length(unique(int$sample[grepl("TCGA",int$sample)]))
length(unique(int$sample[grepl("HTMCP",int$sample)]))

# number of integrated samples
length(unique(int$sample[int$type != "no_detected_integration"]))
length(unique(int$sample[grepl("TCGA",int$sample) & int$type != "no_detected_integration"]))
length(unique(int$sample[grepl("HTMCP",int$sample) & int$type != "no_detected_integration"]))

# number of different types of events
nevents <- length(int$sample[int$type != "no_detected_integration"])
nsites <- nrow(sum[sum$integration_type != "no_detected_integration",])
sum(int$type == "multi-breakpoint_integration")/nevents
sum(int$type == "deletion_integration")
sum(int$type == "duplication_integration")
sum(int$type == "ecDNA_integration")
sum(int$type == "translocation_integration")
sum(int$type == "repeat_integration")
sum(int$nsites > 1)/nevents

# types of events by HPV type
summary(factor(int$type[int$HPV.type == "HPV16"]))
summary(factor(int$type[int$HPV.type == "HPV18"]))
summary(factor(int$type[int$HPV.type == "HPV45"]))

## -------------------------------------------------------------------------
## Source tables and supp tables
## -------------------------------------------------------------------------

# integration types
write.table(nAll, file = "sourceTables/fig1e.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## integration types by HPV type
write.table(ntypes, file = "sourceTables/fig1f.txt", col.names = T, row.names = F, sep = "\t", quote = F)

# integration site locations
write.table(sum, "tables/intSitesSummaryAll.txt", quote = F, sep = "\t", col.names = T, row.names = F)

# integration event types
write.table(int, "tables/integrationTypes.txt", quote = F, sep = "\t", col.names = T, row.names = F)

## -------------------------------------------------------------------------
## Figures - Fig 2
## -------------------------------------------------------------------------

summaryTwo <- summary1[summary1$n == 2,]
summaryTwo$dist <- NA

for (i in summaryTwo$library) {
  for (j in summaryTwo$event) {
    breaks <- sum$pos[sum$library == i & sum$event == j]
    dist <- max(breaks) - min(breaks)
    summaryTwo$dist[summaryTwo$library == i & summaryTwo$event == j] <- dist
  }
}
summaryTwo$intType <- intOnly$type[match(paste0(summaryTwo$library, summaryTwo$event), paste0(intOnly$sample, intOnly$event))]
summaryTwo_sub <- summaryTwo[summaryTwo$intType %in% c("ecDNA_integration", "duplication_integration", "deletion_integration"),]
summaryTwo_sub$intType <- gsub("_integration", "", summaryTwo_sub$intType)

# stat test
kruskal.test(dist ~ intType,
             data = summaryTwo_sub)
dunnTest(dist ~ intType,
         data=summaryTwo_sub,
         method="bh") 

# n vals
table(summaryTwo_sub$intType)

## figure
p3 <- ggplot(summaryTwo_sub, aes(x = intType, y = dist, fill = intType)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.2, size = 3, alpha=0.5) +
  theme_bw() + 
  scale_fill_manual(values = c("#DD4A48", "#C0D8C0", "#f5eedc")) +
  scale_y_log10() +
  labs(x = "Two breakpoint integration types", y = "log10(distance between breakpoints)", fill = NULL, colour = NULL) +
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(), 
        axis.ticks.y = element_line(),
        axis.line = element_line(),
        legend.position = "none")+
  stat_compare_means(method = "kruskal")
p3
ggsave(filename = "figure2/twoBreakBoxplot.pdf", plot = p3, height = 5, width = 5, units = "in")

## -------------------------------------------------------------------------
## genic region of two-break events
## -------------------------------------------------------------------------

dir1 <- "/path/to/htmcp/call_integration/output/"
dir2 <- "/path/to/tcga/call_integration/output/"

files1 <- grep("genic_test/event_location.txt",list.files(dir1, recursive = T), value = T)
files2 <- grep("genic_test/event_location.txt",list.files(dir2, recursive = T), value = T)
name1 <- gsub("intType/", "", files1)
name1 <- gsub("/genic_test/event_location.txt", "", name1)
name2 <- gsub("intType/", "", files2)
name2 <- gsub("/genic_test/event_location.txt", "", name2)

reg1 <- NULL
for (i in 1:length(files1)) {
  reg1[[i]] <- read.delim(paste0(dir1, files1[i]), header = F)
}
names(reg1) <- name1
reg1 <- bind_rows(reg1, .id = "event")


reg2 <- NULL
for (i in 1:length(files2)) {
  reg2[[i]] <- read.delim(paste0(dir2, files2[i]), header = F)
}
names(reg2) <- name2
reg2 <- bind_rows(reg2, .id = "event")

reg <- rbind(reg1, reg2)
reg <- reg %>% separate(event, c("sample", "event"), sep = "/")

# put together with the summaryTwo
summaryTwo_sub$region <- reg$V1[match(paste0(summaryTwo_sub$library, summaryTwo_sub$event), paste0(reg$sample, reg$event))]

## figure
p4 <- ggplot(summaryTwo_sub, aes(x= region,  group=intType)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count", colour = "black") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "% of events", fill=NULL, x = "genomic region") +
    scale_fill_manual(values = c("#f5cac3", "#84a59d", "#f28482")) +
  theme_minimal()+
  facet_grid(~intType) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.y = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 12, hjust = 1, angle = 45),
        axis.title = element_text(colour = "black", size = 12, face = "bold"),
        strip.text = element_text(colour = "black", size = 12, face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(), 
        axis.ticks.y = element_line(),
        axis.line = element_line(),
        legend.position = "none")
ggsave(filename = "figurePanels/figure2/twoBreakEventRegion.pdf", plot = p4, height = 4, width = 7, units = "in")

## p value test between ecDNA vs. deletions/duplications
test <- summaryTwo_sub
test$intType <- ifelse(test$intType != "ecDNA", "deldup", "ecDNA")
mat <- test %>%
  group_by(intType, region) %>%
  summarise(n = n())

# n vals
table(summaryTwo_sub$intType)

# average ecDNA size
median(summaryTwo_sub$dist[summaryTwo_sub$intType == "ecDNA"])

mat <- as.data.frame(spread(mat, key = intType, value = n))
rownames(mat) <- mat$region
mat <- as.matrix(mat[,-1])
chisq.test(mat)

## save table
write.table(summaryTwo_sub, file = "tables/twoBreakIntegrationCharacteristics.txt", quote = F, sep = "\t", col.names = T, row.names = F)

