.libPaths("/path/to/libraries")

## ---------------------------------------------------------------------------
## Stats of Nanopore Runs
## Vanessa Porter, Oct. 2021
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggsci))
suppressMessages(library(tidyr))

## ---------------------------------------------------------------------------
## Load Data
## ---------------------------------------------------------------------------

htmcp <- read.csv("/path/to/tables/NCI_SAIC_HIV_Cervical-Jan-20-2022T10_56_38-nanopore.csv")
tcga <- read.csv("/path/to/tables/Janet Rader-Jan-20-2022T10_57_57-nanopore.csv")
mix <- read.csv("/path/to/tables/Marco_Research-Jan-24-2022T16_39_47-nanopore.csv")

htmcp <- rbind(htmcp, mix[grep("F", mix$lib),])
tcga <- rbind(tcga, mix[grep("D", mix$lib),])

htmcp$cohort <- "HTMCP"
tcga$cohort <- "TCGA"

data <- rbind(htmcp, tcga)
data <- data[-grep("I", data$lib),]
data$coverage <- data$read_length_sum/3049315783 # length of hg38 genome

dataLong <- melt(data[c("cohort", "lib", "coverage", "read_length_sum",
                        "n50", "read_length_median", 
                        "chimerism_prop", "error_rate_by_qualimap")], id.vars = c("cohort", "lib"))

## ---------------------------------------------------------------------------
## Final cohort list
## ---------------------------------------------------------------------------

cohort <- read.delim("/path/to/sourceTables/fig1a-2.txt", header = T)
htmcp_libs <- read.delim("/path/to/htmcp/call_integration/sample_txtfiles/HTMCP_libs_info.txt", header = F)
tcga_libs <- read.delim("/path/to/tcga/call_integration/sample_names/tcga_sequenced_libraries.txt", header = F)
libs <- rbind(htmcp_libs[,2:3], tcga_libs[,2:3])
libs$V3 <- gsub("MaM_485", "TCGA-C5-A2LX",libs$V3)
libs$V3 <- gsub("HTMCP-03-06-02235", "HTMCP-03-06-02210",libs$V3)
libs$V3 <- gsub("MaM_486", "TCGA-C5-A2LY",libs$V3)
libs_used <- libs[libs$V3 %in% cohort$library,]

## ---------------------------------------------------------------------------
## Get the duplicate QC stats
## ---------------------------------------------------------------------------

dups <- data$lib[duplicated(data$lib)]

mergeQCfiles <- list.files("/path/to/htmcp/nanopore_qc_stats/merged_data", full.names = T)
mergeQCnames <- gsub("/path/to/htmcp/nanopore_qc_stats/merged_data/|.csv", "", mergeQCfiles)

mergeQC <- lapply(mergeQCfiles, read.csv, header = TRUE)
names(mergeQC) <- mergeQCnames
mergeQC <- bind_rows(mergeQC, .id = "library")
mergeQC <- mergeQC[mergeQC$X_key == "pass",]
mergeQC$cohort <- data$cohort[match(mergeQC$library, data$lib)]
mergeQC$cohort[mergeQC$library %in% c("D77540-D81468", "D77541-D81469")] <- "HTMCP"
mergeQC$coverage <- mergeQC$sum / 3088269832
mergeQC$chimerism_prop <- NA
mergeQC$error_rate_by_qualimap <- NA

mergeQCsub <- mergeQC[c("cohort", "library", "coverage", "sum",
          "n50", "median", "chimerism_prop", "error_rate_by_qualimap")]
colnames(mergeQCsub) <- c("cohort", "lib", "coverage", "read_length_sum",
                          "n50", "read_length_median", 
                          "chimerism_prop", "error_rate_by_qualimap")

mergeDataLong <- melt(mergeQCsub, id.vars = c("cohort", "lib"))

# remove the duplicate data
datasub <- data[!data$lib %in% c(dups,"D77540", "D81468", "D77541", "D81469"),c("cohort", "lib", "coverage", "read_length_sum",
                                      "n50", "read_length_median", 
                                      "chimerism_prop", "error_rate_by_qualimap")]

datasub <- datasub[datasub$lib %in% libs_used$V2,]

# merge together
alldata <- rbind(datasub,mergeQCsub)
datasubLong <- melt(alldata, id.vars = c("cohort", "lib"))

## ---------------------------------------------------------------------------
## Stats
## ---------------------------------------------------------------------------

# coverage
mean(alldataLongSub$value[alldataLongSub$variable == "read_length_sum" & alldataLongSub$cohort == "HTMCP"])/1000000000
range(alldataLongSub$value[alldataLongSub$variable == "read_length_sum" & alldataLongSub$cohort == "HTMCP"])/1000000000
median(alldataLongSub$value[alldataLongSub$variable == "coverage" & alldataLongSub$cohort == "HTMCP"])

mean(alldataLongSub$value[alldataLongSub$variable == "read_length_sum" & alldataLongSub$cohort == "TCGA"])/1000000000
range(alldataLongSub$value[alldataLongSub$variable == "read_length_sum" & alldataLongSub$cohort == "TCGA"])/1000000000
median(alldataLongSub$value[alldataLongSub$variable == "coverage" & alldataLongSub$cohort == "TCGA"])

# N50
median(alldataLongSub$value[alldataLongSub$variable == "n50" & alldataLongSub$cohort == "HTMCP"])/1000
range(alldataLongSub$value[alldataLongSub$variable == "n50" & alldataLongSub$cohort == "HTMCP"])/1000

median(alldataLongSub$value[alldataLongSub$variable == "n50" & alldataLongSub$cohort == "TCGA"])/1000
range(alldataLongSub$value[alldataLongSub$variable == "n50" & alldataLongSub$cohort == "TCGA"])/1000

## ---------------------------------------------------------------------------
## Make figures
## ---------------------------------------------------------------------------

ggplot(data, aes(y = error_rate_by_qualimap, x = cohort, fill = cohort)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.2, size =2, alpha=0.5) +
  scale_fill_manual(values = c("#DD4A48", "#C0D8C0")) +
  theme_bw() + 
  labs(x = NULL, y = "Base Error Rate", fill = NULL, colour = NULL) +
  theme(axis.text = element_text(colour = "black", size = 13),
        axis.title = element_text(colour = "black", size = 14),
        legend.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(), 
        axis.ticks.y = element_line(),
        axis.line = element_line(),
        legend.position = "none")

p <- ggplot(alldataLongSub, aes(y = value, x = cohort, fill = cohort)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height = 0, width = 0.2, size =2, alpha=0.5) +
  facet_wrap(~ variable, nrow = 1, scales = "free") +
  scale_fill_manual(values = c("#DD4A48", "#C0D8C0")) +
  theme_bw() + 
  labs(x = NULL, y = NULL, fill = NULL, colour = NULL) +
  theme(axis.text.y = element_text(colour = "black", size = 13),
        axis.text.x = element_blank(),
        axis.title = element_text(colour = "black", size = 14),
        legend.text = element_text(size = 12, colour = "black"),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black", size = 0.75),
        axis.ticks.y = element_line(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "#F5EEDC", colour = "black", size = 0.75))

ggsave(p, file = "figurePanels/figure1/suppFig1.pdf", width = 16, height = 2.6)

## ---------------------------------------------------------------------------
## Save table
## ---------------------------------------------------------------------------

write.table(alldata, file = "tables/qcStats.txt", quote = F, col.names = T, row.names = F, sep = "\t")



