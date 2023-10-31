.libPaths("/path/to/libraries")

## ---------------------------------------------------------------------------
## Description of the cohort for ONT sequencing (Figure 1A)
## Vanessa Porter, June 2022
## ---------------------------------------------------------------------------

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyr))
suppressMessages(library(cowplot))

## ---------------------------------------------------------------------------
## Load in relevant files 
## ---------------------------------------------------------------------------

htmcp <- read.delim("/path/to/Supplementary_table_1.txt", header = T)
tcga <- read.csv("/path/to/tcga_metadata_table.csv", header = T)
meta <- read.delim("/path/to/tables/tcga_htmcp_covariates.txt")

tcgaRan <- read.delim("/path/to/ran_libraries.txt", header = F)
htmcpRan <- read.delim("/path/to/ran_libraries.txt", header = F)

# put all the samples together
all <- rbind(htmcpRan, tcgaRan)

## ---------------------------------------------------------------------------
## Combine with the meta data
## ---------------------------------------------------------------------------

# get the HTMCP only meta data
htmcpMeta <- htmcp[htmcp$Patient %in% all$V1,]

# subset the TCGA meta data
tcgaMeta <- tcga[tcga$TCGA %in% all$V1, ]

# load in all the events files
sum <- read.delim("/path/to/Supplementary_Table2.txt", header = T)
sumHTMCP <- sum[grep("HTMCP", sum$library),]
sumTCGA <- sum[grep("TCGA", sum$library),]

# get the HPV types
sumHTMCP <- sumHTMCP[!duplicated(sumHTMCP[,c(1,4)]),c(1,4)]
sumTCGA <- sumTCGA[!duplicated(sumTCGA[,c(1,4)]),c(1,4)]

## ---------------------------------------------------------------------------
## Fix the meta columns - HTMCP
## ---------------------------------------------------------------------------

# binarize linear values
htmcpMeta$Age <- ifelse(htmcpMeta$Age < 45, "Under45", 
                   ifelse(htmcpMeta$Age >= 45 & htmcpMeta$Age <= 65, "A45to65", "Over65"))
htmcpMeta$Max.APOBEC.score <- ifelse(htmcpMeta$Max.APOBEC.score < 0.2, "Under0.2", 
                                ifelse(htmcpMeta$Max.APOBEC.score <= 0.4 & htmcpMeta$Max.APOBEC.score >= 0.2, "A0.2to0.4","Over0.4"))
htmcpMeta$HRD.score <- ifelse(htmcpMeta$HRD.score < 10, "Under10", 
                         ifelse(htmcpMeta$HRD.score <= 30 & htmcpMeta$HRD.score >= 10, "A10to30","Over30"))
htmcpMeta$Ploidy <- ifelse(htmcpMeta$Ploidy == 1, "one", 
                      ifelse(htmcpMeta$Ploidy == 2, "two",
                             ifelse(htmcpMeta$Ploidy == 3, "three","four")))

# fix Stage values
htmcpMeta$Stage <- gsub(" ", "_", htmcpMeta$Stage)
htmcpMeta$Stage <- gsub("B[1-2].*", "", htmcpMeta$Stage)
htmcpMeta$Stage <- gsub("A[1-2].*", "", htmcpMeta$Stage)
htmcpMeta$Stage <- gsub("B|A", "", htmcpMeta$Stage)
htmcpMeta$Stage <- factor(htmcpMeta$Stage, levels = c("Stage_I", "Stage_II","Stage_III", "Stage_IV" ))

# check if the HPV type in the meta matches the integration
x <- htmcpMeta[order(as.character(htmcpMeta$Patient)),c(1,12)]
rownames(x) <- 1:nrow(x)
y <- sumHTMCP[order(as.character(sumHTMCP$library)),]

all(x$HPV.type == y$HPVchr) # true

# long format 
htmcpMeta <- htmcpMeta[,c("Patient","HIV.status", "Stage", "Grade","Final.histology", 
             "Age", "HPV.type", "HPV.clade", "Max.APOBEC.score", "HRD.score", "Ploidy")]
anno <- melt(htmcpMeta, id = "Patient")
colnames(anno) <- c("Patient", "row", "value")
anno$row <- factor(anno$row, levels = c("Ploidy","HRD.score", "Max.APOBEC.score","Estimated.Tumour.purity", 
                                        "Age", "Grade", "Stage", "Final.histology", "HIV.status", 
                                        "num.of.sites","HPV.clade","HPV.type"))

## ---------------------------------------------------------------------------
## Fix the meta columns - TCGA
## ---------------------------------------------------------------------------

# subset just the columns of interest 
tcgaMeta <- tcgaMeta[,c("TCGA", "Final.HPV.for.MCWRPPA", "Final.Hclade.for.MCWRPPA", "Histology", "Stage", "Age")]

# match the column names
colnames(tcgaMeta) <- c("Patient","HPV.type", "HPV.clade", "Final.histology", "Stage", "Age")

# binarize linear values
tcgaMeta$Age <- ifelse(tcgaMeta$Age < 45, "Under45", 
                        ifelse(tcgaMeta$Age >= 45 & tcgaMeta$Age <= 65, "A45to65", "Over65"))

# fix Stage values
tcgaMeta$Stage <- gsub("a[1-2].*", "", tcgaMeta$Stage)
tcgaMeta$Stage <- gsub("b[1-2].*", "", tcgaMeta$Stage)
tcgaMeta$Stage <- gsub("a|b", "", tcgaMeta$Stage)
tcgaMeta$Stage <- paste0("Stage_", tcgaMeta$Stage)
tcgaMeta$Stage <- factor(tcgaMeta$Stage, levels = c("Stage_I", "Stage_II","Stage_III", "Stage_IV" ))

# match the histology
tcgaMeta$Final.histology <- gsub("Squamous Cell Carcinoma", "Squamous", tcgaMeta$Final.histology)

# check if the HPV type in the meta matches the integration
x <- tcgaMeta[order(as.character(tcgaMeta$Patient)),c(1,2)]
rownames(x) <- 1:nrow(x)
y <- sumTCGA[order(as.character(sumTCGA$library)),]

all(x$HPV.type == y$HPVchr) # true

# long format 
annoT <- melt(tcgaMeta, id = "Patient")
colnames(annoT) <- c("Patient", "row", "value")
annoT$row <- factor(annoT$row, levels = c("Age", "Stage", "Final.histology", "HPV.clade", "HPV.type"))

## ---------------------------------------------------------------------------
## ONT Integration sites
## ---------------------------------------------------------------------------

# samples with no integration 
noint <- sum$library[sum$integration_type == "no_detected_integration"]


# summarize the integration events
summary1 <- sum[complete.cases(sum),] %>% 
  group_by(library, event) %>%
  summarize(n = n())

summary2 <- sum[complete.cases(sum),] %>% 
  group_by(library) %>%
  summarize(chromosomes = paste(unique(chr), collapse = ","))

summary3 <- sum[complete.cases(sum),] %>% 
  group_by(library, chr) %>%
  summarize(n.chr = n())

summary3 <- summary3 %>% 
  group_by(library) %>%
  summarize(n.chr = n())

summary <- summary1 %>% 
  group_by(library) %>%
  summarize(events = n(), sites = sum(n))

summary$HPV.type <- sum$HPVchr[match(summary$library, sum$library)]
summary$chromosomes <- summary2$chromosomes[match(summary$library, summary2$library)]
summary$n.chromosomes <- summary3$n.chr[match(summary$library, summary3$library)]


# add the unintegrated samples
summary <- rbind(summary, data.frame(library = noint, events = rep(0, length(noint)), sites = rep(0, length(noint)),
                                     HPV.type = rep(NA, length(noint)), chromosomes = rep(NA, length(noint)), n.chromosomes = rep(NA, length(noint))))

# split by cohort
tcgaSum <- summary[grep("TCGA", summary$library),]
htmcpSum <- summary[grep("HTMCP", summary$library),]

## ---------------------------------------------------------------------------
## Make the figures - HTMCP 
## ---------------------------------------------------------------------------

meta_colours <- list(values = c(Negative="black", Positive="#FACDBE", 
                                Adenocarcinoma="#9E9AC7", Adenosquamous="#DEC9E2", Squamous="#D2E1E9", Neuroendocrine="#6DA3C2", Undifferentiated="#A9A9A9",
                                G1="#FEE5CE", G2="#EAA973", G3="#F47A14",G4="#b35608",
                                Stage_I="#F4D3AE", Stage_II="#FAAD6C", Stage_III="#F06A22",Stage_IV="#8C2F1C",
                                Under45="white", A45to65="grey", Over65="black",
                                Under10="#ECE7F2", A10to30="#A6BDDB", Over30="#2B8CBE",
                                Under0.2="#FDE0DD", A0.2to0.4="#FA9FB5", Over0.4="#C51B8A", A7="#B55F8F", A9="#253083", Other="#A9A9A9",
                                HPV16="#3953A4", HPV18="#9768ad", HPV45="#CC138C", HPV82="#323636", HPV33="#2ba8be", 
                                HPV52="#0B8DCD", HPV31="#d2ecf9", HPV73="#737474", HPV68="#f3b2d4", HPV97="black", HPV58="#8bafba", HPV59="#ae3030",
                                two="#EFEDF5", three="#BCBDDC", four="#756BB1"))



# make the sample order
htmcpMeta$events <- htmcpSum$events[match(htmcpMeta$Patient, htmcpSum$library)]
htmcpMeta$sites <- htmcpSum$sites[match(htmcpMeta$Patient, htmcpSum$library)]
type <- table(factor(htmcpMeta$HPV.type))
type <- names(type[order(type, decreasing = T)])
htmcpMeta$HPV.type <- factor(htmcpMeta$HPV.type, levels = type)
htmcpMeta$HPV.clade <- factor(htmcpMeta$HPV.clade, levels = c("A9", "A7", "Other"))
order <- htmcpMeta %>%
  arrange(factor(HPV.type, levels = c("HPV16", "HPV18", "HPV45", "HPV58", "HPV31", "HPV82", "HPV52", "HPV68","HPV33", "HPV73")), desc(events), desc(sites)) %>%
  select(Patient)

# factor the sample order
anno$Patient <- factor(anno$Patient, levels = order$Patient)
htmcpSum$library <- factor(htmcpSum$library, levels = order$Patient)

all(levels(htmcpSum$library) == levels(anno$Patient))

# number of integration events
a1 <- ggplot(htmcpSum, aes(x = library, y = events)) +
  geom_bar(stat = "identity", colour = "black", fill = "grey", size = 0.25) +
  theme_minimal() +
  labs(y = "# of events") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_line(),
        legend.position = "top") 

# sample information
b1 <- ggplot(anno, aes(y = row, x = Patient, fill = value)) +
  geom_tile(colour = "black") + 
  theme_minimal()+
  scale_fill_manual(values = meta_colours[["values"]]) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.4, "cm"))

# number of integration breakpoints
c1 <- ggplot(htmcpSum, aes(x = library, y = sites)) +
  geom_bar(stat = "identity", colour = "black", fill = "grey", size = 0.25) +
  theme_minimal() +
  scale_y_reverse() +
  labs(y = "# of breakpoints") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_line(),
        legend.position = "none")

p1 <- plot_grid(a1, b1, c1, align = "v", ncol = 1)
p1
ggsave(filename = "figure1/htmcp_cohort_hpvInt.pdf", plot = p1, height = 6, width = 9.5, units = "in")


## ---------------------------------------------------------------------------
## Make a figure - TCGA
## ---------------------------------------------------------------------------

# make the sample order
tcgaMeta$events <- tcgaSum$events[match(tcgaMeta$Patient, tcgaSum$library)]
tcgaMeta$sites <- tcgaSum$sites[match(tcgaMeta$Patient, tcgaSum$library)]
typeT <- table(factor(tcgaMeta$HPV.type))
typeT <- names(typeT[order(typeT, decreasing = T)])
tcgaMeta$HPV.type <- factor(tcgaMeta$HPV.type, levels = typeT)
tcgaMeta$HPV.clade <- factor(tcgaMeta$HPV.clade, levels = c("A9", "A7", "Other"))
orderT <- tcgaMeta %>%
  arrange(factor(HPV.type, levels = c("HPV16", "HPV18", "HPV45", "HPV33", "HPV58", "HPV31", "HPV82", "HPV52", "HPV68", "HPV73")), desc(events), desc(sites)) %>%
  select(Patient)

# factor the sample order
annoT$Patient <- factor(annoT$Patient, levels = orderT$Patient)
tcgaSum$library <- factor(tcgaSum$library, levels = orderT$Patient)

all(levels(tcgaSum$library) == levels(annoT$Patient))

# number of integration events
a2 <- ggplot(tcgaSum, aes(x = library, y = events)) +
  geom_bar(stat = "identity", colour = "black", fill = "grey", size = 0.25) +
  theme_minimal() +
  labs(y = "# of events") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_line(),
        legend.position = "top") 

# sample information
b2 <- ggplot(annoT, aes(y = row, x = Patient, fill = value)) +
  geom_tile(colour = "black") + 
  theme_minimal()+
  scale_fill_manual(values = meta_colours[["values"]]) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.4, "cm"))

# number of integration breakpoints
c2 <- ggplot(tcgaSum, aes(x = library, y = sites)) +
  geom_bar(stat = "identity", colour = "black", fill = "grey", size = 0.25) +
  theme_minimal() +
  scale_y_reverse() +
  labs(y = "# of breakpoints") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black", face = "bold"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_line(),
        legend.position = "none")

p2 <- plot_grid(a2, b2, c2, align = "v", ncol = 1)
p2
ggsave(filename = "figure1/tcga_cohort_hpvInt.pdf", plot = p2, height = 6, width = 6, units = "in")

## ---------------------------------------------------------------------------
## N values
## ---------------------------------------------------------------------------

tcgaMetaSub <- tcgaMeta[tcgaMeta$Patient %in% annoT$Patient,]
htmcpMetaSub <- htmcpMeta[htmcpMeta$Patient %in% anno$Patient,]

###
### Meta Values
###

### HTMCP & TCGA 

# HPV type
tapply(c(table(tcgaMetaSub$HPV.type), table(htmcpMetaSub$HPV.type)), names(c(table(tcgaMetaSub$HPV.type), table(htmcpMetaSub$HPV.type))), sum)

# HPV clade
tapply(c(table(tcgaMetaSub$HPV.clade), table(htmcpMetaSub$HPV.clade)), names(c(table(tcgaMetaSub$HPV.clade), table(htmcpMetaSub$HPV.clade))), sum)

# Stage
tapply(c(table(tcgaMetaSub$Stage), table(htmcpMetaSub$Stage)), names(c(table(tcgaMetaSub$Stage), table(htmcpMetaSub$Stage))), sum)

# histology
tapply(c(table(tcgaMetaSub$Final.histology), table(htmcpMetaSub$Final.histology)), names(c(table(tcgaMetaSub$Final.histology), table(htmcpMetaSub$Final.histology))), sum)

# Age
tapply(c(table(tcgaMetaSub$Age), table(htmcpMetaSub$Age)), names(c(table(tcgaMetaSub$Age), table(htmcpMetaSub$Age))), sum)

### HTMCP ONLY

# HIV Status
table(htmcpMetaSub$HIV.status)

# Grade
table(htmcpMetaSub$Grade)

# HRD Score
table(htmcpMetaSub$HRD.score)

# HRD Score
table(htmcpMetaSub$Max.APOBEC.score)

# Ploidy
table(htmcpMetaSub$Ploidy)

###
### Integration
###

# number of integration sites
sum(summary$sites)

# number of integration events
sum(summary$events)

# number of samples with integration detected
nrow(summary[summary$sites > 0,])

# number of samples with episome HPV
nrow(summary[summary$sites == 0,])

## ---------------------------------------------------------------------------
## Supp Tables
## ---------------------------------------------------------------------------

# meta table changes
htmcpMetaSave <- htmcpMeta
htmcpMetaSave$Age <- htmcp$Age[match(htmcpMetaSave$Patient, htmcp$Patient)]
htmcpMetaSave$Max.APOBEC.score <- htmcp$Max.APOBEC.score[match(htmcpMetaSave$Patient, htmcp$Patient)]
htmcpMetaSave$HRD.score <- htmcp$HRD.score[match(htmcpMetaSave$Patient, htmcp$Patient)]
htmcpMetaSave$Ploidy <- htmcp$Ploidy[match(htmcpMetaSave$Patient, htmcp$Patient)]
htmcpMetaSave$HPV.integrated.chromosomes <- summary$chromosomes[match(htmcpMetaSave$Patient,summary$library)]

tcgaMetaSave <- tcgaMeta
tcgaMetaSave$Age <- tcga$Age[match(tcgaMetaSave$Patient, tcga$TCGA)]
tcgaMetaSave$HPV.integrated.chromosomes <- summary$chromosomes[match(tcgaMetaSave$Patient,summary$library)]
tcgaMetaSave$HIV.status <- NA
tcgaMetaSave$Grade <- NA
tcgaMetaSave$Max.APOBEC.score <- NA
tcgaMetaSave$HRD.score <- NA
tcgaMetaSave$Ploidy <- NA

# put together
tcgaMetaSave <- tcgaMetaSave[,colnames(htmcpMetaSave)]
metaSave <- rbind(htmcpMetaSave,tcgaMetaSave)

# save
write.table(metaSave, "tables/metaDataAllSamples.txt", quote = F, sep = "\t", col.names = T, row.names = F)

## ---------------------------------------------------------------------------
## Source table
## ---------------------------------------------------------------------------

# annotation
annoT$cohort <- "TCGA"
anno$cohort <- "HTMCP"
annoAll <- rbind(anno, annoT)
write.table(annoAll, file = "sourceTables/fig1a-1.txt", col.names = T, sep = "\t", row.names = F, quote = F)

# integration sites/events
tcgaSum$cohort <- "TCGA"
htmcpSum$cohort <- "HTMCP"
sumAll <- rbind(tcgaSum[,c(1:3, 7)], htmcpSum[,c(1:3, 7)])
write.table(sumAll, file = "sourceTables/fig1a-2.txt", col.names = T, sep = "\t", row.names = F, quote = F)







