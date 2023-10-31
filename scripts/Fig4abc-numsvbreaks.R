.libPaths("/path/to/libraries")

## --------------------------------------------
# Load libraries
## --------------------------------------------

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(RColorBrewer)
library(RColorBrewer)
library(cowplot)

## --------------------------------------------
# Load files
## --------------------------------------------

sum <- read.delim("/path/to/tables/integrationTypes.txt")
meta <- read.delim("/path/to/tables/metaDataAllSamples.txt")
sum$event.id <- paste0(sum$sample, '/', sum$event)

htmcp_path <- "/path/to/htmcp/call_integration/output"
tcga_path <- "/path/to/tcga/call_integration/output"

### HPV METHYLATION

# Find files with HPV methylation by read
Files1 <- grep("SVsubset.bedpe", list.files(htmcp_path, recursive = TRUE, full.names = TRUE), value = TRUE)
Files2 <- grep("SVsubset.bedpe", list.files(tcga_path, recursive = TRUE, full.names = TRUE), value = TRUE)
Files <- c(Files1, Files2)

# Extract library names from file paths
event <- gsub(paste0(htmcp_path, "/|", tcga_path, "/|/events|_SVsubset.bedpe"), "", Files)

# Import data
bedpe <- lapply(Files, read.delim, header = FALSE, sep = "\t")
names(bedpe) <- event
bedpe <- dplyr::bind_rows(bedpe, .id = "event.id")

# annotated HPV breakpoints
mb <- sum[sum$type == "multi-breakpoint_integration",]
mb$HPV.clade <- meta$HPV.clade[match(mb$sample, meta$Patient)]

## --------------------------------------------
# SUMMARIZE THE NUMBER OF SV BREAKS IN A EVENT
## --------------------------------------------

svTypes <- bedpe %>% 
    mutate(SV_pair = case_when(
        grepl("chr",V1) & grepl("chr",V4) ~ "human",
        grepl("HPV",V1) & grepl("chr",V4) | grepl("chr",V1) & grepl("HPV",V4) ~ "human_HPV",
        grepl("HPV",V1) & grepl("HPV",V4) ~ "HPV"
    )
)

svTypes$integration_type <- sum$type[match(svTypes$event.id, sum$event.id)]
svTypes$HPV.type <- sum$HPV.type[match(svTypes$event.id, sum$event.id)]
svTypes$HPV.clade <- meta$HPV.clade[match(svTypes$sample, meta$Patient)]
svTypes$sample <- sum$sample[match(svTypes$event.id, sum$event.id)]

# remove the HPV related SVs and replace with the breakpoint calls
nbreaksHuman <- svTypes %>% dplyr::filter(SV_pair == "human" & integration_type == "multi-breakpoint_integration") %>%
  dplyr::group_by(event.id, SV_pair) %>%
  dplyr::summarise(n_breaks = n())

# add the breakpoint calls
nBreaksHPV <- data.frame(event.id = mb$event.id, SV_pair = "human_HPV", n_breaks = mb$nsites)
nbreaks <- rbind(nbreaksHuman, nBreaksHPV)


nbreaksAll <- nbreaks %>%
  dplyr::group_by(event.id) %>%
  dplyr::summarise(n_breaks = sum(n_breaks))

order <- nbreaksAll$event.id[order(nbreaksAll$n_breaks, decreasing = T)]

nbreaks$event.id <- factor(nbreaks$event.id, levels = order)

nbreaksAnno <- gather(mb, variable, value, c("HPV.type", "HPV.clade", "sample"), factor_key=TRUE) 
nbreaksAnno$event.id <- factor(nbreaksAnno$event.id, order)
nbreaksAnno$variable <- factor(nbreaksAnno$variable, levels = c("HPV.clade", "HPV.type", "sample"))

nbreaks$HPV.type <- sum$HPV.type[match(as.character(nbreaks$event.id), sum$event.id)]

## --------------------------------------------
# PLOT
## --------------------------------------------

ann_colors <- list(HPV.clade = c(A7="#B55F8F", A9="#253083", Other="#A9A9A9"),
                   HPV.type = c(HPV16="#3953A4", HPV18="#9768ad", HPV45="#CC138C", HPV82="#369797",
                                HPV52="#0B8DCD", HPV31="#d2ecf9", HPV73="#737474", HPV68="#f3b2d4", HPV97="black", HPV58="#8bafba", HPV59="#ae3030"))
p1 <- nbreaks %>% 
    filter(HPV.type %in% c("HPV16", "HPV18") & SV_pair == "human_HPV") %>%
    ggplot(aes(x = HPV.type, y = n_breaks, colour = HPV.type)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(height = 0, width = 0.2, size =3, alpha=0.5) +
    theme_minimal() +
    xlab("HPV type") + 
    ylab("# of HPV breakpoint in event") +
    scale_colour_manual(values = c(ann_colors[["HPV.type"]])) +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size=12, face = "bold", colour = "black"), 
          axis.ticks.y = element_line(),
          axis.line = element_line(), 
          legend.position = "none") +
    stat_compare_means(method = "wilcox")

ggsave(p1, filename = "figurePanels/figure4/nbreaks_HPVtype.pdf", height = 2.7, width = 2.7, units = "in")


mean(nbreaks$n_breaks[nbreaks$HPV.type == "HPV16" & nbreaks$SV_pair == "human_HPV"])
mean(nbreaks$n_breaks[nbreaks$HPV.type == "HPV18" & nbreaks$SV_pair == "human_HPV"])

## --------------------------------------------
# PLOT
## --------------------------------------------
# find unique colours for the samples
n <- length(unique(nbreaksAnno$value[nbreaksAnno$variable == "sample"]))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col = sample(col_vector, n)
names(col) <- unique(mb$sample)

meta_colours <- list(values = c(A7="#B55F8F", A9="#253083", Other="#A9A9A9",
                                HPV16="#3953A4", HPV18="#9768ad", HPV45="#CC138C", HPV82="#323636", HPV33="#2ba8be", 
                                HPV52="#0B8DCD", HPV31="#d2ecf9", HPV73="#737474", HPV68="#f3b2d4", HPV97="black", HPV58="#8bafba", HPV59="#ae3030", col))

# sample information
p1 <- ggplot(nbreaksAnno, aes(y = variable, x = event.id, fill = value)) +
    geom_tile(colour = "black") + 
    theme_minimal()+
    scale_fill_manual(values = meta_colours[["values"]]) +
    theme(axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size = 12, colour = "black"),
          legend.text = element_text(size = 12, colour = "black"),
          legend.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          legend.key.size = unit(0.4, "cm"))

p2 <- ggplot(nbreaks, aes(y = n_breaks, x = event.id, fill = SV_pair)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(y = "# of breakpoints in the event", x = NULL) +
    scale_fill_manual(values = c("#457b9d","#e63946")) +
    theme(axis.line.x.bottom = element_line(),
          axis.line.y.left = element_line(),
          axis.ticks.y = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")

plot <- plot_grid(p2, p1, rel_heights = c(3,1), ncol = 1, align = "v", axis = "l")

ggsave(plot, filename = "figurePanels/figure4/numSVBreakpointsEvent.pdf", height = 4, width = 7, units = "in")

## --------------------------------------------
# CORRELATION PLOT
## --------------------------------------------

nbreaks2 <- spread(nbreaks, SV_pair, n_breaks)

plot2 <- ggscatter(nbreaks2, x = "human_HPV", y = "human",
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "#457b9d", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE) + stat_cor(method = "spearman") +
    labs(x = "# of HPV-human breakpoints", y = "# of human SV breakpoints")
    
ggsave(plot2, filename = "figurePanels/figure4/hpvSVCorrelation.pdf", height = 3.3, width = 3.3, units = "in")


# save table
write.table(nbreaks2,"tables/multi-breakpoint_svbreaks.txt", sep = "\t", col.names = T, row.names = F, quote = F)
