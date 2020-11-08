# Execute if packages are not installed on PC
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("DESeq2", version = "3.8")
# BiocManager::install("apeglm", version = "3.8")

source("/Users/leonfodoulian/scData/FixSizeAndSave.R")

require(DESeq2)
require(apeglm)
require(ggplot2)
require(extrafont)
require(ggtext)

setwd(dir = "/Users/leonfodoulian/scData/SARS_CoV_2_anosmia/Human_MOE_RE_RNAseq/")

# Read .rds file
# counts <- readRDS(file = "human_moe_and_re_bulk_rnaseq.rds")
counts <- readRDS(file = "human_moe_and_re_bulk_rnaseq.updated_annotation.rds")

# Transform .rds table to wide format; use only counts
counts.wf <- transform(
  tidyr::spread(data = counts[c("Sample", "gene_name", "Counts")],
                key = "Sample",
                value = "Counts",
                fill = 0,
                drop = FALSE),
  row.names = gene_name,
  gene_name = NULL
)

# Convert counts.wf to integer mode
counts.wf <- as.matrix(x = counts.wf)
mode(x = counts.wf) <- "integer"

# Prepare sample info table for DESeq2
samples <- colnames(x = counts.wf)
sample.info <- data.frame(samples = samples,
                          tissues = ifelse(test = endsWith(x = samples, suffix = "moe"),
                                           yes = "moe",
                                           no = "re"),
                          patients = gsub(pattern = "_.*",
                                          replacement = "",
                                          x = samples),
                          row.names = samples,
                          stringsAsFactors = TRUE)

# Filter un-expressed genes
is.expressed <- rowSums(x = counts.wf) > 0
counts.wf <- counts.wf[is.expressed, ]

# Create DESeqDataSet from count matrix
dds <- DESeqDataSetFromMatrix(countData = counts.wf,
                              colData = sample.info,
                              design = ~ tissues)
dds$tissues <- factor(x = dds$tissues,
                      levels = c("re","moe"))

# Perform nbinomWaldTest() on the DESeqDataSet by testing for tissues only
dds <- estimateSizeFactors(object = dds)
dds <- estimateDispersions(object = dds)
dds <- nbinomWaldTest(object = dds,
                      betaPrior = FALSE)

# Save dispersion plot
pdf(
  # file = "Human_MOE_RE_RNAseq_DESeq2_apeglm_dispersion_plot_20200417.pdf",
  file = "Human_MOE_RE_RNAseq_updated_annotation_DESeq2_apeglm_dispersion_plot_20200505.pdf",
  width = 5,
  height = 5
) 
plotDispEsts(object = dds)
dev.off()

# Perform log fold change shrinkage using apeglm
res.apeglm <- lfcShrink(dds = dds,
                        coef = "tissues_moe_vs_re",
                        type = "apeglm")
dds.res.apeglm <- as.data.frame(x = res.apeglm)

write.csv(x = dds.res.apeglm,
          # file = "Human_MOE_RE_RNAseq_DESeq2_apeglm_20200417.csv",
          file = "Human_MOE_RE_RNAseq_updated_annotation_DESeq2_apeglm_20200505.csv",
          quote = FALSE,
          row.names = TRUE)

# Figure 1, Panel A
# Volcano plot
dds.res.apeglm$gene.name <- rownames(x = dds.res.apeglm)

dds.res.apeglm$gene.type <- ifelse(test = grepl(pattern = "^OR[0-9]",
                                                x = dds.res.apeglm$gene.name),
                                   yes = "Olfactory receptor genes",
                                   no = "Other")
dds.res.apeglm$gene.type[dds.res.apeglm$gene.name == "OR2A1-AS1"] <- "Other"
dds.res.apeglm$gene.type <- ifelse(test = dds.res.apeglm$gene.name %in% c("OMP", "CNGA2", "ANO2", "ERMN"),
                                   yes = "Olfactory marker genes",
                                   no = dds.res.apeglm$gene.type)
dds.res.apeglm$gene.type[dds.res.apeglm$gene.name == "ACE2"] <- "ACE2"
dds.res.apeglm$gene.type[dds.res.apeglm$gene.name == "TMPRSS2"] <- "TMPRSS2"
dds.res.apeglm$gene.type <- factor(x = dds.res.apeglm$gene.type,
                                   levels = c("ACE2", "TMPRSS2", "Olfactory marker genes", "Olfactory receptor genes", "Other"))

dds.res.apeglm <- dds.res.apeglm[order(dds.res.apeglm$gene.type, decreasing = TRUE),]

min.x.val <- floor(x = min(dds.res.apeglm$log2FoldChange,
                           na.rm = TRUE))
if (min.x.val %% 2 != 0) {
  min.x.val <- min.x.val - 1
}
max.x.val <- ceiling(x = max(dds.res.apeglm$log2FoldChange,
                             na.rm = TRUE))
if (max.x.val %% 2 != 0) {
  max.x.val <- max.x.val + 1
}
x.breaks <- seq(from = min.x.val,
                to = max.x.val,
                by = 2)

max.y.val <- pretty(x = max(-log10(x = dds.res.apeglm$pvalue),
                            na.rm = TRUE),
                    n = 2)[2]
y.breaks <-  c(0, seq(from = 10,
                      to = max.y.val,
                      by = 10))

volcano.plot <- ggplot(data = dds.res.apeglm[!is.na(x = dds.res.apeglm$padj),],
                       mapping = aes(x = log2FoldChange,
                                     y = -log10(x = padj),
                                     colour = gene.type)) +
  geom_point(size = 1.5,
             shape = 19,
             stroke = 0,
             alpha = 0.75) +
  geom_segment(data = data.frame(p.val.threshold = -log10(x = 0.05),
                                 label = "adjusted p-value = 0.05"),
               mapping = aes(x = -Inf,
                             xend = Inf,
                             y = p.val.threshold,
                             yend = p.val.threshold,
                             linetype = label),
               colour = "red",
               size = 0.25,
               inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = dds.res.apeglm[!dds.res.apeglm$gene.type %in% c("Olfactory receptor genes", "Other"),],
                           mapping = aes(label = gene.name,
                                         segment.colour = gene.type),
                           size = 6 * (1/72 * 25.4),
                           family = "Arial",
                           fontface = "italic",
                           min.segment.length = unit(x = 0,
                                                     units = "mm"),
                           force = TRUE,
                           direction = "both",
                           max.overlaps = Inf,
                           seed = 2,
                           force_pull = 100, # repelled labels version, arbitrarily chosen
                           # box.padding = unit(x = 4.3, units = "mm"), # repelled labels version of 20200505
                           box.padding = unit(x = 2.5, units = "mm"), # repelled labels version of 20200512
                           show.legend = FALSE) +
  scale_x_continuous(breaks = x.breaks,
                     expand = c(0,0)) +
  scale_y_continuous(breaks = y.breaks,
                     expand = c(0,0)) +
  scale_colour_manual(values = c("ACE2" = "violetred1",
                                 "TMPRSS2" = "lightslateblue",
                                 "Olfactory marker genes" = "#064E40",
                                 "Olfactory receptor genes" = "limegreen",
                                 "Other" = "lightgrey"),
                      breaks = c("ACE2",
                                 "TMPRSS2",
                                 "Olfactory marker genes",
                                 "Olfactory receptor genes"),
                      # labels = c("ACE2" = "ACE2",
                      #            "TMPRSS2" = "TMPRSS2",
                      #            "Olfactory marker genes" = "Olfactory marker genes",
                      #            "Olfactory receptor genes" = "Olfactory receptor genes"), # without ggtext
                      labels = c("ACE2" = "<span style='color:violetred1'>ACE2</span>",
                                 "TMPRSS2" = "<span style='color:lightslateblue'>TMPRSS2</span>",
                                 "Olfactory marker genes" = "Olfactory marker genes",
                                 "Olfactory receptor genes" = "Olfactory receptor genes"), # with ggtext
                      aesthetics = c("colour", "segment.colour")) +
  scale_linetype_manual(values = c("adjusted p-value = 0.05" = "dashed")) +
  expand_limits(x = c(min.x.val, max.x.val),
                y = c(0, max.y.val)) +
  labs(x = expression(paste(log[2], " fold change (MOE vs RE)", sep = "")),
       y = expression(paste("-", log[10], "(adjusted p-value)", sep = ""))) +
  guides(colour = guide_legend(order = 1,
                               override.aes = list(size = 1.5,
                                                   alpha = 1)),
         linetype = guide_legend(order = 2,
                                 label.theme = element_text(colour = "red",
                                                            size = 6,
                                                            family = "Arial",
                                                            face = "plain",
                                                            margin = margin(t = 0,
                                                                            r = 0,
                                                                            b = 0,
                                                                            l = -1,
                                                                            unit = "mm")))) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(axis.title = element_text(colour = "black",
                                  size = 7,
                                  family = "Arial"),
        axis.text = element_text(colour = "black",
                                 size = 6,
                                 family = "Arial"),
        axis.line = element_line(colour = "black", 
                                 size = rel(0.5),
                                 lineend = "square"),
        axis.ticks = element_line(colour = "black", 
                                  size = rel(0.5),
                                  lineend = "square"),
        axis.ticks.length = unit(x = 1,
                                 units = "mm"),
        legend.title = element_blank(),
        # legend.text = element_text(colour = "black",
        #                            size = 6,
        #                            family = "Arial",
        #                            face = "italic",
        #                            margin = margin(t = 0,
        #                                            r = 0,
        #                                            b = 0,
        #                                            l = -1,
        #                                            unit = "mm")), # without ggtext
        legend.text = element_markdown(colour = "black",
                                       size = 6,
                                       family = "Arial",
                                       face = "italic",
                                       margin = margin(t = 0,
                                                       r = 0,
                                                       b = 0,
                                                       l = -1,
                                                       unit = "mm")), # with ggtext
        legend.key.size = unit(x = 3,
                               units = "mm"),
        legend.spacing = unit(x = -4,
                              units = "mm"),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        legend.margin = margin(t = 0,
                               r = 0,
                               b = 4,
                               l = 2,
                               unit = "mm"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        plot.title = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(t = 0,
                             r = 0,
                             b = 0,
                             l = 0,
                             unit = "mm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

FixSizeAndSave(plot = volcano.plot,
               # filename = "Human_MOE_RE_RNAseq_DESeq2_apeglm_volcanoplot_20200417.pdf",
               # filename = "Human_MOE_RE_RNAseq_DESeq2_apeglm_volcanoplot_20200417_repelled_labels.pdf",
               # filename = "Human_MOE_RE_RNAseq_updated_annotation_DESeq2_apeglm_volcanoplot_20200505_repelled_labels.pdf",
               # filename = "Human_MOE_RE_RNAseq_updated_annotation_DESeq2_apeglm_volcanoplot_20200512_repelled_labels.pdf", # without ggtext
               filename = "Human_MOE_RE_RNAseq_updated_annotation_DESeq2_apeglm_volcanoplot_20200512_repelled_labels_ggtext.pdf", # with ggtext
               is.ggassemble = FALSE,
               # panel.width = 6, # 20200505
               # panel.height = 4, # 20200505
               panel.width = 4.89, # 20200512
               panel.height = 3.43, # 20200512
               unit.use = "cm",
               margin = 0.5,
               use.ggsave = TRUE,
               useDingbats = FALSE)
