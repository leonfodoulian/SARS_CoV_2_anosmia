source("/Users/leonfodoulian/scData/PlotExpressionOntSNE.R")
source("/Users/leonfodoulian/scData/PlotClustersOntSNE.R")
source("/Users/leonfodoulian/scData/FixSizeAndSave.R")

require(Matrix)
require(ggplot2)
require(RColorBrewer)
require(extrafont)
require(ggtext)

setwd(dir = "/Users/leonfodoulian/scData/SARS_CoV_2_anosmia/Allen_Brain_scRNAseq/Human_Multiple_Cortical_Areas_SMART-seq/")

if (!dir.exists(paths = "./Figures_for_paper")) {
  dir.create(path = "./Figures_for_paper")
}

# Read count matrix file
counts.matrix <- read.table(file = "./raw_data/matrix.csv",
                            header = TRUE,
                            sep = ",",
                            row.names = 1,
                            check.names = FALSE)

# Transform count matrix to sparse matrix
if (!inherits(x = counts.matrix, what = "dgCMatrix")) {
  counts.matrix <- as(object = as.matrix(x = counts.matrix),
                      Class = "dgCMatrix")
}

# Transpose count matrix
counts.matrix <- t(x = counts.matrix)

# Normalize counts.matrix by library size
norm.counts.matrix <- t(x = t(x = counts.matrix) / colSums(x = counts.matrix)) * 1e4

rm(counts.matrix) # save memory on the laptop

# Save normalized count matrix
saveRDS(object = norm.counts.matrix,
        file = "Human_Multiple_Cortical_Areas_SMART-seq_normalized_counts.rds")
# Read normalized count matrix file
norm.counts.matrix <- readRDS(file = "Human_Multiple_Cortical_Areas_SMART-seq_normalized_counts.rds")

# Read meta.data table
meta.data <- read.csv(file = "./raw_data/metadata.csv",
                      header = TRUE,
                      stringsAsFactors = FALSE)
rownames(x = meta.data) <- meta.data$sample_name

# Read tsne table
tsne.data <- read.table(file = "./raw_data/tsne.csv",
                        header = TRUE,
                        sep = ",",
                        stringsAsFactors = FALSE)
rownames(x = tsne.data) <- tsne.data$sample_name

# Subset norm.counts.matrix and keep only clustered cells
norm.counts.matrix <- norm.counts.matrix[,tsne.data$sample_name]

# Subset meta.data and keep only clustered cells
meta.data.sub <- meta.data[tsne.data$sample_name,]


# Figure 6, Panel A
# Plot Clusters on t-SNE

cluster.data <- data.frame(Cluster = meta.data.sub$subclass_label,
                           row.names = meta.data.sub$sample_name)

subclass.cols <- grep(pattern = "^subclass",
                      x = colnames(x = meta.data.sub))
class.cols <- grep(pattern = "^class",
                   x = colnames(x = meta.data.sub))
cluster.details <- meta.data.sub[!duplicated(x = meta.data.sub$subclass_label), c(subclass.cols, class.cols)]
rownames(x = cluster.details) <- cluster.details$subclass_label
cluster.details$class_new_order <- cluster.details$class_order
cluster.details$class_new_order[cluster.details$class_new_order == 2] <- 1
cluster.details$class_new_order[cluster.details$class_new_order == 4] <- 2
cluster.details <- cluster.details[order(cluster.details$class_new_order, cluster.details$subclass_order),]
cluster.details$subclass_new_order <- 1:nrow(x = cluster.details)

cluster.colors <- cluster.details$subclass_color
names(x = cluster.colors) <- cluster.details$subclass_label
cluster.breaks <- cluster.details$subclass_label
cluster.labels <- cluster.details$subclass_label
names(x = cluster.labels) <- cluster.labels

plot.cluster.tsne <- PlotClustersOntSNE(
  tsne.data = tsne.data[,c("tsne_1", "tsne_2")],
  cluster.data = cluster.data,
  cellnames.as.rownames.in.tsne.data = TRUE,
  cellnames.as.rownames.in.cluster.data = TRUE,
  tsne_1_name = "tsne_1",
  tsne_2_name = "tsne_2",
  cluster.ident.name = "Cluster",
  cluster.colors = cluster.colors,
  cluster.breaks = cluster.breaks,
  cluster.labels = cluster.labels,
  add.border.to.points = FALSE,
  point.size = 0.5,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  alpha.use = 0.75,
  legend.title.face = "plain",
  font.family = "Arial",
  import.font = FALSE,
  fix.aspect.ratio = FALSE,
  legend.point.size = 2,
  axis.title.size = 6,
  legend.title.size = 7,
  legend.text.size = 6,
  legend.text.space = unit(x = 3, units = "mm"),
  axis.line.size = rel(0.5),
  legend.ncol = 2,
  arrow.size = rel(2),
  range.scale = 0.3
)

FixSizeAndSave(plot = plot.cluster.tsne,
               # filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_tSNE_subclass_plot_20200521_figure6_panelA.pdf",
               filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_tSNE_subclass_plot_20200521_ggtext_figure6_panelA.pdf",
               is.ggassemble = FALSE,
               panel.width = 5,
               panel.height = 5,
               unit.use = "cm",
               margin = 0.1,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 6, Panel B
# Plot RBFOX3 marker gene on t-SNE

gradient.color <- c(viridis::cividis(2)[2], rev(x = viridis::magma(11))[-c(1:2)])

data.use.tsne.coexpr <- norm.counts.matrix[c("RBFOX3", "ACE2", "TMPRSS2"),]

ACE2.TMPRSS2 <- Matrix::Matrix(data = ifelse(test = data.use.tsne.coexpr["ACE2",, drop = FALSE] * data.use.tsne.coexpr["TMPRSS2",, drop = FALSE] > 0,
                                             yes = 1,
                                             no = 0),
                               nrow = 1,
                               dimnames = list("ACE2/TMPRSS2", colnames(x = data.use.tsne.coexpr)))

data.use.tsne.coexpr <- rbind(data.use.tsne.coexpr, ACE2.TMPRSS2)


plot.rbfox3.expression <- PlotExpressionOntSNE(
  data.use = data.use.tsne.coexpr,
  genes.use = rownames(x = data.use.tsne.coexpr)[c(1,1)], # size do not match between cowplot and patchwork; trick to match size
  tsne.data = tsne.data[,c("tsne_1", "tsne_2")],
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "tsne_1",
  tsne_2_name = "tsne_2",
  is.log.transformed = FALSE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 0.5,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = gradient.color,
  plot.title.size = 7,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title = "norm. expr",
  legend.title.size = 6,
  legend.text.size = 6,
  legend.bar.width = 2,
  legend.bar.height = 0.3,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 2,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.rbfox3.expression,
               # filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_tSNE_RBFOX3_expression_plot_20200521_figure6_panelB.pdf",
               filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_tSNE_RBFOX3_expression_plot_20200521_ggtext_figure6_panelB.pdf",
               is.ggassemble = FALSE, # trick to match size with previous plots
               panel.width = 4.8,
               panel.height = 4.8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 6, Panel C
# Plot ACE2 and TMPRSS2 marker genes on t-SNE

plot.markers.expression <- PlotExpressionOntSNE(
  data.use = data.use.tsne.coexpr,
  genes.use = rownames(x = data.use.tsne.coexpr)[2:3],
  ggtext.gene.names = c("ACE2" = "<span style='color:violetred1'>ACE2</span>",
                        "TMPRSS2" = "<span style='color:lightslateblue'>TMPRSS2</span>"), # with ggtext
  tsne.data = tsne.data[,c("tsne_1", "tsne_2")],
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "tsne_1",
  tsne_2_name = "tsne_2",
  is.log.transformed = FALSE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 0.5,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = gradient.color,
  plot.title.size = 7,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title = "norm. expr",
  legend.title.size = 6,
  legend.text.size = 6,
  legend.bar.width = 2,
  legend.bar.height = 0.3,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 2,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.markers.expression,
               # filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_tSNE_ACE2_TMPRSS2_expression_plot_20200521_figure6_panelC.pdf", # without ggtext
               filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_tSNE_ACE2_TMPRSS2_expression_plot_20200521_ggtext_figure6_panelC.pdf", # with ggtext
               is.ggassemble = TRUE,
               panel.width = 4.8,
               panel.height = 4.8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 6, Panel D
# Plot coexpression of ACE2 with TMPRSS2 on t-SNE

plot.markers.coexpression <- PlotExpressionOntSNE(
  data.use = data.use.tsne.coexpr,
  genes.use = rownames(x = data.use.tsne.coexpr)[c(4,4)], # size do not match between cowplot and patchwork; trick to match size
  ggtext.gene.names = c("ACE2/TMPRSS2" = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>",
                        "ACE2/TMPRSS2" = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>"), # with ggtext
  tsne.data = tsne.data[,c("tsne_1", "tsne_2")],
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "tsne_1",
  tsne_2_name = "tsne_2",
  is.log.transformed = FALSE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 0.5,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = gradient.color,
  plot.title.size = 7,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title = "norm. expr",
  legend.title.size = 6,
  legend.text.size = 6,
  legend.bar.width = 2,
  legend.bar.height = 0.3,
  legend.position = "none",
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 2,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.markers.coexpression,
               # filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_tSNE_ACE2_TMPRSS2_coexpressing_cells_plot_20200521_figure6_panelD.pdf", # without ggtext
               filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_tSNE_ACE2_TMPRSS2_coexpressing_cells_plot_20200521_ggtext_figure6_panelD.pdf", # with ggtext
               is.ggassemble = FALSE, # trick to match size with previous plots
               panel.width = 4.8,
               panel.height = 4.8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 6, Panel E
# Plot average expression of ACE2 and TMPRSS2 per cluster

expr.data <- as.data.frame(x = as.matrix(x = t(x = norm.counts.matrix[c("ACE2", "TMPRSS2"),])))
expr.data$sample_name <- rownames(x = expr.data)
expr.data <- merge(x = expr.data,
                   y = meta.data.sub,
                   by = "sample_name",
                   sort = FALSE)
expr.data.list <- split(x = expr.data,
                        f = expr.data$subclass_label)
expr.data.list <- lapply(X = expr.data.list,
                         FUN = function(expr.data.sub) {
                           data.out <- data.frame(ACE2.mean = mean(x = expr.data.sub$ACE2),
                                                  TMPRSS2.mean = mean(x = expr.data.sub$TMPRSS2),
                                                  subclass_label = unique(x = expr.data.sub$subclass_label))
                           return(data.out)
                         })
expr.data.mean <- dplyr::bind_rows(expr.data.list,
                                   .id = NULL)
max.ace2 <- round(x = max(expr.data.mean$ACE2.mean),
                  digits = 3)
breaks.ace2 <- c(0, round(x = max.ace2/2, digits = 3), max.ace2)
max.tmprss2 <- round(x = max(expr.data.mean$TMPRSS2.mean),
                     digits = 2)
breaks.tmprss2 <- c(0, round(x = max.tmprss2/2, digits = 3), max.tmprss2)

labels.to.plot <- c("L5/6 NP",
                    "Oligodendrocyte",
                    "Astrocyte",
                    "Endothelial",
                    "L5 ET",
                    "OPC",
                    "Microglia",
                    "L4 IT")

mean.expr.plot <- ggplot(data = expr.data.mean,
                         mapping = aes(x = ACE2.mean,
                                       y = TMPRSS2.mean,
                                       color = subclass_label)) +
  geom_point(size = 1) +
  ggrepel::geom_text_repel(data = expr.data.mean[expr.data.mean$subclass_label %in% labels.to.plot,],
                           mapping = aes(label = subclass_label,
                                         segment.colour = subclass_label),
                           size = 6 * (1/72 * 25.4),
                           family = "Arial",
                           min.segment.length = unit(x = 0,
                                                     units = "mm"),
                           nudge_y = 0.01 * max.tmprss2,
                           force = TRUE,
                           direction = "both",
                           max.overlaps = Inf,
                           seed = 0,
                           show.legend = FALSE) +
  scale_x_continuous(limits = c(0, max.ace2),
                     breaks = breaks.ace2,
                     labels = breaks.ace2,
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, max.tmprss2),
                     breaks = breaks.tmprss2,
                     labels = breaks.tmprss2,
                     expand = c(0,0)) +
  scale_color_manual(values = cluster.colors,
                     aesthetics = c("colour", "segment.colour"),
                     name = NULL) +
  # labs(x = expression(paste(italic("ACE2"), " (mean expression)", sep = "")),
  #      y = expression(paste(italic("TMPRSS2"), " (mean expression)", sep = ""))) + # without ggtext
  labs(x = "<span style='color:violetred1'>*ACE2*</span> (mean expression)",
       y = "<span style='color:lightslateblue'>*TMPRSS2*</span> (mean expression)") + # with ggtext
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    # axis.title = ggplot2::element_text(colour = "black",
    #                                        size = 7,
    #                                        family = "Arial"), # without ggtext
    axis.title.x = ggtext::element_markdown(colour = "black",
                                            size = 7,
                                            family = "Arial"), # with ggtext
    axis.title.y = ggtext::element_markdown(colour = "black",
                                            size = 7,
                                            family = "Arial"), # with ggtext
    axis.text = ggplot2::element_text(colour = "black",
                                      size = 6,
                                      family = "Arial"),
    axis.line = ggplot2::element_line(colour = "black", 
                                      size = rel(0.5),
                                      lineend = "square"),
    axis.ticks = ggplot2::element_line(colour = "black", 
                                       size = rel(0.5),
                                       lineend = "square"),
    axis.ticks.length = ggplot2::unit(x = 1, units = "mm"),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(colour = "black",
                                        size = 6,
                                        family = "Arial"),
    legend.key.size = ggplot2::unit(x = 3,
                                    units = "mm"),
    legend.position = "none",
    panel.grid = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.box.background = element_blank(),
    panel.border = ggplot2::element_blank())

FixSizeAndSave(plot = mean.expr.plot,
               # filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_ACE2_TMPRSS2_mean_expression_per_subclass_20200521_figure6_panelE.pdf", # without ggtext
               filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_ACE2_TMPRSS2_mean_expression_per_subclass_20200521_ggtext_figure6_panelE.pdf", # with ggtext
               is.ggassemble = FALSE,
               panel.width = 4,
               panel.height = 4,
               unit.use = "cm",
               margin = 0.5,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 6, Panel F
# Plot the number of cells expressing ACE2 and TMPRSS2 in all subclasses

percentage.data <- as.data.frame(x = as.matrix(x = t(x = data.use.tsne.coexpr[c("ACE2", "TMPRSS2", "ACE2/TMPRSS2"),])))
percentage.data$sample_name <- rownames(x = percentage.data)
percentage.data <- merge(x = percentage.data,
                         y = meta.data.sub,
                         by = "sample_name",
                         sort = FALSE)
percentage.data.list <- split(x = percentage.data,
                              f = percentage.data$subclass_label)
percentage.data.list <- lapply(X = percentage.data.list,
                               FUN = function(percentage.data.sub) {
                                 out.data <- data.frame(subclass_label = unique(x = percentage.data.sub$subclass_label),
                                                        cluster.size = nrow(x = percentage.data.sub),
                                                        ACE2.number = sum(percentage.data.sub$ACE2 > 0),
                                                        ACE2.percentage = sum(percentage.data.sub$ACE2 > 0) / nrow(x = percentage.data.sub) * 100,
                                                        TMPRSS2.number = sum(percentage.data.sub$TMPRSS2 > 0),
                                                        TMPRSS2.percentage = sum(percentage.data.sub$TMPRSS2 > 0) / nrow(x = percentage.data.sub) * 100,
                                                        ACE2.TMPRSS2.number = sum(percentage.data.sub$`ACE2/TMPRSS2` > 0),
                                                        ACE2.TMPRSS2.percentage = sum(percentage.data.sub$`ACE2/TMPRSS2` > 0) / nrow(x = percentage.data.sub) * 100)
                                 return(out.data)
                               })
percentage.data <- dplyr::bind_rows(percentage.data.list,
                                    .id = NULL)
percentage.data$coexpressing <- percentage.data$ACE2.TMPRSS2.number > 0
percentage.data$subclass_label <- factor(x = percentage.data$subclass_label,
                                         levels = cluster.breaks)

max.y.ace2.number <- pretty(x = max(percentage.data$ACE2.number),
                            n = 2)[2]
breaks.ace2.number <- seq(from = 0,
                          to = max.y.ace2.number,
                          by = 100)
labels.ace2.number <- ifelse(test = breaks.ace2.number %in% c(0, max.y.ace2.number),
                             yes = breaks.ace2.number,
                             no = "")

max.y.tmprss2.number <- pretty(x = max(percentage.data$TMPRSS2.number),
                               n = 2)[2]
breaks.tmprss2.number <- seq(from = 0,
                             to = max.y.tmprss2.number,
                             by = 100)
labels.tmprss2.number <- ifelse(test = breaks.tmprss2.number %in% c(0, max.y.tmprss2.number),
                                yes = breaks.tmprss2.number,
                                no = "")

max.y.coexpr.number <- pretty(x = max(percentage.data$ACE2.TMPRSS2.number),
                              n = 2)[2]
breaks.coexpr.number <- seq(from = 0,
                            to = max.y.coexpr.number,
                            by = 10)
labels.coexpr.number <- ifelse(test = breaks.coexpr.number %in% c(0, max.y.coexpr.number),
                               yes = breaks.coexpr.number,
                               no = "")

plot.ace2.number <- ggplot(data = percentage.data,
                           mapping = aes(x = subclass_label,
                                         y = ACE2.number)) +
  geom_bar(stat = "identity",
           colour = NA,
           fill = "violetred1",
           width = 0.5) +
  scale_y_continuous(breaks = breaks.ace2.number,
                     labels = labels.ace2.number,
                     expand = c(0,0)) +
  expand_limits(y = c(0, max.y.ace2.number)) +
  labs(title = "ACE2") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black",
                                   size = 6,
                                   family = "Arial"),
        axis.line = element_line(colour = "black", 
                                 size = rel(0.5),
                                 lineend = "square"),
        axis.ticks = element_line(colour = "black", 
                                  size = rel(0.5),
                                  lineend = "square"),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(x = 1,
                                 units = "mm"),
        legend.position = "none",
        plot.title = element_text(colour = "violetred1",
                                  size = 8,
                                  family = "Arial",
                                  face = "italic",
                                  hjust = 0.5,
                                  vjust = 0),
        plot.background = element_blank(),
        plot.margin = margin(t = 1,
                             r = 1,
                             b = 1,
                             l = 1,
                             unit = "mm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

plot.tmprss2.number <- ggplot(data = percentage.data,
                              mapping = aes(x = subclass_label,
                                            y = TMPRSS2.number)) +
  geom_bar(stat = "identity",
           colour = NA,
           fill = "lightslateblue",
           width = 0.5) +
  scale_y_continuous(breaks = breaks.tmprss2.number,
                     labels = labels.tmprss2.number,
                     expand = c(0,0)) +
  expand_limits(y = c(0, max.y.tmprss2.number)) +
  labs(title = "TMPRSS2",
       y = "Expressing cells (#)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black",
                                    size = 7,
                                    family = "Arial"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black",
                                   size = 6,
                                   family = "Arial"),
        axis.line = element_line(colour = "black", 
                                 size = rel(0.5),
                                 lineend = "square"),
        axis.ticks = element_line(colour = "black", 
                                  size = rel(0.5),
                                  lineend = "square"),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(x = 1,
                                 units = "mm"),
        legend.position = "none",
        plot.title = element_text(colour = "lightslateblue",
                                  size = 8,
                                  family = "Arial",
                                  face = "italic",
                                  hjust = 0.5,
                                  vjust = 0),
        plot.background = element_blank(),
        plot.margin = margin(t = 1,
                             r = 1,
                             b = 1,
                             l = 1,
                             unit = "mm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

plot.coexpr.number <- ggplot(data = percentage.data,
                             mapping = aes(x = subclass_label,
                                           y = ACE2.TMPRSS2.number,
                                           fill = subclass_label)) +
  geom_bar(stat = "identity",
           colour = NA,
           width = 0.5) +
  scale_y_continuous(breaks = breaks.coexpr.number,
                     labels = labels.coexpr.number,
                     expand = c(0,0)) +
  scale_fill_manual(values = cluster.colors[levels(x = percentage.data$subclass_label)]) +
  expand_limits(y = c(0, max.y.coexpr.number)) +
  # labs(title = "ACE2/TMPRSS2") + # without ggtext
  labs(title = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>") + # with ggtext
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = cluster.colors[levels(x = percentage.data$subclass_label)],
                                   size = 6,
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1,
                                   family = "Arial"),
        axis.text.y = element_text(colour = "black",
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
        legend.position = "none",
        # plot.title = element_text(colour = "black",
        #                           size = 8,
        #                           family = "Arial",
        #                           face = "italic",
        #                           hjust = 0.5,
        #                           vjust = 0), # without ggtext
        plot.title = element_markdown(colour = "black",
                                      size = 8,
                                      family = "Arial",
                                      face = "italic",
                                      hjust = 0.5,
                                      vjust = 0), # with ggtext
        plot.background = element_blank(),
        plot.margin = margin(t = 1,
                             r = 1,
                             b = 1,
                             l = 1,
                             unit = "mm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

plot.comb.number <- patchwork::wrap_plots(plot.ace2.number,
                                          plot.tmprss2.number,
                                          plot.coexpr.number,
                                          ncol = 1)

FixSizeAndSave(plot = plot.comb.number,
               # filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_ACE2_TMPRSS2_number_of_expressing_cells_all_subclasses_plot_20200521_figure6_panelF.pdf", # without ggtext
               filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_ACE2_TMPRSS2_number_of_expressing_cells_all_subclasses_plot_20200521_ggtext_figure6_panelF.pdf", # with ggtext
               is.ggassemble = TRUE,
               panel.width = 9.5,
               panel.height = 1.5,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 6, Panel G
# Plot the percentage of cells expressing ACE2 and TMPRSS2 in coexpressing subclasses

percentage.data.coexpressing <- percentage.data[percentage.data$coexpressing,]
percentage.data.coexpressing$subclass_label <- droplevels(x = percentage.data.coexpressing$subclass_label)

max.y.ace2.percentage <- ceiling(x = max(percentage.data.coexpressing$ACE2.percentage))
breaks.ace2.percentage <- seq(from = 0,
                              to = max.y.ace2.percentage,
                              by = 1)
labels.ace2.percentage <- ifelse(test = breaks.ace2.percentage %in% c(0, max.y.ace2.percentage),
                                 yes = breaks.ace2.percentage,
                                 no = "")

max.y.tmprss2.percentage <- ceiling(x = max(percentage.data.coexpressing$TMPRSS2.percentage))
breaks.tmprss2.percentage <- seq(from = 0,
                                 to = max.y.tmprss2.percentage,
                                 by = 1)
labels.tmprss2.percentage <- ifelse(test = breaks.tmprss2.percentage %in% c(0, max.y.tmprss2.percentage),
                                    yes = breaks.tmprss2.percentage,
                                    no = "")

max.y.coexpr.percentage <- ceiling(x = max(percentage.data.coexpressing$ACE2.TMPRSS2.percentage))
breaks.coexpr.percentage <- seq(from = 0,
                                to = max.y.coexpr.percentage,
                                by = 0.1)
labels.coexpr.percentage <- ifelse(test = breaks.coexpr.percentage %in% c(0, max.y.coexpr.percentage),
                                   yes = breaks.coexpr.percentage,
                                   no = "")

plot.ace2.percentage <- ggplot(data = percentage.data.coexpressing,
                               mapping = aes(x = subclass_label,
                                             y = ACE2.percentage)) +
  geom_bar(stat = "identity",
           colour = NA,
           fill = "violetred1",
           width = 0.5) +
  scale_y_continuous(breaks = breaks.ace2.percentage,
                     labels = labels.ace2.percentage,
                     expand = c(0,0)) +
  expand_limits(y = c(0, max.y.ace2.percentage)) +
  labs(title = "ACE2") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black",
                                   size = 6,
                                   family = "Arial"),
        axis.line = element_line(colour = "black", 
                                 size = rel(0.5),
                                 lineend = "square"),
        axis.ticks = element_line(colour = "black", 
                                  size = rel(0.5),
                                  lineend = "square"),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(x = 1,
                                 units = "mm"),
        legend.position = "none",
        plot.title = element_text(colour = "violetred1",
                                  size = 8,
                                  family = "Arial",
                                  face = "italic",
                                  hjust = 0.5,
                                  vjust = 0),
        plot.background = element_blank(),
        plot.margin = margin(t = 1,
                             r = 1,
                             b = 1,
                             l = 1,
                             unit = "mm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

plot.tmprss2.percentage <- ggplot(data = percentage.data.coexpressing,
                                  mapping = aes(x = subclass_label,
                                                y = TMPRSS2.percentage)) +
  geom_bar(stat = "identity",
           colour = NA,
           fill = "lightslateblue",
           width = 0.5) +
  scale_y_continuous(breaks = breaks.tmprss2.percentage,
                     labels = labels.tmprss2.percentage,
                     expand = c(0,0)) +
  expand_limits(y = c(0, max.y.tmprss2.percentage)) +
  labs(title = "TMPRSS2",
       y = "Expressing cells (%)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black",
                                    size = 7,
                                    family = "Arial"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black",
                                   size = 6,
                                   family = "Arial"),
        axis.line = element_line(colour = "black", 
                                 size = rel(0.5),
                                 lineend = "square"),
        axis.ticks = element_line(colour = "black", 
                                  size = rel(0.5),
                                  lineend = "square"),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(x = 1,
                                 units = "mm"),
        legend.position = "none",
        plot.title = element_text(colour = "lightslateblue",
                                  size = 8,
                                  family = "Arial",
                                  face = "italic",
                                  hjust = 0.5,
                                  vjust = 0),
        plot.background = element_blank(),
        plot.margin = margin(t = 1,
                             r = 1,
                             b = 1,
                             l = 1,
                             unit = "mm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

plot.coexpr.percentage <- ggplot(data = percentage.data.coexpressing,
                                 mapping = aes(x = subclass_label,
                                               y = ACE2.TMPRSS2.percentage,
                                               fill = subclass_label)) +
  geom_bar(stat = "identity",
           colour = NA,
           width = 0.5) +
  scale_y_continuous(breaks = breaks.coexpr.percentage,
                     labels = labels.coexpr.percentage,
                     expand = c(0,0)) +
  scale_fill_manual(values = cluster.colors[levels(x = percentage.data.coexpressing$subclass_label)]) +
  expand_limits(y = c(0, max.y.coexpr.percentage)) +
  # labs(title = "ACE2/TMPRSS2") + # without ggtext
  labs(title = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>") + # with ggtext
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = cluster.colors[levels(x = percentage.data.coexpressing$subclass_label)],
                                   size = 6,
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1,
                                   family = "Arial"),
        axis.text.y = element_text(colour = "black",
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
        legend.position = "none",
        # plot.title = element_text(colour = "black",
        #                           size = 8,
        #                           family = "Arial",
        #                           face = "italic",
        #                           hjust = 0.5,
        #                           vjust = 0), # without ggtext
        plot.title = element_markdown(colour = "black",
                                      size = 8,
                                      family = "Arial",
                                      face = "italic",
                                      hjust = 0.5,
                                      vjust = 0), # with ggtext
        plot.background = element_blank(),
        plot.margin = margin(t = 1,
                             r = 1,
                             b = 1,
                             l = 1,
                             unit = "mm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

plot.comb.percentage <- patchwork::wrap_plots(plot.ace2.percentage,
                                              plot.tmprss2.percentage,
                                              plot.coexpr.percentage,
                                              ncol = 1)

FixSizeAndSave(plot = plot.comb.percentage,
               # filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_ACE2_TMPRSS2_percentage_of_expressing_cells_coexpressing_subclasses_plot_20200521_figure6_panelG.pdf", # without ggtext
               filename = "./Figures_for_paper/Human_Multiple_Cortical_Areas_SMART-seq_ACE2_TMPRSS2_percentage_of_expressing_cells_coexpressing_subclasses_plot_20200521_ggtext_figure6_panelG.pdf", # with ggtext
               is.ggassemble = TRUE,
               panel.width = 7,
               panel.height = 1.5,
               unit.use = "cm",
               margin = 1,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Save percentage.data per subclass as .csv table for paper
out.percentage.subclass  <- percentage.data
colnames(x = out.percentage.subclass) <- gsub(pattern = "\\.",
                                              replacement = "_",
                                              x = colnames(x = out.percentage.subclass))
colnames(x = out.percentage.subclass) <- gsub(pattern = "label",
                                              replacement = "name",
                                              x = colnames(x = out.percentage.subclass))
colnames(x = out.percentage.subclass) <- gsub(pattern = "cluster",
                                              replacement = "subclass",
                                              x = colnames(x = out.percentage.subclass))
colnames(x = out.percentage.subclass) <- gsub(pattern = "ACE2_TMPRSS2",
                                              replacement = "ACE2/TMPRSS2",
                                              x = colnames(x = out.percentage.subclass))
out.percentage.subclass <- out.percentage.subclass[order(out.percentage.subclass$subclass_name),]

write.csv(x = out.percentage.subclass,
          file = "Human_Multiple_Cortical_Areas_SMART-seq_ACE2_TMPRSS2_expressing_cells_percentage_per_subclass.csv",
          quote = TRUE,
          row.names = FALSE)


# Compute the number and percentage of cells expressing ACE2 and TMPRSS2 in all clusters
percentage.data.c <- as.data.frame(x = as.matrix(x = t(x = data.use.tsne.coexpr[c("ACE2", "TMPRSS2", "ACE2/TMPRSS2"),])))
percentage.data.c$sample_name <- rownames(x = percentage.data.c)
percentage.data.c <- merge(x = percentage.data.c,
                           y = meta.data.sub,
                           by = "sample_name",
                           sort = FALSE)
percentage.data.c.list <- split(x = percentage.data.c,
                                f = percentage.data.c$cluster_label)
percentage.data.c.list <- lapply(X = percentage.data.c.list,
                                 FUN = function(percentage.data.sub) {
                                   out.data <- data.frame(cluster_label = unique(x = percentage.data.sub$cluster_label),
                                                          cluster.size = nrow(x = percentage.data.sub),
                                                          ACE2.number = sum(percentage.data.sub$ACE2 > 0),
                                                          ACE2.percentage = sum(percentage.data.sub$ACE2 > 0) / nrow(x = percentage.data.sub) * 100,
                                                          TMPRSS2.number = sum(percentage.data.sub$TMPRSS2 > 0),
                                                          TMPRSS2.percentage = sum(percentage.data.sub$TMPRSS2 > 0) / nrow(x = percentage.data.sub) * 100,
                                                          ACE2.TMPRSS2.number = sum(percentage.data.sub$`ACE2/TMPRSS2` > 0),
                                                          ACE2.TMPRSS2.percentage = sum(percentage.data.sub$`ACE2/TMPRSS2` > 0) / nrow(x = percentage.data.sub) * 100)
                                   return(out.data)
                                 })
percentage.data.c <- dplyr::bind_rows(percentage.data.c.list,
                                      .id = NULL)
percentage.data.c$coexpressing <- percentage.data.c$ACE2.TMPRSS2.number > 0
c.s.corresp <- meta.data.sub[!duplicated(x = meta.data.sub$cluster_label), c("cluster_label", "subclass_label")]
percentage.data.c <- merge(x = percentage.data.c,
                           y = c.s.corresp,
                           by = "cluster_label",
                           sort = FALSE)
cols <- c(grep(pattern = "\\_",
               x = colnames(x = percentage.data.c),
               value = TRUE),
          grep(pattern = "\\.|coexpressing",
               x = colnames(x = percentage.data.c),
               value = TRUE))
percentage.data.c <- percentage.data.c[,cols]
percentage.data.c$subclass_label <- factor(x = percentage.data.c$subclass_label,
                                           levels = cluster.breaks)
percentage.data.c <- percentage.data.c[order(percentage.data.c$subclass_label, percentage.data.c$cluster_label),]


# Save percentage.data.c per cluster as .csv table for paper
colnames(x = percentage.data.c) <- gsub(pattern = "\\.",
                                        replacement = "_",
                                        x = colnames(x = percentage.data.c))
colnames(x = percentage.data.c) <- gsub(pattern = "label",
                                        replacement = "name",
                                        x = colnames(x = percentage.data.c))
colnames(x = percentage.data.c) <- gsub(pattern = "ACE2_TMPRSS2",
                                        replacement = "ACE2/TMPRSS2",
                                        x = colnames(x = percentage.data.c))

write.csv(x = percentage.data.c,
          file = "Human_Multiple_Cortical_Areas_SMART-seq_ACE2_TMPRSS2_expressing_cells_percentage_per_cluster.csv",
          quote = TRUE,
          row.names = FALSE)







table.c.s <- as.data.frame(x = table(cluster_label = meta.data$cluster_label, 
                                     subclass_label = meta.data$subclass_label))
table.c.s <- table.c.s[table.c.s$Freq != 0,]

write.csv(x = table.c.s,
          file = "Human_Multiple_Cortical_Areas_SMART-seq_cluster_subclass_correspondence.csv",
          quote = FALSE,
          row.names = FALSE)
