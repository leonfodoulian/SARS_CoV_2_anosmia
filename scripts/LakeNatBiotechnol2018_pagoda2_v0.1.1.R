source("/Users/leonfodoulian/scData/PlotExpressionOntSNE.R")
source("/Users/leonfodoulian/scData/PlotClustersOntSNE.R")
source("/Users/leonfodoulian/scData/FixSizeAndSave.R")

require(Matrix)
require(pagoda2) # from github
require(uwot)
require(ggplot2)
require(RColorBrewer)
require(extrafont)

setwd(dir = "/Users/leonfodoulian/scData/SARS_CoV_2_anosmia/LakeNatBiotechnol2018/")

if (!dir.exists(paths = "./Figures_for_paper")) {
  dir.create(path = "./Figures_for_paper")
}

set.seed(seed = 315)

# Read meta.data table from LakeNatBiotechnol2018-ST2
meta.data <- read.csv(file = "LakeNatBiotechnol2018-ST2.csv",
                      header = TRUE,
                      check.names = FALSE,
                      stringsAsFactors = FALSE)
colnames(x = meta.data) <- gsub(pattern = " ",
                                replacement = "_",
                                x = colnames(x = meta.data))
meta.data$cell.names <- paste(meta.data$Identity,
                              meta.data$`Sample_Names_(Library_Barcode)`,
                              sep = "_")

# Rename cell clusters with their full names
cluster.labels <- data.frame(Identity = sort(x = unique(x = meta.data$Identity)),
                             stringsAsFactors = FALSE)
cluster.labels$cluster.labels <- ifelse(test = grepl(pattern = "^Ex",
                                                     x = cluster.labels$Identity),
                                        yes = paste("Excitatory neurons",
                                                    cluster.labels$Identity,
                                                    sep = " "),
                                        no = cluster.labels$Identity)
cluster.labels$cluster.labels <- ifelse(test = grepl(pattern = "^In",
                                                     x = cluster.labels$Identity),
                                        yes = paste("Inhibitory neurons",
                                                    cluster.labels$Identity,
                                                    sep = " "),
                                        no = cluster.labels$cluster.labels)
cluster.labels$cluster.labels <- ifelse(test = grepl(pattern = "^Purk",
                                                     x = cluster.labels$Identity),
                                        yes = paste("Purkinje neurons",
                                                    cluster.labels$Identity,
                                                    sep = " "),
                                        no = cluster.labels$cluster.labels)
cluster.labels$cluster.labels[cluster.labels$Identity == "Ast"] <- "Astrocytes"
cluster.labels$cluster.labels[cluster.labels$Identity == "Ast_Cer"] <- "Cerebellar astrocytes"
cluster.labels$cluster.labels[cluster.labels$Identity == "End"] <- "Endothelial cells"
cluster.labels$cluster.labels[cluster.labels$Identity == "Gran"] <- "Cerebellar granule cells"
cluster.labels$cluster.labels[cluster.labels$Identity == "Mic"] <- "Microglia"
cluster.labels$cluster.labels[cluster.labels$Identity == "Oli"] <- "Oligodendrocytes"
cluster.labels$cluster.labels[cluster.labels$Identity == "OPC"] <- "Oligodendrocyte precursor cells"
cluster.labels$cluster.labels[cluster.labels$Identity == "OPC_Cer"] <- "Cerebellar oligodendrocyte precursor cells"
cluster.labels$cluster.labels[cluster.labels$Identity == "Per"] <- "Pericytes"

# Merge cluster.labels to meta.data table from LakeNatBiotechnol2018-ST2
meta.data <- merge(x = meta.data,
                   y = cluster.labels,
                   by = "Identity",
                   sort = FALSE)

# Define broader cluster definitions
meta.data$broad.labels <- ifelse(test = grepl(pattern = "^Ex",
                                              x = meta.data$cluster.labels),
                                 yes = "Excitatory neurons",
                                 no = meta.data$cluster.labels)
meta.data$broad.labels <- ifelse(test = grepl(pattern = "^In",
                                              x = meta.data$cluster.labels),
                                 yes = "Inhibitory neurons",
                                 no = meta.data$broad.labels)
meta.data$broad.labels <- ifelse(test = grepl(pattern = "^Purk",
                                              x = meta.data$cluster.labels),
                                 yes = "Purkinje neurons",
                                 no = meta.data$broad.labels)

# Get names of count matrix files
files <- list.files(path = "GSE97930/",
                    pattern = "Count_Matrix",
                    full.names = TRUE)
names(x = files) <- gsub(pattern = ".*GSE97930_|_snDrop-seq.*",
                         replacement = "",
                         x = files)

# Read count matrix files and store in a list
counts.list <- lapply(X = files,
                      FUN = read.table,
                      header = TRUE,
                      row.names = 1)

# Transform count matrices to sparse matrix
counts.list <- lapply(X = counts.list,
                      FUN = function(count.matrix) {
                        if (!inherits(x = count.matrix, what = "dgCMatrix")) {
                          count.matrix <- as(object = as.matrix(x = count.matrix),
                                             Class = "dgCMatrix")
                        }
                        return(count.matrix)
                      })

# Define for each cell its brain region of origin
brain.regions <- dplyr::bind_rows(mapply(count.matrix = counts.list,
                                         brain.region = names(x = counts.list),
                                         FUN = function(count.matrix,
                                                        brain.region) {
                                           data.frame(cell.names = colnames(x = count.matrix),
                                                      brain.region = brain.region,
                                                      stringsAsFactors = FALSE)
                                         },
                                         SIMPLIFY = FALSE,
                                         USE.NAMES = TRUE),
                                  .id = NULL)

# Merge brain region information to meta.data table from LakeNatBiotechnol2018-ST2
meta.data <- merge(x = meta.data,
                   y = brain.regions,
                   by = "cell.names",
                   sort = FALSE)

# Extract common genes in all 3 datasets
common.genes <- sort(x = Reduce(f = function(x,y) { intersect(x = x,
                                                              y = y) },
                                lapply(X = counts.list,
                                       FUN = rownames)))

# Subset count matrices and keep only common genes
counts.list.sub <- lapply(X = counts.list,
                          FUN = function(count.matrix,
                                         common.genes) {
                            count.matrix <- count.matrix[common.genes,]
                            return(count.matrix)
                          },
                          common.genes = common.genes)

# Create a single matrix for subsequent analysis using Seurat
counts.matrix <- Reduce(f = function(x, y) { Matrix::cbind2(x = x,
                                                            y = y)},
                        x = counts.list.sub)

# Normalize counts.matrix by library size
norm.counts.matrix <- t(x = t(x = counts.matrix) / colSums(x = counts.matrix)) * 1e4

rm(counts.list) # save memory on the laptop
rm(counts.list.sub) # save memory on the laptop

# Create a vector of identities for batch correction
# batch <- brain.regions$brain.region
# names(x = batch) <- brain.regions$cell.names
batch <- unlist(x = lapply(X = strsplit(x = colnames(x = counts.matrix),
                                        split = "_"),
                           FUN = function(split.cols) {
                             split.cols[(length(x = split.cols) - 1)]
                           }))
names(x = batch) <- colnames(x = counts.matrix)

if (!all(names(x = batch) == colnames(x = counts.matrix))) {
  cat("Cells order does not match between 'batch' and 'counts.matrix': reordering 'batch'")
  batch <- batch[colnames(x = counts.matrix)]
}

# Create Pagoda2 object and correct for batch
obj <- Pagoda2$new(x = counts.matrix,
                   n.cores = 4,
                   trim = 10,
                   batch = batch)

# Adjust the variance of genes
obj$adjustVariance(gam.k = 10,
                   alpha = 0.05,
                   plot = TRUE,
                   use.raw.variance = (obj$modelType == "raw"),
                   use.unadjusted.pvals = FALSE,
                   do.par = TRUE,
                   max.adjusted.variance = 1000,
                   min.adjusted.variance = 0.001,
                   cells = NULL,
                   verbose = TRUE,
                   min.gene.cells = 0,
                   persist = TRUE,
                   n.cores = obj$n.cores)

# Calculate the first 150 PCs using 2000 overdispersed genes
obj$calculatePcaReduction(nPcs = 150,
                          type = "counts",
                          name = "PCA",
                          use.odgenes = TRUE,
                          n.odgenes = 2.e3,
                          odgenes = NULL,
                          center = TRUE,
                          cells = NULL,
                          fastpath = TRUE,
                          maxit = 10000,
                          verbose = TRUE,
                          var.scale = TRUE)

# Generate an embedding with UMAP on the basis of the PCA reduction
obj$getEmbedding(type = "PCA",
                 embeddingType = "UMAP",
                 name = NULL,
                 dims = 2,
                 M = 1,
                 gamma = 1/1,
                 perplexity = 30,
                 sgd_batches = NULL,
                 diffusion.steps = 0,
                 diffusion.power = 0.5,
                 # distance = "pearson",
                 distance = "cosine",
                 # n.cores = obj$n.cores,
                 n.cores = 1,
                 # n.sgd.cores = obj$n.cores,
                 n.sgd.cores = 1,
                 n_neighbors = 15, # uwot::umap() argument
                 min_dist = 0.1) # uwot::umap() argument

colnames(x = obj$embeddings$PCA$UMAP) <- c("UMAP_1", "UMAP_2")

rds.file <- list(obj = obj,
                 norm.counts.matrix = norm.counts.matrix,
                 meta.data = meta.data,
                 batch = batch)

saveRDS(object = rds.file,
        file = "LakeNatBiotechnol2018_2000odg_150pc_pagoda2_v0.1.1.rds")
rds.file <- readRDS(file = "LakeNatBiotechnol2018_2000odg_150pc_pagoda2_v0.1.1.rds")

obj <- rds.file$obj
norm.counts.matrix <- rds.file$norm.counts.matrix
meta.data <- rds.file$meta.data

# Figure X, Panel A
# Plot Clusters on UMAP

all.clusters <- sort(x = unique(x = meta.data$cluster.labels))

excitatory.clusters <- grepl(pattern = "Excitatory",
                             x = all.clusters)
excitatory.colors <- colorRampPalette(colors = c("royalblue4", "deepskyblue1"))(sum(excitatory.clusters))
names(x = excitatory.colors) <- all.clusters[excitatory.clusters]

inhibitory.clusters <- grepl(pattern = "Inhibitory",
                             x = all.clusters)
inhibitory.colors <- colorRampPalette(colors = c("deeppink4", "hotpink1"))(sum(inhibitory.clusters))
names(x = inhibitory.colors) <- all.clusters[inhibitory.clusters]

purkinje.colors <- c("Purkinje neurons Purk1" = "darkgreen",
                     "Purkinje neurons Purk2" = "limegreen")

granule.colors <- c("Cerebellar granule cells" = "darkcyan")

oligo.colors <- c("Oligodendrocytes" = "#B15928",
                  "Oligodendrocyte precursor cells" = "#D95F02",
                  "Cerebellar oligodendrocyte precursor cells" = "#E6AB02")

astrocyte.colors <- c("Astrocytes" = "mediumorchid1",
                      "Cerebellar astrocytes" = "mediumpurple1")

microglia.colors <- c("Microglia" = "#CAB2D6")

blood.cell.colors <- c("Endothelial cells" = "#666666",
                       "Pericytes" = "#B3B3B3")

all.colors <- c(excitatory.colors,
                inhibitory.colors,
                purkinje.colors,
                granule.colors,
                oligo.colors,
                astrocyte.colors,
                microglia.colors,
                blood.cell.colors)

broad.clusters <- c(sort(x = grep(pattern = "Excitatory",
                                  x = names(x = all.colors),
                                  value = TRUE))[1],
                    sort(x = grep(pattern = "Inhibitory",
                                  x = names(x = all.colors),
                                  value = TRUE))[1],
                    sort(x = grep(pattern = "Purkinje",
                                  x = names(x = all.colors),
                                  value = TRUE))[1],
                    names(x = all.colors)[!grepl(pattern = "Excitatory|Inhibitory|Purkinje",
                                                 x = names(x = all.colors))])
broad.colors <- all.colors[broad.clusters]
names(x = broad.colors) <- gsub(pattern = "\\sEx{0-9}.*|\\sIn1{0-9}.*|\\sPurk{0-9}.*",
                                replacement = "",
                                x = names(x = broad.colors))

cluster.breaks <- names(x = all.colors)
cluster.labels  <- cluster.breaks
names(x = cluster.labels) <- cluster.labels

cluster.data <- data.frame(Cluster = meta.data$cluster.labels,
                           row.names = meta.data$cell.names)
cluster.data <- cluster.data[rownames(x = obj$embeddings$PCA$UMAP),,drop = FALSE]

plot.cluster.umap <- PlotClustersOntSNE(
  tsne.data = obj$embeddings$PCA$UMAP,
  cluster.data = cluster.data,
  cellnames.as.rownames.in.tsne.data = TRUE,
  cellnames.as.rownames.in.cluster.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  cluster.ident.name = "Cluster",
  cluster.colors = all.colors,
  cluster.breaks = cluster.breaks,
  cluster.labels = cluster.labels,
  add.border.to.points = FALSE,
  point.size = 0.5,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  alpha.use = 0.5,
  font.family = "Arial",
  import.font = FALSE,
  fix.aspect.ratio = FALSE,
  legend.point.size = 2,
  axis.title.size = 6,
  legend.title.size = 7,
  legend.text.size = 6,
  legend.text.space = unit(x = 3, units = "mm"),
  axis.line.size = rel(0.5),
  arrow.size = rel(2),
  range.scale = 0.285
)

plot.cluster.umap <- plot.cluster.umap +
  labs(x = "UMAP 1 (a.u.)",
       y = "UMAP 2 (a.u.)",
       colour = "Cluster IDs (# cells)") +
  theme(legend.title = element_text(colour = "black",
                                    size = 7,
                                    family = "Arial",
                                    face = "plain"))

FixSizeAndSave(plot = plot.cluster.umap,
               # filename = "./Figures_for_paper/LakeNatBiotechnol2018_UMAP_clusters_plot_20200429_figureX_panelA.pdf",
               filename = "./Figures_for_paper/LakeNatBiotechnol2018_UMAP_clusters_plot_20200429_ggtext_suppfigure6_panelA.pdf",
               is.ggassemble = FALSE,
               panel.width = 5,
               panel.height = 5,
               unit.use = "cm",
               margin = 1,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure X, Panel B
# Plot ACE2 and TMPRSS2 marker genes on UMAP

gradient.color <- c(viridis::cividis(2)[2], rev(x = viridis::magma(11))[-c(1:2)])

data.use.umap.coexpr <- norm.counts.matrix[c("ACE2", "TMPRSS2"),]

ACE2.TMPRSS2 <- Matrix::Matrix(data = ifelse(test = data.use.umap.coexpr["ACE2",, drop = FALSE] * data.use.umap.coexpr["TMPRSS2",, drop = FALSE] > 0,
                                             yes = 1,
                                             no = 0),
                               nrow = 1,
                               dimnames = list("ACE2/TMPRSS2", colnames(x = data.use.umap.coexpr)))

data.use.umap.coexpr <- rbind(data.use.umap.coexpr, ACE2.TMPRSS2)

plot.markers.expression <- PlotExpressionOntSNE(
  data.use = data.use.umap.coexpr,
  genes.use = rownames(x = data.use.umap.coexpr)[1:2],
  ggtext.gene.names = c("ACE2" = "<span style='color:violetred1'>ACE2</span>",
                        "TMPRSS2" = "<span style='color:lightslateblue'>TMPRSS2</span>"), # with ggtext
  tsne.data = obj$embeddings$PCA$UMAP,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
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
               # filename = "./Figures_for_paper/LakeNatBiotechnol2018_UMAP_ACE2_TMPRSS2_expression_plot_20200429_figureX_panelB.pdf", # without ggtext
               filename = "./Figures_for_paper/LakeNatBiotechnol2018_UMAP_ACE2_TMPRSS2_expression_plot_20200429_ggtext_suppfigure6_panelB.pdf", # with ggtext
               is.ggassemble = TRUE,
               panel.width = 4.8,
               panel.height = 4.8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure X, Panel C
# Plot coexpression of ACE2 with TMPRSS2 on UMAP

plot.markers.coexpression <- PlotExpressionOntSNE(
  data.use = data.use.umap.coexpr,
  genes.use = rownames(x = data.use.umap.coexpr)[c(3,3)], # size do not match between cowplot and patchwork; trick to match size
  ggtext.gene.names = c("ACE2/TMPRSS2" = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>",
                        "ACE2/TMPRSS2" = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>"), # with ggtext
  tsne.data = obj$embeddings$PCA$UMAP,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
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
               # filename = "./Figures_for_paper/LakeNatBiotechnol2018_UMAP_ACE2_TMPRSS2_coexpressing_cells_plot_20200429_figureX_panelC.pdf", # without ggtext
               filename = "./Figures_for_paper/LakeNatBiotechnol2018_UMAP_ACE2_TMPRSS2_coexpressing_cells_plot_20200429_ggtext_suppfigure6_panelC.pdf", # with ggtext
               is.ggassemble = FALSE, # trick to match size with previous plots
               panel.width = 4.8,
               panel.height = 4.8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure X, Panel D
# Plot the number of cells expressing ACE2 and TMPRSS2

percentage.data <- as.data.frame(x = as.matrix(x = t(x = data.use.umap.coexpr)))
percentage.data$cell.names <- rownames(x = percentage.data)
percentage.data <- merge(x = percentage.data,
                         y = meta.data,
                         by = "cell.names",
                         sort = FALSE)
percentage.data.list <- split(x = percentage.data,
                              f = percentage.data$broad.labels)
percentage.data.list <- lapply(X = percentage.data.list,
                               FUN = function(percentage.data.sub) {
                                 out.data <- data.frame(broad.labels = unique(x = percentage.data.sub$broad.labels),
                                                        cluster.size = nrow(x = percentage.data.sub),
                                                        ACE2.number = sum(percentage.data.sub$ACE2 > 0),
                                                        ACE2.percentage = sum(percentage.data.sub$ACE2 > 0) / nrow(x = percentage.data.sub) * 100,
                                                        TMPRSS2.number = sum(percentage.data.sub$TMPRSS2 > 0),
                                                        TMPRSS2.percentage = sum(percentage.data.sub$TMPRSS2 > 0) / nrow(x = percentage.data.sub) * 100)
                                 return(out.data)
                               })
percentage.data <- dplyr::bind_rows(percentage.data.list,
                                    .id = NULL)
percentage.data$broad.labels <- factor(x = percentage.data$broad.labels,
                                       levels = names(x = broad.colors))

max.y.ace2.number <- max(percentage.data$ACE2.number)
breaks.ace2.number <- seq(from = 0,
                          to = max.y.ace2.number,
                          by = 2)
labels.ace2.number <- ifelse(test = breaks.ace2.number %in% c(0, max.y.ace2.number),
                             yes = breaks.ace2.number,
                             no = "")

max.y.tmprss2.number <- max(percentage.data$TMPRSS2.number)
breaks.tmprss2.number <- seq(from = 0,
                             to = max.y.tmprss2.number,
                             by = 5)
labels.tmprss2.number <- ifelse(test = breaks.tmprss2.number %in% c(0, max.y.tmprss2.number),
                                yes = breaks.tmprss2.number,
                                no = "")

plot.ace2.number <- ggplot(data = percentage.data,
                           mapping = aes(x = broad.labels,
                                         y = ACE2.number)) +
  geom_bar(stat = "identity",
           colour = NA,
           fill = "violetred1",
           width = 0.5) +
  scale_y_continuous(breaks = breaks.ace2.number,
                     labels = labels.ace2.number,
                     expand = c(0,0)) +
  expand_limits(y = c(0, max.y.ace2.number)) +
  labs(title = "ACE2",
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
                              mapping = aes(x = broad.labels,
                                            y = TMPRSS2.number)) +
  geom_bar(stat = "identity",
           colour = NA,
           fill = "lightslateblue",
           width = 0.5) +
  scale_y_continuous(breaks = breaks.tmprss2.number,
                     labels = labels.tmprss2.number,
                     expand = c(0,0)) +
  expand_limits(y = c(0, max.y.tmprss2.number)) +
  labs(title = "TMPRSS2") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = broad.colors[levels(x = percentage.data$broad.labels)],
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

plot.comb.number <- patchwork::wrap_plots(plot.ace2.number,
                                          plot.tmprss2.number,
                                          ncol = 1)

FixSizeAndSave(plot = plot.comb.number,
               # filename = "./Figures_for_paper/LakeNatBiotechnol2018_ACE2_TMPRSS2_number_of_expressing_cells_broad_clusters_plot_20200430_figureX_panelD.pdf",
               filename = "./Figures_for_paper/LakeNatBiotechnol2018_ACE2_TMPRSS2_number_of_expressing_cells_broad_clusters_plot_20200430_ggtext_suppfigure6_panelD.pdf",
               is.ggassemble = TRUE,
               panel.width = 6,
               panel.height = 1.5,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure X, Panel E
# Plot the percentage of cells expressing ACE2 and TMPRSS2

max.y.ace2.percentage <- ceiling(x = max(percentage.data$ACE2.percentage) * 10) / 10
breaks.ace2.percentage <- seq(from = 0,
                              to = max.y.ace2.percentage,
                              by = 0.1)
labels.ace2.percentage <- ifelse(test = breaks.ace2.percentage %in% c(0, max.y.ace2.percentage),
                                 yes = breaks.ace2.percentage,
                                 no = "")

max.y.tmprss2.percentage <- ceiling(x = max(percentage.data$TMPRSS2.percentage) * 10) / 10
breaks.tmprss2.percentage <- seq(from = 0,
                                 to = max.y.tmprss2.percentage,
                                 by = 0.1)
labels.tmprss2.percentage <- ifelse(test = breaks.tmprss2.percentage %in% c(0, max.y.tmprss2.percentage),
                                    yes = breaks.tmprss2.percentage,
                                    no = "")

plot.ace2.percentage <- ggplot(data = percentage.data,
                               mapping = aes(x = broad.labels,
                                             y = ACE2.percentage)) +
  geom_bar(stat = "identity",
           colour = NA,
           fill = "violetred1",
           width = 0.5) +
  scale_y_continuous(breaks = breaks.ace2.percentage,
                     labels = labels.ace2.percentage,
                     expand = c(0,0)) +
  expand_limits(y = c(0, max.y.ace2.percentage)) +
  labs(title = "ACE2",
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

plot.tmprss2.percentage <- ggplot(data = percentage.data,
                                  mapping = aes(x = broad.labels,
                                                y = TMPRSS2.percentage)) +
  geom_bar(stat = "identity",
           colour = NA,
           fill = "lightslateblue",
           width = 0.5) +
  scale_y_continuous(breaks = breaks.tmprss2.percentage,
                     labels = labels.tmprss2.percentage,
                     expand = c(0,0)) +
  expand_limits(y = c(0, max.y.tmprss2.percentage)) +
  labs(title = "TMPRSS2") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = broad.colors[levels(x = percentage.data$broad.labels)],
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

plot.comb.percentage <- patchwork::wrap_plots(plot.ace2.percentage,
                                              plot.tmprss2.percentage,
                                              ncol = 1)

FixSizeAndSave(plot = plot.comb.percentage,
               # filename = "./Figures_for_paper/LakeNatBiotechnol2018_ACE2_TMPRSS2_percentage_of_expressing_cells_broad_clusters_plot_20200430_figureX_panelE.pdf",
               filename = "./Figures_for_paper/LakeNatBiotechnol2018_ACE2_TMPRSS2_percentage_of_expressing_cells_broad_clusters_plot_20200430_ggtext_suppfigure6_panelE.pdf",
               is.ggassemble = TRUE,
               panel.width = 6,
               panel.height = 1.5,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Save percentage.data as .csv table for paper
out.percentage  <- percentage.data
out.percentage$ACE2.TMPRSS2.number <- sum(data.use.umap.coexpr["ACE2/TMPRSS2",])
out.percentage$ACE2.TMPRSS2.percentage <- sum(data.use.umap.coexpr["ACE2/TMPRSS2",])
out.percentage$coexpressing <- out.percentage$ACE2.TMPRSS2.number > 0
colnames(x = out.percentage) <- gsub(pattern = "\\.",
                                     replacement = "_",
                                     x = colnames(x = out.percentage))
colnames(x = out.percentage) <- gsub(pattern = "broad_labels",
                                     replacement = "broad_cluster_name",
                                     x = colnames(x = out.percentage))
colnames(x = out.percentage) <- gsub(pattern = "cluster_size",
                                     replacement = "broad_cluster_size",
                                     x = colnames(x = out.percentage))
colnames(x = out.percentage) <- gsub(pattern = "ACE2_TMPRSS2",
                                     replacement = "ACE2/TMPRSS2",
                                     x = colnames(x = out.percentage))
out.percentage <- out.percentage[order(out.percentage$broad_cluster_name),]

write.csv(x = out.percentage,
          file = "LakeNatBiotechnol2018_ACE2_TMPRSS2_expressing_cells_percentage_per_broad_cluster.csv",
          quote = TRUE,
          row.names = FALSE)
