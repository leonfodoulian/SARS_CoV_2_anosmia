source("/Users/leonfodoulian/scData/PlotExpressionOntSNE.R")
source("/Users/leonfodoulian/scData/PlotClustersOntSNE.R")
source("/Users/leonfodoulian/scData/FixSizeAndSave.R")

require(Matrix)
# require(Seurat, lib.loc = "/Library/Frameworks/R.framework/Resources/library/Seurat_v3.0.0/")
require(Seurat, lib.loc = "/Library/Frameworks/R.framework/Resources/library/Seurat_v3.1.4/")
packageVersion(pkg = "Seurat")
require(ggplot2)
require(RColorBrewer)
require(extrafont)
require(ggtext)

setwd(dir = "/Users/leonfodoulian/scData/SARS_CoV_2_anosmia/DuranteNatNeurosci2020/")

if (!dir.exists(paths = "./Figures_for_paper")) {
  dir.create(path = "./Figures_for_paper")
}

## Read meta.data table from Durante
## Rename 'Respiratory Columnar Cells' to 'Respiratory Early Secretory Cells'
## and reduce all cell type words except the first one to lower case
meta.data.table <- read.csv(file = "20190720.OE.4samp.standard.integration.FINAL.meta.data.patient.csv",
                            header = TRUE,
                            row.names = 1,
                            stringsAsFactors = FALSE)
colnames(x = meta.data.table) <- paste(colnames(x = meta.data.table),
                                       "Durante",
                                       sep = "_")
meta.data.table$cell.names_Durante <- rownames(x = meta.data.table)
meta.data.table$orig.ident <- gsub(pattern = "P",
                                   replacement = "Patient",
                                   x = meta.data.table$orig.ident_Durante)
meta.data.table$cell.names <- gsub(pattern = "\\-.*",
                                   replacement = "",
                                   x = meta.data.table$cell.names_Durante)
meta.data.table$cell.names <- paste(meta.data.table$orig.ident,
                                    meta.data.table$cell.names,
                                    sep = "_")
rownames(x = meta.data.table) <- meta.data.table$cell.names
meta.data.table$CellType <- ifelse(test = meta.data.table$CellType_Durante %in% "Respiratory Columnar Cells",
                                   yes = "Respiratory Early Secretory Cells",
                                   no = meta.data.table$CellType_Durante)
meta.data.table$CellType <- unlist(x = lapply(X = strsplit(x = meta.data.table$CellType,
                                                           split = "\\s"),
                                              FUN = function(full.string) {
                                                sub.string <- ifelse(test = grepl(pattern = "\\+",
                                                                                  x = full.string[-1]),
                                                                     yes = full.string[-1],
                                                                     no = tolower(x = full.string[-1]))
                                                new.string <- paste(c(full.string[1],
                                                                      sub.string),
                                                                    collapse = " ")
                                                return(new.string)
                                              }))
meta.data.table$CellType <- factor(x = meta.data.table$CellType,
                                   levels = sort(x = unique(x = meta.data.table$CellType)))

# Get names of tar files and extract them
tar.files <- list.files(path = ".",
                        pattern = ".tar$",
                        full.names = FALSE,
                        recursive = FALSE)
# Define directories where files should be extracted
# One directory per Patient (cleaner and easy to trace back)
ex.dirs <- gsub(pattern = ".tar$",
                replacement = "",
                x = tar.files)
invisible(x = mapply(tarfile = tar.files,
                     exdir = ex.dirs,
                     FUN = untar))
# List all extracted files
old.file.names <- list.files(path = ".",
                             pattern = ".gz$",
                             recursive = TRUE)
# Compute checksum on all extracted files
checksum.old <- tools::md5sum(files = old.file.names)

# Seurat::Read10X expects 3 files:
# 1: barcodes.tsv.gz
# 2: features.tsv.gz
# 3: matrix.mtx.gz
# Rename all files as expected so that Seurat is able to read them
new.file.names <- gsub(pattern = "/GSM(\\d+)_Patient(\\d+)_",
                       replacement = "/",
                       x = old.file.names)
file.rename(from = old.file.names,
            to = new.file.names)
# Compute checksum on renamed files
checksum.new <- tools::md5sum(files = new.file.names)
# Check if all checksums match
all(checksum.old == checksum.new)

# Get names of directories containing 10X matrices
# There are 4 directories, one for each Patient
dirs <- list.dirs(path = ".",
                  full.names = FALSE,
                  recursive = FALSE)
dirs <- dirs[grepl(pattern = "Patient{1-9}", x = dirs)]
# Give Patient IDs as names to directories so that when Seurat::Read10X 
# reads the files it prefixes the cell barcode names with the Patient ID
# This will facilitate to trace each cell back to its Patient of origin
names(x = dirs) <- gsub(pattern = ".*\\_",
                        replacement = "",
                        x = dirs)
# Transform 'dirs' to list and add a name to the vector in the list because 
# otherwise lapply omits the vector names when passing them to Seurat::Read10X
dirs <- as.list(x = dirs)
dirs <- mapply(file = dirs,
               file.name = names(x = dirs),
               FUN = function(file,
                              file.name) {
                 names(x = file) <- file.name
                 return(file)
               },
               SIMPLIFY = FALSE,
               USE.NAMES = TRUE)

# Read all 10X files using Seurat::Read10X
matrix.list <- lapply(X = dirs,
                      FUN = Read10X,
                      gene.column = 2,
                      unique.features = TRUE)

## Create a single matrix for subsequent analysis using Seurat
matrix <- Reduce(f = function(x, y) { Matrix::cbind2(x = x,
                                                     y = y)},
                 x = matrix.list)

# # Create Seurat objects without any filters
# obj.list <- lapply(X = matrix.list,
#                    FUN = CreateSeuratObject,
#                    project = "HumanMOE",
#                    assay = "RNA",
#                    min.cells = 0,
#                    min.features = 0,
#                    names.field = 1,
#                    names.delim = "_",
#                    meta.data = NULL)
## Create Seurat object by removing genes not expressed in at least 3 cells
obj.raw <- CreateSeuratObject(counts = matrix,
                              project = "HumanMOE",
                              assay = "RNA",
                              min.cells = 3,
                              min.features = 0,
                              names.field = 1,
                              names.delim = "_",
                              meta.data = NULL)

rm(matrix.list) # save memory on the laptop
rm(matrix) ## save memory on the laptop

# # Add mitochondrial count percentages to each Seurat object
# obj.list <- lapply(X = obj.list,
#                    FUN = function(seurat.object) {
#                      seurat.object[["percent.mt"]] <- PercentageFeatureSet(object = seurat.object,
#                                                                            pattern = "^MT-",
#                                                                            assay = "RNA")
#                      return(seurat.object)
#                    })
## Add mitochondrial count percentages to Seurat object
obj.raw[["percent.mt"]] <- PercentageFeatureSet(object = obj.raw,
                                                pattern = "^MT-",
                                                features = NULL,
                                                col.name = NULL,
                                                assay = "RNA")

# # Get initial number of cells pre filtering
# n.cells.pre.filtering <- sum(unlist(x = lapply(X = obj.list,
#                                                FUN = ncol)))
## Get initial number of cells pre filtering
n.cells.pre.filtering <- ncol(x = obj.raw)

# # Filter cells using the criteria described in DuranteNatNeurosci2020
# obj.list.filtered <- lapply(X = obj.list,
#                             FUN = function(seurat.object) {
#                               seurat.object <- subset(x = seurat.object,
#                                                       subset = nCount_RNA > 400 & nFeature_RNA >= 100  & nFeature_RNA <= 8000 & percent.mt < 10)
#                               return(seurat.object)
#                             })
## Filter cells using the criteria described in DuranteNatNeurosci2020
obj.filtered <- subset(x = obj.raw,
                       subset = nCount_RNA > 400 & nCount_RNA < Inf & nFeature_RNA > 100  & nFeature_RNA < 6000 & percent.mt > -Inf & percent.mt < 10)

# # Get number of cells post filtering
# # This filtering (which is supposed to be identical) results in 29619 cells, 
# # almost a 1000 more than what was reported in the paper (28726)
# # Note that nFeature_RNA >= 300 results in 28721 cells (close to what they reported)
# (n.cells.post.filtering <- sum(unlist(x = lapply(X = obj.list.filtered,
#                                                  FUN = ncol))))
## Get number of cells post filtering
## This filtering matches the one reported in the paper
(n.cells.post.filtering <- ncol(x = obj.filtered))

# rm(obj.list) # save memory on the laptop
rm(obj.raw) # save memory on the laptop

## Split Seurat object into a list of objects, one per patient
obj.list.filtered <- SplitObject(object = obj.filtered,
                                 split.by = "orig.ident")

# Normalize data using the standard library normalization and log-transformation method
# Weird that they didn't use the new SCTransform method
obj.list.filtered <- lapply(X = obj.list.filtered,
                            FUN = NormalizeData,
                            assay = "RNA",
                            normalization.method = "LogNormalize",
                            scale.factor = 10000,
                            margin = 1,
                            verbose = TRUE)

# Find variable features using the criteria described in DuranteNatNeurosci2020
obj.list.filtered <- lapply(X = obj.list.filtered,
                            FUN = FindVariableFeatures,
                            assay = "RNA",
                            selection.method = "vst",
                            loess.span = 0.3,
                            clip.max = "auto",
                            # mean.function = FastExpMean, # does not work if defined
                            # dispersion.function = FastLogVMR, # does not work if defined
                            num.bin = 20,
                            binning.method = "equal_width", 
                            nfeatures = 5000,
                            mean.cutoff = c(0.1, 8),
                            dispersion.cutoff = c(1, Inf),
                            verbose = TRUE)

# Find integration anchors to integrate all datasets
obj.anchors <- FindIntegrationAnchors(object.list = obj.list.filtered,
                                      assay = rep(x = "RNA",
                                                  times = length(x = obj.list.filtered)),
                                      reference = NULL, ## v3.1.4
                                      anchor.features = 5000,
                                      scale = TRUE,
                                      normalization.method = "LogNormalize", ## v3.1.4
                                      sct.clip.range = NULL, ## v3.1.4
                                      reduction = "cca", ## v3.1.4
                                      l2.norm = TRUE,
                                      dims = 1:30,
                                      k.anchor = 5,
                                      k.filter = 200,
                                      k.score = 30,
                                      max.features = 200,
                                      nn.method = "rann", ## v3.1.4
                                      eps = 0,
                                      verbose = TRUE)

# Integrate all data into a single Seurat object
obj <- IntegrateData(anchorset = obj.anchors,
                     new.assay.name = "Integrated",
                     normalization.method = "LogNormalize", ## v3.1.4
                     features = NULL,
                     features.to.integrate = NULL,
                     dims = 1:30,
                     k.weight = 100,
                     weight.reduction = NULL,
                     sd.weight = 1,
                     sample.tree = NULL,
                     preserve.order = FALSE,
                     do.cpp = TRUE,
                     eps = 0,
                     verbose = TRUE)

rm(obj.list.filtered) # save memory on the laptop
rm(obj.anchors) # save memory on the laptop

# Set default assay to 'integrated' (previous default assay was 'RNA')
DefaultAssay(object = obj) <- "Integrated"

# Scale data for dimensionality reduction
obj <- ScaleData(object = obj,
                 features = NULL,
                 assay = "Integrated",
                 vars.to.regress = NULL,
                 split.by = NULL, ## v3.1.4
                 model.use = "linear",
                 use.umi = FALSE,
                 do.scale = TRUE,
                 do.center = TRUE,
                 scale.max = 10,
                 block.size = 1000,
                 min.cells.to.block = 3000,
                 verbose = TRUE)

# Compute the first 30 PCs
obj <- RunPCA(object = obj,
              assay = "Integrated",
              features = NULL,
              npcs = 100,
              rev.pca = FALSE,
              weight.by.var = TRUE,
              verbose = TRUE,
              ndims.print = 1:5,
              nfeatures.print = 30,
              reduction.name = "pca",
              reduction.key = "PC_",
              seed.use = 42)

# Compute UMAP
obj <- RunUMAP(object = obj,
               dims = 1:30,
               reduction = "pca",
               features = NULL,
               graph = NULL,
               assay = "Integrated",
               umap.method = "uwot", ## v3.1.4
               n.neighbors = 30L,
               n.components = 2L,
               # metric = "correlation", ## v3.0.0
               metric = "cosine", ## v3.1.4
               n.epochs = NULL,
               learning.rate = 1,
               min.dist = 0.3,
               spread = 1,
               set.op.mix.ratio = 1,
               local.connectivity = 1L,
               repulsion.strength = 1,
               negative.sample.rate = 5L,
               a = NULL,
               b = NULL,
               uwot.sgd = FALSE, ## v3.1.4
               seed.use = 42L,
               metric.kwds = NULL,
               angular.rp.forest = FALSE,
               verbose = TRUE,
               reduction.name = "umap",
               reduction.key = "UMAP_")

DimPlot(object = obj,
        reduction = "umap")

# Compute SNN graph
obj <- FindNeighbors(object = obj,
                     reduction = "pca",
                     dims = 1:30,
                     assay = "Integrated",
                     features = NULL,
                     k.param = 20,
                     compute.SNN = TRUE,
                     prune.SNN = 1/15,
                     nn.method = "rann", ## v3.1.4
                     annoy.metric = "euclidean", ## v3.1.4
                     nn.eps = 0,
                     verbose = TRUE,
                     force.recalc = FALSE,
                     do.plot = FALSE,
                     graph.name = NULL)

# Cluster cells
# They report using a resolution of 1.8 and having 26 clusters
# I get 50 (47 if patients analysed individually) clusters using those parameters
# Using a resolution of 0.4 results in 26 clusters
# but does not match exactly what they have
obj <- FindClusters(object = obj,
                    graph.name = NULL,
                    modularity.fxn = 1,
                    initial.membership = NULL,
                    weights = NULL,
                    node.sizes = NULL,
                    # resolution = 0.4,
                    resolution = 1.8,
                    method = "matrix", ## v3.1.4
                    algorithm = 1,
                    n.start = 10,
                    n.iter = 10,
                    random.seed = 0,
                    group.singletons = TRUE, ## v3.1.4
                    temp.file.location = NULL,
                    edge.file.name = NULL,
                    verbose = TRUE)

DimPlot(object = obj,
        reduction = "umap",
        group.by = "ident")

obj@meta.data$cell.names <- rownames(x = obj@meta.data)
obj@meta.data <- merge(x = obj@meta.data,
                       y = meta.data.table,
                       by = c("cell.names", "orig.ident"),
                       sort = FALSE)
rownames(x = obj@meta.data) <- obj@meta.data$cell.names
if (!all(rownames(x = obj@meta.data) == colnames(x = obj))) {
  cat("Cells order does not match between 'obj@meta.data' and the other slots in 'obj': reordering 'obj@meta.data'")
  obj@meta.data <- obj@meta.data[colnames(x = obj),]
}

table(obj@meta.data$CellType, obj@meta.data$seurat_clusters)

Idents(object = obj) <- "CellType"

# Check that idents were assigned correctly
all(names(x = obj@active.ident) == rownames(x = obj@meta.data))
all(obj@active.ident ==  obj@meta.data$CellType)

saveRDS(object = obj,
        # file = "DuranteNatNeurosci2020_seurat_clustering_as_reported_in_paper_seurat_v3.0.0.rds")
        file = "DuranteNatNeurosci2020_seurat_clustering_as_reported_in_paper_seurat_v3.1.4.rds")
# obj <- readRDS(file = "DuranteNatNeurosci2020_seurat_clustering_as_reported_in_paper_seurat_v3.0.0.rds")
obj <- readRDS(file = "DuranteNatNeurosci2020_seurat_clustering_as_reported_in_paper_seurat_v3.1.4.rds")


# Figure 2, Panel A
# Plot violin plots of ACE2 and TMPRSS2 marker genes for all cells gated on ACE2

data.use <- obj@assays$RNA@data
data.use <- data.use[rownames(x = data.use) %in% c("ACE2", "TMPRSS2"), , drop = FALSE]
data.use <- as.matrix(x = data.use)
ACE2pos <- colnames(x = data.use)[data.use["ACE2",] > 0]
TMPRSS2pos <- colnames(x = data.use)[data.use["TMPRSS2",] > 0]
data.use <- reshape2::melt(data = data.use,
                           varnames = c("Gene", "Cell"),
                           value.name = "Count")
data.use$ACE2pos <- ifelse(test = data.use$Cell %in% ACE2pos,
                           yes = "ACE2 positive",
                           no = "ACE2 negative")
data.use$ACE2pos <- factor(x = data.use$ACE2pos,
                           levels = c("ACE2 negative", "ACE2 positive"))
data.use$TMPRSS2pos <- ifelse(test = data.use$Cell %in% TMPRSS2pos,
                              yes = "TMPRSS2 positive",
                              no = "TMPRSS2 negative")
data.use$TMPRSS2pos <- factor(x = data.use$TMPRSS2pos,
                              levels = c("TMPRSS2 negative", "TMPRSS2 positive"))
data.use$Count <- expm1(x = data.use$Count)

axis.text.x.label.ace2 <- c("ACE2 positive" = expression(paste(italic("ACE2"),
                                                               " positive",
                                                               sep = "")),
                            "ACE2 negative" = expression(paste(italic("ACE2"),
                                                               " negative",
                                                               sep = "")))
axis.text.x.label.tmprss2 <- c("TMPRSS2 positive" = expression(paste(italic("TMPRSS2"),
                                                                     " positive",
                                                                     sep = "")),
                               "TMPRSS2 negative" = expression(paste(italic("TMPRSS2"),
                                                                     " negative",
                                                                     sep = "")))

vln.plot.ace2.all <- ggplot(data = data.use,
                            mapping = aes(x = ACE2pos,
                                          y = Count + 1,
                                          fill = Gene,
                                          colour = Gene)) +
  geom_violin(position = position_dodge(width = 1),
              scale = "width",
              colour = "black",
              size = 0.1) +
  # geom_point(position = position_jitterdodge(jitter.width = 0.3,
  #                                            jitter.height = 0,
  #                                            dodge.width = 1),
  #            colour = "black",
  #            size = 0.1,
  #            show.legend = FALSE) +
  scale_x_discrete(label = axis.text.x.label.ace2) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1, 6, 31),
                     labels = c(0, 5, 30),
                     expand = c(0,0)) +
  expand_limits(y = c(1, 31)) +
  # scale_fill_manual(values = c("ACE2" = "violetred1",
  #                              "TMPRSS2" = "lightslateblue")) + # without ggtext
  scale_fill_manual(values = c("ACE2" = "violetred1",
                               "TMPRSS2" = "lightslateblue"),
                    labels = c("ACE2" = "<span style='color:violetred1'>ACE2</span>",
                               "TMPRSS2" = "<span style='color:lightslateblue'>TMPRSS2</span>")) + # with ggtext
  guides(fill = guide_legend(override.aes = list(colour = NA),
                             keywidth = 0.75,
                             keyheight = 0.75)) +
  facet_wrap(facets = ~ ACE2pos,
             scales = "free_x") +
  labs(title = "All cells",
       y = expression(paste("Normalized counts (log"[10], " scale)", sep = ""))) +
  theme_classic() +
  theme(axis.title.x = ggplot2::element_blank(),
        # axis.text.x = ggplot2::element_text(colour = "black",
        #                                     size = 7,
        #                                     family = "Arial",
        #                                     angle = 45,
        #                                     hjust = 1,
        #                                     vjust = 1),
        axis.title.y = ggplot2::element_text(colour = "black",
                                             size = 7,
                                             family = "Arial"),
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
        # legend.text = ggplot2::element_text(colour = "black",
        #                                     size = 6,
        #                                     family = "Arial",
        #                                     face = "italic"), # without ggtext
        legend.text = ggtext::element_markdown(colour = "black",
                                               size = 6,
                                               family = "Arial",
                                               face = "italic"), # with ggtext
        plot.title = ggplot2::element_text(colour = "black",
                                           size = 8,
                                           family = "Arial",
                                           # face = "bold",
                                           hjust = 0.5,
                                           vjust = 0),
        panel.spacing = unit(x = 2, units = "mm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.box.background = element_blank(),
        panel.border = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

FixSizeAndSave(plot = vln.plot.ace2.all,
               # filename = "DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_20200325.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_all_cells_20200330_figure2_panelA.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_all_cells_20200429_figure2_panelA.pdf", # without ggtext
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_all_cells_20200429_ggtext_figure2_panelA.pdf", # with ggtext
               is.ggassemble = FALSE,
               # panel.width = 4,
               # panel.height = 4,
               panel.width = 1.6,
               panel.height = 4,
               unit.use = "cm",
               margin = 0.5,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 2, Panel B
# Plot Clusters on UMAP

colors.of.interest <- c("#1F78B4", "#FF7F00", "#66C2A5", "#E7298A")
types.of.interest <- c("Respiratory ciliated cells",
                       "Respiratory early secretory cells",
                       "Respiratory horizontal basal cells",
                       "Sustentacular cells")
abbreviation.of.interest <- c("Resp. ciliated",
                              "Resp. early secretory",
                              "Resp. HBC",
                              "Sustentacular")
names(x = abbreviation.of.interest) <- types.of.interest
colors.to.remove <- c("#FFFF99", "#666666")
all.clusters <- unique(x = as.character(x = obj@active.ident))

n.clusters <- length(x = all.clusters)
# colors <- scales::hue_pal()(n.clusters)
qual.colors <- RColorBrewer::brewer.pal.info
qual.colors <- qual.colors[qual.colors$category == "qual" & qual.colors$colorblind == TRUE,]
colors <- rev(x = unique(x = unlist(x = mapply(FUN = brewer.pal,
                                               n = qual.colors$maxcolors,
                                               name = rownames(x = qual.colors)))))
colors <- colors[!colors %in% colors.to.remove]
colors <- c(colors.of.interest, colors[!colors %in% colors.of.interest])
types <- c(types.of.interest, all.clusters[!all.clusters %in% types.of.interest])
names(x = colors) <- types

plot.cluster.umap <- PlotClustersOntSNE(
  tsne.data = obj@reductions$umap@cell.embeddings,
  cluster.data = data.frame(Cluster = obj@active.ident),
  cellnames.as.rownames.in.tsne.data = TRUE,
  cellnames.as.rownames.in.cluster.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  cluster.ident.name = "Cluster",
  cluster.colors = colors,
  add.border.to.points = FALSE,
  point.size = 0.5,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  alpha.use = 0.75,
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

plot.data.list <- split(x = plot.cluster.umap$data, f = plot.cluster.umap$data$Cluster)
centroid.data.list <- lapply(X = plot.data.list,
                             FUN = function(cluster.umap) {
                               centroid.data <- data.frame(Cluster.label = unique(x = cluster.umap$Cluster),
                                                           UMAP_1 = median(x = cluster.umap$UMAP_1),
                                                           UMAP_2 = median(x = cluster.umap$UMAP_2))
                               return(centroid.data)
                             })
centroid.data <- dplyr::bind_rows(centroid.data.list)
centroid.data <- centroid.data[centroid.data$Cluster.label %in% types.of.interest,]
centroid.data$Cluster.label.abbrev <- abbreviation.of.interest[as.character(x = centroid.data$Cluster.label)]
centroid.data$Cluster.label.ggtext <- paste("<span style='color:",
                                            colors[as.character(x = centroid.data$Cluster.label)],
                                            "'>",
                                            centroid.data$Cluster.label.abbrev,
                                            "</span>",
                                            sep = "")

plot.cluster.umap <- plot.cluster.umap +
  # geom_text(data = centroid.data,
  #           mapping = aes(x = UMAP_1,
  #                         y = UMAP_2,
  #                         label = Cluster.label.abbrev),
  #           colour = "black",
  #           size = 6 * (1/72 * 25.4),
  #           family = "Arial",
  #           show.legend = FALSE,
  #           inherit.aes = FALSE) + # without ggtext
  geom_richtext(data = centroid.data,
                mapping = aes(x = UMAP_1,
                              y = UMAP_2,
                              label = Cluster.label.ggtext),
                colour = "black",
                size = 6 * (1/72 * 25.4),
                family = "Arial",
                show.legend = FALSE,
                inherit.aes = FALSE,
                fill = NA,
                label.color = NA,
                label.padding = grid::unit(x = rep(x = 0,
                                                   times = nrow(x = centroid.data)),
                                           units = "pt")) + # with ggtext
  labs(x = "UMAP 1 (a.u.)",
       y = "UMAP 2 (a.u.)",
       colour = "Cluster IDs (# cells)") +
  theme(legend.title = element_text(colour = "black",
                                    size = 7,
                                    family = "Arial",
                                    face = "plain"))

FixSizeAndSave(plot = plot.cluster.umap,
               # filename = "DuranteNatNeurosci2020_UMAP_clusters_plot_20200325.pdf",
               # filename = "DuranteNatNeurosci2020_UMAP_clusters_plot_20200329_version2.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_clusters_plot_20200330_figure2_panelB.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_clusters_plot_20200429_figure2_panelB.pdf", # without ggtext
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_clusters_plot_20200429_ggtext_figure2_panelB.pdf", # with ggtext
               is.ggassemble = FALSE,
               panel.width = 5,
               panel.height = 5,
               unit.use = "cm",
               margin = 0.1,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 2, Panel C
# Plot ACE2 and TMPRSS2 marker genes on UMAP

gradient.color <- c(viridis::cividis(2)[2], rev(x = viridis::magma(11))[-c(1:2)])

data.use.umap.coexpr <- expm1(x = obj@assays$RNA@data[c("ACE2", "TMPRSS2", "ERMN", "GSTA2"),])

ACE2.TMPRSS2 <- Matrix::Matrix(data = ifelse(test = data.use.umap.coexpr["ACE2",, drop = FALSE] * data.use.umap.coexpr["TMPRSS2",, drop = FALSE] > 0,
                                             yes = 1,
                                             no = 0),
                               nrow = 1,
                               dimnames = list("ACE2/TMPRSS2", colnames(x = data.use.umap.coexpr)))

data.use.umap.coexpr <- rbind(data.use.umap.coexpr, ACE2.TMPRSS2)


plot.markers.expression <- PlotExpressionOntSNE(
  data.use = data.use.umap.coexpr,
  genes.use = rownames(x = data.use.umap.coexpr)[1:4],
  ggtext.gene.names = c("ACE2" = "<span style='color:violetred1'>ACE2</span>",
                        "TMPRSS2" = "<span style='color:lightslateblue'>TMPRSS2</span>",
                        "ERMN" = "ERMN",
                        "GSTA2" = "GSTA2"), # with ggtext
  tsne.data = obj@reductions$umap@cell.embeddings,
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
  ncol.use = 4,
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
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_ERMN_GSTA2_expression_plot_20200330_figure2_panelC.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_ERMN_GSTA2_expression_plot_20200429_figure2_panelC.pdf", # without ggtext
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_ERMN_GSTA2_expression_plot_20200429_ggtext_figure2_panelC.pdf", # with ggtext
               is.ggassemble = TRUE,
               panel.width = 4.8,
               panel.height = 4.8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 2, Panel D
# Plot coexpression of ACE2 with TMPRSS2 on UMAP

plot.markers.coexpression <- PlotExpressionOntSNE(
  data.use = data.use.umap.coexpr,
  genes.use = rownames(x = data.use.umap.coexpr)[c(5,5)], # size do not match between cowplot and patchwork; trick to match size
  ggtext.gene.names = c("ACE2/TMPRSS2" = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>",
                        "ACE2/TMPRSS2" = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>"), # with ggtext
  tsne.data = obj@reductions$umap@cell.embeddings,
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
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_coexpressing_cells_plot_20200330_figure2_panelD.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_coexpressing_cells_plot_20200429_figure2_panelD.pdf", # without ggtext
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_coexpressing_cells_plot_20200429_ggtext_figure2_panelD.pdf", # with ggtext
               is.ggassemble = FALSE, # trick to match size with previous plots
               panel.width = 4.8,
               panel.height = 4.8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 2, Panel F
# Plot violin plots of ACE2 and TMPRSS2 marker genes for sustentacular cells gated on ACE2

data.use <- merge(x = data.use,
                  y = obj@meta.data,
                  by.x = "Cell",
                  by.y = "cell.names")

data.use.sustentacular <- data.use[data.use$CellType == "Sustentacular cells",]

vln.plot.ace2.sustentacular <- ggplot(data = data.use.sustentacular,
                                      mapping = aes(x = ACE2pos,
                                                    y = Count + 1,
                                                    fill = Gene,
                                                    colour = Gene)) +
  geom_violin(position = position_dodge(width = 1),
              scale = "width",
              colour = "black",
              size = 0.1) +
  # geom_point(position = position_jitterdodge(jitter.width = 0.3,
  #                                            jitter.height = 0,
  #                                            dodge.width = 1),
  #            colour = "black",
  #            size = 0.1,
  #            show.legend = FALSE) +
  scale_x_discrete(label = axis.text.x.label.ace2) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1, 6, 31),
                     labels = c(0, 5, 30),
                     expand = c(0,0)) +
  expand_limits(y = c(1, 31)) +
  scale_fill_manual(values = c("ACE2" = "violetred1",
                               "TMPRSS2" = "lightslateblue")) +
  guides(fill = guide_legend(override.aes = list(colour = NA),
                             keywidth = 0.75,
                             keyheight = 0.75)) +
  facet_wrap(facets = ~ ACE2pos,
             scales = "free_x") +
  labs(title = "Sustentacular cells",
       y = expression(paste("Normalized counts (log"[10], " scale)", sep = ""))) +
  theme_classic() +
  theme(axis.title.x = ggplot2::element_blank(),
        # axis.text.x = ggplot2::element_text(colour = "black",
        #                                     size = 7,
        #                                     family = "Arial",
        #                                     angle = 45,
        #                                     hjust = 1,
        #                                     vjust = 1),
        axis.title.y = ggplot2::element_text(colour = "black",
                                             size = 7,
                                             family = "Arial"),
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
                                            family = "Arial",
                                            face = "italic"),
        plot.title = ggplot2::element_text(colour = "black",
                                           size = 8,
                                           family = "Arial",
                                           # face = "bold",
                                           hjust = 0.5,
                                           vjust = 0),
        legend.position = "none",
        panel.spacing = unit(x = 2, units = "mm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.box.background = element_blank(),
        panel.border = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

FixSizeAndSave(plot = vln.plot.ace2.sustentacular,
               # filename = "DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_20200325.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_sustentacular_cells_20200330_figure2_panelE.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_sustentacular_cells_20200429_figure2_panelF.pdf",
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_sustentacular_cells_20200429_ggtext_figure2_panelF.pdf",
               is.ggassemble = FALSE,
               # panel.width = 4,
               # panel.height = 4,
               panel.width = 1.6,
               panel.height = 4,
               unit.use = "cm",
               margin = 0.5,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 2, Panel G
# Plot violin plots of ACE2 and TMPRSS2 marker genes for respiratory ciliated cells gated on ACE2

data.use.ciliated <- data.use[data.use$CellType == "Respiratory ciliated cells",]

vln.plot.ace2.ciliated <- ggplot(data = data.use.ciliated,
                                 mapping = aes(x = ACE2pos,
                                               y = Count + 1,
                                               fill = Gene,
                                               colour = Gene)) +
  geom_violin(position = position_dodge(width = 1),
              scale = "width",
              colour = "black",
              size = 0.1) +
  # geom_point(position = position_jitterdodge(jitter.width = 0.3,
  #                                            jitter.height = 0,
  #                                            dodge.width = 1),
  #            colour = "black",
  #            size = 0.1,
  #            show.legend = FALSE) +
  scale_x_discrete(label = axis.text.x.label.ace2) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1, 6, 31),
                     labels = c(0, 5, 30),
                     expand = c(0,0)) +
  expand_limits(y = c(1, 31)) +
  # scale_fill_manual(values = c("ACE2" = "violetred1",
  #                              "TMPRSS2" = "lightslateblue")) + # without ggtext
  scale_fill_manual(values = c("ACE2" = "violetred1",
                               "TMPRSS2" = "lightslateblue"),
                    labels = c("ACE2" = "<span style='color:violetred1'>ACE2</span>",
                               "TMPRSS2" = "<span style='color:lightslateblue'>TMPRSS2</span>")) + # with ggtext
  guides(fill = guide_legend(override.aes = list(colour = NA),
                             keywidth = 0.75,
                             keyheight = 0.75)) +
  facet_wrap(facets = ~ ACE2pos,
             scales = "free_x") +
  labs(title = "Respiratory ciliated cells") +
  theme_classic() +
  theme(axis.title = ggplot2::element_blank(),
        # axis.text.x = ggplot2::element_text(colour = "black",
        #                                     size = 7,
        #                                     family = "Arial",
        #                                     angle = 45,
        #                                     hjust = 1,
        #                                     vjust = 1),
        # axis.title.y = ggplot2::element_text(colour = "black",
        #                                      size = 7,
        #                                      family = "Arial"),
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
        # legend.text = ggplot2::element_text(colour = "black",
        #                                     size = 6,
        #                                     family = "Arial",
        #                                     face = "italic"), # without ggtext
        legend.text = ggtext::element_markdown(colour = "black",
                                               size = 6,
                                               family = "Arial",
                                               face = "italic"), # with ggtext
        plot.title = ggplot2::element_text(colour = "black",
                                           size = 8,
                                           family = "Arial",
                                           # face = "bold",
                                           hjust = 0.5,
                                           vjust = 0),
        panel.spacing = unit(x = 2, units = "mm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.box.background = element_blank(),
        panel.border = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

FixSizeAndSave(plot = vln.plot.ace2.ciliated,
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_respiratory_ciliated_cells_20200330_figure2_panelF.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_respiratory_ciliated_cells_20200429_figure2_panelG.pdf", # without ggtext
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_ACE2_respiratory_ciliated_cells_20200429_ggtext_figure2_panelG.pdf", # with ggtext
               is.ggassemble = FALSE,
               # panel.width = 4,
               # panel.height = 4,
               panel.width = 1.6,
               panel.height = 4,
               unit.use = "cm",
               margin = 0.5,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 2, Panel E
# Plot average expression of ACE2 and TMPRSS2 per cluster

expr.data <- as.data.frame(x = t(x = expm1(x = obj@assays$RNA@data[c("ACE2", "TMPRSS2"),])))
expr.data$cell.names <- rownames(x = expr.data)
expr.data <- merge(x = expr.data,
                   y = obj@meta.data,
                   by = "cell.names")
expr.data.list <- split(x = expr.data,
                        f = expr.data$CellType)
expr.data.list <- lapply(X = expr.data.list,
                         FUN = function(expr.data.sub) {
                           data.out <- data.frame(ACE2.mean = mean(x = expr.data.sub$ACE2),
                                                  TMPRSS2.mean = mean(x = expr.data.sub$TMPRSS2),
                                                  CellType = unique(x = expr.data.sub$CellType))
                           return(data.out)
                         })
expr.data.mean <- dplyr::bind_rows(expr.data.list,
                                   .id = NULL)
max.ace2 <- round(x = max(expr.data.mean$ACE2.mean),
                  digits = 3)
breaks.ace2 <- c(0, round(x = max.ace2/2, digits = 2), max.ace2)
max.tmprss2 <- round(x = max(expr.data.mean$TMPRSS2.mean),
                     digits = 2)
breaks.tmprss2 <- c(0, round(x = max.tmprss2/2, digits = 2), max.tmprss2)


mean.expr.plot <- ggplot(data = expr.data.mean,
                         mapping = aes(x = ACE2.mean,
                                       y = TMPRSS2.mean,
                                       color = CellType)) +
  geom_point(size = 1) +
  ggrepel::geom_text_repel(data = expr.data.mean[expr.data.mean$CellType %in% types.of.interest,],
                           mapping = aes(label = CellType,
                                         segment.colour = CellType),
                           size = 6 * (1/72 * 25.4),
                           family = "Arial",
                           min.segment.length = unit(x = 0,
                                                     units = "mm"),
                           nudge_y = 0.01 * max.tmprss2,
                           force = TRUE,
                           direction = "both",
                           max.overlaps = Inf,
                           show.legend = FALSE) +
  scale_x_continuous(limits = c(0, max.ace2),
                     breaks = breaks.ace2,
                     labels = breaks.ace2,
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, max.tmprss2),
                     breaks = breaks.tmprss2,
                     labels = breaks.tmprss2,
                     expand = c(0,0)) +
  scale_color_manual(values = colors,
                     aesthetics = c("colour", "segment.colour"),
                     name = NULL) +
  # labs(x = expression(italic(ACE2)~"(mean expression)"),
  #      y = expression(italic(TMPRSS2)~"(mean expression)")) +
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
               # filename = "DuranteNatNeurosci2020_ACE2_TMPRSS2_mean_expression_per_cluster_20200329_version2.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_mean_expression_per_cluster_20200330_figure2_panelG.pdf",
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_mean_expression_per_cluster_20200429_figure2_panelE.pdf", # without ggtext
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_mean_expression_per_cluster_20200429_ggtext_figure2_panelE.pdf", # with ggtext
               is.ggassemble = FALSE,
               panel.width = 4,
               panel.height = 4,
               unit.use = "cm",
               margin = 0.5,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 2, Panel H
# Plot the number of cells expressing ACE2 and TMPRSS2 in all clusters

percentage.data <- as.data.frame(x = as.matrix(x = t(x = data.use.umap.coexpr[c("ACE2", "TMPRSS2", "ACE2/TMPRSS2"),])))
percentage.data$cell.names <- rownames(x = percentage.data)
percentage.data <- merge(x = percentage.data,
                         y = obj@meta.data,
                         by = "cell.names",
                         sort = FALSE)
percentage.data.list <- split(x = percentage.data,
                              f = percentage.data$CellType)
percentage.data.list <- lapply(X = percentage.data.list,
                               FUN = function(percentage.data.sub) {
                                 out.data <- data.frame(CellType = unique(x = percentage.data.sub$CellType),
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

max.y.ace2.number <- pretty(x = max(percentage.data$ACE2.number),
                            n = 2)[2]
breaks.ace2.number <- seq(from = 0,
                          to = max.y.ace2.number,
                          by = 10)
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
                           mapping = aes(x = CellType,
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
                              mapping = aes(x = CellType,
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
                             mapping = aes(x = CellType,
                                           y = ACE2.TMPRSS2.number,
                                           fill = CellType)) +
  geom_bar(stat = "identity",
           colour = NA,
           width = 0.5) +
  scale_y_continuous(breaks = breaks.coexpr.number,
                     labels = labels.coexpr.number,
                     expand = c(0,0)) +
  scale_fill_manual(values = colors[levels(x = percentage.data$CellType)]) +
  expand_limits(y = c(0, max.y.coexpr.number)) +
  # labs(title = "ACE2/TMPRSS2") + # without ggtext
  labs(title = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>") + # with ggtext
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = colors[levels(x = percentage.data$CellType)],
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
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_number_of_expressing_cells_all_clusters_plot_20200429_figure2_panelH.pdf", # without ggtext
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_number_of_expressing_cells_all_clusters_plot_20200429_ggtext_figure2_panelH.pdf", # with ggtext
               is.ggassemble = TRUE,
               panel.width = 13,
               panel.height = 1.5,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Figure 2, Panel I
# Plot the percentage of cells expressing ACE2 and TMPRSS2 in coexpressing clusters

percentage.data.coexpressing <- percentage.data[percentage.data$coexpressing,]
percentage.data.coexpressing$CellType <- droplevels(x = percentage.data.coexpressing$CellType)

max.y.ace2.percentage <- ceiling(x = max(percentage.data.coexpressing$ACE2.percentage))
breaks.ace2.percentage <- seq(from = 0,
                              to = max.y.ace2.percentage,
                              by = 1)
labels.ace2.percentage <- ifelse(test = breaks.ace2.percentage %in% c(0, max.y.ace2.percentage),
                                 yes = breaks.ace2.percentage,
                                 no = "")

max.y.tmprss2.percentage <- pretty(x = max(percentage.data.coexpressing$TMPRSS2.percentage),
                                   n = 2)[2]
breaks.tmprss2.percentage <- seq(from = 0,
                                 to = max.y.tmprss2.percentage,
                                 by = 10)
labels.tmprss2.percentage <- ifelse(test = breaks.tmprss2.percentage %in% c(0, max.y.tmprss2.percentage),
                                    yes = breaks.tmprss2.percentage,
                                    no = "")

max.y.coexpr.percentage <- ceiling(x = max(percentage.data.coexpressing$ACE2.TMPRSS2.percentage))
breaks.coexpr.percentage <- seq(from = 0,
                                to = max.y.coexpr.percentage,
                                by = 1)
labels.coexpr.percentage <- ifelse(test = breaks.coexpr.percentage %in% c(0, max.y.coexpr.percentage),
                                   yes = breaks.coexpr.percentage,
                                   no = "")

plot.ace2.percentage <- ggplot(data = percentage.data.coexpressing,
                               mapping = aes(x = CellType,
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
                                  mapping = aes(x = CellType,
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
                                 mapping = aes(x = CellType,
                                               y = ACE2.TMPRSS2.percentage,
                                               fill = CellType)) +
  geom_bar(stat = "identity",
           colour = NA,
           width = 0.5) +
  scale_y_continuous(breaks = breaks.coexpr.percentage,
                     labels = labels.coexpr.percentage,
                     expand = c(0,0)) +
  scale_fill_manual(values = colors[levels(x = percentage.data.coexpressing$CellType)]) +
  expand_limits(y = c(0, max.y.coexpr.percentage)) +
  # labs(title = "ACE2/TMPRSS2") + # without ggtext
  labs(title = "<span style='color:violetred1'>*ACE2*</span>/<span style='color:lightslateblue'>*TMPRSS2*</span>") + # with ggtext
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = colors[levels(x = percentage.data.coexpressing$CellType)],
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
               # filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_percentage_of_expressing_cells_coexpressing_clusters_plot_20200429_figure2_panelI.pdf", # without ggtext
               filename = "./Figures_for_paper/DuranteNatNeurosci2020_ACE2_TMPRSS2_percentage_of_expressing_cells_coexpressing_clusters_plot_20200429_ggtext_figure2_panelI.pdf", # with ggtext
               is.ggassemble = TRUE,
               panel.width = 4,
               panel.height = 1.5,
               unit.use = "cm",
               margin = 1,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Save percentage.data as .csv table for paper
out.percentage  <- percentage.data
colnames(x = out.percentage) <- gsub(pattern = "\\.",
                                     replacement = "_",
                                     x = colnames(x = out.percentage))
colnames(x = out.percentage) <- ifelse(test = grepl(pattern = "CellType",
                                                    x = colnames(x = out.percentage)),
                                       yes = "cluster_name",
                                       no = colnames(x = out.percentage))
colnames(x = out.percentage) <- gsub(pattern = "ACE2_TMPRSS2",
                                     replacement = "ACE2/TMPRSS2",
                                     x = colnames(x = out.percentage))

write.csv(x = out.percentage,
          file = "DuranteNatNeurosci2020_ACE2_TMPRSS2_expressing_cells_percentage_per_cluster.csv",
          quote = TRUE,
          row.names = FALSE)






# Not used for the graph
# Plot violin plots of ACE2 and TMPRSS2 marker genes gated on TMPRSS2

vln.plot.tmprss2.all <- ggplot(data = data.use,
                               mapping = aes(x = TMPRSS2pos,
                                             y = Count + 1,
                                             fill = Gene,
                                             colour = Gene)) +
  geom_violin(position = position_dodge(width = 1),
              scale = "width",
              colour = "black",
              size = 0.1) +
  # geom_point(position = position_jitterdodge(jitter.width = 0.3,
  #                                            jitter.height = 0,
  #                                            dodge.width = 1),
  #            colour = "black",
  #            size = 0.1,
  #            show.legend = FALSE) +
  scale_x_discrete(label = axis.text.x.label.tmprss2) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1, 6, 31),
                     labels = c(0, 5, 30),
                     expand = c(0,0)) +
  expand_limits(y = c(1, 31)) +
  scale_fill_manual(values = c("ACE2" = "violetred1",
                               "TMPRSS2" = "lightslateblue")) +
  guides(fill = guide_legend(override.aes = list(colour = NA),
                             keywidth = 0.75,
                             keyheight = 0.75)) +
  facet_wrap(facets = ~ TMPRSS2pos,
             scales = "free_x") +
  labs(y = expression("Normalized counts"~(log[10]~scale))) +
  theme_classic() +
  theme(axis.title.x = ggplot2::element_blank(),
        # axis.text.x = ggplot2::element_text(colour = "black",
        #                                     size = 7,
        #                                     family = "Arial",
        #                                     angle = 45,
        #                                     hjust = 1,
        #                                     vjust = 1),
        axis.title.y = ggplot2::element_text(colour = "black",
                                             size = 7,
                                             family = "Arial"),
        # axis.text.y = ggplot2::element_text(colour = "black",
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
                                            family = "Arial",
                                            face = "italic"),
        legend.position = "top",
        panel.spacing = unit(x = 2, units = "mm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.box.background = element_blank(),
        panel.border = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

FixSizeAndSave(plot = vln.plot.tmprss2.all,
               filename = "DuranteNatNeurosci2020_ACE2_TMPRSS2_violin_plot_gated_on_TMPRSS2_all_cells_20200330.pdf",
               is.ggassemble = FALSE,
               panel.width = 1.6,
               panel.height = 3.6,
               unit.use = "cm",
               margin = 0.5,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Plot Patients on UMAP

orig.data <- obj@meta.data[,"orig.ident", drop = FALSE]
orig.data$orig.ident <- factor(x = orig.data$orig.ident,
                               levels = rev(x = sort(x = unique(x = orig.data$orig.ident))))

plot.patient.umap <- PlotClustersOntSNE(
  tsne.data = obj@reductions$umap@cell.embeddings,
  cluster.data = obj@meta.data[,"orig.ident", drop = FALSE],
  cellnames.as.rownames.in.tsne.data = TRUE,
  cellnames.as.rownames.in.cluster.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  cluster.ident.name = "orig.ident",
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  alpha.use = 0.75,
  font.family = "Arial",
  import.font = FALSE,
  fix.aspect.ratio = FALSE,
  legend.point.size = 2,
  axis.title.size = 8,
  legend.title.size = 8,
  legend.text.size = 7,
  legend.text.space = unit(x = 3, units = "mm"),
  axis.line.size = rel(0.5),
  arrow.size = rel(1),
  range.scale = 0.2
)

plot.patient.umap <- plot.patient.umap +
  labs(x = "UMAP 1 (a.u.)",
       y = "UMAP 2 (a.u.)",
       colour = "Patient IDs (# cells)")

plot.patient.umap.1.top <- plot.patient.umap
plot.patient.umap.1.top$data <- plot.patient.umap.1.top$data[order(plot.patient.umap.1.top$data$orig.ident, decreasing = TRUE),]
plot.patient.umap.4.top <- plot.patient.umap +
  theme(legend.position = "none")
plot.patient.umap.4.top$data <- plot.patient.umap.4.top$data[order(plot.patient.umap.4.top$data$orig.ident, decreasing = FALSE),]

FixSizeAndSave(plot = plot.patient.umap.1.top,
               filename = "DuranteNatNeurosci2020_UMAP_patients_plot_patient_1_top_20200325.pdf",
               is.ggassemble = FALSE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)

FixSizeAndSave(plot = plot.patient.umap.4.top,
               filename = "DuranteNatNeurosci2020_UMAP_patients_plot_patient_4_top_20200325.pdf",
               is.ggassemble = FALSE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Plot ACE2 marker gene on UMAP

plot.ace2.expression <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = "ACE2",
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = TRUE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.ace2.expression,
               filename = "DuranteNatNeurosci2020_UMAP_ACE2_expression_plot_20200325.pdf",
               is.ggassemble = FALSE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Combine both previous plots
caption <- 
  "Figure 1. \t Re-analysis of the data reported in Durante et al. (2020) Nat. Neurosci.
(A) \t \t \t \t \t \t \t \t \t \t UMAP plot displaying the various clusters identified using Seurat version 3.0.0.
(B) \t \t \t \t \t \t \t \t \t \t ACE2 normalized UMI expression displayed on the UMAP plot.
(C-D) \t \t \t \t \t UMAP plot displaying the patient of origin for each cell. 
\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t (C) \t \t Patient 1 cells are plotted on top.
\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t (D) \t \t Patient 4 cells are plotted on top.
"
comb.plot <- patchwork::wrap_plots(plot.cluster.umap,
                                   expression.plot,
                                   plot.patient.umap.1.top,
                                   plot.patient.umap.4.top,
                                   ncol = 2) +
  plot_annotation(caption = caption,
                  tag_levels = "A",
                  theme = theme(plot.caption = element_text(family = "Arial",
                                                            colour = "black",
                                                            size = 12,
                                                            hjust = 0,
                                                            lineheight = unit(x = 1,
                                                                              units = "mm"),
                                                            margin = margin(t = 10,
                                                                            r = 0,
                                                                            b = 0,
                                                                            l = 0,
                                                                            unit = "mm")))) & 
  theme(text = element_text(family = "Arial",
                            colour = "black",
                            size = 12))

FixSizeAndSave(plot = comb.plot,
               filename = "DuranteNatNeurosci2020_UMAP_combined_clusters_patients_ACE2_expression_plots_20200325.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# TMPRSS2 expression
# Plot ACE2 and TMPRSS2 marker genes on UMAP

plot.ace2.tmprss2.expression <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = c("ACE2", "TMPRSS2"),
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
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

FixSizeAndSave(plot = plot.ace2.tmprss2.expression,
               filename = "DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_expression_plot_20200325.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# TMPRSS2 and ERMN expression
# Check double and triple positive cells with ACE2

data.use.umap.exp <- expm1(x = obj@assays$RNA@data[c("ACE2", "TMPRSS2", "ERMN"),])
# ACE2.TMPRSS2 <- data.use.umap.exp["ACE2",, drop = FALSE] * data.use.umap.exp["TMPRSS2",, drop = FALSE]
# rownames(x = ACE2.TMPRSS2) <- "ACE2/TMPRSS2"
# ACE2.ERMN <- data.use.umap.exp["ACE2",, drop = FALSE] * data.use.umap.exp["ERMN",, drop = FALSE]
# rownames(x = ACE2.ERMN) <- "ACE2/ERMN"
# ACE2.TMPRSS2.ERMN <- data.use.umap.exp["ACE2",, drop = FALSE] * data.use.umap.exp["TMPRSS2",, drop = FALSE] * data.use.umap.exp["ERMN",, drop = FALSE]
# rownames(x = ACE2.TMPRSS2.ERMN) <- "ACE2/TMPRSS2/ERMN"
ACE2.TMPRSS2 <- Matrix::Matrix(data = ifelse(test = data.use.umap.exp["ACE2",, drop = FALSE] * data.use.umap.exp["TMPRSS2",, drop = FALSE] > 0,
                                             yes = 1,
                                             no = 0),
                               nrow = 1,
                               dimnames = list("ACE2/TMPRSS2", colnames(x = data.use.umap.exp)))
ACE2.ERMN <- Matrix::Matrix(data = ifelse(test = data.use.umap.exp["ACE2",, drop = FALSE] * data.use.umap.exp["ERMN",, drop = FALSE] > 0,
                                          yes = 1,
                                          no = 0),
                            nrow = 1,
                            dimnames = list("ACE2/ERMN", colnames(x = data.use.umap.exp)))
ACE2.TMPRSS2.ERMN <- Matrix::Matrix(data = ifelse(test = data.use.umap.exp["ACE2",, drop = FALSE] * data.use.umap.exp["TMPRSS2",, drop = FALSE] * data.use.umap.exp["ERMN",, drop = FALSE] > 0,
                                                  yes = 1,
                                                  no = 0),
                                    nrow = 1,
                                    dimnames = list("ACE2/TMPRSS2/ERMN", colnames(x = data.use.umap.exp)))

data.use.umap.exp <- rbind(data.use.umap.exp, ACE2.TMPRSS2, ACE2.ERMN, ACE2.TMPRSS2.ERMN)

plot.ace2.tmprss2.ermn.expression <- PlotExpressionOntSNE(
  data.use = data.use.umap.exp,
  genes.use = rownames(x = data.use.umap.exp)[1:3],
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = FALSE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.ace2.tmprss2.ermn.expression,
               # filename = "DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_ERMN_expression_plot_20200325.pdf",
               filename = "DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_ERMN_expression_plot_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)

plot.ace2.tmprss2.ermn.coexpression <- PlotExpressionOntSNE(
  data.use = data.use.umap.exp,
  genes.use = rownames(x = data.use.umap.exp)[4:6],
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = FALSE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = "none",
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.ace2.tmprss2.ermn.coexpression,
               filename = "DuranteNatNeurosci2020_UMAP_ACE2_TMPRSS2_ERMN_coexpressing_cells_plot_20200325.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Check sustentacular cells markers on the UMAP plot

sustentacular.cells.genes <- c("SOX2", "CYP1A2", "NOTCH2", # From FletcherCellStemCell2017
                               "CYP2A13", "CYP2J2", "GPX6", "ERMN") # From DuranteNatNeurosci2020

plot.sustentacular.cells <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = sustentacular.cells.genes,
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.sustentacular.cells,
               filename = "DuranteNatNeurosci2020_UMAP_sustentacular_cells_gene_expression_plot_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


sustentacular.cells.genes.2 <- c("SLC2A3", "EPAS1", "SIX1") # From SaraivaSciRep2015 (first 2) and SaraivaSciAdv2019 (last one)

plot.sustentacular.cells.2 <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = sustentacular.cells.genes.2,
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.sustentacular.cells.2,
               filename = "DuranteNatNeurosci2020_UMAP_sustentacular_cells_gene_expression_plot_2_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


sustentacular.cells.genes.3 <- c("KITLG", paste0("REEP", 1:6), "SOX2", "PAX6", "HES1", "HES5", "SIX1", "SMARCC1") # From a weird publication that Ivan gave us (KITLG = Steel; SMARCC1 = BAF155)

plot.sustentacular.cells.3 <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = sustentacular.cells.genes.3,
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 5,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.sustentacular.cells.3,
               filename = "DuranteNatNeurosci2020_UMAP_sustentacular_cells_gene_expression_plot_3_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Check Bowman's gland markers on the UMAP plot

bowmans.gland.genes <- c("SOX9", "SOX10", "GPX3", "MUC5AC", "MUC5B")

plot.bowmans.gland <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = bowmans.gland.genes,
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.bowmans.gland,
               filename = "DuranteNatNeurosci2020_UMAP_Bowmans_gland_gene_expression_plot_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


bowmans.gland.genes.2 <- c("ASCL3", "KRT18") # From a weird publication that Ivan gave us (KRT18 = Ck18)

plot.bowmans.gland.2 <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = bowmans.gland.genes.2,
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
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

FixSizeAndSave(plot = plot.bowmans.gland.2,
               filename = "DuranteNatNeurosci2020_UMAP_Bowmans_gland_gene_expression_plot_2_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Check respiratory ciliated cells markers on the UMAP plot

resp.ciliated.cells.genes <- c("CAPS", "C9orf24", "C20orf85", "GSTA2", "TMEM190", "AL357093.2")

plot.resp.ciliated.cells <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = resp.ciliated.cells.genes,
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.resp.ciliated.cells,
               # filename = "DuranteNatNeurosci2020_UMAP_respiratory_ciliated_cells_gene_expression_plot_20200325.pdf",
               filename = "DuranteNatNeurosci2020_UMAP_respiratory_ciliated_cells_gene_expression_plot_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Check respiratory secretory cells markers on the UMAP plot

resp.secretory.cells.genes <- c("SCGB1A1", "STATH", "C6orf58", "BPIFA1", "DMBT1", "PRB3", "LTF", "ZG16B", "PIGR")

plot.resp.secretory.cells <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = resp.secretory.cells.genes,
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.resp.secretory.cells,
               # filename = "DuranteNatNeurosci2020_UMAP_respiratory_secretory_cells_gene_expression_plot_20200325.pdf",
               filename = "DuranteNatNeurosci2020_UMAP_respiratory_secretory_cells_gene_expression_plot_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Check aquaporin genes experession on the UMAP plot

plot.aquaporin <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = c("AQP3", "AQP4", "AQP5"),
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 3,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.aquaporin,
               filename = "DuranteNatNeurosci2020_UMAP_aquaporin_genes_expression_plot_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Check aquaporin genes and MUC5AC experession on the UMAP plot

plot.aquaporin.muc5ac <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = c("AQP3", "AQP4", "AQP5", "MUC5AC"),
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
  legend.position = c(1,0),
  plot.colourbar.ticks = FALSE,
  colourbar.ticks.colour = "black",
  colourbar.frame.colour = NA,
  fix.aspect.ratio = FALSE,
  add.border = FALSE,
  ncol.use = 4,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.aquaporin.muc5ac,
               filename = "DuranteNatNeurosci2020_UMAP_aquaporin_MUC5AC_genes_expression_plot_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Check MUC5AC and MUC15 experession on the UMAP plot

plot.muc5ac.muc15 <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = c("MUC5AC", "MUC15"),
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  add.border.to.points = FALSE,
  point.size = 1,
  border.stroke = 0.01,
  border.colour = "#C3C3C3",
  color.gradient.use = rev(x = viridis::magma(10)),
  plot.title.size = 8,
  plot.title.position = "top",
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  legend.title.size = 7,
  legend.text.size = 7,
  legend.bar.width = 2.5,
  legend.bar.height = 0.4,
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

FixSizeAndSave(plot = plot.muc5ac.muc15,
               filename = "DuranteNatNeurosci2020_UMAP_MUC5AC_MUC15_expression_plot_20200329.pdf",
               is.ggassemble = TRUE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)


# Check BSG experession on the UMAP plot

plot.bsg <- PlotExpressionOntSNE(
  data.use = obj@assays$RNA@data,
  genes.use = "BSG",
  tsne.data = obj@reductions$umap@cell.embeddings,
  cellnames.as.rownames.in.tsne.data = TRUE,
  tsne_1_name = "UMAP_1",
  tsne_2_name = "UMAP_2",
  is.log.transformed = TRUE,
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
  ncol.use = 1,
  panel.margin = 0.1,
  legend.margin = 0.1,
  round.to.ceiling = TRUE,
  verbose = TRUE,
  legend.title.hjust = 0.5,
  legend.title.vjust = -1,
  legend.label.hjust = 0.5,
  legend.label.vjust = 3
)

FixSizeAndSave(plot = plot.bsg,
               filename = "DuranteNatNeurosci2020_UMAP_BSG_expression_plot_20200407.pdf",
               is.ggassemble = FALSE,
               panel.width = 8,
               panel.height = 8,
               unit.use = "cm",
               margin = 0.01,
               use.ggsave = TRUE,
               useDingbats = FALSE)

