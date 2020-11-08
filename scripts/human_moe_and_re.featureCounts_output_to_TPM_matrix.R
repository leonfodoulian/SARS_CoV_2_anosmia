require(stringr)

setwd("D:/Projects/ace2_expr/human_moe_and_re/")

source(file = "leonlib/cat.verbose.R")
source(file = "leonlib/PasteNames.R")
source(file = "leonlib/ReadGTF.R")
source(file = "leonlib/ReadfeatureCountsTable.R")
source(file = "leonlib/ComputeTPM.R")

pattern <- ".featureCounts_counts.updated_annotation.txt$"
path <- "."
gtf.file <- "E:/Annotations/Homo_sapiens.GRCh38.99.no_AC097625.2.gtf"

dirs <- c("moe", "re")
names(dirs) <- c("main olfactory", "respiratory")

files <- lapply(X = dirs, FUN = function(dir) {
  files <- list.files(path = dir,
                      pattern = pattern,
                      full.names = TRUE,
                      recursive = TRUE)
  names(x = files) <- str_match(pattern = "patient[0-9.]+_[a-z]+", string = files)
  return(files)
})

gtf.data <- ReadGTF(gtf.file = gtf.file,
                    columns.to.keep = c("gene_id", "gene_name"),
                    types.keep = "gene",
                    check.gene.name.dups = TRUE,
                    verbose = TRUE)


counts.list <- mapply(file = files,
                      condition = names(x = files),
                      FUN = function(file,
                                     condition,
                                     gtf.data) {
                        
                        counts.list <- mapply(FUN = ReadfeatureCountsTable,
                                              featureCounts.file = file,
                                              sample.name = names(x = file),
                                              SIMPLIFY = FALSE,
                                              USE.NAMES = TRUE)
                        
                        counts.list <- lapply(X = counts.list,
                                              FUN = ComputeTPM,
                                              counts.column = "Counts",
                                              length.column = "Length",
                                              verbose = TRUE)
                        
                        counts <- data.table::data.table(dplyr::bind_rows(counts.list,
                                                                          .id = NULL))
                        
                        counts <- merge(x = counts,
                                        y = gtf.data,
                                        by.x = "Geneid",
                                        by.y = "gene_id",
                                        all = TRUE,
                                        allow.cartesian = TRUE)
                        
                        cat.verbose(x = "Summing gene counts and length for genes with duplicated names (but not IDs)",
                                    verbose = TRUE)
                        counts <- counts[, list(Counts = sum(Counts,
                                                             na.rm = TRUE),
                                                TPM = sum(TPM,
                                                          na.rm = TRUE),
                                                Length = mean(x = Length,
                                                              na.rm = TRUE)),
                                         by = c("Sample", "gene_name")]
                        
                        counts <- data.frame(counts)
                        
                        return(counts)
                        
                        },
                      MoreArgs = list(gtf.data = gtf.data),
                      SIMPLIFY = FALSE,
                      USE.NAMES = TRUE)

counts <- dplyr::bind_rows(counts.list,
                           .id = NULL)
counts$tissue <- ifelse(endsWith(counts$Sample, "moe"), "olfactory", 
                        ifelse(endsWith(counts$Sample, "re"), "respiratory",
                               "unsure"))

saveRDS(object = counts, file = "Human_OE_RE_RNAseq.quantifications.rds")
