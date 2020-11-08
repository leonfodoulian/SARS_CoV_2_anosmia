require(ggplot2)
require(stringr)
require(patchwork)

counts <- readRDS(file = "Human_OE_RE_RNAseq.quantifications.rds")

markers <- read.csv("markers_list.txt", sep = "\t", stringsAsFactors = F)

cell.types <- markers$cluster
names(cell.types) <- markers$gene.name
counts.subset <- counts[counts$gene_name %in% markers$gene.name,]
counts.subset$cell.type <- counts.subset$gene_name
counts.subset$patient <- str_match(counts.subset$Sample, "patient[0-9.]+")

plot.list <- mapply(
       gene.name=markers$gene.name,
       cell.type=markers$cluster,
       i=seq_along(markers$cluster),
       FUN = function(gene.name, cell.type, i, d) {
         d <- d[d$gene_name == gene.name,]
         print(nrow(d))
         if (nrow(d) != 0) {
           
           plot <- ggplot() +
             geom_segment(data=aggregate(data=d, TPM+1 ~ tissue, FUN = mean),
                          aes(x=ifelse(tissue=="olfactory", 1, 2)-0.25, 
                              xend=ifelse(tissue=="olfactory", 1, 2)+0.25, 
                              y = `TPM + 1`, yend=`TPM + 1`)) +
             geom_point(data=d, aes(x=tissue, y=TPM+1, color=patient)) +
             scale_y_log10(limits=c(1, NA),
                           expand=expand_scale(mult = c(0, 0.3))) +
             scale_x_discrete(labels=c("olfactory"="olf.", "respiratory"="resp."),
                                       expand=c(1/3, 1/3)) +
             ggtitle(gene.name, cell.type) +
             theme_classic(base_family = "ArialMT",
                           base_size = 8) +
             theme(
               legend.position = "none",
               panel.grid = element_blank(),
               plot.title = element_text(size = 8, face = "italic", hjust = 0.5),
               panel.background = element_blank(),
               panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
               axis.line = element_blank(),
               plot.subtitle = element_text(size = 6, hjust = 0.5),
               axis.title.y = element_text(size = 7),
               axis.ticks = element_line(color = "black", size = 0.5),
               axis.text = element_text(color = "black"),
               axis.title.x = element_blank()
             ) +
             ylab("TPM+1 (log scale)")
             
           if (i != 1) {
             plot <- plot + theme(axis.title.y = element_blank()) 
           }
           
           if (i %in% 1:5) {
             plot <- plot + theme(axis.ticks.x = element_blank(),
                                  axis.text.x = element_blank())
           }
           
           return(plot) }
         else {
           return(NA)
         }
       },
       MoreArgs = list(d = counts.subset),
       SIMPLIFY = FALSE,
       USE.NAMES = TRUE
       )

final.plot <- wrap_plots(plot.list[!is.na(plot.list)], ncol = 5)

ggsave(filename = "human_bulk_oere.pdf", 
       plot = final.plot,
       device = "pdf", 
       units = "cm",
       width = 11.8, 
       height = 6, 
       useDingbats=FALSE)
