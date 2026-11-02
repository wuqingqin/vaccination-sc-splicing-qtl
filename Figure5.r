plot_root <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/PLOT_GS"
plot_out <- file.path(plot_root, "plots", "Figure5")
dir.create(plot_out, recursive = TRUE, showWarnings = FALSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggforce)
library(locuscomparer)


# Plot ldsc enrichment.
ldsc_enrichment_data <- readRDS(file.path(plot_root, "plot_data", "Figure5_ldsc_enrichment.rds"))
ggplot() +
  geom_point(data = ldsc_enrichment_data,
             aes(x = cell_type, y = trait,
                 size = Enrichment,
                 fill = -log10(Enrichment_p)),
             shape = 21,
             color = "black",
             stroke = 0.5) +
ggtitle("Significant GWAS heritability enrichment across cell types")+
  labs(
    x = NULL,
    y = NULL
  )+scale_size(range = c(3, 8))+
  scale_fill_gradient( low = "#4A90E2",
                       high = "#E94E77") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12,angle = 90, hjust = 1,colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),legend.position = "bottom",
        panel.grid = element_blank(),plot.title = element_text(
          hjust = 0.5
        ))

ggsave(file.path(plot_out, "GWAS_IAS_heatmap_final.pdf"), width = 10, height = 8)

# Plot colocalization genes.
coloc_gene_data <- readRDS(file.path(plot_root, "plot_data", "Figure5_colocalization_genes.rds"))
p <- ggplot() +

  geom_arc_bar(
    data = coloc_gene_data$single_disease,
    aes(
      x0 = y1,
      y0 = x1,
      r0 = 0,
      r = 0.4,
      start = 2*pi*ymin,
      end = 2*pi*ymax,
      fill = category
    ),
    color = NA
  ) +

  geom_arc_bar(
    data = coloc_gene_data$shared_diseases,
    aes(
      x0 = y1,
      y0 = x1,
      r0 = 0,
      r = 0.4,
      start = 2*pi*ymin,
      end = 2*pi*ymax,
      fill = category
    ),
    color = "black",
    size = 0.25
  ) +

  coord_fixed() +

  scale_x_continuous(
    breaks = unique(coloc_gene_data$plot_positions$y1),
    labels = unique(coloc_gene_data$plot_positions$gid)
  ) +

  scale_y_continuous(
    breaks = unique(coloc_gene_data$plot_positions$x1),
    labels = unique(coloc_gene_data$plot_positions$cell_type)
  ) +

  theme_bw() +

  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 16,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 16,
      colour = "black"
    ),
    panel.grid = element_blank(),
    legend.text = element_text(
      size = 16,
      colour = "black"
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(
      hjust = 0.5
    )
  ) +

  labs(
    title = "Cell type–specific GWAS colocalization of genes",
    fill = "GWAS category"
  )

p
ggsave(file.path(plot_out, "GWAS_AS_COLOC_final_gene_lunwen_1023.pdf"),width = 16,height = 9)

# Plot GWAS and splicing QTL locus comparisons.
locuscompare_input <- file.path(plot_root, "plot_data", "Figure5_locuscompare_files")
locus_gwas_files <- list.files(locuscompare_input, pattern = "_gwas\\.txt$", full.names = TRUE)
locus_qtl_files <- list.files(locuscompare_input, pattern = "_qtl\\.txt$", full.names = TRUE)
locus_prefix <- function(filename) sub("_(gwas|qtl)\\.txt$", "", basename(filename))
for (locus_name in intersect(vapply(locus_gwas_files, locus_prefix, ""), vapply(locus_qtl_files, locus_prefix, ""))) {
  gwas_file <- locus_gwas_files[vapply(locus_gwas_files, locus_prefix, "") == locus_name][1]
  qtl_file <- locus_qtl_files[vapply(locus_qtl_files, locus_prefix, "") == locus_name][1]
  locus_tokens <- strsplit(locus_name, "_", fixed=TRUE)[[1]]
  if (length(locus_tokens) < 3) next
  pdf(file.path(plot_out, paste0(locus_name, "_locuscompare.pdf")), width = 7, height = 5)
  print(locuscomparer::locuscompare(in_fn1=gwas_file, in_fn2=qtl_file, combine=TRUE, genome="hg38", population="EAS", marker_col1="SNP", pval_col1="P", title1=locus_tokens[1], marker_col2="SNP", pval_col2="pval_gi", title2=paste(locus_tokens[-(1:2)], collapse="_")))
  dev.off()
}

# Plot colocalization counts.
coloc_counts_data <- readRDS(file.path(plot_root, "plot_data", "Figure5_colocalization_counts.rds"))
ggplot(
  coloc_counts_data,
  aes(
    x = cell_type,
    y = cat,
    size = n,
    fill = ratio
  )
) +
  geom_point(
    shape = 23,
    color = "black",
    stroke = 0.5,
    alpha = 0.9
  ) +
  scale_size(range = c(2, 12)) +
  scale_fill_gradient(
    low = "#deebf7",
    high = "#08519c",
    name = "Ratio (%)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=12,angle = 90, hjust = 1,colour = "black"),
    axis.text.y = element_text(size=12,colour = "black"),legend.position = "bottom",
    panel.grid.major = element_line(color = "grey90"),plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Count",
    title = "Distribution of GWAS-sQTL colocalized signals"
  )

ggsave(file.path(plot_out, "GWAS_AS_COLOC_final.pdf"), width = 12, height = 10)
