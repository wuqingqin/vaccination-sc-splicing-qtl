plot_root <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/PLOT_GS"
plot_out <- file.path(plot_root, "plots", "Figure4")
dir.create(plot_out, recursive = TRUE, showWarnings = FALSE)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(patchwork)
library(ggrepel)
library(VennDiagram)
library(grid)


# Plot iqtl event counts.
iqtl_counts_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_iqtl_event_counts.rds"))
ggplot(
  iqtl_counts_data,
  aes(x = reorder(cell_type, -row_count), y = row_count)
) +
  geom_col(fill ="#3C5488FF", alpha = 0.8) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab(NULL) +
  ylab("Number of significant iQTL-associated sEvents count") +
  ggtitle("Significant iQTL-associated sEvents count per cell type")

ggsave(file.path(plot_out, "cell_type_SIQTL_events.pdf"),height = 5,width = 10)

# Plot iqtl event sharing.
iqtl_sharing_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_iqtl_event_sharing.rds"))
ggplot(
  iqtl_sharing_data,
  aes(
    x = n_celltype,
    y = n_event
  )
) +
  geom_segment(
    aes(
      x = n_celltype,
      xend = n_celltype,
      y = 1,
      yend = n_event
    ),
    color = "grey70"
  ) +
  geom_point(
    shape=21,
    fill = "#3C5488FF",
    size = 5
  ) +

  ggtitle("Sharing pattern of iQTL event across cell types")+
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey50"
  )+labs(
    x = "Number of shared cell types",
    y = "Number of iQTL event"
  )+
  scale_y_log10() +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()

  )
ggsave(file.path(plot_out, "event\u2013lead SNP pairs_share.pdf"),height =4,width = 3)

# Plot iqtl go enrichment.
iqtl_go_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_iqtl_go_enrichment.rds"))
ggplot(iqtl_go_data, aes(x = -log10(pvalue), y = reorder(Description, -log10(pvalue)))) +
  geom_col(aes(fill = -log10(pvalue)), width = 0.7, alpha = 0.8) +

  scale_fill_gradient(low = "#deebf7", high = "#E94E77") +

  geom_text(
    aes(label = Description,
        hjust = ifelse(-log10(pvalue) > 1, 0, -0.05),
        x = ifelse(-log10(pvalue) > 1, 0.2, NA),
        color = ifelse(-log10(pvalue) > 1, "black", "black")),
    size = 5,
    lineheight = 1.2
  ) +ggtitle("Top 15 of GO Biological Process")+
  scale_color_identity() +
  labs(x = "-log10(pvalue)", y = "", fill = "-log10(pvalue)") +
  theme_bw() +
  theme(
    legend.text = element_text(
      size = 12,
      colour = "black"
    ),
    legend.position  = "right",
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 12,
      colour = "black"
    ),
    axis.text.x = element_text(size=12,colour = "black"),
    axis.title.x = element_text(size=12,colour = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()  )

ggsave(file.path(plot_out, "iqtl_bp_df_list_unique_top15.pdf"), width = 10, height = 10)

# Plot genotype-stratified PSI panels.
genotype_panel_plots <- readRDS(file.path(plot_root, "plot_data", "Figure4_genotype_panel_plots.rds"))
pdf(file.path(plot_out, "merged_plots_2.pdf"), width = 7.5, height = 5)
for (page_start in seq(1, length(genotype_panel_plots), by = 12)) {
  page_end <- min(page_start + 11, length(genotype_panel_plots))
  print(patchwork::wrap_plots(genotype_panel_plots[page_start:page_end], ncol = 3))
}
dev.off()

# Plot sf enrichment.
iqtl_sf_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_sf_enrichment.rds"))
ggplot(iqtl_sf_data,
       aes(x = Feature,
           y = celltype)) +

  geom_point(
    aes(color = logOR),size=5
  ) +

  geom_text(
    aes(label = signif),
    size = 4,
    vjust = 0.8
  ) +

  scale_color_gradient2(
    low = "#4A90E2",
    mid = "white",
    high = "#E94E77",
    midpoint = 0,
    name = "logOR"
  ) +

  theme_bw() +

  labs(
    x = NULL,
    y = NULL,

    title  =  "Enrichment of iQTL variants in SF binding sites"
  ) + coord_flip()+

  theme(
    axis.text.x = element_text(angle = 90,size = 12,hjust =1, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank()

  )
ggsave(file.path(plot_out, "iQTL_variants_in_SF_binding_sites_1.pdf"),height = 5,width = 8)

# Plot atac module enrichment.
iqtl_atac_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_atac_module_enrichment.rds"))
p <- ggplot(
  iqtl_atac_data,
  aes(
    x = Feature,
    y = celltype,
    fill = logOR
  )
) +
  geom_tile(color = "white") +

  scale_fill_gradient(
    low = "#FCE8EE",
    high = "#E94E77",
    name = "logOR"
  ) +

  theme_bw() +

  labs(
    x = NULL,
    y = NULL,
    title = "Enrichment of ATAC peak across cell types"
  ) +

  coord_flip() +

  theme(
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    legend.text = element_text(
      size = 12,
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 12,
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 1
    ),
    plot.title = element_text(
      hjust = 0.5
    )
  )

print(p)

ggsave(file.path(plot_out, "iQTL_variants_in_ATACpeak_some.pdf"),height = 5,width = 8)

# Plot finemap event counts.
finemapping_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_finemap_event_counts.rds"))
ggplot(finemapping_data, aes(x = cell_type, y = n, fill = cs_group)) +
  geom_bar(stat = "identity",alpha=0.8) +
  scale_fill_manual(values = c("#4A90E2","#E94E77"))+
  theme_bw() +
  labs(
    x = "",
    y = "Number of sQTL events",
    fill = "Independent signals (CS)",
    title = "Distribution of independent signals across cell types"
  ) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(size=12,color = "black",angle = 90,hjust = 1),
    axis.text = element_text(size=12,color = "black"),
    plot.title = element_text(color = "black", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(plot_out, "splicing_qtls_final.pdf"), width = 6, height = 5)

# Plot iqtl cell abundance.
iqtl_cell_number_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_iqtl_cell_abundance.rds"))
p <- ggplot(
  iqtl_cell_number_data$normalized_iqtl_counts,
  aes(
    x = cell_type,
    y = Value,
    color = Type,
    group = Type
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +

  scale_color_manual(
    values = c(
      Cell_number = "grey50",
      iQTL_events = "#3C5488FF"
    ),
    labels = c(
      "Cell number",
      "iQTL-associated events"
    )
  ) +geom_text_repel(
    data = iqtl_cell_number_data$label_data,
    aes(
      x = cell_type,
      y = iQTL_events,
      label = cell_type
    ),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf,
    inherit.aes = FALSE
  )+

  labs(
    x = "",
    y = "Normalized value",
    color = "",
    title = "Comparison of cell abundance and iQTL-associated sEvents count"
  ) +

  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 12,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    axis.title = element_text(
      size = 15,
      colour = "black"
    ),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(
      size = 15,
      hjust = 0.5
    )
  )

p

ggsave(
  file.path(plot_out, "isevent_cellnum_relationship.pdf"),
  p,
  width = 7,
  height = 4.8
)

# Plot iqtl sensitivity.
iqtl_sensitivity_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_iqtl_sensitivity.rds"))
p <- ggplot(
  iqtl_sensitivity_data,
  aes(
    x = b_gi,
    y = delta_beta
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.4
  ) +
  geom_point(
    aes(fill = category),
    shape = 21,
    size = 2.5,
    colour = "black",
    stroke = 0.25,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c(
      "Consistent & nominal P < 0.05" = "#3C5488FF",
      "Consistent only" = "#F39B7FFF",
      "Not directionally consistent" = "grey75"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Concordance of iQTL and delta PSI-based effects",
    x = "Original iQTL interaction effect",
    y = expression("Delta PSI-based genotype effect"),
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(colour = "black", size = 12),
    plot.subtitle = element_text(colour = "black", size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(
  file.path(plot_out, "iQTL_deltaPSI_quadrant.pdf"),
  plot = p,
  width = 5,
  height = 5
)

# Plot the overlap of differential splicing and iQTL-associated events.
das_iqtl_overlap_counts <- readRDS(file.path(plot_root, "plot_data", "Figure4_DAS_iQTL_Venn_counts.rds"))
n_das <- das_iqtl_overlap_counts$n_das
n_iqtl <- das_iqtl_overlap_counts$n_iqtl
n_overlap <- das_iqtl_overlap_counts$n_overlap
venn.plot <- draw.pairwise.venn(
  area1 = n_das,
  area2 = n_iqtl,
  cross.area = n_overlap,
  category = c("DAS events", "iQTL-associated sEvents"),
  fill = c("#3C5488FF", "#F39B7FFF"),
  alpha = c(0.6, 0.6),
  scaled = TRUE,
  euler.d = TRUE,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = c(0.06, 0.06),
  lty = "blank",
  ind = FALSE
)
draw_venn <- function() {

  grid.newpage()

  pushViewport(
    viewport(
      x = 0.5,
      y = 0.44,
      width = 0.88,
      height = 0.78
    )
  )

  grid.draw(venn.plot)

  popViewport()

  grid.text(
    "Overlap between DASEs and iQTL-associated sEvents",
    x = 0.5,
    y = 0.95,
    gp = gpar(
      fontsize = 13,
      fontface = "bold",
      col = "black"
    )
  )
}

pdf(file.path(plot_out, "DAS_iQTL_Venn.pdf"), width = 7, height = 7)
draw_venn()
dev.off()

# Plot matched depth concordance.
matched_depth_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_matched_depth_concordance.rds"))
p <- ggplot(
  matched_depth_data,
  aes(
    x = Mean_Pearson_r,
    y = celltype
  )
) +
  geom_point(
    size = 3,
    color = "#3C5488FF",
    alpha = 0.8
  ) +geom_text(
    aes(label = sprintf("%.3f", Mean_Pearson_r)),
    nudge_x = 0.015,
    size = 4,
    color = "black"
  ) +
  labs(
    title = "Robustness of DAS effect sizes to matched-depth downsampling",
    x = "Pearson's r",
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(color = "black")
  )

p

ggsave(
  file.path(plot_out, "reads_depth_sensitivity.pdf"),
  p,
  height = 10,
  width = 6
)

# Plot genotype-stratified delta PSI validation.
delta_psi_genotype_plot <- readRDS(file.path(plot_root, "plot_data", "Figure4_delta_PSI_genotype_plot.rds"))
ggsave(file.path(plot_out, "delta_PSI_validation_PRMT5_ZNF226_IFI6.pdf"), plot = delta_psi_genotype_plot, width = 12, height = 4, device = cairo_pdf)

# Plot iqtl adt heatmap.
iqtl_adt_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_iqtl_adt_heatmap.rds"))
p_all <- ggplot(
  iqtl_adt_data$correlation_table,
  aes(
    x = event_label,
    y = marker_label,
    fill = rho
  )
) +

  geom_tile(
    color = "white",
    linewidth = 0.4,
    width = 0.95,
    height = 0.95
  ) +

  geom_text(
    aes(label = sig_label),
    size = 5,
    color = "black"
  ) +

  facet_grid(
    category ~ .,
    scales = "free_y",
    space = "free_y",
    drop = TRUE
  ) +

  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -iqtl_adt_data$color_limit,
      iqtl_adt_data$color_limit
    ),
    name = "Spearman's r",
    na.value = "grey90"
  ) +

  scale_x_discrete(
    drop = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL,
    title = "Correlations between ΔPSI and ΔADT"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    panel.grid = element_blank(),

    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),

    strip.background = element_blank(),

    strip.text.y = element_text(
      size = 11,
      color = "black",
      angle = 90
    ),

    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),

    axis.title = element_text(
      color = "black"
    ),

    legend.title = element_text(
      size = 10,
      color = "black"
    ),

    legend.text = element_text(
      size = 9,
      color = "black"
    ),

    plot.caption = element_text(
      size = 9,
      color = "black",
      hjust = 0
    )
  )

print(p_all)

ggsave(
  filename = file.path(plot_out, 
    file.path(plot_out, "Three_events_ADT_function_heatmap.pdf")
  ),
  plot = p_all,
  width = 5,
  height = iqtl_adt_data$figure_height,
  limitsize = FALSE
)

# Plot iqtl antibody heatmap.
iqtl_antibody_data <- readRDS(file.path(plot_root, "plot_data", "Figure4_iqtl_antibody_heatmap.rds"))
ggplot(iqtl_antibody_data, aes(x = label, y = phenotype)) +
  geom_point(aes(fill = rho, size = abs(rho)), shape = 21, stroke = 1) +
  geom_text(aes(label = sig), color = "black", size = 7, vjust = 0.8) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                        midpoint = 0, limits = c(-0.3, 0.3)) +
  scale_size_continuous(range = c(3,10), limits = c(0,0.3)) +
  labs(title = "Correlations between delat PSI and delta antibody",
    x = NULL,
       y = NULL,
       color = "Spearman rho",
       size = "|rho|") +
  theme_bw() +
  theme(plot.title = element_text(size = 13, color = "black", hjust = 0.5),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        panel.grid = element_blank())

ggsave(file.path(plot_out, "three_event_antibody.pdf"),height = 5,width = 5)
