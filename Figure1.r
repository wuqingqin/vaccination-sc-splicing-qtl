plot_root <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/PLOT_GS"
plot_out <- file.path(plot_root, "plots", "Figure1")
dir.create(plot_out, recursive = TRUE, showWarnings = FALSE)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(patchwork)
library(ggrepel)


# Plot as event distribution.
as_event_distribution_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_as_event_distribution.rds"))
ggplot(as_event_distribution_data, aes(x = celltype2, y = count,fill = time)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +scale_fill_manual(values = c("#3C5488FF","#F39B7FFF"))+
  labs(
    x = NULL,
    y = "Number of AS events per sample"
  ) +ggtitle("Distribution of AS events per sample")+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(plot_out, "all_celltypes_events_heatmap.pdf"),height = 5,width = 10)

# Plot as events vs cells.
as_event_cell_correlation_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_as_events_vs_cells.rds"))
ggplot(as_event_cell_correlation_data$event_counts, aes(n_cells, count, colour = time)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_bw()+scale_colour_manual(values = c("#3C5488FF","#F39B7FFF"))+
  geom_text(
    data = as_event_cell_correlation_data$correlation_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.5
  ) +
  labs(
    title = "AS event detection increases with pseudobulk cell number",
    x = "Number of cells",
    y = "AS events detected"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text( colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(plot_out, "all_celltypes_coreelation_ASevents.pdf"),height = 5,width = 5)

# Plot dase counts.
dase_counts_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_dase_counts.rds"))
ggplot(dase_counts_data, aes(Freq, reorder(Var1, Freq))) +
  geom_segment(
    aes(x = 0, xend = Freq, y = Var1, yend = Var1),
    color = "grey70"
  ) +
  geom_point(size = 5, shape = 21, fill = "lightblue", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  labs(
    x = "Number of DAS events",
    y = NULL,
    title = "Distribution of differential AS events across cell types"
  ) +
  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      color = "black",
      family = "serif",
      size = 12
    ),

    axis.text.y = element_text(
      color = "black",
      family = "serif",
      size = 12
    ),

    axis.title.x = element_text(
      color = "black",
      family = "serif",
      size = 12
    ),

    axis.title.y = element_text(
      color = "black",
      family = "serif",
      size = 12
    ),

    plot.subtitle = element_text(
      family = "serif",
      size = 12
    )
  )
ggsave(file.path(plot_out, "all_celltypes_DASES_plot_0916.pdf"),height = 4,width = 8)

# Plot dase cell sharing.
dase_sharing_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_dase_cell_sharing.rds"))
p_share <- ggplot(
  dase_sharing_data,
  aes(x = factor(n_celltype), y = n_events)
) +
  geom_col(
    width = 0.75,
    fill = "#8499C4",
    color = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = n_events),
    vjust = -0.35,
    size = 4,
    family = "serif"
  ) +
  theme_bw() +
  labs(
    x = "Number of cell types",
    y = "Number of DAS events",
    title = "Cell type sharing of DAS events"
  ) +
  theme(
    panel.grid = element_blank(),

    axis.text = element_text(
      color = "black",
      family = "serif",
      size = 12
    ),

    axis.title = element_text(
      color = "black",
      family = "serif",
      size = 12
    ),

    plot.title = element_text(
      family = "serif",
      size = 12,
      hjust = 0.5
    )
  )

p_share
ggsave(file.path(plot_out, "all_celltypes_DASES_share_0916.pdf"),height = 4,width = 4)

# Plot ZAP70 and SSH2 splice-junction coverage.
splice_coverage_plots <- readRDS(file.path(plot_root, "plot_data", "Figure1_splice_coverage_objects.rds"))
ggsave(file.path(plot_out, "ZAP70_coverage_final.pdf"), plot = splice_coverage_plots$ZAP70, width = 10, height = 5, useDingbats = FALSE)
ggsave(file.path(plot_out, "SSH2_coverage_final.pdf"), plot = splice_coverage_plots$SSH2, width = 10, height = 6.5, useDingbats = FALSE)

# Plot sex concordance.
sex_concordance_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_sex_concordance.rds"))
ggplot(sex_concordance_data, aes(x = recorded, y = inferred, color = qc_flag)) +
  geom_jitter(width = 0.05, height = 0.05, size = 2.5) +

  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +

  scale_x_continuous(breaks = c(0,1), labels = c("Female", "Male")) +
  scale_y_continuous(breaks = c(0,1), labels = c("Female", "Male")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62"))+
  labs(
    x = "Recorded Gender",
    y = "Genotype-inferred Sex",
    color = "QC Status",
    title = "Sample Sex Concordance Check"
  ) +
  theme_bw()+theme(axis.text.x = element_text(size = 12, colour = "black"),
                   axis.text.y = element_text(size = 12, colour = "black"),
                   axis.title.y = element_text(size = 12),
                   legend.text = element_text(size = 10),
                   legend.title = element_text(size = 10),
                   plot.title = element_text(hjust = 0.5),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.box = "horizontal")

ggsave(file.path(plot_out, "sex_QC.pdf"),width = 5,height = 2.5)

# Plot heterozygosity.
heterozygosity_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_heterozygosity.rds"))
ggplot(heterozygosity_data, aes(x = het_rate)) +
  geom_histogram(
    bins = 40,
    fill = "#FC8D62",
    color = "black"
  ) +
  geom_vline(xintercept = mu, linewidth = 1) +
  geom_vline(
    xintercept = c(lower, upper),
    linetype = "dashed",
    color = "red",
    linewidth = 1
  ) +
  theme_classic() +
  labs(
    title = "Sample heterozygosity",
    subtitle = paste0(
      "Mean = ", round(mu, 4),
      ", ±3 SD = [", round(lower, 4), ", ", round(upper, 4), "]"
    ),
    x = "Heterozygosity rate",
    y = "Number of samples"
  )+theme(axis.text.x = element_text(size = 12, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.y = element_text(size = 12))
ggsave(file.path(plot_out, "het_QC.pdf"),width = 8,height = 5)

# Plot participant demographics.
participant_demographics_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_participant_demographics.rds"))
p_gender <- ggplot(participant_demographics_data, aes(x = gender, fill = gender)) +
  geom_bar(width = 0.65) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.4,
    size = 5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of individuals") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

p_age <- ggplot(participant_demographics_data, aes(x = age)) +
  geom_histogram(
    binwidth = 5,
    boundary = 0,
    color = "white",
    fill = "#4DBBD5"
  ) +
  labs(x = "Age (years)", y = "Number of individuals") +
  theme_classic(base_size = 14)

p_gender + p_age

ggsave(
  file.path(plot_out, "age_gender_distribution.pdf"),
  p_gender + p_age,
  width = 9,
  height = 4
)

# Plot RNA and protein marker expression.
rna_adt_marker_plots <- readRDS(file.path(plot_root, "plot_data", "Figure1_RNA_ADT_marker_plots.rds"))
ggsave(file.path(plot_out, "AS_DOTPLOT.pdf"), plot = rna_adt_marker_plots$RNA + rna_adt_marker_plots$ADT, width = 10, height = 6)

# Plot celltype counts.
celltype_counts_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_celltype_counts.rds"))
ggplot(
  celltype_counts_data,
  aes(
    x = Var2,
    y = Freq,
    fill = Var1
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "before" = "#3C5488FF",
      "after" = "#F39B7FFF"
    )
  ) +
  labs(
    title = "Number of cells identified for each cell type",
    x = NULL,
    y = "Number of cells",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 16,
      colour = "black"
    ),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 16, colour = "black"),
    legend.position = "top",
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
ggsave(file.path(plot_out, "celltype_AS_DOTPLOT.pdf"),height = 5,width=12)

# Plot as gene distribution.
as_gene_distribution_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_as_gene_distribution.rds"))
ggplot(as_gene_distribution_data, aes(x = celltype2, y = n_AS_genes,fill = time)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +scale_fill_manual(values = c("#3C5488FF","#F39B7FFF"))+
  labs(
    x = NULL,
    y = "Number of AS genes per sample"
  ) +ggtitle("Distribution of AS genes per sample")+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(plot_out, "all_celltypes_heatmap.pdf"),height = 5,width = 10)

# Plot as genes vs cells.
as_gene_cell_correlation_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_as_genes_vs_cells.rds"))
ggplot(as_gene_cell_correlation_data$as_gene_counts, aes(n_cells, n_AS_genes, colour = time)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_bw()+scale_colour_manual(values = c("#3C5488FF","#F39B7FFF"))+
geom_text(
  data = as_gene_cell_correlation_data$correlation_labels,
  aes(x = x, y = y, label = label),
  inherit.aes = FALSE,
  hjust = 1.1,
  vjust = 1.5
) +
  labs(
    title = "AS gene detection increases with pseudobulk cell number",
    x = "Number of cells",
    y = "AS genes detected"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text( colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(plot_out, "all_celltypes_coreelation.pdf"),height = 5,width = 5)

# Plot dase adt heatmap.
dase_adt_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_dase_adt_heatmap.rds"))
p_all <- ggplot(
  dase_adt_data$correlation_table,
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
    aes(
      label = sig_label
    ),
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
      -dase_adt_data$color_limit,
      dase_adt_data$color_limit
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
    title = "Correlations between ΔPSI and ΔADT",
    caption = dase_adt_data$caption
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
    file.path(plot_out, "DAS_ADT_function_heatmap.pdf")
  ),
  plot = p_all,
  width = 5,
  height = dase_adt_data$figure_height,
  limitsize = FALSE
)

# Plot dase antibody heatmap.
dase_antibody_data <- readRDS(file.path(plot_root, "plot_data", "Figure1_dase_antibody_heatmap.rds"))
ggplot(dase_antibody_data, aes(x = label, y = phenotype)) +
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
