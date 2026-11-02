plot_root <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/PLOT_GS"
plot_out <- file.path(plot_root, "plots", "Figure2")
dir.create(plot_out, recursive = TRUE, showWarnings = FALSE)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(patchwork)
library(ggrepel)


# Plot significant sqtl counts.
sqtl_counts_data <- readRDS(file.path(plot_root, "plot_data", "Figure2_significant_sqtl_counts.rds"))
ggplot(sqtl_counts_data, aes(x = reorder(cell_type, -count), y = count, fill = time)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), alpha = 0.8) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8), vjust = -0.4,  size = 3.5) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("#3C5488FF","#F39B7FFF","#D9C0A0"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab(NULL) +
  ylab("Number of Significant Sevent") +
  ggtitle("Splicing events with significant sQTLs across cell types")
ggsave(file.path(plot_out, "Significant_Sevents_Cell_Type.pdf"), width = 13, height = 7)

# Plot sqtl lead snp counts.
sqtl_lead_snp_data <- readRDS(file.path(plot_root, "plot_data", "Figure2_sqtl_lead_snp_counts.rds"))
ggplot(sqtl_lead_snp_data,
       aes(x = cell_type,
           y = proportion,
           fill = class)) +
  geom_col(
    width = 0.8,
    colour = "black",
    linewidth = 0.3
  ) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(
    "Pre-vacc. only" = "#3C5488FF",
    "Post-vacc. only" = "#F39B7FFF",
    "Pre and Post-vacc." = "#91D1C2FF"
  )) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL,
    title = "Composition of shared and time-specific event across cell types"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.8
    ),
    axis.text.x = element_text(
      size = 14,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    legend.text = element_text(
      size = 14,
      colour = "black"
    ),
    legend.position  = "bottom",
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      colour = "black"
    )
  )

ggsave(file.path(plot_out, "Significant_Sevents_leadsnp_Cell_Type.pdf"), width = 8, height = 8)

# Plot sqtl event sharing.
sqtl_sharing_data <- readRDS(file.path(plot_root, "plot_data", "Figure2_sqtl_event_sharing.rds"))
ggplot(sqtl_sharing_data,
       aes(x = n_celltype,
           y = log10(n_event_snp),
           color = group_label,
           group = group_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 4) +
  scale_x_continuous(
    breaks = seq(
      min(sqtl_sharing_data$n_celltype),
      max(sqtl_sharing_data$n_celltype),
      by = 2
    )
  ) +labs(
    title = "Cell type sharing patterns of Sevent",
    x = "Number of shared cell types",
    y = expression(log[10]*"(Number of Sevent)")
  )+
  scale_color_manual(values = c(
    "Pre-vacc. only" = "#3C5488FF",
    "Post-vacc. only" = "#F39B7FFF",
    "Pre and Post-vacc." = "#91D1C2FF"
  )) +
  theme_bw()+ theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.8
    ),
    axis.text.x = element_text(
      size = 14,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    legend.text = element_text(
      size = 14,
      colour = "black"
    ),
    legend.position  = "bottom",
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      colour = "black"
    )
  )

ggsave(file.path(plot_out, "Significant_Sevents_leadsnp_Cell_Type_share.pdf"), width = 5, height = 6)

# Plot sqtl go enrichment.
sqtl_go_data <- readRDS(file.path(plot_root, "plot_data", "Figure2_sqtl_go_enrichment.rds"))
ggplot(
  sqtl_go_data,
  aes(
    x=class,
    y=Description,
    size=Count,
    color=-log10(p.adjust)
  )
)+
  geom_point()+ggtitle("Top 5 of GO Biological Process")+
  scale_color_gradient(low = "#4A90E2", high = "#E94E77")+
  theme_bw()+scale_size(range = c(5, 12)) +
  labs(
    x=NULL,
    y=NULL,
    color="-log10(adj.P)",
    size="Gene count"
  )+
  theme_bw()+ theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.8
    ),
    axis.text.x = element_text(
      size = 14,angle = 90,hjust=1,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    legend.text = element_text(
      size = 14,
      colour = "black"
    ),
    legend.position  = "right",
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      colour = "black"
    )
  )
ggsave(file.path(plot_out, "Significant_Sevents_leadsnp_Cell_Type_share_BP_plot.pdf"), width = 10, height = 10)

# Plot sqtl effect concordance.
sqtl_effects_data <- readRDS(file.path(plot_root, "plot_data", "Figure2_sqtl_effect_concordance.rds"))
ggplot(sqtl_effects_data$matched_effects,
       aes(x = e1_beta, y = e2_beta)) +

  geom_point(alpha = 0.4, size = 1) +

  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "black") +

  geom_smooth(method = "lm", se = FALSE,
              color = "red", linewidth = 0.8) +

  theme_classic() +

  labs(
    x = "Effect size (e1)",
    y = "Effect size (e2)",
    title = "Conservation and reprogramming of sQTL effects across immune stimulation"
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  annotate(
    "text",
    x = min(sqtl_effects_data$matched_effects$e1_beta, na.rm = TRUE),
    y = max(sqtl_effects_data$matched_effects$e2_beta, na.rm = TRUE),
    hjust = 0,
    vjust = 1,
    size = 4,
    label = paste0(
      "Spearman r = ", round(sqtl_effects_data$spearman_r, 3), "\n",
      "Slope = ", round(sqtl_effects_data$slope, 3), "\n",
      "Intercept = ", round(sqtl_effects_data$intercept, 3), "\n",
      "Direction discordance = ", round(sqtl_effects_data$discordance_rate * 100, 2), "%"
    )
  )+
  theme_bw(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(plot_out, "shift_of_sQTL_effect_size_1.pdf"), width = 7, height = 8)

# Plot genotype-stratified PSI panels.
genotype_panel_plots <- readRDS(file.path(plot_root, "plot_data", "Figure2_genotype_panel_plots.rds"))
pdf(file.path(plot_out, "merged_plots_1.pdf"), width = 15, height = 5)
for (page_start in seq(1, length(genotype_panel_plots), by = 12)) {
  page_end <- min(page_start + 11, length(genotype_panel_plots))
  print(patchwork::wrap_plots(genotype_panel_plots[page_start:page_end], ncol = 3))
}
dev.off()

# Plot celltype normalized comparison.
sqtl_cell_number_data <- readRDS(file.path(plot_root, "plot_data", "Figure2_celltype_normalized_comparison.rds"))
ggplot(sqtl_cell_number_data$normalized_celltype_data,
       aes(x = cell_type,
           y = Value,
           color = Type,
           group = Type)) +

  geom_line(linewidth = 1.2) +
  geom_text_repel(
    data = sqtl_cell_number_data$label_data,
    aes(
      x = cell_type,
      y = DAS_events,
      label = cell_type
    ),
    inherit.aes = FALSE,
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  geom_point(size = 3) +

  facet_wrap(~time, scales = "free_x") +

  scale_color_manual(
    values = c(
      Cell_number = "grey50",
      DAS_events = "#3C5488FF"
    ),
    labels = c(
      "Cell number",
      "Differential AS events"
    )
  ) +

  labs(
    x = "",
    y = "Normalized value",
    color = "",
    title = "Differential AS events are not completely explained by cell abundance"
  ) +

  theme_bw() +

  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 10,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    axis.title = element_text(
      size = 14,
      colour = "black"
    ),
    strip.text = element_text(
      size = 14,
      face = "bold"
    ),
    plot.title = element_text(
      size = 15,
      hjust = 0.5
    ),
    legend.position = "top"
  )
ggsave(file.path(plot_out, "Significant_Sevents_leadsnp_Cell_Type_cc.pdf"), width = 8, height = 8)

# Plot IFI6 sensitivity analysis excluding the rare homozygous genotype.
ifi6_excluding_tt_plot <- readRDS(file.path(plot_root, "plot_data", "Figure2_ifi6_excluding_TT_plot.rds"))
ggsave(file.path(plot_out, "IFI6_rs117237130_exclude_TT_e1_e2_iQTL.pdf"), plot = ifi6_excluding_tt_plot, width = 12, height = 4, device = cairo_pdf)
