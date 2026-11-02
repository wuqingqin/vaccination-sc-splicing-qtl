plot_root <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/PLOT_GS"
plot_out <- file.path(plot_root, "plots", "Figure3")
dir.create(plot_out, recursive = TRUE, showWarnings = FALSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(ComplexHeatmap)
library(circlize)
library(grid)


# Plot sf enrichment categories.
sf_enrichment_categories_data <- readRDS(file.path(plot_root, "plot_data", "Figure3_sf_enrichment_categories.rds"))
ggplot(sf_enrichment_categories_data, aes(x = celltype, y = count, fill = type)) +
  geom_col() +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  ggtitle("Number of significantly enriched splicing factors") +
  labs(
    x = NULL,
    y = "Number of significant SFs",
    fill = "Experiment"
  ) +
  scale_fill_manual(
    values = c(
      "Pre-vacc. only" = "#3C5488FF",
      "Pre- and Post-vacc." = "#91D1C2FF",
      "Post-vacc. only" = "#F39B7FFF"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      size = 12,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )
ggsave(file.path(plot_out, "iQTL_variants_in_SF_binding_sites_gg_1005.pdf"),height = 5,width = 10)

# Plot sf representative.
sf_examples_data <- readRDS(file.path(plot_root, "plot_data", "Figure3_sf_representative.rds"))
ggplot(sf_examples_data, aes(x = celltype, y = Feature)) +
  geom_point(aes(size = logOR, color = time)) +
  scale_color_manual(values = c("Pre-vacc." = "#3C5488FF", "Post-vacc." = "#F39B7FFF")) +
  scale_size_continuous(range = c(1, 8), name = "Log2OR") +
  labs(x = NULL, y = "Splicing factor", color = "Time point",
       title = "Enrichment scores of significant splicing factors") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())

ggsave(file.path(plot_out, "iQTL_variants_in_SF_binding_sites_gg_example.pdf"),height = 5,width = 4)

# Plot sf composition pies.
sf_pie_data <- readRDS(file.path(plot_root, "plot_data", "Figure3_sf_composition_pies.rds"))
ggplot() +

  geom_arc_bar(
    data = sf_pie_data$single_sf,
    aes(
      x0 = x,
      y0 = y,
      r0 = 0,
      r = 0.4,
      start = 2*pi*ymin,
      end = 2*pi*ymax,
      fill = Feature
    ),
    color = NA
  ) +

  geom_arc_bar(
    data = sf_pie_data$shared_sf,
    aes(
      x0 = x,
      y0 = y,
      r0 = 0,
      r = 0.4,
      start = 2*pi*ymin,
      end = 2*pi*ymax,
      fill = Feature
    ),
    color = "black",
    linewidth = 0.25
  ) +

  coord_fixed() +

  scale_x_continuous(
    breaks = unique(sf_pie_data$plot_positions$x),
    labels = unique(sf_pie_data$plot_positions$condition)
  ) +

  scale_y_continuous(
    breaks = unique(sf_pie_data$plot_positions$y),
    labels = unique(sf_pie_data$plot_positions$celltype)
  ) +

  theme_bw() +

  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size=14,
      colour="black"
    ),
    axis.text.y = element_text(
      size=14,
      colour="black"
    ),
    axis.title = element_blank(),
    legend.position="bottom"
  ) +

  labs(
    fill="Feature"
  )

ggsave(file.path(plot_out, "est_all_plot_some_ball.pdf"),width =4,height = 12)

# Plot sf estimates ci.
sf_ci_data <- readRDS(file.path(plot_root, "plot_data", "Figure3_sf_estimates_ci.rds"))
p <- ggplot(sf_ci_data, aes(x = Feature, y = logOR,
                         group = condition_label, color = condition_label)) +
  geom_line(size = 1) +
  geom_point(size = 2.5) +

  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  facet_wrap(~celltype, scales = "free_x") +
  scale_color_manual(values = c("Pre-vacc." = "#3C5488FF",
                                "Post-vacc." = "#F39B7FFF")) +
  labs(x = "Feature", y = "logOR (Enrichment score)", color = "Condition") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size =14),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 14), legend.title = element_text(size = 14)
  )

print(p)

ggsave(file.path(plot_out, "est_all_plot_some.pdf"),width =12,height = 6)

# Plot gwas triangular heatmap.
gwas_triangular_heatmap_data <- readRDS(file.path(plot_root, "plot_data", "Figure3_gwas_triangular_heatmap.rds"))
mat_e1 <- gwas_triangular_heatmap_data$mat_e1
mat_e2 <- gwas_triangular_heatmap_data$mat_e2
col_fun_e1 <- gwas_triangular_heatmap_data$col_fun_e1
col_fun_e2 <- gwas_triangular_heatmap_data$col_fun_e2
pdf(file.path(plot_out, "GWAS_AS_heatmap_final.pdf"), width = 12, height = 10)
ComplexHeatmap::draw(gwas_triangular_heatmap_data$heatmap, annotation_legend_list = gwas_triangular_heatmap_data$legends, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

# Plot sf triangular heatmap.
sf_triangular_heatmap_data <- readRDS(file.path(plot_root, "plot_data", "Figure3_sf_triangular_heatmap.rds"))
mat_e1 <- sf_triangular_heatmap_data$mat_e1
mat_e2 <- sf_triangular_heatmap_data$mat_e2
col_fun_e1 <- sf_triangular_heatmap_data$col_fun_e1
col_fun_e2 <- sf_triangular_heatmap_data$col_fun_e2
pdf(file.path(plot_out, "RBP_AS_heatmap_final.pdf"), width = 12, height = 10)
ComplexHeatmap::draw(sf_triangular_heatmap_data$heatmap, annotation_legend_list = sf_triangular_heatmap_data$legends, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()

# Plot sf enrichment heatmap.
sf_heatmap_data <- readRDS(file.path(plot_root, "plot_data", "Figure3_sf_enrichment_heatmap.rds"))
ggplot(sf_heatmap_data,
       aes(x = Feature,
           y = celltype,
           fill = logOR)) +
  geom_tile(color = "grey80", linewidth = 0.3) +
  facet_wrap(~ condition) +
  scale_fill_gradient2(
    low = "#4A90E2",
    mid = "white",
    high = "#E94E77",
    midpoint = 0
  ) +labs(
    x = NULL,
    y = NULL,
    fill = "log2(OR)",
    title = "Significantly enriched sQTL variants in accessible chromatin regions"
  )+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5,size=12),

    axis.text.y = element_text(size=12,colour = "black"),
    axis.text.x = element_text(
      colour = "black",size=12,
      angle = 45,
      hjust = 1
    ),

    panel.grid = element_blank(),

    strip.background = element_blank(),
    strip.text = element_text(
      colour = "black",
      face = "plain"
    )
  )
ggsave(file.path(plot_out, "est_all_plot.pdf"),width =12,height = 8)
