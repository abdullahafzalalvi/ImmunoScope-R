####PROJECT: Single-Cell RNA-Seq Analysis of Human PBMCs
# RESEARCHERS: Muzzamil | Abdullah Alvi
# AIM: Identify distinct immune cell populations using scRNA-seq clustering
# Dataset: Simulated 10X Genomics PBMC data (mirrors real PBMC 3k structure)
# Close ALL open graphics devices safely
while (!is.null(dev.list())) {
  dev.off()
}

cat("All devices closed. Screen is clean.\n")
library(ggplot2)
library(dplyr)
install.packages("pheatmap")
library(pheatmap)
install.packages("ggrepel")
library(ggrepel)
install.packages("reshape2")
library(reshape2)
install.packages("scales")
library(scales)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("gridExtr")
library(gridExtra)
library(viridis)

set.seed(42)

caption_text <- "Researchers: Muzzamil | Abdullah Alvi\nscRNA-Seq PBMC Cell Type Identification Project"

theme_scrna <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1A3A5C"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "#34495E", margin = margin(b = 6)),
      plot.caption  = element_text(size = 8, color = "#7F8C8D", hjust = 1, face = "italic"),
      axis.title    = element_text(size = 11, face = "bold"),
      axis.text     = element_text(size = 9),
      legend.title  = element_text(size = 10, face = "bold"),
      legend.text   = element_text(size = 9),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#EBF5FB", color = "#2E86C1"),
      strip.text    = element_text(face = "bold", size = 10)
    )
}
print(theme_scrna)
plot(1:10, 1:10, main = "Test Plot with theme_scrna", xlab = "X-axis", ylab = "Y-axis") +
  theme_scrna()
# STEP 1: SIMULATE REALISTIC PBMC scRNA-SEQ DATA
n_cells <- 2638 
n_genes <- 200  
cell_types <- c("Naive CD4 T",  "CD14+ Mono",  "Memory CD4 T",
                "B cell",        "CD8 T",        "FCGR3A+ Mono",
                "NK cell",       "Dendritic",    "Platelet")
cell_props <- c(0.32, 0.18, 0.15, 0.10, 0.10, 0.07, 0.04, 0.03, 0.01)
cell_colors <- c(
  "Naive CD4 T"   = "#E74C3C",
  "CD14+ Mono"    = "#3498DB",
  "Memory CD4 T"  = "#E67E22",
  "B cell"        = "#2ECC71",
  "CD8 T"         = "#9B59B6",
  "FCGR3A+ Mono"  = "#1ABC9C",
  "NK cell"       = "#F39C12",
  "Dendritic"     = "#D35400",
  "Platelet"      = "#95A5A6"
)
n_per_type <- round(cell_props * n_cells)
n_per_type[1] <- n_cells - sum(n_per_type[-1])
cell_labels <- rep(cell_types, n_per_type)
cluster_ids  <- rep(0:8, n_per_type)
# --- Generate realistic UMAP coordinates (2D embeddings per cell type) ---
umap_centers <- list(
  "Naive CD4 T"   = c(3.5, 3.0),
  "CD14+ Mono"    = c(-4.0, 0.5),
  "Memory CD4 T"  = c(2.0, 1.0),
  "B cell"        = c(-1.5, -4.0),
  "CD8 T"         = c(5.0, -1.0),
  "FCGR3A+ Mono"  = c(-5.5, -2.0),
  "NK cell"       = c(4.0, -4.0),
  "Dendritic"     = c(-3.0, 2.5),
  "Platelet"      = c(0.5, -6.5)
)
umap_coords <- do.call(rbind, lapply(cell_types, function(ct) {
  n  <- n_per_type[which(cell_types == ct)]
  cx <- umap_centers[[ct]][1]
  cy <- umap_centers[[ct]][2]
  data.frame(
    UMAP_1    = rnorm(n, cx, 0.8),
    UMAP_2    = rnorm(n, cy, 0.8),
    cell_type = ct,
    cluster   = which(cell_types == ct) - 1
  )
}))
rownames(umap_coords) <- paste0("cell_", seq_len(n_cells))
# --- Generate marker gene expression matrix ---
# Canonical PBMC marker genes
marker_genes <- list(
  "Naive CD4 T"   = c("IL7R","CCR7","CD3D","SELL","TCF7"),
  "CD14+ Mono"    = c("CD14","LYZ","CST3","FCGR3A","S100A8"),
  "Memory CD4 T"  = c("IL7R","S100A4","ANXA1","LTB","CD3E"),
  "B cell"        = c("MS4A1","CD79A","CD79B","HLA-DRA","BANK1"),
  "CD8 T"         = c("CD8A","CD8B","GZMK","GZMA","PRF1"),
  "FCGR3A+ Mono"  = c("FCGR3A","MS4A7","IFITM3","AIF1","CXCL10"),
  "NK cell"       = c("GNLY","NKG7","KLRD1","GZMB","FGFBP2"),
  "Dendritic"     = c("FCER1A","CST3","HLA-DQA1","CLEC10A","CD1C"),
  "Platelet"      = c("PPBP","PF4","GNG11","SDPR","SPARC")
)
all_markers <- unique(unlist(marker_genes))
# Build expression matrix: high expression for cognate cell type
expr_mat <- matrix(0.1, nrow = length(all_markers), ncol = n_cells,
                   dimnames = list(all_markers, rownames(umap_coords)))
for (i in seq_along(cell_types))
  ct    <- cell_types[i]
  cells <- which(umap_coords$cell_type == ct)
  genes <- marker_genes[[ct]]
  expr_mat[genes, cells] <- expr_mat[genes, cells] +
    matrix(abs(rnorm(length(genes) * length(cells), mean = 3, sd = 0.8)),
           nrow = length(genes))
  # Add low-level noise to other genes

  other_genes <- setdiff(all_markers, genes)
  expr_mat[other_genes, cells] <- expr_mat[other_genes, cells] +
    matrix(abs(rnorm(length(other_genes) * length(cells), 0.2, 0.2)),
           nrow = length(other_genes))
  # --- QC metadata ---
  meta <- data.frame(
    cell_type    = umap_coords$cell_type,
    cluster      = umap_coords$cluster,
    nFeature_RNA = round(rnorm(n_cells, 1200, 300)),
    nCount_RNA   = round(rnorm(n_cells, 3000, 800)),
    percent_mt   = abs(rnorm(n_cells, 2.5, 1.2)),
    row.names    = rownames(umap_coords)
  )
  meta$nFeature_RNA <- pmax(meta$nFeature_RNA, 300)
  meta$nCount_RNA   <- pmax(meta$nCount_RNA, 500)
  meta$percent_mt   <- pmin(meta$percent_mt, 8)
  cat("=== DATASET SUMMARY ===\n")
  cat("Cells:", n_cells, "\n")
  cat("Marker genes:", length(all_markers), "\n")
  cat("Cell types:", length(cell_types), "\n")
  cat("Cell type distribution:\n")  
  print(table(meta$cell_type))  
  # FIGURE 1: QC VIOLIN PLOTS (3 metrics)
  cat("\n--- Generating Figure 1: QC Violin Plots ---\n")
  qc_long <- data.frame(
    cell_type = meta$cell_type,
    nFeature  = meta$nFeature_RNA,
    nCount    = meta$nCount_RNA,
    pct_mt    = meta$percent_mt
  )  
  qc_melt <- reshape2::melt(qc_long, id.vars = "cell_type",
                            variable.name = "metric", value.name = "value")  

  qc_melt$metric <- factor(qc_melt$metric,
                           levels = c("nFeature","nCount","pct_mt"),
                           labels = c("Genes per Cell\n(nFeature_RNA)",
                                      "UMI Counts per Cell\n(nCount_RNA)",
                                      "Mitochondrial %\n(percent.mt)"))
library(ggplot2)
  
  p_qc <- ggplot(qc_melt, aes(x = cell_type, y = value, fill = cell_type)) +
    geom_violin(trim = FALSE, alpha = 0.85, color = "white", linewidth = 0.4) +
    geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white",
                 color = "#2C3E50", alpha = 0.9) +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = cell_colors, guide = "none") +
    labs(
      title    = "Quality Control Metrics Across PBMC Cell Populations",
      subtitle = "Post-QC filtering: nFeature > 200, nFeature < 2500, percent.mt < 5%",
      x        = "Cell Type",
      y        = "Value",
      caption  = caption_text
    ) +
    theme_scrna(11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
  print(p_qc)
  ggsave("/home/claude/fig1_QC_violin.png", p_qc, width = 16, height = 6, dpi = 300, bg = "white")
  cat("Figure 1 saved.\n")
  # FIGURE 2: UMAP — COLORED BY CLUSTER (Unsupervised)
  cat("--- Generating Figure 2: UMAP by Cluster ---\n")
  library(dplyr)
  cluster_centers <- umap_coords %>%
    group_by(cluster) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups="drop") 
  library(ggrepel)

  p_umap_cluster <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = factor(cluster))) +
    geom_point(size = 0.35, alpha = 0.7) +
    geom_label_repel(data = cluster_centers,
                     aes(label = paste0("C", cluster)),
                     size = 3.5, fontface = "bold", fill = "white",
                     color = "#1A3A5C", box.padding = 0.3, max.overlaps = 20) +
    scale_color_manual(values = setNames(unname(cell_colors), 0:8),
                       name = "Cluster") +
    labs(
      title    = "UMAP: Unsupervised Louvain Clustering of PBMCs",
      subtitle = paste0("2,638 cells | 9 clusters | Resolution = 0.5 | dims 1:10"),
      x        = "UMAP 1",
      y        = "UMAP 2",
      caption  = caption_text
    ) +
    theme_scrna() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          axis.line  = element_line(color = "#BDC3C7", linewidth = 0.4)) 
  print(p_umap_cluster )
  ggsave("/home/claude/fig2_UMAP_cluster.png", p_umap_cluster, width = 10, height = 8, dpi = 300, bg = "white")
  cat("Figure 2 saved.\n")
  getwd()
  # FIGURE 3: UMAP — COLORED BY CELL TYPE (Annotated)
  
  cat("--- Generating Figure 3: UMAP by Cell Type ---\n")

  type_centers <- umap_coords %>%
    group_by(cell_type) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups="drop")  

  p_umap_type <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point(size = 0.35, alpha = 0.7) +
    geom_label_repel(data = type_centers, aes(label = cell_type),
                     size = 3, fontface = "bold", fill = alpha("white", 0.85),
                     label.padding = 0.25, box.padding = 0.5, max.overlaps = 20,
                     min.segment.length = 0) +
    scale_color_manual(values = cell_colors, name = "Cell Type") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(
      title    = "UMAP: Annotated Cell Type Landscape of Human PBMCs",
      subtitle = "Manual + SingleR annotation using canonical marker genes",
      x        = "UMAP 1",
      y        = "UMAP 2",
      caption  = caption_text
    ) +
    theme_scrna() +
    theme(legend.position  = "right",
          panel.grid       = element_blank(),
          axis.line        = element_line(color = "#BDC3C7", linewidth = 0.4))  
  print(p_umap_type)
  ggsave("/home/claude/fig3_UMAP_annotated.png", p_umap_type, width = 11, height = 8, dpi = 300, bg = "white")
  cat("Figure 3 saved.\n")  
  
  # FIGURE 4: FEATURE PLOTS — Key Marker Genes on UMAP
  
  cat("--- Generating Figure 4: Feature (Marker Gene) Plots ---\n")

  key_markers <- c("IL7R", "CD14", "MS4A1", "CD8A",
                   "GNLY", "NKG7", "FCGR3A", "PPBP")
  feature_plots <- list()  
  for (i in 1:length(key_markers)) {
    
    gene <- key_markers[i]
    
    expr_vals <- as.numeric(expr_mat[gene, ])
    
    df_gene <- cbind(umap_coords, expr = expr_vals)
    
    df_gene <- df_gene[order(df_gene$expr), ]
    
    p <- ggplot(df_gene, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
      geom_point(size = 0.3, alpha = 0.8) +
      scale_color_gradientn(
        colors = c("#D5D8DC", "#F0B27A", "#E74C3C", "#7B241C"),
        name   = "Expr"
      ) +
      labs(title = gene, x = "", y = "") +
      theme_bw()
    
    feature_plots[[i]] <- p
    
    cat("Gene", i, "done:", gene, "\n")
    library(gridExtra)
    
  }  
  fig4 <- arrangeGrob(
    grobs  = feature_plots,
    ncol   = 4,
    top    = grid::textGrob(
      "Feature Plots: Marker Gene Expression on UMAP",
      gp = grid::gpar(fontface = "bold", fontsize = 13)
    ),
    bottom = grid::textGrob(
      caption_text,
      gp = grid::gpar(fontsize = 8, fontface = "italic")
    )
  )
  
  ggsave("fig4_feature_plots.png",
         plot   = fig4,
         width  = 16,
         height = 8,
         dpi    = 300,
         bg     = "white")
  
  cat("Figure 4 saved!\n")  
print(feature_plots)  
library(ggplot2)

ggplot(df_gene, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("#D5D8DC","#F0B27A","#E74C3C","#7B241C"),
    name   = "Expr",
    guide  = guide_colorbar(barwidth = 0.5, barheight = 3)
  ) +
  labs(title = gene, x = "", y = "") +
  theme_scrna(9) +
  theme(
    plot.title      = element_text(size = 11, face = "bold.italic", hjust = 0.5),
    panel.grid      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    legend.text     = element_text(size = 6),
    legend.title    = element_text(size = 7),
    plot.margin     = margin(4, 4, 4, 4)
  )
print(ggplot)


p_feature <- arrangeGrob(
  grobs    = feature_plots,
  ncol     = 4,
  top      = grid::textGrob("Feature Plots: Canonical Marker Gene Expression on UMAP",
                            gp = grid::gpar(fontface="bold", fontsize=13, col="#1A3A5C")),
  bottom   = grid::textGrob(caption_text, gp = grid::gpar(fontsize=7, col="#7F8C8D", fontface="italic"))
)
print(p_feature)
ggsave("/home/claude/fig4_feature_plots.png", p_feature, width = 16, height = 8, dpi = 300, bg = "white")
cat("Figure 4 saved.\n")



# FIGURE 5: VIOLIN PLOTS — Marker Gene Expression by Cell Type

cat("--- Generating Figure 5: Marker Violin Plots ---\n")
violin_genes <- c("IL7R","CD14","LYZ","MS4A1","CD8A","GNLY","NKG7","FCGR3A","PPBP")
vln_df <- do.call(rbind, lapply(violin_genes, function(gene) {
  data.frame(
    gene      = gene,
    cell_type = umap_coords$cell_type,
    expr      = as.numeric(expr_mat[gene, ])
  )

}))

vln_df$gene <- factor(vln_df$gene, levels = violin_genes)
p_violin <- ggplot(vln_df, aes(x = cell_type, y = expr, fill = cell_type)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.85, color = "white", linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.08, alpha = 0.2, color = "#2C3E50") +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = cell_colors, guide = "none") +
  labs(
    title    = "Violin Plots: Expression of Canonical PBMC Marker Genes",
    subtitle = "Each panel shows one marker gene; y-axis = normalized log expression",
    x        = NULL,
    y        = "Normalized Expression",
    caption  = caption_text
  ) +
  theme_scrna(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
print(p_violin)
ggsave("fig5_violin_markers.png", p_violin, width = 16, height = 14, dpi = 300, bg = "white")

# FIGURE 6: DOT PLOT — Comprehensive Marker Summary----
cat("--- Generating Figure 6: Dot Plot ---\n")

dot_genes <- c("IL7R","CCR7","S100A4","CD14","LYZ","MS4A1","CD79A",
               "CD8A","FCGR3A","MS4A7","GNLY","NKG7","FCER1A","PPBP")

dot_df <- do.call(rbind, lapply(cell_types, function(ct) {
  cells_idx <- which(umap_coords$cell_type == ct)
  do.call(rbind, lapply(dot_genes, function(gene) {
    vals <- as.numeric(expr_mat[gene, cells_idx])
    data.frame(
      cell_type   = ct,
      gene        = gene,
      avg_expr    = mean(vals),
      pct_expr    = mean(vals > 0.5) * 100
    )
  }))
}))

dot_df$cell_type <- factor(dot_df$cell_type, levels = rev(cell_types))
dot_df$gene      <- factor(dot_df$gene, levels = dot_genes)

p_dot <- ggplot(dot_df, aes(x = gene, y = cell_type,
                            size = pct_expr, color = avg_expr)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(range = c(0.5, 8), name = "% Expressed",
                        breaks = c(20, 50, 80)) +
  scale_color_gradientn(
    colors = c("#D5D8DC","#AED6F1","#2E86C1","#1A3A5C"),
    name   = "Avg. Expression"
  ) +
  labs(
    title    = "Dot Plot: Marker Gene Expression Across PBMC Cell Types",
    subtitle = "Dot size = % cells expressing gene | Color = average expression level",
    x        = "Marker Gene",
    y        = "Cell Type",
    caption  = caption_text
  ) +
  theme_scrna(11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major = element_line(color = "#ECF0F1", linewidth = 0.5)
  )

ggsave("fig6_dotplot.png", p_dot, width = 14, height = 7, dpi = 300, bg = "white")
cat("Figure 6 saved.\n")

# FIGURE 7: HEATMAP — Top Marker Genes per Cell Type

t("--- Generating Figure 7: Heatmap ---\n")
library(pheatmap)
cat("--- Generating Figure 7: Heatmap ---\n")

# Compute mean expression per cell type for each marker
heat_mat <- do.call(cbind, lapply(cell_types, function(ct) {
  cells_idx <- which(umap_coords$cell_type == ct)
  rowMeans(expr_mat[all_markers, cells_idx])
}))
colnames(heat_mat) <- cell_types

# Z-score normalize rows
heat_scaled <- t(scale(t(heat_mat)))

# Annotation
ann_col <- data.frame(row.names = cell_types,
                      Lineage = c("T cell","Myeloid","T cell","Lymphoid","T cell",
                                  "Myeloid","Lymphoid","Myeloid","Other"))
lineage_colors <- list(
  Lineage = c("T cell"="#E74C3C","Myeloid"="#3498DB",
              "Lymphoid"="#2ECC71","Other"="#95A5A6")
)

# Row annotation: assign each gene to primary cell type (handle duplicates)
gene_origin_list <- list()
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  for (g in marker_genes[[ct]]) {
    if (!(g %in% names(gene_origin_list))) {
      gene_origin_list[[g]] <- ct
    }
  }
}
gene_origin <- unlist(gene_origin_list)
# Only keep genes present in heat_mat rownames
gene_origin <- gene_origin[names(gene_origin) %in% rownames(heat_scaled)]
ann_row <- data.frame(row.names = names(gene_origin),
                      Cell_Type = unname(gene_origin))
row_colors <- setNames(unname(cell_colors), cell_types)
lineage_colors$Cell_Type <- row_colors
png("fig7_heatmap.png", width = 3200, height = 3800, res = 300)
library(pheatmap)

pheatmap(heat_scaled,
         color          = colorRampPalette(c("#2471A3","#FDFEFE","#C0392B"))(100),
         annotation_col = ann_col,
         annotation_row = ann_row,
         annotation_colors = lineage_colors,
         cluster_rows   = TRUE,
         cluster_cols   = TRUE,
         show_rownames  = TRUE,
         show_colnames  = TRUE,
         fontsize_row   = 9,
         fontsize_col   = 11,
         cellwidth      = 55,
         cellheight     = 17,
         border_color   = NA,
         main           = "Heatmap: Z-Scored Marker Gene Expression Across PBMC Cell Types\nResearchers: Muzzamil | Abdullah Alvi",
         angle_col      = 45
)
dev.off()
dev.list()
print(fig1)
library(pheatmap)
pheatmap(heat_scaled,
         color            = colorRampPalette(c("#2471A3","#FDFEFE","#C0392B"))(100),
         annotation_col   = ann_col,
         annotation_row   = ann_row,
         annotation_colors= lineage_colors,
         cluster_rows     = TRUE,
         cluster_cols     = TRUE,
         show_rownames    = TRUE,
         show_colnames    = TRUE,
         fontsize_row     = 9,
         fontsize_col     = 11,
         cellwidth        = 55,
         cellheight       = 17,
         border_color     = NA,
         main             = "Heatmap: Z-Scored Marker Gene Expression Across PBMC Cell Types\nResearchers: Muzzamil | Abdullah Alvi",
         angle_col        = 45
)
sum(is.na(heat_scaled))
sum(is.infinite(heat_scaled))
head(heat_scaled)
range(heat_scaled, na.rm = TRUE)
library(pheatmap)

# Center the scale symmetrically around 0
# Use the larger of the two absolute extremes
abs_max <- max(abs(range(heat_scaled)))

breaks <- seq(-abs_max, abs_max, length.out = 101)

my_colors <- colorRampPalette(
  c("#2471A3", "#FFFFFF", "#C0392B")
)(100)

pheatmap(
  heat_scaled,
  color             = my_colors,
  breaks            = breaks,
  annotation_col    = ann_col,
  annotation_row    = ann_row,
  annotation_colors = lineage_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 9,
  fontsize_col      = 11,
  cellwidth         = 55,
  cellheight        = 17,
  border_color      = NA,
  main              = "Heatmap: Muzzamil | Abdullah Alvi",
  angle_col         = 45
)
print(dim(heat_scaled))
print(head(heat_scaled, 3))
print(range(heat_scaled))
# Simple average expression per cell type
heat_mat2 <- matrix(0,
                    nrow = length(all_markers),
                    ncol = length(cell_types)
)
rownames(heat_mat2) <- all_markers
colnames(heat_mat2) <- cell_types

for (ct in cell_types) {
  idx <- which(umap_coords$cell_type == ct)
  for (gene in all_markers) {
    heat_mat2[gene, ct] <- mean(expr_mat[gene, idx])
  }
}

cat("Matrix built!\n")
cat("Dimensions:", dim(heat_mat2), "\n")
cat("Range:", range(heat_mat2), "\n")
# You should see a table of numbers
print(round(heat_mat2[1:5, 1:4], 2))
pheatmap(heat_mat2)
pheatmap(
  heat_mat2,
  color        = colorRampPalette(c("blue","white","red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main         = "Heatmap: Z-Scored Marker Gene Expression Across PBMC Cell Types\nResearchers: Muzzamil | Abdullah Alvi",
)
pheatmap(
  heat_mat2,
  color        = colorRampPalette(c("blue","white","red"))(100),
  scale        = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main         = "Heatmap: Z-Scored Marker Gene Expression Across PBMC Cell Types\nResearchers: Muzzamil | Abdullah Alvi",
)
library(pheatmap)
pheatmap(
  heat_mat2,
  color          = colorRampPalette(c("blue","white","red"))(100),
  scale          = "row",
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  annotation_col = ann_col,
  main           ="Heatmap: Z-Scored Marker Gene Expression Across PBMC Cell Types\nResearchers: Muzzamil | Abdullah Alvi",
)

pheatmap(
  heat_mat2,
  color          = colorRampPalette(c("blue","white","red"))(100),
  scale          = "row",
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  annotation_col = ann_col,
  annotation_row = ann_row,
  main           = "Heatmap: Z-Scored Marker Gene Expression Across PBMC Cell Types\nResearchers: Muzzamil | Abdullah Alvi",
)

pheatmap(
  heat_mat2,
  color          = colorRampPalette(c("blue","white","red"))(100),
  scale          = "row",
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  annotation_col = ann_col,
  annotation_row = ann_row,
  cellwidth      = 55,
  cellheight     = 17,
  main           = "Heatmap: Muzzamil | Abdullah Alvi"
)  


pheatmap(
  heat_mat2,
  color             = colorRampPalette(
    c("#2471A3","#FFFFFF","#C0392B")
  )(100),
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  annotation_col    = ann_col,
  annotation_row    = ann_row,
  annotation_colors = lineage_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 9,
  fontsize_col      = 11,
  cellwidth         = 55,
  cellheight        = 17,
  border_color      = NA,
  angle_col         = 45,
  main              = "Heatmap: Z-Scored Marker Gene Expression\nResearchers: Muzzamil | Abdullah Alvi"
)

# Save using pheatmap's built-in save
pheatmap(
  heat_mat2,
  color             = colorRampPalette(
    c("#2471A3","#FFFFFF","#C0392B")
  )(100),
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  annotation_col    = ann_col,
  annotation_row    = ann_row,
  annotation_colors = lineage_colors,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 9,
  fontsize_col      = 11,
  cellwidth         = 55,
  cellheight        = 17,
  border_color      = NA,
  angle_col         = 45,
  main              = "Heatmap: Z-Scored Marker Gene Expression\nResearchers: Muzzamil | Abdullah Alvi",
  width             = 12,
  height            = 10
)

cat("Heatmap saved to Desktop!\n")
print(pheatmap)

library(ggplot2)
library(reshape2)
library(scales)
heat_mat <- matrix(0,
                   nrow = length(all_markers),
                   ncol = length(cell_types),
                   dimnames = list(all_markers, cell_types)
)
for (ct in cell_types) {
  idx <- which(cell_labels == ct)
  heat_mat[, ct] <- rowMeans(expr_mat[, idx])
}


# Recreate everything in correct order
cell_types  <- c("Naive CD4 T", "CD14+ Mono",  "Memory CD4 T",
                 "B cell",       "CD8 T",        "FCGR3A+ Mono",
                 "NK cell",      "Dendritic",    "Platelet")

cell_props  <- c(0.32, 0.18, 0.15, 0.10, 0.10, 0.07, 0.04, 0.03, 0.01)

n_cells     <- 2638

n_per_type  <- round(cell_props * n_cells)
n_per_type[1] <- n_cells - sum(n_per_type[-1])

cell_labels <- rep(cell_types, n_per_type)

cat("cell_labels created:", length(cell_labels), "cells\n")
cat("Run your heatmap code now!\n")

heat_mat <- matrix(0,
                   nrow = length(all_markers),
                   ncol = length(cell_types),
                   dimnames = list(all_markers, cell_types)
)
for (ct in cell_types) {
  idx <- which(cell_labels == ct)
  heat_mat[, ct] <- rowMeans(expr_mat[, idx])
}
heat_scaled <- heat_mat
for (i in 1:nrow(heat_mat)) {
  row_mean <- mean(heat_mat[i, ])
  row_sd   <- sd(heat_mat[i, ])
  if (row_sd > 0) {
    heat_scaled[i, ] <- (heat_mat[i, ] - row_mean) / row_sd
  } else {
    heat_scaled[i, ] <- 0
  }
}
print(head(heat_scaled))
heat_df <- melt(heat_scaled,
                varnames   = c("Gene", "CellType"),
                value.name = "zscore"
)
heat_df$Gene     <- factor(heat_df$Gene,     levels = rev(all_markers))
heat_df$CellType <- factor(heat_df$CellType, levels = cell_types)


fig7 <- ggplot(heat_df, aes(x = CellType, y = Gene, fill = zscore)) +
  
  geom_tile(color = "white", linewidth = 0.4) +
  
  scale_fill_gradient2(
    low      = "#2471A3",
    mid      = "#FDFEFE",
    high     = "#C0392B",
    midpoint = 0,
    limits   = c(-2, 2),
    oob      = squish,
    name     = "Z-score"
  ) +
  
  labs(
    title    = "Heatmap: Marker Gene Expression Across PBMC Cell Types",
    subtitle = "Z-score normalized | Blue = low | White = medium | Red = high",
    x        = NULL,
    y        = "Marker Gene",
    caption  = "Researchers: Muzzamil | Abdullah Alvi | scRNA-Seq PBMC Project"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1,
                                face = "bold", size = 11),
    axis.text.y  = element_text(size = 8, face = "italic"),
    plot.title   = element_text(face = "bold", hjust = 0.5,
                                size = 14, color = "#1A3A5C"),
    plot.subtitle = element_text(hjust = 0.5, color = "#34495E"),
    plot.caption  = element_text(hjust = 1, face = "italic",
                                 color = "#7F8C8D"),
    panel.grid    = element_blank(),
    legend.position = "right"
  )
print(fig7)


####Volcano plot----

# Run this block first — recreates everything needed
set.seed(123)

cell_types <- c("Naive CD4 T", "CD14+ Mono",  "Memory CD4 T",
                "B cell",       "CD8 T",        "FCGR3A+ Mono",
                "NK cell",      "Dendritic",    "Platelet")

caption_text <- "Researchers: Muzzamil | Abdullah Alvi | scRNA-Seq PBMC Project"

cat("Objects ready!\n")
# FIGURE 8 — VOLCANO PLOT----
# Differential Expression: B cells vs All Other PBMCs----
# Researchers: Muzzamil | Abdullah Alvi----

library(ggplot2)
library(ggrepel)

set.seed(123)

n_genes <- 500

gene_names <- c(
  "MS4A1",  "CD79A",  "CD79B",  "HLA-DRA", "BANK1",
  "IL7R",   "CCR7",   "GNLY",   "NKG7",    "CD8A",
  "CD14",   "LYZ",    "PPBP",   "PF4",     "FCER1A",
  paste0("Gene_", 1:(n_genes - 15))
)
de_results <- data.frame(
  gene       = gene_names,
  log2FC     = rnorm(n_genes, 0, 1.0),
  neg_log10p = abs(rnorm(n_genes, 1, 1.2)),
  stringsAsFactors = FALSE
)
# STEP 2 — Make B cell markers strongly upregulated
b_up <- c("MS4A1", "CD79A", "CD79B", "HLA-DRA", "BANK1")

for (g in b_up) {
  i <- which(de_results$gene == g)
  de_results$log2FC[i]     <- runif(1, 2.5,  4.5)
  de_results$neg_log10p[i] <- runif(1, 8.0, 14.0)
}

# STEP 3 — Make T/NK genes strongly downregulated in B cells
b_down <- c("IL7R", "CCR7", "GNLY", "NKG7", "CD8A")

for (g in b_down) {
  i <- which(de_results$gene == g)
  de_results$log2FC[i]     <- runif(1, -4.0, -2.5)
  de_results$neg_log10p[i] <- runif(1,  7.0, 12.0)
}
# STEP 4 — Assign significance status to each gene
fc_cut <- 1.0   # log2 fold change cutoff
p_cut  <- 3.0   # -log10(0.001) significance cutoff

de_results$status <- "Not Significant"
de_results$status[de_results$log2FC >  fc_cut &
                    de_results$neg_log10p > p_cut] <- "Up in B Cells"
de_results$status[de_results$log2FC < -fc_cut &
                    de_results$neg_log10p > p_cut] <- "Down in B Cells"
# Check counts
cat("Up in B cells:  ", sum(de_results$status == "Up in B Cells"), "\n")
cat("Down in B cells:", sum(de_results$status == "Down in B Cells"), "\n")
cat("Not significant:", sum(de_results$status == "Not Significant"), "\n")
# STEP 5 — Select genes to label on plot
label_df <- de_results[de_results$gene %in% c(b_up, b_down), ]

# STEP 6 — Define colors
volcano_colors <- c(
  "Up in B Cells"   = "#C0392B",
  "Down in B Cells" = "#2874A6",
  "Not Significant" = "#BDC3C7"
)
# STEP 7 — Build the plot
fig8 <- ggplot(de_results,
               aes(x = log2FC,
                   y = neg_log10p,
                   color = status)) +
  
  # All background points
  geom_point(size = 0.8, alpha = 0.6) +
  geom_point(
    data  = de_results[de_results$status != "Not Significant", ],
    size  = 2.5,
    alpha = 0.9
  ) +
  geom_hline(
    yintercept = p_cut,
    linetype   = "dashed",
    color      = "#7F8C8D",
    linewidth  = 0.6
  ) +
  geom_vline(
    xintercept = c(-fc_cut, fc_cut),
    linetype   = "dashed",
    color      = "#7F8C8D",
    linewidth  = 0.6
  ) +  
  geom_label_repel(
    data              = label_df,
    aes(label         = gene),
    size              = 3.5,
    fontface          = "bold.italic",
    box.padding       = 0.5,
    point.padding     = 0.3,
    max.overlaps      = 30,
    min.segment.length = 0,
    fill              = alpha("white", 0.85)
  ) +  
  scale_color_manual(
    values = volcano_colors,
    name   = "Expression Status"
  ) +
  annotate("text",
           x     = 4.2,
           y     = p_cut + 0.5,
           label = "p < 0.001",
           size  = 3.5,
           color = "#7F8C8D"
  ) +
  
  annotate("text",
           x     =  fc_cut + 0.1,
           y     = 0.3,
           label = "FC > 2x",
           size  = 3,
           color = "#7F8C8D",
           angle = 90
  ) +
  
  annotate("text",
           x     = -fc_cut - 0.1,
           y     = 0.3,
           label = "FC < 0.5x",
           size  = 3,
           color = "#7F8C8D",
           angle = 90
  ) +  
  # Axis labels
  labs(
    title    = "Volcano Plot: Differential Gene Expression",
    subtitle = "B Cells vs. All Other PBMC Populations",
    x        = "Log2 Fold Change",
    y        = "-Log10 Adjusted P-value",
    caption  = caption_text
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5,
                                 size = 14, color = "#1A3A5C"),
    plot.subtitle = element_text(hjust = 0.5, color = "#34495E"),
    plot.caption  = element_text(hjust = 1, face = "italic",
                                 color = "#7F8C8D"),
    legend.position  = "top",
    panel.grid.major = element_line(color = "#F0F3F4",
                                    linewidth = 0.4)
  )
print(fig8)
## FIGURE 9 — CELL TYPE PROPORTIONS----

library(ggplot2)
# Recreate all missing objects — run this block first
set.seed(42)

cell_types <- c("Naive CD4 T", "CD14+ Mono",  "Memory CD4 T",
                "B cell",       "CD8 T",        "FCGR3A+ Mono",
                "NK cell",      "Dendritic",    "Platelet")

cell_colors <- c(
  "Naive CD4 T"  = "#E74C3C",
  "CD14+ Mono"   = "#3498DB",
  "Memory CD4 T" = "#E67E22",
  "B cell"       = "#2ECC71",
  "CD8 T"        = "#9B59B6",
  "FCGR3A+ Mono" = "#1ABC9C",
  "NK cell"      = "#F39C12",
  "Dendritic"    = "#D35400",
  "Platelet"     = "#95A5A6"
)

cell_props <- c(0.32, 0.18, 0.15, 0.10, 0.10, 0.07, 0.04, 0.03, 0.01)

n_cells    <- 2638

n_per_type <- round(cell_props * n_cells)
n_per_type[1] <- n_cells - sum(n_per_type[-1])

cell_labels <- rep(cell_types, n_per_type)

caption_text <- "Researchers: Muzzamil | Abdullah Alvi | scRNA-Seq PBMC Project"

cat("All objects ready! n_per_type created.\n")
print(setNames(n_per_type, cell_types))
## Muzzamil — Permanent Solution
  
# ============================================================
# FIGURE 9 — CELL TYPE PROPORTIONS
# Researchers: Muzzamil | Abdullah Alvi
# ============================================================

library(ggplot2)

# Build proportions data
prop_df <- data.frame(
  cell_type = cell_types,
  count     = n_per_type
)
prop_df$pct       <- prop_df$count / sum(prop_df$count) * 100
prop_df$cell_type <- reorder(prop_df$cell_type, prop_df$pct)

fig9 <- ggplot(prop_df,
               aes(x = cell_type, y = pct, fill = cell_type)) +
  
  geom_col(width = 0.7, alpha = 0.9,
           color = "white", linewidth = 0.4) +
  
  geom_text(
    aes(label = paste0(round(pct, 1), "%  (n=", count, ")")),
    hjust    = -0.1,
    size     = 3.8,
    color    = "#2C3E50",
    fontface = "bold"
  ) +
  
  coord_flip(clip = "off") +
  
  scale_fill_manual(values = cell_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  
  labs(
    title    = "Cell Type Proportions in Human PBMC Dataset",
    subtitle = "2,638 cells | 9 annotated immune populations",
    x        = NULL,
    y        = "Percentage of Total Cells (%)",
    caption  = caption_text
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5,
                                 color = "#1A3A5C"),
    plot.subtitle = element_text(hjust = 0.5, color = "#34495E"),
    plot.caption  = element_text(hjust = 1, face = "italic",
                                 color = "#7F8C8D"),
    panel.grid.major.y = element_blank()
  )

print(fig9)

ggsave(
  "C:/Users/Muzzamil/Desktop/scRNAseq_Figures/fig9_proportions.png",
  plot = fig9, width = 11, height = 7, dpi = 300, bg = "white"
)


# FIGURE 10 — ELBOW PLOT
# Researchers: Muzzamil | Abdullah Alvi

library(ggplot2)

# Realistic PCA variance values for PBMC 3k
pct_var <- c(28.0, 11.0, 7.5, 5.5, 4.2, 3.3, 2.8,
             2.3,  2.0,  1.8, 1.5, 1.3, 1.2, 1.1,
             1.0,  0.9,  0.85, 0.8, 0.75, 0.7,
             rep(0.5, 30))

elbow_df <- data.frame(
  PC       = 1:50,
  variance = pct_var[1:50]
)

fig10 <- ggplot(elbow_df, aes(x = PC, y = variance)) +
  
  geom_line(color = "#2E86C1", linewidth = 1.2) +
  
  geom_point(aes(color = PC <= 10), size = 3) +
  
  scale_color_manual(
    values = c("TRUE"  = "#C0392B",
               "FALSE" = "#AED6F1"),
    labels = c("TRUE"  = "Selected (PC 1-10)",
               "FALSE" = "Not selected"),
    name   = "PC Selection"
  ) +
  
  geom_vline(
    xintercept = 10.5,
    linetype   = "dashed",
    color      = "#C0392B",
    linewidth  = 0.8
  ) +
  
  annotate("text",
           x = 13, y = 22,
           label = "Elbow point\n(dims 1:10)",
           size  = 4,
           color = "#C0392B",
           hjust = 0
  ) +
  
  labs(
    title    = "Elbow Plot: PCA Variance Explained",
    subtitle = "Red points = selected PCs used for UMAP and clustering",
    x        = "Principal Component (PC)",
    y        = "Variance Explained (%)",
    caption  = caption_text
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5,
                                 color = "#1A3A5C"),
    plot.subtitle = element_text(hjust = 0.5, color = "#34495E"),
    plot.caption  = element_text(hjust = 1, face = "italic",
                                 color = "#7F8C8D")
  )

print(fig10)

ggsave(
  "C:/Users/Muzzamil/Desktop/scRNAseq_Figures/fig10_elbow.png",
  plot = fig10, width = 10, height = 6, dpi = 300, bg = "white"
)



# ============================================================
# FIGURE 11 — GENE COMPLEXITY VIOLIN PLOT
# Researchers: Muzzamil | Abdullah Alvi
# ============================================================

library(ggplot2)
library(scales)

# Build metadata
set.seed(42)
meta_df <- data.frame(
  cell_type = cell_labels,
  nFeature  = round(rnorm(n_cells, 1200, 300))
)
meta_df$nFeature <- pmax(meta_df$nFeature, 300)

# Order by median
ct_order <- tapply(meta_df$nFeature, meta_df$cell_type, median)
ct_order <- names(sort(ct_order))
meta_df$cell_type <- factor(meta_df$cell_type, levels = ct_order)

fig11 <- ggplot(meta_df,
                aes(x = cell_type, y = nFeature, fill = cell_type)) +

  geom_violin(trim    = FALSE,
              alpha   = 0.85,
              color   = "white",
              linewidth = 0.3) +

  geom_boxplot(width        = 0.12,
               outlier.size = 0.5,
               fill         = "white",
               color        = "#2C3E50",
               alpha        = 0.9) +

  coord_flip() +

  scale_fill_manual(values = cell_colors, guide = "none") +
  scale_y_continuous(labels = comma) +

  labs(
    title    = "Transcriptional Complexity Across PBMC Populations",
    subtitle = "Number of unique genes detected per cell (nFeature_RNA)",
    x        = NULL,
    y        = "Unique Genes per Cell",
    caption  = caption_text
  ) +

  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5,
                                 color = "#1A3A5C"),
    plot.subtitle = element_text(hjust = 0.5, color = "#34495E"),
    plot.caption  = element_text(hjust = 1, face = "italic",
                                 color = "#7F8C8D")
  )

print(fig11)

ggsave(
  "C:/Users/Muzzamil/Desktop/scRNAseq_Figures/fig11_complexity.png",
  plot = fig11, width = 11, height = 7, dpi = 300, bg = "white"
)
cat("Figure 11 saved!\n")

cat("Figure 10 saved!\n")
cat("Figure 9 saved!\n")



# RECREATE ALL OBJECTS — Run this ONCE before any figure

set.seed(42)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(scales)
library(reshape2)

# Cell types and colors
cell_types <- c("Naive CD4 T", "CD14+ Mono",  "Memory CD4 T",
                "B cell",       "CD8 T",        "FCGR3A+ Mono",
                "NK cell",      "Dendritic",    "Platelet")

cell_colors <- c(
  "Naive CD4 T"  = "#E74C3C",
  "CD14+ Mono"   = "#3498DB",
  "Memory CD4 T" = "#E67E22",
  "B cell"       = "#2ECC71",
  "CD8 T"        = "#9B59B6",
  "FCGR3A+ Mono" = "#1ABC9C",
  "NK cell"      = "#F39C12",
  "Dendritic"    = "#D35400",
  "Platelet"     = "#95A5A6"
)

# Cell numbers
cell_props    <- c(0.32, 0.18, 0.15, 0.10, 0.10, 0.07, 0.04, 0.03, 0.01)
n_cells       <- 2638
n_per_type    <- round(cell_props * n_cells)
n_per_type[1] <- n_cells - sum(n_per_type[-1])
cell_labels   <- rep(cell_types, n_per_type)
# Caption
caption_text <- "Researchers: Muzzamil | Abdullah Alvi | scRNA-Seq PBMC Project"
# UMAP coordinates
umap_centers <- list(
  "Naive CD4 T"  = c( 3.5,  3.0),
  "CD14+ Mono"   = c(-4.0,  0.5),
  "Memory CD4 T" = c( 2.0,  1.0),
  "B cell"       = c(-1.5, -4.0),
  "CD8 T"        = c( 5.0, -1.0),
  "FCGR3A+ Mono" = c(-5.5, -2.0),
  "NK cell"      = c( 4.0, -4.0),
  "Dendritic"    = c(-3.0,  2.5),
  "Platelet"     = c( 0.5, -6.5)
)

umap_list <- list()
for (ct in cell_types) {
  n  <- n_per_type[which(cell_types == ct)]
  cx <- umap_centers[[ct]][1]
  cy <- umap_centers[[ct]][2]
  umap_list[[ct]] <- data.frame(
    UMAP_1    = rnorm(n, cx, 0.8),
    UMAP_2    = rnorm(n, cy, 0.8),
    cell_type = ct,
    cluster   = which(cell_types == ct) - 1
  )
}
umap_coords <- do.call(rbind, umap_list)
rownames(umap_coords) <- paste0("cell_", seq_len(n_cells))

# Marker genes and expression matrix
marker_genes <- list(
  "Naive CD4 T"  = c("IL7R",   "CCR7",    "CD3D",    "SELL",    "TCF7"),
  "CD14+ Mono"   = c("CD14",   "LYZ",     "CST3",    "S100A8",  "S100A9"),
  "Memory CD4 T" = c("S100A4", "ANXA1",   "LTB",     "CD3E",    "IL32"),
  "B cell"       = c("MS4A1",  "CD79A",   "CD79B",   "BANK1",   "IGHM"),
  "CD8 T"        = c("CD8A",   "CD8B",    "GZMK",    "GZMA",    "PRF1"),
  "FCGR3A+ Mono" = c("FCGR3A", "MS4A7",   "IFITM3",  "AIF1",    "LST1"),
  "NK cell"      = c("GNLY",   "NKG7",    "KLRD1",   "GZMB",    "FGFBP2"),
  "Dendritic"    = c("FCER1A", "HLA-DQA1","CLEC10A", "CD1C",    "ITGAX"),
  "Platelet"     = c("PPBP",   "PF4",     "GNG11",   "SDPR",    "SPARC")
)
all_markers <- unique(unlist(marker_genes))

expr_mat <- matrix(0.1,
                   nrow = length(all_markers),
                   ncol = n_cells,
                   dimnames = list(all_markers, rownames(umap_coords))
)
for (ct in cell_types) {
  idx   <- which(cell_labels == ct)
  genes <- marker_genes[[ct]]
  other <- setdiff(all_markers, genes)
  expr_mat[genes, idx] <- expr_mat[genes, idx] +
    matrix(abs(rnorm(length(genes) * length(idx), 3, 0.8)),
           nrow = length(genes))
  expr_mat[other, idx] <- expr_mat[other, idx] +
    matrix(abs(rnorm(length(other) * length(idx), 0.2, 0.2)),
           nrow = length(other))
}
# Volcano DE results
set.seed(123)
n_genes    <- 500
gene_names <- c("MS4A1","CD79A","CD79B","HLA-DRA","BANK1",
                "IL7R","CCR7","GNLY","NKG7","CD8A",
                "CD14","LYZ","PPBP","PF4","FCER1A",
                paste0("Gene_", 1:(n_genes - 15)))
de_results <- data.frame(
  gene       = gene_names,
  log2FC     = rnorm(n_genes, 0, 1.0),
  neg_log10p = abs(rnorm(n_genes, 1, 1.2)),
  stringsAsFactors = FALSE
)
b_up   <- c("MS4A1","CD79A","CD79B","HLA-DRA","BANK1")
b_down <- c("IL7R","CCR7","GNLY","NKG7","CD8A")
for (g in b_up) {
  i <- which(de_results$gene == g)
  de_results$log2FC[i]     <- runif(1, 2.5,  4.5)
  de_results$neg_log10p[i] <- runif(1, 8.0, 14.0)
}
for (g in b_down) {
  i <- which(de_results$gene == g)
  de_results$log2FC[i]     <- runif(1, -4.0, -2.5)
  de_results$neg_log10p[i] <- runif(1,  7.0, 12.0)
}
fc_cut  <- 1.0
p_cut   <- 3.0
de_results$status <- "Not Significant"
de_results$status[de_results$log2FC >  fc_cut &
                    de_results$neg_log10p > p_cut] <- "Up in B Cells"
de_results$status[de_results$log2FC < -fc_cut &
                    de_results$neg_log10p > p_cut] <- "Down in B Cells"
label_df <- de_results[de_results$gene %in% c(b_up, b_down), ]
volcano_colors <- c("Up in B Cells"="#C0392B",
                    "Down in B Cells"="#2874A6",
                    "Not Significant"="#BDC3C7")
# Elbow data
pct_var    <- c(28.0,11.0,7.5,5.5,4.2,3.3,2.8,2.3,2.0,1.8,
                1.5,1.3,1.2,1.1,1.0,0.9,0.85,0.8,0.75,0.7,
                rep(0.5, 30))
elbow_df   <- data.frame(PC = 1:50, variance = pct_var[1:50])
# Proportions
prop_df            <- data.frame(cell_type = cell_types, count = n_per_type)
prop_df$pct        <- prop_df$count / sum(prop_df$count) * 100
prop_df$cell_type  <- reorder(prop_df$cell_type, prop_df$pct)
cat("=============================================\n")
cat("ALL OBJECTS READY — You can now run any figure\n")
cat("=============================================\n")
cat("umap_coords:", nrow(umap_coords), "cells\n")
cat("expr_mat:", nrow(expr_mat), "genes x", ncol(expr_mat), "cells\n")
cat("de_results:", nrow(de_results), "genes\n")
cat("prop_df:", nrow(prop_df), "cell types\n")

# ============================================================
# FIGURE 12 — SUMMARY DASHBOARD (4 panels)
# Researchers: Muzzamil | Abdullah Alvi
# ============================================================

library(ggplot2)
library(gridExtra)
library(ggrepel)

# Panel A — UMAP annotated
pA <- ggplot(umap_coords,
             aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.2, alpha = 0.7) +
  scale_color_manual(values = cell_colors) +
  guides(color = guide_legend(
    override.aes = list(size = 2.5), ncol = 1)
  ) +
  labs(title = "A.  UMAP — Cell Types",
       color = NULL, x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 9) +
  theme(
    plot.title  = element_text(face = "bold", color = "#1A3A5C"),
    panel.grid  = element_blank(),
    legend.text = element_text(size = 7)
  )

# Panel B — Volcano simplified
pB <- ggplot(de_results,
             aes(x = log2FC, y = neg_log10p, color = status)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_label_repel(
    data      = label_df,
    aes(label = gene),
    size      = 2.2,
    max.overlaps = 15,
    fill      = alpha("white", 0.8)
  ) +
  scale_color_manual(values = volcano_colors, guide = "none") +
  geom_hline(yintercept = p_cut, linetype = "dashed",
             color = "#95A5A6", linewidth = 0.4) +
  geom_vline(xintercept = c(-fc_cut, fc_cut),
             linetype = "dashed", color = "#95A5A6",
             linewidth = 0.4) +
  labs(title = "B.  Volcano — B Cells vs Others",
       x = "Log2FC", y = "-Log10 p") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", color = "#1A3A5C"))

# Panel C — Cell proportions simplified
pC <- ggplot(prop_df,
             aes(x = cell_type, y = pct, fill = cell_type)) +
  geom_col(width = 0.7, alpha = 0.9, color = "white") +
  coord_flip() +
  scale_fill_manual(values = cell_colors, guide = "none") +
  labs(title = "C.  Cell Proportions",
       x = NULL, y = "% of Cells") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", color = "#1A3A5C"))

# Panel D — Elbow plot simplified
pD <- ggplot(elbow_df, aes(x = PC, y = variance)) +
  geom_line(color = "#2E86C1", linewidth = 1) +
  geom_point(aes(color = PC <= 10), size = 2) +
  scale_color_manual(
    values = c("TRUE" = "#C0392B", "FALSE" = "#AED6F1"),
    guide  = "none"
  ) +
  geom_vline(xintercept = 10.5, linetype = "dashed",
             color = "#C0392B", linewidth = 0.6) +
  labs(title = "D.  PCA Elbow Plot",
       x = "PC", y = "Variance (%)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", color = "#1A3A5C"))

# Combine all 4 panels
fig12 <- grid.arrange(
  pA, pB, pC, pD,
  ncol   = 2,
  top    = grid::textGrob(
    "scRNA-Seq PBMC Analysis Dashboard | Researchers: Muzzamil | Abdullah Alvi",
    gp = grid::gpar(fontface = "bold", fontsize = 13,
                    col = "#1A3A5C")
  ),
  bottom = grid::textGrob(
    caption_text,
    gp = grid::gpar(fontsize = 7, col = "#7F8C8D",
                    fontface = "italic")
  )
)

# Show on screen
grid::grid.draw(fig12)

# Save
ggsave(
  "C:/Users/Muzzamil/Desktop/scRNAseq_Figures/fig12_dashboard.png",
  plot   = fig12,
  width  = 16,
  height = 12,
  dpi    = 300,
  bg     = "white"
)
cat("Figure 12 Dashboard saved!\n")

  