## PCA plot  
    library(DESeq2)
    library(ggplot2)

    vsd <- vst(dds, blind = FALSE) # PCA plot dùng “dds”

    pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 4, box.padding = 0.4, point.padding = 0.3,
                            max.overlaps = Inf, segment.color = "grey50") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA plot") +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title = element_text(face = "bold"),
        legend.position = "right",        # 👈 legend ra ngoài
        legend.title = element_blank(),
        legend.text = element_text(size = 12)
    )
    ggsave("/Volumes/Extreme SSD/Thesis/RNAseq_Test/PCA_plot_final.png",
        plot = pca_plot, width = 8, height = 6, dpi = 300)

    ggsave("/Volumes/Extreme SSD/Thesis/RNAseq_Test/PCA_plot_final.pdf",
        plot = pca_plot, width = 8, height = 6)

## Vocalno plot
    library(ggplot2) 
    library(dplyr)
    library(ggrepel)

    # dùng file res01 chỉ loại NA theo padj
    res01<- res_df[!is.na(res_df$padj), ] 
    res01$SYMBOL <- mapIds(org.Hs.eg.db,
                        keys = res01$gene_id_clean,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = function(x) paste(unique(x), collapse = ";")) # giữ lại tất cả annotation

    res01$status <- "NotSig"
    res01$status[res01$log2FoldChange >= 1 & res01$pvalue < 0.05] <- "Up"
    res01$status[res01$log2FoldChange <= -1 & res01$pvalue < 0.05] <- "Down"

    top_genes <- bind_rows(Top10_up, Top10_down) # top 10 gene up và 10 gene down

    # Vẽ volcano plot + chú giải top 10 gene up và 10 gene down
    p <- ggplot(res01, aes(x=log2FoldChange, y=-log10(pvalue))) +
    geom_point(aes(color=status), alpha=0.4, size=1.2) +
    scale_color_manual(values=c("Up"="red","Down"="blue","NotSig"="grey")) +

    geom_vline(xintercept=c(-1,1), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +

    geom_text_repel(data=top_genes,
        aes(label=SYMBOL),
        size=2,
        fontface="bold",
        box.padding=0.5,
        point.padding=0.3,
        segment.color="black",
        segment.size=0.4,
        max.overlaps = Inf
    ) +

    scale_x_continuous(limits=c(-25,25), breaks=seq(-25,25,5)) +

    theme_minimal(base_size=13) +
    labs(title="Volcano Plot: Cancer vs Normal",
        subtitle="Top 10 Up & Down genes by |log2FC| and padj",
        x="log2 Fold Change",
        y="-log10(p-value)",
        color="Group")

    p
    ggsave("/Volumes/Extreme SSD/Thesis/RNAseq_Test/volcano_res01.pdf",
        width = 7, height = 6)

##Dispersion estimate
    dds <- DESeq(dds)
    #Vẽ dispersion plot
    plotDispEsts(dds)
    #lưu file 
    pdf("/Volumes/Extreme SSD/Thesis/RNAseq_Test/dispersion_plot.pdf",
        width = 7, height = 6)
    plotDispEsts(dds)
    dev.off()


## Heatmap Plot
# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# Load kết quả DEseq2
res <- results(dds)

#1.Lấy top genes có significant differential expression

# Sắp xếp theo log2FoldChange để lấy top 50 gene up/down 
res03 <- res03[!is.na(res03$padj), ]
res03 <- res03[order(res03$log2FoldChange, decreasing = TRUE), ]
top_up <- res03 %>% 
  filter(padj < 0.05, log2FoldChange > 1) %>% 
  arrange(desc(log2FoldChange)) %>% 
  head(50)
top_down <- res03 %>% 
  filter(padj < 0.05, log2FoldChange < -1) %>% 
  arrange(desc(log2FoldChange)) %>% 
  tail(50)
res_top <- rbind(top_up, top_down) #top 50 gene up, 50 gene down
rownames(res_top) <- res_top$gene_id
top_genes <- res_top
gene_list <- rownames(top_genes)

#2. Chuẩn bị ma trận expression cho heatmap
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
mat_top <- mat[rownames(res_top), ]
# Z-score normalization (standardize theo hàng)
heatmap_data_scaled <- t(scale(t(mat_top)))

#3. Tạo annotation cho samples
# Tạo annotation dataframe cho samples
sample_annotation <- data.frame(
  Condition = colData(dds)$condition,  # thay đổi theo tên column thực tế
  row.names = colnames(heatmap_data_scaled)
)
# Tạo annotation cho genes (log2FC và adjusted p-value)
gene_annotation <- data.frame(
  logFC   = top_genes[gene_list, "log2FoldChange"],
  AveExpr = top_genes[gene_list, "baseMean"],
  SYMBOL  = top_genes[gene_list, "SYMBOL"],
  row.names = gene_list   # 👈 ENSG ID
)
# Sắp xếp gene theo logFC giảm dần
gene_annotation_sorted <- gene_annotation[order(-gene_annotation$logFC), ]
# Đồng bộ lại với Z-score
heatmap_data_scaled_sorted <- heatmap_data_scaled[rownames(gene_annotation_sorted), ]

#4. Tạo heatmap với pheatmap
library(ComplexHeatmap)
library(circlize)
# Tạo color function 
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
# Column annotation
col_ha <- HeatmapAnnotation(
  Condition = sample_annotation$Condition,
  col = list(Condition = c("normal" = "#FF9999", "cancer" = "#9999FF")),
  annotation_name_gp = gpar(fontsize = 8)   # 👈 không bold
)
# Row annotation  
row_ha_sorted <- rowAnnotation(
  logFC = gene_annotation_sorted$logFC,
  AveExpr = anno_simple(
    log10(gene_annotation_sorted$AveExpr + 1),
    col = colorRamp2(
      range(log10(gene_annotation_sorted$AveExpr + 1), na.rm = TRUE),
      c("white", "darkgreen")
    )
  ),
  col = list(logFC = colorRamp2(c(-15, 0, 25), c("blue", "white", "red"))),
  annotation_name_gp = gpar(fontsize = 8)   
)
# Tạo main heatmap
ht <- Heatmap(
  heatmap_data_scaled_sorted,
  name = "Z-score",
  col = col_fun,
  top_annotation = col_ha,
  right_annotation = row_ha_sorted,
  cluster_rows = FALSE,   
  cluster_columns = TRUE,
  show_row_names = TRUE,
row_labels = gene_annotation_sorted$SYMBOL,
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
column_title = "Heatmap Gene expression",   # 👈 TIÊU ĐỀ TRÊN
  row_title = "Top 100 DEGs ordered by log2FC"         # 👈 TIÊU ĐỀ BÊN TRÁI
)
draw(ht,
     annotation_legend_list = list(
       Legend(
         title = "Mean expression",
         at = range(log10(gene_annotation_sorted$AveExpr + 1), na.rm = TRUE),
         labels = c("Low", "High"),
         col_fun = colorRamp2(
           range(log10(gene_annotation_sorted$AveExpr + 1), na.rm = TRUE),
           c("white", "darkgreen")
         )
       )
     )
)
# Lưu hình
pdf("/Volumes/Extreme SSD/Thesis/RNAseq_Test/ heatmap_complex.pdf", width = 12, height = 10)
draw(ht,
     annotation_legend_list = list(
       Legend(
         title = "Mean expression",
         at = range(log10(gene_annotation_sorted$AveExpr + 1), na.rm = TRUE),
         labels = c("Low", "High"),
         col_fun = colorRamp2(
           range(log10(gene_annotation_sorted$AveExpr + 1), na.rm = TRUE),
           c("white", "darkgreen")
         )
       )
     )
)
dev.off()


