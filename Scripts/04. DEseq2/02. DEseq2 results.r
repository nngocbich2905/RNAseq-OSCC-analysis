#Data Visualization
##DEG _ DEseq2 results
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)

    #clean GeneID version
    res_df$gene_id_clean <- sub("\\..*$", "", res_df$gene_id)

    # Loại NA theo padj, Lọc padj<0.05
    res02<- res_df[!is.na(res_df$padj) &
                        res_df$padj < 0.05, ]
    # map ID từ ENSEMBL SANG GENE SYMBOL
    library(org.Hs.eg.db)
    library(AnnotationDbi)

    res02$SYMBOL <- mapIds(org.Hs.eg.db,
                        keys = res02$gene_id_clean,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = function(x) paste(unique(x), collapse = ";")) # giữ lại tất cả annotation vì có 1 ENSEMBL map được nhiều SYMBOL

    sum(res02$SYMBOL == "NA") # hàng có symbol = “NA”
    [1] 942
    res02_clean <- res02[!is.na(res02$SYMBOL) & res02$SYMBOL != "NA", ] # xoá hang symbol = NA
    res03 <-res02_clean
    sum(res03$SYMBOL == "NA") 
    [1] 0
    res03<- res03[order(res03 $log2FoldChange
    , decreasing = TRUE), ]

    write.csv(res03, file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/DESeq2_results_cancer_normal_clean.csv", row.names = FALSE) # có loại NA symbol, loại NA padj, lọc padj < 0.05
    # Số lượng gene upregulation và downregulation
    # lọc file up gene và down gene
    res03_df <- as.data.frame(res03)

    up_genes03 <- res03_df[res03_df$padj < 0.05 &
                        res03_df$log2FoldChange > 1, ]
    up_genes03 <- up_genes03[order(up_genes03 $log2FoldChange
    , decreasing = TRUE), ]

    down_genes03 <- res03_df[res03_df$padj < 0.05 &
                        res03_df$log2FoldChange < -1, ]
    down_genes03 <- down_genes03[order(down_genes03 $log2FoldChange
    , decreasing = TRUE), ]

    cat("Up genes:", nrow(up_genes03), "\n") #Up genes: 3773

    cat("Down genes:", nrow(down_genes03), "\n") #Down genes: 5562

    write.csv(up_genes03, file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/DESeq2_up_genes03.csv", row.names = FALSE)
    write.csv(down_genes03, file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/DESeq2_down_genes03.csv", row.names = FALSE)

    # Top 10 gene up_regulated và 10 gene down_regulated
    library(dplyr)

    top10_up <- up_genes03 %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice(1:10)
    top10_up[, c("SYMBOL", "log2FoldChange", "padj")]
    write.csv(top10_up,
            file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/top10_up_by_LFC.csv",
            row.names = FALSE)


    top10_down <- down_genes03 %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice(1:10)
    top10_down[, c("SYMBOL", "log2FoldChange", "padj")]
    write.csv(top10_down,
            file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/top10_down_by_LFC.csv",
            row.names = FALSE)

    # Lọc gene marker ("MGAT1", "BAG1", "MARCKS", "OGFR", "GIGYF1", "FKBP8")
    gene_list <- c("MGAT1", "BAG1", "MARCKS", "OGFR", "GIGYF1", "FKBP8")

    genes_in_res03 <- res03 %>%
    filter(SYMBOL %in% gene_list) # Lọc xem trong res03 có những gene này không

    genes_in_res03[, c("SYMBOL", "log2FoldChange", "padj")]

    genes_sorted <- genes_in_res03[, c("SYMBOL", "log2FoldChange", "padj")]

    setdiff(gene_list, res03$SYMBOL) # ktra gene không có 

    write.csv(genes_in_res03,
            file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/genes_interest_in_res03.csv",
            row.names = FALSE)

