# Functional analysis
## GSEA
    # Cài package
    R
install.packages("BiocManager")
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}



    library(clusterProfiler)
    library(org.Hs.eg.db)
    # 0. reading in data from deseq2
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$gene_id_clean <- sub("\\..*", "", res_df$gene) #loại ENSG version

    # 1. tạo gene_list cho GO
    gene_list <- res_df$log2FoldChange
    names(gene_list) <- res_df$gene_id_clean # gán tên gene_id cho vector gene list
    gene_list<-na.omit(gene_list) # Loại LFC có giá trị NA
    write.csv(geneList_df , 
          file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/genelist_GSEA-2.csv", 
          row.names = FALSE)

        # Gene_list có nhiều SYMBOL trùng nhau -> thực hiện gộp SYMBOL thành 1 giá trị và chỉ lấy giá trị có LFC lớn nhất. Tức là, lấy 1 giá trị LFC lớn nhất cho mỗi gene
        res_df2_unique <- res_df2 %>%
        filter(!is.na(log2FoldChange)) %>%     # 👈 loại LFC=NA trước
        group_by(SYMBOL) %>%
        summarise(
        log2FoldChange = max(log2FoldChange, na.rm = TRUE),
        .groups = "drop"
        )
        write.csv(res_df2_unique , 
          file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/ref_df2_gop_genesymbol_2.csv", 
          row.names = FALSE)


        geneList <- res_df2_unique$log2FoldChange
        names(geneList) <- res_df2_unique$SYMBOL
        geneList <- sort(geneList, decreasing = TRUE)  #geneLIst gồm 2 cột: LFC và Symbol

        geneList_df <- data.frame(
        Gene_symbol = names(geneList),
        log2FoldChange = geneList
        )
        write.csv(geneList_df , 
                file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/genelist_GSEA_gop_genesymbol_2.csv", 
                row.names = FALSE)
        any(duplicated(names(geneList))) # Ktra lại còn dup ko -> FALSE -> OKE

    

    # 2. Chạy GSEA - GO 
    gse_bp <- gseGO(geneList=gene_list, 
                ont ="BP",  # Biological Process
                keyType = "ENSEMBL", 
                minGSSize = 3, # loại những gene set < 3 gene, tránh pathway quá nhỏ, gây nhiễu
                maxGSSize = 800, # loại những gene set > 800 gene, tránh pathway quá lớn, quá chung
                pvalueCutoff = 0.05, # cutoff của pvalue = 0.05
                verbose = TRUE, # Không in log dài dòng khi chạy
                OrgDb = org.Hs.eg.db) # Organism là Human

    ## Kết quả cảnh báo: There are ties in the preranked stats (44,29%)
    ## Tránh cảnh báo trùng giá trị LFC này
    set.seed(123)
    gene_list2 <- gene_list + rnorm(length(gene_list), 0, 1e-6) # cộng thêm giá trị rất rất nhỏ cho LFC
    gene_list2 <- sort(gene_list2, decreasing = TRUE)

    ## chạy lại GSEA_BP
        gsea_bp <- gseGO(
        geneList     = gene_list2,
        ont          = "BP", # Biological Process
        keyType      = "ENSEMBL",
        OrgDb        = org.Hs.eg.db,
        minGSSize    = 3,
        maxGSSize    = 800,
        pvalueCutoff = 0.05,
        verbose      = TRUE
        ) 
        gse_bp_sorted <- gsea_bp@result[order(gsea_bp@result$NES, decreasing = TRUE), ] # Sắp xếp theo thứ tự giảm dần của NES
        head(gse_bp_sorted[, c("ID","Description","NES","p.adjust")])
        write.csv(gse_bp_sorted, 
                file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/gsea_go_bp.csv", 
                row.names = FALSE) # Lưu file

    ## Chạy GSEA_MF
        gsea_mf <- gseGO(
        geneList     = gene_list2,
        ont          = "MF", #Molecular function
        keyType      = "ENSEMBL",
        OrgDb        = org.Hs.eg.db,
        minGSSize    = 3,
        maxGSSize    = 800,
        pvalueCutoff = 0.05,
        verbose      = TRUE
        )
        gse_res_sorted <- gsea_mf@result[order(gsea_mf@result$NES, decreasing = TRUE), ]
        head(gse_res_sorted[, c("ID","Description","NES","p.adjust")])
        write.csv(gse_res_sorted, 
                file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/gsea_go_mf.csv", 
                row.names = FALSE)

    ## chạy GSEA_CC
        gsea_cc <- gseGO(
        geneList     = gene_list2,
        ont          = "CC", #Cellular component
        keyType      = "ENSEMBL",
        OrgDb        = org.Hs.eg.db,
        minGSSize    = 3,
        maxGSSize    = 800,
        pvalueCutoff = 0.05,
        verbose      = TRUE
        )
        gse_cc_sorted <- gsea_cc@result[order(gsea_cc@result$NES, decreasing = TRUE), ]
        head(gse_cc_sorted[, c("ID","Description","NES","p.adjust")])
        write.csv(gse_cc_sorted, 
                file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/gsea_go_cc.csv", 
                row.names = FALSE)

   
    # 3. Tạo gene_list cho KEGG
    ids<-bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db) # chuyển gene id dạng ENSEMBL dang ENTREZID
    dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),] # chỉ map 1 ENSEMBL = 1 ENTREZID tránh lỗi 1 ENSEMBL map với nhiều ENTREZID
    df2 = res_df[res_df$gene_id_clean %in% dedup_ids$ENSEMBL,] #lọc data cho những gene map thành công
    df2$Y = dedup_ids$ENTREZID # thêm cột ENTREZID

    kegg_gene_list <- df2$log2FoldChange # tạo vector kegg_gene_list theo LFC
    names(kegg_gene_list) <- df2$Y # tên hàng là ENTREZID
    kegg_gene_list<-na.omit(kegg_gene_list) # loại NA 
    kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE) # sort data

    #gộp các ENTREZID trùng -> giữ LFC lớn nhất
    kegg_df <- data.frame(
    ENTREZID = names(kegg_gene_list),
    log2FC   = kegg_gene_list   
    )
    kegg_df_unique <- kegg_df %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarise(log2FC = max(log2FC, na.rm = TRUE)) %>%
    dplyr::ungroup() 

    kegg_gene_list2 <- kegg_df_unique$log2FC
    names(kegg_gene_list2) <- kegg_df_unique$ENTREZID
    kegg_gene_list2 <- sort(kegg_gene_list2, decreasing = TRUE)

    kegg_df_2 <- data.frame(
    ENTREZID = names(kegg_gene_list2),
    log2FC   = kegg_gene_list2
    )
    write.csv(kegg_df_2, 
            file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/kegg_gene_list.csv", 
            row.names = FALSE)

    #4. Chạy GSEA_KEGG
    kk2 <- gseKEGG(
    geneList     = kegg_gene_list2,
    organism     = "hsa",
    minGSSize    = 3,
    maxGSSize    = 800,
    pvalueCutoff = 0.05,
    keyType      = "ncbi-geneid"
    )
    ## Kết quả cảnh báo: There are ties in the preranked stats (26,52%)
    ## Tránh cảnh báo trùng giá trị LFC này
    set.seed(123)
    kegg_gene_list2 <- kegg_gene_list2 + rnorm(length(kegg_gene_list2), 0, 1e-6)
    kegg_gene_list2 <- sort(kegg_gene_list2, decreasing = TRUE)

    ## chạy lại GSEA KEGG
    ## lưu file 
    gse_kegg_sorted <- kk2@result[order(kk2@result$NES, decreasing = TRUE), ]
    head(gse_kegg_sorted[, c("ID","Description","NES","p.adjust")])
    write.csv(gse_kegg_sorted, 
            file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/gsea_kegg.csv", 
            row.names = FALSE)










