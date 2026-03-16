# Differential Expression Analysis
## DEseq2 
### 1. Tạo metadata (sample and condition)
samples <-c("cancer_164","cancer_170","cancer_260","normal_1","normal_2","normal_3")
conditions <- ifelse(grepl("cancer", samples), "cancer", "normal")
df <- data.frame(sample=samples, condition=conditions)
print(df)
write.table(df, file="samples.txt", sep="\t", quote=FALSE, row.names=FALSE) #lưu file metadata
### 2. Tạo count matrix: tách cột các file abundance.tsv và gộp lại thành file count mới
python 3 
import pandas as pd
#### Đọc file abundance.tsv và tách cột 
    RNA_170 = pd.read_csv('/Volumes/Extreme SSD/Thesis/RNAseq_Test/Kallisto/170/abundance.tsv',sep = '\t')
    RNA_260 = pd.read_csv('/Volumes/Extreme SSD/Thesis/RNAseq_Test/Kallisto/260/abundance.tsv',sep = '\t')
    RNA_N1=pd.read_csv('/Volumes/Extreme SSD/Thesis/RNAseq_Test/Kallisto/N1/abundance.tsv',sep = '\t')
    RNA_N2=pd.read_csv('/Volumes/Extreme SSD/Thesis/RNAseq_Test/Kallisto/N2/abundance.tsv',sep = '\t')
    RNA_N3=pd.read_csv('/Volumes/Extreme SSD/Thesis/RNAseq_Test/Kallisto/N3/abundance.tsv',sep = '\t')
    RNA_164=pd.read_csv('/Volumes/Extreme SSD/Thesis/RNAseq_Test/Kallisto/167/abundance.tsv',sep = '\t’)
#### Tạo cột 'gene_id' và 'gene name'. dữ liệu được tách ra từ cột 'target_id' của fiel abundance.tsv
    RNA_164['gene id'] = [i.split('|')[1] for i in RNA_164['target_id']]
    RNA_164['gene name'] = [i.split('|')[5] for i in RNA_164['target_id']]
    RNA_170['gene id'] = [i.split('|')[1] for i in RNA_170['target_id']]
    RNA_170['gene name'] = [i.split('|')[5] for i in RNA_170['target_id']]
    RNA_260['gene id'] = [i.split('|')[1] for i in RNA_260['target_id']]
    RNA_260['gene name'] = [i.split('|')[5] for i in RNA_260['target_id']]
    RNA_N1['gene id'] = [i.split('|')[1] for i in RNA_N1['target_id']]
    RNA_N1['gene name'] = [i.split('|')[5] for i in RNA_N1['target_id']]
    RNA_N2['gene id'] = [i.split('|')[1] for i in RNA_N2['target_id']]
    RNA_N2['gene name'] = [i.split('|')[5] for i in RNA_N2['target_id']]
    RNA_N3['gene id'] = [i.split('|')[1] for i in RNA_N3['target_id']]
    RNA_N3['gene name'] = [i.split('|')[5] for i in RNA_N3['target_id']]
#### Tạo data mới chỉ gồm các cột geneid, gene name, est_count và TPM
    RNA_164 = RNA_164[['gene id','gene name','est_counts', 'tpm']]
    RNA_170 = RNA_170[['gene id','gene name','est_counts', 'tpm']]
    RNA_260 = RNA_260[['gene id','gene name','est_counts', 'tpm']]
    RNA_N1 = RNA_N1[['gene id','gene name','est_counts', 'tpm']]
    RNA_N2 = RNA_N2[['gene id','gene name','est_counts', 'tpm']]
    RNA_N3 = RNA_N3[['gene id','gene name','est_counts', 'tpm’]]
#### Tạo thêm cột sample và gán điều kiện của mẫu 
    RNA_164['sample']='cancer_164'
    RNA_170['sample']='cancer_170'
    RNA_260['sample']='cancer_260'
    RNA_N1['sample']='normal_1'
    RNA_N2['sample']='normal_2'
    RNA_N3['sample']='normal_3'
#### Tạo count matrix
    data = pd.concat([RNA_164,RNA_170,RNA_260,RNA_N1,RNA_N2,RNA_N3], axis=0) #gộp các mẫu thành 1 data
    Data
    counts = data.groupby(['gene id', 'sample'])['est_counts'].sum().unstack() # tạo matrix 
    counts = counts.round().astype(int) # chỉ lấy số nguyên của est_counts
    counts.to_csv("count_matrix.csv") # lưu file count_matrix

### 3. Thực hiện chạy DEseq2
conda create -n deseq2_env -c bioconda -c conda-forge r-base bioconductor-deseq2 bioconductor-tximport
conda activate deseq2_env # cài pagekage DEseq2
R
library(DESeq2)
count_matrix <- read.table(
  "/Volumes/Extreme SSD/Thesis/RNAseq_Test/counts_matrix_integer.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
samples <- read.table(
  "/Volumes/Extreme SSD/Thesis/RNAseq_Test/samples.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = TRUE
)
ncol(count_matrix)
nrow(samples) #kiểm tra số cột của count_matrix bằng với số hàng của sample

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = samples,
                              design = ~ condition) 
dds <- DESeq(dds)   #chạy lệnh DESeq2
res <- results(dds, contrast = c("condition", "cancer", "normal")) #lấy kết quả so sánh giữa mẫu cancer và normal
res <- res[order(res$padj),] # sắp xếp thứ tự thep padj
write.csv(as.data.frame(res),
          file = "/Volumes/Extreme SSD/Thesis/RNAseq_Test/DESeq2_results_cancer_vs_normal.csv") # lưu file thành tệp csv




