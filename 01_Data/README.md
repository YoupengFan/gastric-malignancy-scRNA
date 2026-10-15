# 01_Data — 外部数据文件

本目录存放管线引用的外部数据。**仓库只随带 `GSE150290_genes_hg19.tsv`（798K）**，大文件不上传，由复现者按下方来源自行下载后放入。

> 04 双种子的种子基因集（`04_Bulk_DEG/data/step1_top50_*`）已随仓库。正常复现无需下载任何大文件；仅当需要从头重推导种子基因时才需 `tcga_RSEM_gene_tpm.gz` / `TCGA.STAD.sampleMap_HiSeqV2.gz`。

| 文件/目录 | 大小 | 随仓库 | 来源 |
|------|------|:---:|------|
| `GSE150290_genes_hg19.tsv` | 798K | ✅ | ENSEMBL→symbol 映射（32,738 行，取自 GSE150290 的 genes.tsv），04 双种子 DEG 的基因名映射 |
| `raw/` | 1.3 GB（tar）/ 解压后 ~5.6G | ❌ | GEO [GSE150290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150290) 的 `GSE150290_RAW.tar`（52 样本，GSM4546300–GSM4546351，Kim et al. 2022 NPJ Precis Oncol），按 01_QC/README 约定的逐样本目录结构解压 |
| `tcga_RSEM_gene_tpm.gz` | 707M | ❌ | UCSC Xena Toil hub：`https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz`（04 主种子输入） |
| `TCGA.STAD.sampleMap_HiSeqV2.gz` | 24M | ❌ | UCSC Xena：`https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.STAD.sampleMap_HiSeqV2.gz`（04 测试种子输入） |
