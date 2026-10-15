# 04_Bulk_DEG — TCGA-STAD 配对 bulk DEG（函数1）

## 双种子结构

| 种子 | 数据源 | 配方 | 角色 |
|------|------|------|------|
| **主种子（主线）** | Toil hub `tcga_RSEM_gene_tpm` | plain limma（TPM 非计数 → 无 voom），log2(TPM+1) 变换，TPM>1 in ≥20% 过滤，33 配对 | 论文指定数据源，主流程使用 |
| **测试种子（敏感性）** | Xena sampleMap `TCGA.STAD.sampleMap_HiSeqV2.gz` | limma-voom（round(2^x−1) → counts>1 in ≥20% → TMM → voom → eBayes(trend, robust)），32 配对 | 迭代聚类种子稳健性对照 |

## 主种子数据源核实与 Table S1 复现

- **数据源**：论文方法正文原引 = UCSC Xena Toil hub `tcga_RSEM_gene_tpm`
  （Toil 流水线重算的基因级 RSEM TPM，值 = log2(TPM+0.001)）。STAD
  子集 = Toil 文件自身全部 STAD 前缀列 = **33 配对**——比"HiSeqV2 样本
  交集"口径多 TCGA-BR-7703-11 与 TCGA-BR-8060-11 两个 normal 样本；
  HiSeqV2 矩阵仅用于提供 STAD 患者前缀集合
- **输入变换 log2(TPM+1)**：文件值 log2(TPM+0.001) 直接输入会把地板
  基因（TPM<1 的低表达 lncRNA）的 logFC 放大 ~1.6×（CST1 9.118 vs
  论文 5.781）；log2(TPM+1) 输入下 CST1 5.781 / PGA3 −3.965 与论文
  Table S1 逐位一致，top50 重叠 94/100
- **低表达过滤 TPM>1 in ≥20% 配对样本（n≈15,800）**：论文方法正文
  未写明，由 Table S1 证据推定——无过滤 31,439 基因时 GKN1/GKN2/PGA3
  adj.p 0.010-0.014 恰好越线；过滤后 ~0.005-0.007 与论文表
  0.0055/0.0067/0.0079 吻合，且论文表内无地板 lncRNA 霸榜
- **DEG 口径 = padj<0.01，无 FC 截断**（论文方法正文只写 adj.p<0.01；
  论文 supp Table 2 迭代内 DEG 含 |FC|<1 显著基因佐证）；每方向按
  |logFC| 取前 50。火山图 ±1 logFC 虚线 = 论文 supp007 Fig A 画法（阈值线 +
  着色边界），不参与基因筛选
- **Ensembl → symbol**：本地映射（`01_Data/GSE150290_genes_hg19.tsv`，
  GSE150290 原始 10x 目录的 hg19 注释副本，32,738 基因，离线，随仓库
  分发）；未映射基因丢弃；重名 symbol 按平均表达取最大的一行

## 测试种子的纠错史（诚实交代）

早期（R 管线搭建期）误将 UCSC Xena sampleMap 的 HiSeqV2 矩阵
（值 = log2(RSEM+1)）当作论文数据源，由此产出第一批种子基因。
后面对照论文 Table S1 逐位核验以及method描述，确认论文实际使用
Toil hub 的 `tcga_RSEM_gene_tpm`
（CST1 5.781 / PGA3 −3.965 与论文表逐位一致），主线随即切换
到 Toil + plain limma（见 scripts/step1_bulk_deg.R）。

HiSeqV2 种子**保留为测试种子**（testseed），而不是丢弃。

## 运行

```bash
export PATH="$HOME/softwares/bin:$PATH"
# 主线（主种子）
Rscript 04_Bulk_DEG/scripts/step1_bulk_deg.R
# 测试种子（敏感性对照，输出 data/testseed/）
Rscript 04_Bulk_DEG/scripts/sensitivity/step1_bulk_deg_testseed.R
```

## 数据文件

```
data/
├── top50_tumor_limma.txt                # 主种子（主线输入，33 配对 Toil）
├── top50_normal_limma.txt
├── deg_limma.csv                        # 全基因 DEG 表
├── bulk_deg_params.json                 # 参数记录
├── paper_S1/                            # 论文 supp Table S1 提取（核验参考，非管线产物）
│   ├── paper_S1_top50_tumor.txt
│   └── paper_S1_top50_normal.txt
└── testseed/
    ├── top50_tumor_limma.txt            # 测试种子（32 配对 HiSeqV2，六次 testseed Run 输入）
    ├── top50_normal_limma.txt
    └── bulk_deg_testseed_params.json

figures/
├── fig_bulk_deg.png            # 主种子火山图（Toil + plain limma）
└── fig_bulk_deg_testseed.png   # 测试种子火山图（HiSeqV2 + limma-voom）
                                 # 画法同款，两张可作补充材料（数据源×检验链对照）
```
