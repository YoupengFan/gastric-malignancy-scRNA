# 10_R_pipeline — R 复现：Zhang et al. 2021 Gut 迭代 K-means 恶性细胞分类

> 本目录是 issue #26 的答复代码：在 **R/Seurat** 生态下忠实复现
> *Zhang M, et al. Dissecting transcriptional heterogeneity in primary gastric
> adenocarcinoma by single cell RNA sequencing. Gut 2021;70(3):464–475* 的
> "malignant vs non-malignant epithelium" 迭代分类方法，并做批判性检验。

## 数据源声明

原论文数据（GSA HRA000051）为申请制，不向学生/公众直接开放。本复现改用
**GSE150290**（Kim et al. 2022 *NPJ Precis Oncol*，韩国 NCC，52 样本，10x v2）
——同样是胃癌 paired tumor/adjacent 单细胞数据。**好的方法应该能在另一批
数据上复现**，这也是本复现的意义所在。

## 与原论文参数对照

| 步骤 | Zhang 2021（supp001 p3–5） | 本复现 | 说明 |
|------|------|------|------|
| QC 下限 | UMI≥1001, genes≥501, MT<50% | **一致** | 忠实保留 |
| QC 上限 | **UMI≤6000** | **UMI≤15000** | ⚠️ 放宽，见下"批判点" |
| 看家基因过滤 | mean<2（curated list） | **省略** | 论文未公开列表，无法忠实复现 |
| 归一化 | total-UMI 归一 + scale 回归总 UMI | log1p(×M/lib)，M=median | 时代共识替代 |
| HVG | 0.125–3 均值 + 方差>0.5 | vst 3000→回归→残差方差重排→2000 | 时代共识 + 技术噪音黑名单（MT-/核糖体/MALAT1） |
| 降维/聚类 | top20 PCs + tSNE + DBClustDimension | 50 PCs + UMAP + Leiden r=0.5 | 时代共识（DBSCAN 系密度聚类在 2021 后已非主流） |
| 批次矫正 | 无 | **Harmony（sample）** | 时代共识；原论文 9 患者 5 簇无矫正（见批判节） |
| 初始 DEG | TCGA Xena 配对 limma，top50/方向，adj.p<0.01 | **一致**（tcga_RSEM_gene_tpm + plain limma + log2(TPM+1) + 低表达过滤） | 数据集 ID 与论文方法正文一致（Xena Toil hub RSEM TPM）；TPM 非计数 → plain limma（无 voom）；输入变换与过滤按论文 Table S1 逐位校准：log2(TPM+1) 下 CST1 5.781 / PGA3 −3.965 与论文表逐位一致、top50 重叠 94/100（残留 6 基因 = top50 边界尾部互换 + 1 个本地映射缺失）；Fig 3A 的 \|logFC\|=1 虚线仅为火山图图注（阈值线与上色边界）、不参与筛选（论文 supp Table 2 迭代内 DEG 含 \|FC\|<1 显著基因佐证）；**双种子结构**：早期误用 Xena sampleMap HiSeqV2 的种子保留为测试种子（敏感性对照，纠错史见 04_Bulk_DEG/README.md） |
| 打分 | 模块打分 | AddModuleScore（25 bins × 100 同 bin 对照） | 一致 |
| 二分 | K-means k=2 | k=2, nstart=25 | 一致 |
| 收敛 | 标签一致（原文 90%） | **≥90%** | 一致 |
| 迭代内 DEG | 单细胞级 | **pseudobulk（patient×label 配对）+ Wilcoxon** | 时代共识（单细胞伪重复假阳性是 2020 后公认问题） |
| 恶性率 | 活检 79.7% / 手术 29.4%（S3D） | 按组织报告 | 对照检验 |

## 上皮定义声明（与论文一致）

原论文的"上皮"是**狭义定义**：九大谱系中仅 `epithelium`（10,411 细胞，37.6%，
marker = EPCAM/KRT18/KRT8）参与恶性/非恶性迭代分类；**壁细胞**（parietal cell，
215 细胞，ATP4A）与**内分泌细胞**（endocrine cell，486 细胞，CHGA）在论文中是
独立谱系，**从未被纳入上皮类型**（主文 figure 1C–E 谱系注释 + Figure 2 标题
"Classification of 10 411 epithelial cells"；最终 5,635 恶性 + 4,776 非恶性 =
10,411）。论文九谱系中也没有 chief cell —— PGC+ 主细胞表达 EPCAM/KRT8，
落在论文的 "epithelium" 内。

因此**主流程 05 即此忠实口径 = Epithelial + Chief**，排除 Parietal 与
Enteroendocrine **与论文一致，并非我们遗漏**。六 Runs 扫描中 run2（4 类型，
广义口径）保留为敏感性对照（Python 版 46,554 的 R 对应物），用于展示
"上皮定义不影响收敛主题（IM 轴）"这一结论的稳健性——放宽定义会拉入
GHRL/CHGA 神经内分泌轴（EE 伪影，见结果节）。

## 批判点（QC 6000 UMI 上限）

在本数据（GSE150290）中，UMI 6000 上限会系统性截断大细胞：

- 巨液滴普查：每样本 19–384 个 >15k UMI 液滴，肿瘤侧尾部更重（P17B max 250,637）
- 断崖检验：肿瘤平滑无断崖 vs 对照 30k 处消失（P27A 零个 >15k）
- 组间频率：肿瘤 88.9 vs 对照 12.8 /10 万液滴（p=0.0225）；>100k 液滴 45 vs 0
- DNA 级锚：Kim Supp Table 6 突变上皮细胞集中在深液滴（>20k 才有热点覆盖）

**推论**：原论文的 6000 上限可能恰好在第一步把（偏大的）恶性细胞排除在外，
这可能是其"恶性=肠化生"结论的部分成因。原论文自身数据可能有 6000 附近断崖；
但"在我们数据中 6000 是败笔"是本复现放宽上限的理由。

（其余病理/术语批判见 issue 答复正文，不在本 README 展开。）

## 快速开始

```bash
export PATH="$HOME/softwares/bin:$PATH"
Rscript 01_QC/scripts/qc_filter.R
Rscript 02_Normalize_Batch/scripts/normalize_batch.R
Rscript 03_Cluster_Annotate/scripts/cluster_annotate.R
Rscript 04_Bulk_DEG/scripts/step1_bulk_deg.R
Rscript 05_Iterative_Classify/scripts/step2_iterative_classify.R
```

绘图（按工作流拆分：WSL 做重计算，Windows 原生 R 渲染——ragg 高质量抗锯齿）：

```bash
# ── WSL（计算/聚合）──
export PATH="$HOME/softwares/bin:$PATH"
Rscript 01_QC/scripts/qc_plot_compute.R            # ~12 min, raw 矩阵单遍
Rscript 02_Normalize_Batch/scripts/validate_pcs_compute.R  # ~10 min, 6.8GB RDS

# ── Windows 原生 R 4.5.1（渲染）──
"E:\R-4.5.1\bin\Rscript.exe" 01_QC\scripts\qc_plot_figure.R
"E:\R-4.5.1\bin\Rscript.exe" 02_Normalize_Batch\scripts\validate_pcs_figure.R
"E:\R-4.5.1\bin\Rscript.exe" 02_Normalize_Batch\scripts\viz_library_scaling.R  # 计算轻量, 全在 Windows
```

产出图：

- `01_QC/figures/qc_metrics_overview.png` — 过滤前 QC 四面板（含论文 6,000 上限对照线）
- `02_Normalize_Batch/figures/library_and_cells_overview.png` — 每样本文库大小分布 + 细胞数（scaling_factor_by_sample 因 CP10K 切换于 2026-10-08 删除）
- `02_Normalize_Batch/figures/pc_validation.png` — PC 质量四面板（方差/批次 η²/CC-MT 残差/loading 组成 + Harmony 指标）

两个核心函数（issue 要求的"第一步一个函数、第二三步合成一个函数"）：

- `step1_bulk_deg(bulk_expr, sample_table, n_top=50, padj=0.01, ...)` → `list(tumor_genes, normal_genes)` — TCGA 配对 plain limma DEG（主种子），输入 log2(TPM+1)、低表达过滤 TPM>1 in ≥20%（`04_Bulk_DEG/scripts/step1_bulk_deg.R`）
- `step1_testseed_deg(bulk_expr, sample_table, n_top=50, padj=0.01)` → `list(tumor_genes, normal_genes, deg_table, design_info)` — 测试种子生成器（HiSeqV2 + limma-voom 母本配方，敏感性对照；`04_Bulk_DEG/scripts/sensitivity/step1_bulk_deg_testseed.R`）
- `step2_iterative_classify(counts, cell_meta, tumor_genes, normal_genes, n_top=50, convergence=0.90, max_iter=10, ...)` → `list(labels, scores, iteration_log, final_signatures)` — AddModuleScore → K-means k=2 → pseudobulk 配对 Wilcoxon → 迭代（`05_Iterative_Classify/scripts/step2_iterative_classify.R`）。**`counts` 必须是原始 UMI counts**（未归一化）：函数内部按传入子集自行 log1p(CP10K) 归一化，与 02 的归一化结果相互独立

## 目录结构

```
.
├── README.md
├── 01_QC/                  # UMI 1001–15000, genes≥501, MT<50%（无看家基因约束）
├── 02_Normalize_Batch/     # log1p(CP10K) + CC/MT 回归 + HVG 黑名单 + Harmony
├── 03_Cluster_Annotate/    # Leiden r=0.5 + Hybrid v2（Kim+Zhang 12 类）注释 + 红细胞/死细胞剔除
├── 04_Bulk_DEG/            # ★ 函数1: limma only（TCGA-STAD 配对，双种子：主线 Toil + 测试种子 HiSeqV2）
├── 05_Iterative_Classify/  # ★ 函数2: 90% 收敛, pseudobulk+Wilcoxon
└── 06_Results/             # 汇总结果（待生成）
```

## 环境

- WSL2 Ubuntu 22.04, R 4.5.1, Seurat 5.x, harmony, limma, edgeR（仅 voom 输入）, Matrix
- 原始数据：GSE150290 raw 10x 矩阵（GEO 下载，见 data_reference.json）
- TCGA-STAD：Toil hub tcga_RSEM_gene_tpm（主线种子）+ Xena sampleMap HiSeqV2（测试种子，仓库 01_Data/）

## 结果（管线重跑后回填）

- [x] 01_QC 细胞数: **160,359**（UMI 1001–15000 / genes≥501 / MT<50%）
- [x] 03 注释类型分布: 159,769 cells（剔除 590 红细胞），23 簇 11 类；
      Plasma 27.8%（Python 版 27.4%）、T 16.4、B 11.8、Fibro 5.1 等逐类一致；
      NK_cells 无主导簇缺失
- [x] 05 迭代收敛轨迹 + 恶性率（按组织，忠实口径 Epi+Chief）:
      **4 轮收敛 94.66%**；27,460 恶性 / 22,697 非恶性（54.7%）；
      AdjA 48.4% / TumorB 68.0%；轨迹 20,399 → 27,625 → 29,298 → 27,460
- [x] 六 Runs 敏感性对比（上皮定义 × 基因选择 × 看家黑名单，主种子）：

| run | 上皮定义 | 基因选择 | HK 黑名单 | n | 恶性% | Ctrl% | AdjA% | TumorB% | 轮数 | agree |
|-----|------|------|:---:|-----:|-----:|-----:|-----:|-----:|:---:|-----:|
| run1 | Epi+Chief | logFC | ✗ | 50,157 | 54.7 | 36.4 | 48.4 | 68.0 | 4 | 0.9466 |
| run2 | 4 类型 | logFC | ✗ | 54,638 | 53.2 | 34.9 | 46.1 | 68.2 | 4 | 0.9429 |
| run3 | Epi+Chief | expression | ✓ | 50,157 | 52.6 | 35.7 | 46.5 | 65.1 | 4 | 0.9417 |
| run4 | 仅 Epi | expression | ✓ | 29,213 | 56.2 | 47.0 | 54.5 | 58.7 | 3 | 0.9305 |
| run5 | 仅 Epi | logFC | ✓ | 29,213 | 57.0 | 43.2 | 54.4 | 60.5 | 3 | 0.9293 |
| run6 | Epi+Chief+低置信剔除 | logFC | ✗ | 43,234 | 60.7 | 39.7 | 53.4 | 74.8 | 4 | 0.9584 |

要点（完整版见 05_Iterative_Classify/README.md）：

- **恶性主题全部 = IM 轴**（S100A6/PHGR1/APOA1/FABP1/TFF3/REG4 + 角蛋白）；
  非恶性签名 5 runs 完全一致（PGC/LIPF/LYZ/TFF1/REG3A/TFF2/AGR2/REG1A）
- **Ctrl 恶性率非 0 且高**（34.9–47.0%）：对照组为慢性胃炎，推测含少量 IM
  细胞（慢性胃炎→肠化生是 Correa 级联的常见伴随病变），方法按身份轴一律划为
  "恶性"——"恶性=肠化生混淆"批判论点；run4/5（仅 Epi）最高，
  上皮成分更稀 → 分类更不稳定
- **测试种子对照**（HiSeqV2 早期误用种子，run1 + run6 两种配置）：run1 对照
  恶性率 36.4% → 12.5% ≈ **3× 漂移**——切点无生物锚定（种子换源即整体平移）；
  但收敛主题同为 IM 轴、最终签名 J(non) 0.82、label_agree 0.79——方法缺陷跨
  数据源复现。测试种子的角色 = 检验 bulk 种子集的引导作用（数据集构成敏感性
  由主种子六 Runs 矩阵承担）
- **vs Python 版**（8.9 / 10.0 / 58.4 / 63.4 / 12.8%，run3/4 逃逸到 Chief 轴）：
  R 版无逃逸——pseudobulk 配对 Wilcoxon 迭代 DEG 比单细胞级 DEG 稳定；
  **收敛主题一致 = IM 轴**（方法缺陷跨实现复现）
- **数据版本注记**：早期五 Runs 记录（50,300/44,880/35,015 细胞）为旧数据版本，
  01/02/03 上游重算后恶性率系统性上移（如 4 类型口径 33.4% → 53.2%），本表为
  当前版本
