# Issue #26 — 迭代恶性细胞分类：数据处理概览

目录




> **参考论文**: Zhang M et al. Dissecting transcriptional heterogeneity in primary gastric adenocarcinoma by single cell RNA sequencing. *Gut* 2021;70:464–475. DOI: [10.1136/gutjnl-2019-320368](https://doi.org/10.1136/gutjnl-2019-320368)
> **论文 method：**(/00_Doc/zhang2021_paper/supply_material/gutjnl-2019-320368supp001.pdf)
---

## 2. 框架 · 我们的分析系统设计

基于原论文的方法学框架，我们设计了一个 **双跑道、单汇合** 的分析系统（详见项目 [README.md](../README.md)）：

```
══════════════════════════════════════════════════════════
跑道 A: scRNA-seq 预处理 → 上皮细胞鉴定
══════════════════════════════════════════════════════════

  01_Data/ ──→ 02_Preprocessing/ ──→ 07_Clustering_Annotation/
  GSE150290      QC → 归一化          Marker gene scoring
  52 样本           → Harmony 校正      → 9 种细胞类型
  136K cells        → UMAP + Leiden    → 33,934 Epithelial cells

══════════════════════════════════════════════════════════
跑道 B: Bulk DEG → 种子基因集
══════════════════════════════════════════════════════════

  03_Bulk_DEG/
  TCGA-STAD (415 tumor + 35 normal; paired DEG on 32 matched pairs)
  → edgeR + DESeq2 + limma-voom (三方法 Venn consensus)
  → Top 50 Tumor-up + Top 50 Normal-up

══════════════════════════════════════════════════════════
汇合: 迭代恶性细胞分类
══════════════════════════════════════════════════════════

  08_Iterative_Classification/
  上皮细胞 (33,934) + 种子基因 (50+50)
  → AddModuleScore → K-means (k=2) → Pseudobulk DEG
  → 迭代至收敛 (≥95%) → 最终恶性/非恶性标签
══════════════════════════════════════════════════════════
```

---

## 3. 数据 · 为什么是 GSE150290 而不是 HRA000051

### 3.1 原论文数据不可获取

Zhang 2021 的原始 scRNA-seq 数据存放于 **GSA-Human**（国家基因组科学数据中心），accession **HRA000051**。尽管论文声明"data are available on reasonable request"，实际访问需要经过 **Data Access Committee (DAC) 审核**。

![HRA000051 GSA-Human 页面截图 — 显示受控访问](/00_Doc/zhang2021_paper/屏幕截图%202026-06-27%20193343.png)

*图 3.1a: HRA000051 在 GSA-Human 数据库的访问页面。数据标注为"Controlled Access"，需要向 DAC (Omics_Data) 提交申请方可下载。*

![NGDC 账户注册/申请受阻截图](/00_Doc/zhang2021_paper/屏幕截图%202026-06-27%20193611.png)

*图 3.1b: 数据访问申请界面。NGDC 注册要求机构邮箱及身份验证，学生身份无法独立完成 DAC 审核流程。*

### 3.2 GSE150290

我们选择了 **Kim et al. (2022) *NPJ Precision Oncology*** 的公开数据集 **GSE150290**（韩国 National Cancer Center, PMID: 35087207）：

**原始数据结构**

GSE150290 包含 **52 个 10x Genomics 样本**（GEO 页面概述为 "47 biopsies"，指 47 例癌症患者内镜活检：23 GC 肿瘤 + 24 癌旁正常；另有 5 例非癌对照 Pat25–Pat29，GEO 标记为 "Control"；合计 52）。全部为内镜活检（GEO: "biopsies"），来自 24 例胃癌患者（Pat01–Pat24，其中 Pat21 仅剩 A 样本，配对 B 未通过 QC 故实际肿瘤样本为 23 例）。每个样本提供标准的 10x 输出：barcodes、features、matrix。


> 📎 数据来源：[GSE150290 — NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150290)（Status: Public on Jan 30, 2022，无需注册即可访问元数据及下载 count matrix。）

选择 GSE150290 的优势是多重的：(1) 公开可获取，确保完全可复现；(2) 47 例癌症活检 vs 12 样本，统计功效显著增强；(3) 韩国独立队列为 Zhang 2021 方法提供了跨人群、跨机构的外部验证——这是原论文未曾进行的检验；

---

## 4. 结果 — 跑道 A：从原始数据到上皮细胞

### 4.1 Qc

**QC 阈值：** UMI：1,000 ~ 15,000 ; gene ≥ 500 ; MT ≤ 50%。

> 论文UMI上限为6000，在此数据集上截断大量细胞，因此选择更宽松的15,000  
> **关于线粒体 50% 阈值**：胃上皮细胞代谢活动活跃，胃上皮组织中可能含有大量的线粒体。因此论文采用了更宽松的 50% 标准。

![QC metrics overview](/01_QC/figures/qc_metrics_overview.png)

**Figure 1**  质控变量的分布图 (GSE150290).
**(A)** total UMI counts 直方图.
**(B)** gene numbers 直方图.
**(C)** UMI 计数的Rank图（双对数坐标系）。虚线标示了各计数阈值与相应 barcode 排序位置的交点。
**(D)** 基因数量与 UMI 总计数之间的关系，随机抽取 50,000 个的barcode；图中颜色代表线粒体基因占比（MT%）。


Data: GSE150290, 52 gastric cancer samples (10X Genomics). After filtering (min_genes = 500, min_counts = 1,000, max_counts = 6,000, MT% < 50%), **136,352 high-quality cells** retained.

### 4.2 归一化(CP10K + log1p)

![](/02_Normalize_Batch/figures/library_and_cells_overview.png)
**Figure 2. library and cell的分布图.**  上图为libraries直方图，下图为cells直方图

### 4.3 批次校正：Harmony PCA

多样本 scRNA-seq 数据整合中，技术批次效应会混淆生物学信号。我们使用 **Harmony** (v2.0.5) 在 PCA 空间进行迭代校正，以 `sample_id` 为批次变量。Harmony 在 17 次迭代后收敛，校正后 sample-level R² 从 0.153 降至 0.125（降低 18.3%）。

CC 评分 → HVG 3000（黑名单后精选）→回归 S/G2M/pct_mt（残差，不缩放）→ 残差方差重排 top 2000 → z-score（clip 10）→ PCA 50 → Harmony（key=sample；theta=2,
sigma=0.1, nclust=100）


![UMAP before and after Harmony](/02_Normalize_Batch/figures/pc_validation.png)
**Figure 3.1 Batch effect correction validation.** **(A)** 每细胞相邻细胞 k=20 是属于同样本的比例（post后降低），虚线为随机线**(B)** PC 1-30 样本相关性，矫正后降低

![UMAP before and after Harmony](/02_Normalize_Batch/figures/umap_batch_before_after.png)
**Figure 3.2 批次矫正 UMAP validation.** 抽取patient01、17、21 的 adjacent cancer侧 上皮类型细胞做示意图



Data: GSE150290, 52 gastric cancer samples (10X Genomics). 136,352 cells after QC filtering. Harmony: R package v2.0.5 (Korsunsky et al., 2019 *Nature Methods*).

### 4.4 聚类与细胞类型注释

使用 Harmony 校正后的前 50 PCs， Leiden 聚类（分辨率 r=0.5）。


我们采用原论文的 marker gene scoring 策略——对 9 种经典细胞类型，使用 canonical marker genes 计算每个 cluster 的平均表达得分，取 `argmax` 分配细胞类型（详见 [](../07_Clustering_Annotation/README.md)）。

![Cell type annotation](/03_Cluster_Annotate/figures/qc_umap_clean_nolabel.png)



*图 4.7: 细胞类型注释结果。9 种细胞类型在 UMAP 上呈清晰分区。*

| 细胞类型 | 细胞数 | 比例 |
|:---|:---:|:---:|
| Plasma cells | 38,144 | 28.0% |
| **Epithelial cells** | **33,934** | **24.9%** |
| T cells | 20,033 | 14.7% |
| B cells | 19,720 | 14.5% |
| Fibroblasts | 7,460 | 5.5% |
| NK cells | 6,074 | 4.5% |
| Myeloid | 4,952 | 3.6% |
| Endothelial cells | 3,477 | 2.6% |
| Mast cells | 2,558 | 1.9% |

**33,934 个上皮细胞**（EPCAM⁺, KRT8⁺, KRT18⁺）被提取，作为迭代分类的输入。这一数量与原论文的 10,411 个上皮细胞在数量级上相当，但来自更多样、更多样的样本群。

---

## 5. 结果 — 跑道 B：TCGA-STAD 种子基因集

### 5.1 数据与方法

从 UCSC Xena 下载 **TCGA-STAD** bulk RNA-seq 数据（HiSeqV2, RSEM normalized）。表达矩阵包含 **415 例肿瘤**和 **35 例正常**样本。取其中 **32 对**配对的肿瘤-正常样本进行差异表达分析（仅患者同时有 tumor 和 normal 样本入组），

![](/04_Bulk_DEG/figures/fig_bulk_deg.png)
**种子基因代表性成员**：

| 方向 | 代表基因 | 功能 |
|:---|:---|:---|
| Tumor-up (50) | CST1, HOXC10, CLDN3, CLDN4, KRT7, S100A2, REG4, TFF3, MUC13... | 增殖、紧密连接、肠型粘蛋白、癌胚抗原 |
| Normal-up (50) | S100G, ATP4A, ATP4B, PGC, LIPF, GKN1, GKN2, MUC5AC, TFF1, TFF2... | 胃酸分泌、消化酶、胃黏液保护 |

⚠️ 注意：种子基因来自 bulk 组织 DEG，包含非上皮来源的基因（如免疫细胞和基质细胞的差异表达）。这正是原论文需要迭代步骤的原因——bulk 种子中混杂的非上皮信号将在 sc-level 迭代中被逐步过滤。

---

## 6. 结果 — 汇合：迭代恶性细胞分类

### 6.1 算法实现

我们在 R 中完整实现了迭代分类框架（[`08_Iterative_Classification/scripts/step2_iterative_classify.R`](../08_Iterative_Classification/scripts/step2_iterative_classify.R)，~1,226 行），核心循环为：

```
种子基因 (50 Tumor-up + 50 Normal-up)
       │
       ▼
  AddModuleScore (25 bins, 100 ctrl per gene)
       │
       ▼
  2D scores × 33,934 epithelial cells
       │
       ▼
  K-means (k=2, nstart=25) ─→ "推定恶性" / "推定非恶性"
       │
       ▼
  Pseudobulk 聚合 (patient × condition)
       │
       ▼
     limma(paired design ~ patient + condition)
       │
       ▼
  CC 保护: 97 基因黑名单 + Spearman |ρ| > 0.25 过滤
       │
       ▼
  重新打分 → 检查收敛 (标签一致性 ≥ 95%)
       │         │
       │    ┌────┘
       │    ▼
       │  未收敛: 继续迭代
       │
       ▼
     收敛: 输出最终标签 + 基因签名
```

### 6.2 收敛轨迹

迭代在 **4 轮后收敛于 98.0% 一致性**：

| 轮次 | 恶性细胞 | 非恶性 | mal_score (恶性组) | mal_score (非恶性组) | 一致性 | 收敛 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 (种子) | 8,298 (24.5%) | 25,636 | 0.125 | 0.001 | — | ❌ |
| 1 | 22,602 (66.6%) | 11,332 | 0.074 | −0.114 | 55.6% | ❌ |
| 2 | 3,764 (11.1%) | 30,170 | 1.699 | −0.202 | 44.4% | ❌ |
| 3 | **3,480 (10.3%)** | **30,454** | **2.564** | −0.282 | **98.0%** | ✅ |

![Convergence trajectory](../08_Iterative_Classification/figures/convergence_trajectory.png)

*图 6.1: 收敛轨迹。左：恶性评分在两组间的分化（信号放大 20.5×，从 0.125 → 2.564）；中：恶性细胞比例变化（Iter 1 的 66.6% 为种子基因噪音所致的过度膨胀）；右：标签一致性（Iter 3 达到 98.0% 收敛）。*


### 6.3 最终 DEG 结果

最终 pseudobulk 三方法 DEG 鉴定出 **50 个恶性基因**和 **50 个非恶性基因**（adj.p < 0.01, 三方法 consensus ≥ 2/3）。

![Volcano plot: Final DEG](../08_Iterative_Classification/figures/volcano_deg_final.png)

*图 6.2: 最终迭代 DEG 火山图。x 轴为 log₂ fold change (malignant vs non-malignant)，y 轴为 −log₁₀(adjusted p-value)。红色点：恶性方向 (Top 50)；蓝色点：非恶性方向 (Top 50)。关键基因标注。所有 Top 50 基因均通过 CC 保护过滤（0 个细胞周期基因混入）。*

### 6.4 基因签名热图

![Gene heatmap: Top 30](../08_Iterative_Classification/figures/gene_heatmap_top30.png)

*图 6.3: Top 30 恶性 + Top 30 非恶性基因在 33,934 个上皮细胞中的表达热图。细胞按恶性评分排序，左侧为非恶性，右侧为恶性。恶性基因（红色条）在右侧高表达；非恶性基因（蓝色条）在左侧高表达。基因签名纯净度极高——无一免疫标记 (CD8A, CD19, CD163)，无一基质标记 (COL1A1, ACTA2)，全部为上皮/肠化生特异性基因。*

---

## 7. 评价：一个"自洽但无锚"的算法

<!-- Claude 起草底稿（2026-10-15），待你手写定稿——请按 diff 审改 -->

评价分三层：先承认方法与论文各自的价值（公正），再指出三重缺陷（核心），最后给出这件事作为"反面教材"的教学落点。

### 7.1 先承认两件事

**方法的机械价值**：k=2 强制二分 + 迭代 DEG，本质是一台"找数据中最大表达分歧轴"的自动化探针——实现简单、收敛快、无需先验。作为粗筛工具（"这张图谱里最大的转录分裂在哪"）它是合格的；我们六 Runs + 双种子全部稳定收敛（3-4 轮，agreement 91-96%）证明其机械稳定性。

**论文的数据价值**：Zhang 2021 是首张胃早癌/癌前 scRNA 图谱，其 chief→SPEM 轨迹发现真实且重要；置于 2019-21 语境，单细胞方法本身即新颖性，发表在 Gut 有历史合理性。批判不针对数据与发现，针对方法与结论的对应关系。

### 7.2 缺陷一：无生物学锚定——收敛到最大方差轴，而非恶性轴

算法的"恶性"与"非恶性"只是 k=2 在模块分数空间画的一条线。迭代 DEG 是正反馈：每一轮的切点都在强化上一轮的对比，最终收敛到数据方差最大的轴。在我们的复现中，这个轴是肠化生轴（IM axis），不是恶性轴：

- 六 Runs（上皮口径 / 基因选择 / 黑名单 / 低置信排除四种扰动）× 双种子（TCGA Toil 与 HiSeqV2 两套 bulk DEG），收敛主题**全部**是 IM 轴；
- 非恶性方向签名跨种子高度一致（Jaccard 0.82-1.00）；
- 测试种子第 1 轮签名即 4.7× 富集 IM 基因——种子自带 IM 偏倚，迭代不纠正反而放大。

算法没有任何环节引入"恶性是什么"的生物学知识；种子基因集只是初值。**收敛方向由数据方差结构决定，与疾病无关。**

### 7.3 缺陷二：切点无锚——输出是配置的函数，不是数据的属性

因为切点没有生物学含义，任何配置扰动都会平移切点：

| 扰动 | 恶性率 | 对照恶性率（Pat25-29） |
|------|:---:|:---:|
| 主种子 run1 | 50.5% | 30.6% |
| 仅换种子（HiSeqV2 早期误用） | 31.6% | **6.6%（4.6× 漂移）** |
| 仅加 12 基因管家黑名单 | 27.6% | **3.3%** |
| 仅剔除 14.1% 低置信注释细胞 | 57.6% | 38.4% |

对照组恶性率 3.3-41.6%：Pat25-29 是慢性胃炎患者（无癌），其细胞被判"恶性"的比例随配置任意漂移。低置信剔除实验中，测试种子下保留细胞 22.9% 跨越切点翻转标签（主种子仅 5.0%）——**切点附近的细胞成片翻转**，说明"恶性率"是种子与配置的函数，不是数据固有属性。

### 7.4 缺陷三（核心）：参照自污——迭代参照无纯度约束，且零保险

迭代 DEG 的"非恶性"参照并非固定：第 N 轮的参照 = 第 N−1 轮 k=2 画线剩下的池。算法从未使用临床标签，其"参照是正常细胞"的隐含前提没有任何机制保证。我们量化了主流程 run1 的两池构成：

| 池 | 对照（无癌） | A 癌旁 | B 肿瘤组织 |
|------|:---:|:---:|:---:|
| **非恶性池**（21,804） | 34.5% | 44.5% | **20.9%** |
| 恶性池（22,207） | **15.0%** | 44.9% | 40.2% |

- 非恶性参照池中 **65.5% 的细胞来自癌症患者，20.9% 直接来自肿瘤组织**；四 Run（双种子 × 低置信排除）下该构成几乎不变（20.9-22.4%）——污染是结构性的，不是运气；
- 反向同样失守：恶性池中 15.0% 来自无癌对照样本；
- 算法唯一的"保险"是 90% agreement 收敛阈值，而收敛检查的是**切点是否停止移动，不是参照是否纯净**——我们 91.0% 收敛时参照已 1/5 是肿瘤细胞；
- 论文自身提供了参照被污染的痕迹：其 no-malignant 群体呈 SPEM 演化方向——参照中混入了化生程序；
- 全流程无任何正交真值（inferCNV / 突变 / 独立标记体系）校验两池身份；论文唯一的"验证"是标记物表达——而标记物就是签名基因本身（循环论证）。

**参照冻结是本可选择的保险**：同一数据，我们的 MSDA 管线将参照模具从纯净对照一次性修剪后冻结、永不更新，迭代只发生在检测端。张 2021 选择了"让数据自己重定义对比"——迭代本身没错，**迭代而没有固定纯净支点**才是错。

### 7.5 语义与验证：malignant 标签的偷换

输出端还有一层独立的语义问题：算法输出的两群是"肠化生样 vs 胃分化样"（身份对比），论文将其表述为"malignant vs non-malignant"（疾病对比）。两个论据：① 收敛主题是 IM 轴；② 对照样本（无癌）30.6% 被判"恶性"，与"恶性=癌"的语义直接冲突。标签先行还导致 Fig 4A 的验证循环——用标签定义签名、再用签名解释标签。

### 7.6 教训（这件事作为反面教材的题）

1. **收敛 ≠ 正确**——agreement 是自指一致性（切点稳定），不是外部验证；
2. **参照必须锚定**——迭代可以，支点不能动；参照的正当性不能交给算法自己输出；
3. **语义与算法分离**——"最大的转录分歧" ≠ "恶性 vs 正常"，命名前先问轴是什么；
4. **scRNA 层无法自我救赎**——CNV score 测的是表达偏离（身份信号 100× 淹没 CNV 的 2×），真值锚必须来自 DNA 级或正交先验。

---