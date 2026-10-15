# 02_Normalize_Batch — 归一化 + 批次矫正

GSE150290 QC 矩阵（159,308 保留细胞，已剔除 1,051 Scrublet 异型 doublet）→ log1p(CP10K) 归一化 + CC/MT 回归 + Harmony（batch key = sample，52 个 10x run）→ 归一化 Seurat 对象 + 矫正前后验证图与定量指标。

## 工作流分工：服务器计算 → 轻量数据 → 本地渲染

与 01_QC 同款分工：**重计算留在 Linux/WSL（服务器语义），只有轻量结果数据集被带回本地渲染**。各脚本自包含、跨环境只通过 `data/` 里的 CSV/JSON 交接，不共享内存态。

| 脚本 | 运行环境 | 作用 | 输入 | 产出 | 耗时 |
|------|----------|------|------|------|:---:|
| `normalize_batch.R` | WSL/Linux（服务器） | log1p(CP10K) + CC/MT 回归 + Harmony 全管线 | 01_QC 的 `gse150290_qc_matrix.rds` + `gse150290_qc_meta.csv` + `gse150290_doublet_scores.csv` | `gse150290_norm_harmony.rds`（6.7 GB）、`normalize_batch_params.json`、`normalize_batch.log` | ~14 min |
| `umap_batch_compute.R` | WSL/Linux（服务器） | pre-Harmony UMAP（uwot，seed 42）+ 导出 post UMAP 与 meta | 03 的 `gse150290_annotated.rds` | `umap_batch.csv`（150,490 cells） | ~15 min |
| `umap_batch_figure.R` | Windows/Linux 本地 | 批次矫正 before/after demo 图（上皮谱系 × A 侧 × Pat01/17/21，原数据直出） | `umap_batch.csv` | `umap_batch_before_after.png` | <1 min |
| `validate_pcs_compute.R` | WSL/Linux（服务器） | per-PC 批次信号（η²）+ FNN k=20 混合率 | 03 的 `gse150290_annotated.rds` | `pc_validation.csv`、`pc_validation_summary.json`、`fnn_mixing_percell.csv` | ~40-90 min |
| `validate_pcs_figure.R` | Windows/Linux 本地 | PC 验证图（A batch-mixing 背靠背直方图 + B per-PC η² 柱） | 上表 3 个 CSV/JSON | `pc_validation.png` | <1 min |
| `viz_library_scaling.R` | Windows/Linux 本地 | 逐样本文库大小 + 细胞数总览（与 02 归一化同批细胞） | 01_QC 的 `gse150290_qc_meta.csv` + `gse150290_doublet_scores.csv` | `library_and_cells_overview.png`、`library_summary.json` | <1 min |

- **执行依赖**：`normalize_batch.R` 在 03 的 `cluster_annotate.R` 之前；`umap_batch_compute.R` 与 `validate_pcs_compute.R` 读 03 的 annotated RDS，必须在 03 之后运行。
- **计算端**（Linux）：R ≥ 4.5（Seurat / harmony / FNN / uwot），脚本路径均为仓库根相对路径。
- **渲染端**（本地）：只需 R + ragg / ggplot2 / patchwork，只读轻量 CSV/JSON，可迁移任意操作系统。

## 关键数字（论文图注 / 附表引用）

**批次矫正效果**（pc_validation，10K 口径）：

| 指标 | 矫正前 | 矫正后 |
|------|:---:|:---:|
| mean η²(sample) PC1-30 | 0.1014 | 0.0345 |
| same-sample ratio（FNN k=20） | 0.283 | 0.098（random 1/52 ≈ 0.019） |

## 与论文对照

- Zhang 2021 在 scale 步骤回归 total UMI + tSNE；我们用 CC/MT 回归 + Harmony。
- Harmony batch key = sample（52 个 10x run）；矫正后 η²(sample) 从 0.1014 降至 0.0345。
- Scrublet 异型 doublet 剔除（1,051 / 0.66%）——异型 doublet 是偏离类分析混淆项；刻意不做大细胞/高 UMI 硬过滤（会误杀 Parietal/Chief/恶性大细胞）。

## 复现

01_QC 之后：`normalize_batch.R` → 03 的 `cluster_annotate.R` → `umap_batch_compute.R` + `validate_pcs_compute.R` → 本地跑 3 个渲染脚本。原始数据不随仓库分发（GSE150290 公开下载）。