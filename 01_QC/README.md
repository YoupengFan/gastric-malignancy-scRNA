# 01_QC — 质量控制

GSE150290（52 个 10x v2 样本，38.4M raw barcodes）→ QC 过滤合并矩阵 + doublet 评分 + 四面板 QC 总览图。

**QC 阈值**：UMI ∈ [1,001, 15,000]（下界逐字沿用 Zhang 2021，上界从论文的 6,000 放宽——见批判记录）、genes ≥ 501、MT ≤ 50%。看家基因过滤省略：论文 curated list 未随文公开。

## 工作流分工：服务器计算 → 轻量数据 → 本地渲染

按真实服务器 + 个人电脑工作流设计：**重计算留在 Linux/WSL（服务器语义），只有轻量结果数据集被带回本地渲染**。各脚本自包含、跨环境只通过 `data/` 里的 CSV/JSON 交接，不共享内存态。

| 脚本 | 运行环境 | 作用 | 输入 | 产出 | 耗时 |
|------|----------|------|------|------|:---:|
| `qc_filter.R` | WSL/Linux（服务器） | 读 52 个 raw 10x 矩阵 → QC 过滤 → 合并 | `01_Data/raw/`（原始数据，仓库外） | `gse150290_qc_matrix.rds`（合并 dgCMatrix）、`gse150290_qc_meta.csv`、`qc_stats.csv`、`qc_params.json` | ~15 min |
| `qc_doublets.py` | WSL/Linux（服务器） | 逐样本 Scrublet doublet 评分（checkpoint 断点续跑） | raw 矩阵 + `gse150290_qc_meta.csv` | `gse150290_doublet_scores.csv`、`qc_doublets_summary.json` | ~15 min |
| `qc_plot_compute.R` | WSL/Linux（服务器） | 过滤前指标聚合（38M barcode 单遍聚合 + knee 图抽稀到 ~200k 点） | raw 矩阵 | `qc_plot_hist_total.csv`、`qc_plot_hist_genes.csv`、`qc_plot_knee.csv`、`qc_plot_scatter.csv`、`qc_plot_summary.json` | ~12 min |
| `qc_plot_figure.R` | Windows/Linux 本地 | 四面板 QC 总览图（ragg 渲染，仅消费轻量 CSV/JSON） | 上表 5 个 CSV/JSON | `qc_metrics_overview.png` | <1 min |

- **计算端**（Linux）：需要 R ≥ 4.5（Matrix/data.table/jsonlite）、Python 3 + Scrublet。脚本路径均为仓库根相对路径。
- **渲染端**（本地）：只需 R + ragg/ggplot2/patchwork，路径 `E:/...` 按本机修改。渲染端可以自由迁移到任何操作系统——它只读轻量数据。

## 阈值注释（批判点留档）

- 下界逐字沿用 Zhang 2021 Gut：UMI ≥ 1,001、genes ≥ 501、MT ≤ 50%。
- **上界放宽 6,000 → 15,000**：数据显示 6,000 天花板截断大细胞（无自然悬崖，6k 处系统性误删）。15,000 = P99.97 数据驱动选择。详见 `qc_params.json` 与四面板图中的对照虚线。
- **看家基因过滤省略**：论文提到"curated list"但列表未随文公开（supp001 p5），无法忠实复现。

## 复现

按上表顺序执行三个计算脚本（任一 Linux/WSL 环境）后，在本地执行渲染脚本。原始数据不随仓库分发（GSE150290 公开下载，52 个 GSM4546300–GSM4546351）。
