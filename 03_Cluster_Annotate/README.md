# 03_Cluster_Annotate — Leiden 聚类 + Marker注释

02 的 Harmony 归一化对象（159,308 细胞）→ 红细胞剔除 → harmony 50PC 上 Leiden 聚类 + UMAP → 12 类 marker 投票注释（逐细胞 z-score）→ 注释 QC 表与图件。

## 工作流分工：服务器计算 → 轻量数据 → 本地渲染

与 01/02 同款分工：**重计算留在 Linux/WSL（服务器语义），只有轻量结果数据集被带回本地渲染**。各脚本自包含、跨环境只通过 `data/` 里的 CSV/JSON 交接。

| 脚本 | 运行环境 | 作用 | 输入 | 产出 | 耗时 |
|------|----------|------|------|------|:---:|
| `cluster_annotate.R` | WSL/Linux（服务器） | 红细胞剔除 + Leiden 聚类 + UMAP + Marker 注释 | 02 的 `gse150290_norm_harmony.rds`（6.7 GB） | `gse150290_annotated.rds`、`cell_annotation.csv`、`cluster_annotation.csv`、`cluster_annotate_params.json`、`cluster_annotate.log` | ~8 min |
| `annotation_qc_export.R` | WSL/Linux（服务器） | 从 annotated RDS 算 6 张注释 QC 轻量表 | 本目录 `gse150290_annotated.rds` | `annotation_qc_*.csv` ×6 | ~2.5 min |
| `export_epi_counts.R` | WSL/Linux（服务器） | 导出 05 的上皮谱系极简输入（RAW counts + meta，xz 压缩） | 本目录 `gse150290_annotated.rds` | `gse150290_epi_counts.rds` | ~3 min |
| `annotation_qc_figures.R` | Windows/Linux 本地 | 零重算渲染 9 张注释 QC 图 | 上表 6 张 CSV | `figures/` 下 9 张 PNG | <1 min |

- 补充诊断（非核心流程）：`scripts/contamination_analysis/`——嫌疑污染细胞的患者分布与配对枚举（零重算控制台统计）+ 污染分析 NOTES。
- 执行依赖：`cluster_annotate.R` 在 02 的 `normalize_batch.R` 之后；`annotation_qc_export.R` 读本目录 annotated RDS；`export_epi_counts.R` 在 `cluster_annotate.R` 之后（05 的输入准备，产物随仓库分发）；`annotation_qc_figures.R` 在 export 之后。
- 计算端（Linux）：R ≥ 4.5（Seurat / Matrix / data.table / jsonlite），脚本路径均为仓库根相对路径。
- 渲染端（本地）：只需 R + ragg / ggplot2 / data.table，只读轻量 CSV/JSON，可迁移任意操作系统。

## Marker 集：来源与 hg19 适配（12 类）

**来源细分**：

| 面板 | 来源 |
|------|------|
| 免疫面板（T/B/NK/Myeloid/Mast/Plasma） | Zhang 2021——Kim 2022 免疫细胞仅列 6 类名称（Supp Fig. 2），未公开 marker 基因，免疫注释必须依赖 Zhang |
| 基质面板 | Zhang 2021 基础 + Kim 2022 补充基因（Fibroblasts += MMP2/PDGFRA/MYL9；Endothelial += PLVAP/KDR/PTPRB） |
| 胃上皮亚型 | Kim 2022（Enteroendocrine = CHGA/GAST/PROX1；Chief = PGC）+ **本数据 DEG 发现的 Parietal**（GIF/ATP4B/ATP4A，无文献 marker 表来源：Zhang 9 类注释下 Cls 15 被投为 Fibroblasts 但纯度仅 6.3%，其 per-cluster DEG top 基因恰为壁细胞特征基因 → 重标注并确立面板；Kim 论文未设壁细胞类型）+ 通用胃黏膜 marker（Epithelial 面板、LIPF/PGA3/PGA4、GHRL/SST/PCSK1N） |

**针对本 hg19 参考与数据的适配**（逐样本验证）：

| 面板 | 适配 | 原因 |
|------|------|------|
| T_cells | +CD7/LCK/IL32/CCL5 | 本 hg19 参考缺失 TRAC/TRBC1 → 面板原本只剩 3 基因；重投票实验 VA：NK 误投 2662→1411、CD3D+ 污染 38.8%→27.8% |
| Plasma_cells | 仅 MZB1 + SDC1 | IGHG1/IGKC 不在本 hg19 参考注释（免疫球蛋白位点）；JCHAIN 以别名 IGJ 存在但 85% 全局检出（ambient RNA）；TNFRSF17/DERL3/FKBP11 在 26-39% 全部细胞检出（太泛） |
| Fibroblasts | −FAP | 本数据死基因（det 4%、mean 0.024） |
| Endothelial | −PECAM1 | hg19 参考缺失；7/8 marker 存活 |

## 注释方法要点

- 打分：AddModuleScore 协议（25 bins × 100 同 bin 对照扣除，稀疏安全；对照扣除去掉 ambient 基线——PGC ambient 在 57% 全部细胞检出）
- 逐细胞 z-score（12 类 marker score 均值）→ argmax 投票 = 最终注释
- 低置信 = z-delta（最优 − 次优）< 0.5
- CD3D 谱系救援：NK 面板投票 + CD3D 检出 → T_cells（T 面板失去 TRAC/TRBC1 后 CD3D+ T 被拉向 NK 面板）
- 簇 mode 仅参照列（cluster_annotation + purity），不覆盖逐细胞投票
- 剔除：红细胞（globin 计数占比 > 5%，聚类前完成，簇号对最终细胞集自洽）+ 死细胞簇（簇中位 MT% > 25 且 max 中位 marker z < 1.5）

## 关键数字（论文图注 / 附表引用）

| 指标 | 值 |
|------|------|
| 输入细胞 | 159,308（02 输出，已剔除 1,051 Scrublet 异型 doublet） |
| 红细胞剔除 | 583（0.37%） |
| 最终细胞 | 150,490 / 20 clusters（Leiden r=0.5）/ 12 类型 |
| 低置信（z-delta < 0.5） | 15,940（10.6%） |
| CD3D 救援 | 691/2,120 NK 投票（32.6%）→ T_cells |
| 死细胞簇 | Cls 7 Epithelial（8,235，5.19%；Parietal 高 MT 受保护豁免） |

逐类型计数（逐细胞投票，丰度序）：Plasma_cells 39,326 / Epithelial 24,946 / T_cells 24,721 / Chief 19,065 / B_cells 17,925 / Fibroblasts 7,999 / Myeloid 3,848 / Endothelial 3,706 / Mast_cells 3,152 / Enteroendocrine 2,933 / Parietal 1,455 / NK_cells 1,414。

低置信率分类型跨度大：稀有类型高（NK 29.5%、Parietal 27.6%、Enteroendocrine 22.8%、Chief 19.8%），丰度类型低（Fibroblasts 5.4%、Plasma 6.7%）——投票边际（z-delta）低的天然落在面板薄或面板间竞争强的类型。

## 与论文对照

- Zhang 2021 手动 marker 注释 9 类；我们用 Hybrid （Kim + Zhang 12 类、逐细胞 z-score、低置信量化）。免疫注释必须依赖 Zhang 2021——Kim 2022 未公开免疫 marker。
- 注释 QC 图参照 Zhang 2021 Fig D/E 布局（代表 marker violin = Fig D 风格；UMAP feature 图 = Fig E 风格）。
- Plasma 面板仅两基因（MZB1+SDC1）→ Plasma 簇内 2,798 疑似标记（Chief 2,468）：分泌上皮 ↔ 浆细胞互吸。

## 复现

02 之后：`cluster_annotate.R` → `annotation_qc_export.R` → `export_epi_counts.R` → 本地跑 `annotation_qc_figures.R`。原始数据不随仓库分发（GSE150290 公开下载）。
