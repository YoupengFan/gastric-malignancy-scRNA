# 嫌疑污染细胞：标记标准、构成、机制与患者分布

> 数据源：03_Cluster_Annotate/data/annotation_qc_umap.csv（150,490 cells）
> 推导脚本：annotation_qc_export.R（flag 定义）、annotation_qc_figures.R（图）
> 分析脚本：03_Cluster_Annotate/scripts/contamination_analysis/contam_pairs_enumerate.R、contam_patient_dist.R

## 1. 标记标准：只标跨谱系，同谱系豁免

定义（annotation_qc_export.R L319-334）：

```
suspected_contam := lineage(vote) != lineage(cluster_mode)
```

- **谱系三分**：epi（Epithelial/Chief/Parietal/Enteroendocrine）/ immune
  （T/B/NK/Myeloid/Mast/Plasma）/ stroma（Fibroblasts/Endothelial）
- **豁免**：同谱系混杂（Chief↔Epithelial 等）是分化连续体
  （differentiation continua），互相渗透是正常生物学，不是污染
- **定性条件、无表达量门槛、嫌疑非裁决**（suspect = flag, never a verdict）
- **理由**：跨谱系转录程序互斥。投票与落点簇分属不同谱系时，唯一自洽
  解释 = 外部 RNA 引入（解离破裂泄漏 / doublet）或投影误差；两者下游
  DEG 都应视为嫌疑

验证：同谱系豁免 0%（Chief↔Epithelial 等不标记）；
low_conf 富集 4.8×（42.8% vs 8.9%）；PGC+MZB1 双阳富集 1.9×
（43.3% vs 23.2%）。

## 2. 总体构成：7,549 cells（5.02%）

### 2.1 谱系方向

| 方向 | N | 占比 |
|------|:---:|:---:|
| epi → immune | 5,555 | 73.6% |
| immune → epi | 549 | 7.3% |
| stroma → immune | 580 | 7.7% |
| epi → stroma | 546 | 7.2% |
| immune → stroma | 294 | 3.9% |
| stroma → epi | 25 | 0.3% |

### 2.2 Top 混杂配对

| 配对 | N | 泄漏源 |
|------|:---:|------|
| Chief → Plasma | 2,468 | PGC |
| Chief → B | 998 | PGC |
| Chief → T | 733 | PGC |
| Chief → Endothelial | 243 | PGC |
| Chief → Fibroblasts | 220 | PGC |
| Chief → Myeloid | 217 | PGC |
| Endothelial → Myeloid | 203 | VWF 系 |
| Fibroblasts → Myeloid | 203 | COL1A1 系 |
| Plasma → Epithelial | 166 | MZB1 |
| Enteroendocrine → B | 136 | CHGA |

## 3. 为什么 Chief 污染最高 — 四大泄漏源机制

| 泄漏源 | 嫌疑数 | 占比 | 机制 |
|------|:---:|:---:|------|
| Chief / PGC 系 | 4,882 | 65% | PGC 中位 ~770 counts（纯 Chief log1p 6.65），比典型
marker 高 1-2 个数量级；主细胞是消化酶工厂，解离时最易破裂 →
超比例局部 ambient，对照减法穿透，全方向泄漏（immune 簇为主） |
| Plasma / MZB1 | 330 | 4% | 浆细胞是抗体工厂、解离中最脆弱；落点几乎全在
Epithelial 簇（166/330） |
| B 系 | 104 | 1% | CD79A/B、MS4A1；落点 Epithelial 37/104 |
| Enteroendocrine / CHGA | 375 | 5% | CHGA 是分泌性颗粒蛋白，内分泌细胞天然分泌 —
污染不纯是解离损伤（分泌生理 vs 破裂泄漏并存） |

**Chief 污染高不是 Chief 特殊，而是 PGC 表达量特殊** —
同样破碎率下，表达量高 1-2 个数量级的基因泄漏拷贝数最多，最先穿透
空液滴对照减法。这是 scRNA-seq 局部 ambient 的剂量效应，不是注释
方法缺陷。

## 4. 患者分布：人人都有，头部集中，且头部在对照患者

### 4.1 集中度统计

- 29/29 名患者全部有嫌疑细胞（不存在"仅少数患者被标记"）
- 比率范围 1.00%（Pat20）– 35.90%（Pat29），中位 2.72%，均值 6.05%
- top-3 占 41.1%，top-5 占 55.6%，top-10 占 74.2%
- **Gini = 0.520**（0 = 完全均匀，1 = 一人独占）→ 中等偏集中：
  不是均匀铺开，也非一两人包揽
- corr(比率, 患者细胞数) = **−0.347** — 排除"细胞多所以嫌疑多"的
  伪影，反而小样本患者比率更高

### 4.2 头部患者

| 患者 | 分组 | 嫌疑数 | 比率 | epi→immune 占比 |
|------|------|:---:|:---:|:---:|
| Pat28 | 慢性胃炎对照 | 1,249 | 33.6% | 93.7% |
| Pat24 | 癌配对 | 1,078 | 20.5% | 90.7% |
| Pat29 | 慢性胃炎对照 | 779 | **35.9%** | 95.0% |
| Pat13 | 癌配对 | 705 | 12.3% | 68.9% |
| Pat08 | 癌配对 | 384 | 5.7% | 57.8% |

### 4.3 解读：样本层面因素，非肿瘤生物学

1. **比率最高的恰是两名慢性胃炎对照（Pat28/29）**，且其嫌疑细胞
   90%+ 是同一方向 epi→immune（Chief→Plasma 类 PGC 泄漏）。对照
   患者不存在恶性细胞，高比率更可能来自**样本层面的解离质量 /
   组织状态**：胃炎状态下主细胞更脆弱，或这两个样本组织处理质量
   较差（候选解释，不裁决）。佐证：Pat28 在同型 doublet 分析中也
   反复出现（解离破碎多的样本）。
2. **两种模式并存**：Pat28/29/24 是单方向碾压（91-95% epi→immune =
   PGC 泄漏主导）；Pat13/08 方向较分散（epi→immune 58-69%，
   Pat13 另有 22.8% epi→stroma）— 头部集中不是单一机制。

## 5. 去除嫌疑后的注释质量（干净热图证据）

qc_confusion_heatmap_clean.png（排除 7,549 嫌疑后）：
对角线纯度 Fibroblasts 1.00 / Mast 0.99 / Endothelial 0.99 /
Plasma 0.96 / T 0.95 / Parietal 0.93 / B 0.92 / Epithelial 0.91 /
Myeloid 0.84 / Chief 0.82 / Enteroendocrine 0.59。
残留非对角线全部为同谱系分化连续体（豁免设计使然：Chief 簇内
17.2% Epithelial vote、EE 簇内 35.9% Epithelial vote）。
→ **跨谱系混杂是注释内部不一致的主要来源**；去除后各簇自洽性
大幅提升，同时不牺牲任何细胞类型（嫌疑只 5.02%，且只是标记）。

## 6. 立场

- **标记而非删除**：suspected_contam 是 per-cell 标记列，随全量
  数据一起发布；任何下游可从全量 + 标记重建去除版（一条 filter）
- **嫌疑非裁决**：7,549 细胞未被从任何数据中删除，主分析默认
  全量；去除版仅作敏感性检验
- **透明报告**：两份成分热图（全量红版 / 去除蓝版）+ 诊断 UMAP
  （qc_suspected_contam.png + qc_umap_clean.png）随图发布

## 7. 相关产出

| 文件 | 说明 |
|------|------|
| data/annotation_qc_umap.csv | per-cell 导出，含 suspected_contam 标记列 |
| data/annotation_qc_confusion.csv | 全量 vote×cluster 列联表（红版热图数据） |
| figures/qc_suspected_contam.png | 嫌疑细胞诊断 UMAP（补充材料） |
| figures/qc_umap_clean.png | 去除嫌疑版注释总览 UMAP |
| figures/qc_confusion_heatmap.png | 全量成分热图（红系 #8b0000） |
| figures/qc_confusion_heatmap_clean.png | 去除嫌疑版成分热图（蓝系 #0b3d91） |
| scripts/annotation_qc_export.R | flag 定义（L319-334）+ 导出 |
| scripts/annotation_qc_figures.R | 诊断图与干净热图（疑似污染 UMAP、clean UMAP、干净混淆热图）|
