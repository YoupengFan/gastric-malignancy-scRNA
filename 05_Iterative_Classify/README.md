# 05_Iterative_Classify — 迭代 K-means 恶性/非恶性分类（函数2）

## 主流程结果（忠实口径 = Epithelial + Chief，与论文一致）

| 指标 | 值 |
|------|------|
| 上皮细胞数 | 44,011（Epithelial + Chief，排除 Parietal/Enteroendocrine = 论文口径） |
| 收敛 | **4 轮，agreement 91.0%**（阈值 ≥90%） |
| 轨迹 | 16,347 → 13,115 → 20,281 → 22,207 |
| 恶性 | 22,207 / 44,011 = **50.5%**（AdjA 43.5% / TumorB 66.2%） |
| 对照恶性率 | 30.6%（Pat25-29，10,853 细胞） |

## 六 Runs 敏感性矩阵（主种子专属）

Run 定义 × 结果：

| run | 上皮定义 | 基因选择 | HK 黑名单 | n | 恶性% | Ctrl% | AdjA% | TumorB% | 轮数 | agree |
|-----|------|------|:---:|-----:|-----:|-----:|-----:|-----:|:---:|-----:|
| run1 | Epi+Chief | logFC | ✗ | 44,011 | 50.5 | 30.6 | 43.5 | 66.2 | 4 | 0.9101 |
| run2 | 4 类型 | logFC | ✗ | 48,399 | 51.9 | 35.8 | 45.6 | 66.6 | 4 | 0.9349 |
| run3 | Epi+Chief | expression | ✗ | 44,011 | 52.1 | 35.6 | 46.5 | 64.8 | 4 | 0.9627 |
| run4 | Epi+Chief | logFC | ✓ | 44,011 | 27.6 | 3.3 | 20.4 | 43.9 | 3 | 0.9582 |
| run5 | 仅 Epi | logFC | ✗ | 24,946 | 57.2 | 41.6 | 54.5 | 61.4 | 3 | 0.9220 |
| run6 | Epi+Chief | logFC | ✗ | 37,798 | 57.6 | 38.4 | 51.1 | 71.2 | 4 | 0.9609 |

- run1 = 主流程复现（与主 05 输出逐位一致 = 复现锚）
- run2 = 广义口径敏感性（Python 版 46,554 的 R 对应物）
- run3 = 基因选择轴（expression 排序替代 logFC）
- run4 = 管家黑名单轴（额外排除 12 HK 基因）
- run5 = 上皮口径窄（仅 Epi）
- run6 = 注释置信度敏感性（剔除低置信注释细胞 6,213 = 14.1%）

- **Ctrl 恶性率非 0**（3.3–41.6%）：对照组为慢性胃炎（Pat25-29），
  推测含少量 IM 细胞（慢性胃炎→肠化生是 Correa 级联的常见伴随病变），
  导致方法按身份轴将ctrl划为"恶性"——————"恶性=肠化生混淆"

六 Runs 矩阵只对主种子做——检验数据集构成（上皮定义/基因选择/
低置信排除）对结果的影响：恶性率 27.6–57.6%、对照 3.3–41.6%，
收敛主题一致，结论对数据集条件稳健。

## 测试种子对照（run1 + run6，testseed = HiSeqV2 早期误用种子）

| run | 恶性% | 轮数 | label_agree | J(mal top10) | J(non top10) | AdjA% | TumorB% |
|-----|-----:|:---:|-----:|-----:|-----:|-----:|-----:|
| test_run1 | 31.6 | 3 | 0.791 | 0.429 | 0.818 | 23.4 | 50.0 |
| test_run6 | 57.5 | 4 | 0.983 | 1.000 | 1.000 | 50.6 | 72.0 |

测试种子对照恶性率：test_run1 6.6 / test_run6 37.7（Pat25-29）。

测试种子只保留与论文数据定义符合的 run1（Epi+Chief, logFC）与
run6（低置信注释排除）两种配置。其余 run 与论文工作条件差异大
（数据源不同 × 细胞类型条件不同），不具说明力——测试种子的角色
是检验 **bulk 种子集对结果的引导作用**，数据集构成敏感性由主种子
六 Runs 矩阵承担。

**双种子两个论点**：

1. **切点无生物锚定**：run1 对照恶性率 30.6% → 6.6% ≈ 4.6× 漂移
   （种子换源即整体平移）——"恶性率"由种子锚定，非数据固有属性
2. **跨数据源收敛稳健**：收敛主题同为 IM 轴、最终签名 J(non) 0.82、
   label_agree 0.79——方法缺陷（强制二分收敛于最大方差轴）跨源复现

## 对照 per-patient 拆列（main/run1 vs test_run1）

| patient | n | main/run1 恶性% | test_run1 恶性% |
|------|-----:|-----:|-----:|
| Pat25 | 1,749 | 48.5 | 15.0 |
| Pat26 | 4,889 | 36.1 | 6.6 |
| Pat27 | 545 | 68.3 | 0.7 |
| Pat28 | 2,341 | 10.0 | 5.2 |
| Pat29 | 1,329 | 7.8 | 0.6 |

## paired 患者 A/B 拆侧（main/run1 口径）

| patient | n_A | n_B | A% | B% |
|------|-----:|-----:|-----:|-----:|
| Pat01 | 1,110 | 562 | 97.8 | 94.7 |
| Pat02 | 1,161 | 1,309 | 66.7 | 58.4 |
| Pat03 | 751 | 812 | 83.2 | 64.9 |
| Pat04 | 626 | 428 | 36.4 | 56.3 |
| Pat05 | 1,345 | 453 | 23.3 | 53.2 |
| Pat06 | 1,832 | 672 | 29.0 | 40.8 |
| Pat07 | 603 | 275 | 41.3 | 85.8 |
| Pat08 | 480 | 279 | 20.8 | 58.8 |
| Pat09 | 1,096 | 414 | 24.2 | 76.1 |
| Pat10 | 806 | 446 | 92.3 | 89.7 |
| Pat11 | 493 | 285 | 47.3 | 68.1 |
| Pat12 | 650 | 529 | 22.8 | 83.2 |
| Pat13 | 1,023 | 648 | 35.5 | 95.4 |
| Pat14 | 240 | 44 | 39.2 | 100.0 |
| Pat15 | 990 | 674 | 20.0 | 37.4 |
| Pat16 | 875 | 1,157 | 96.6 | 93.9 |
| Pat17 | 1,451 | 897 | 99.5 | 97.7 |
| Pat18 | 976 | 699 | 25.4 | 50.8 |
| Pat19 | 512 | 864 | 64.3 | 59.7 |
| Pat20 | 228 | 64 | 37.3 | 85.9 |
| Pat22 | 790 | 744 | 49.4 | 68.3 |
| Pat23 | 445 | 98 | 56.4 | 63.3 |
| Pat24 | 760 | 1,127 | 9.6 | 19.3 |

Pat21（unpaired，仅 A 侧）435 细胞 80.2%。

## 重跑方式

**输入依赖**：`03_Cluster_Annotate/data/gse150290_epi_counts.rds`——上皮谱系
（Epithelial/Chief/Enteroendocrine/Parietal）极简子集（RAW counts + meta），
随仓库分发；本地全链复现时由 03 的 `export_epi_counts.R` 从全量注释对象导出。

```bash
export PATH="$HOME/softwares/bin:$PATH"
# 均在仓库根目录执行

# ① 主流程（= run1 配置，主种子）
#    产出：data/ 根目录五件套（labels / scores / iteration_log /
#    final_signatures / iterative_params.json）+ figures/
Rscript 05_Iterative_Classify/scripts/step2_iterative_classify.R

# ② 六 Runs 敏感性扫描（主种子专属）
#    用法：run_iterative_runs.R [run{1..6} | all] [种子目录]
#    缺省 all = 六 Runs 全扫。星形设计：run1 复现锚，run2-6 各变一个条件
#    （run2 上皮口径宽 / run3 基因选择 expression / run4 管家黑名单 /
#    run5 上皮口径窄 / run6 低置信排除）
#    产出：data/run{1..6}_*/ 四件套（labels / iteration_log /
#    final_signatures / run_params.json）；run1 与主流程逐位一致（复现锚）
Rscript 05_Iterative_Classify/scripts/sensitivity/run_iterative_runs.R all

# ③ 测试种子对照（HiSeqV2 早期误用种子，只跑论文定义符合的 run1/run6）
#    第二参数指定种子目录；产出加 testseed_ 前缀
Rscript 05_Iterative_Classify/scripts/sensitivity/run_iterative_runs.R run1 04_Bulk_DEG/data/testseed
Rscript 05_Iterative_Classify/scripts/sensitivity/run_iterative_runs.R run6 04_Bulk_DEG/data/testseed

# ④ 复现锚检验：逐位比对主流程 labels.csv 与 run1 产出（一致 = PASS）
Rscript 05_Iterative_Classify/scripts/sensitivity/check_run1_anchor.R

# ⑤ 结果汇总：9 份 labels.csv 回连注释 → data/run_summary.csv（每 Run 一行：
#    细胞数 / 恶性% / A 侧% / B 侧% / 对照 Pat25-29%）+ data/per_patient_malignant.csv
Rscript 05_Iterative_Classify/scripts/sensitivity/summarize_runs.R

# ⑥ 双种子对比：主 vs 测试种子 run1/run6 → data/seed_comparison.csv
#    （label_agree 逐细胞一致率 + 签名 top10 Jaccard）
Rscript 05_Iterative_Classify/scripts/sensitivity/compare_seeds.R
```

## 数据文件说明

### data/ 根目录（主流程五件套 + 三份汇总）

| 文件 | 来源 | 内容 |
|------|------|------|
| labels.csv | 主流程 | 逐细胞标签（malignant / non_malignant） |
| scores.csv | 主流程 | 逐细胞恶性/非恶性模块评分（K-means k=2 输入） |
| iteration_log.csv | 主流程 | 迭代轨迹：每轮恶性数 + 标签一致率 |
| final_signatures.csv | 主流程 | 最终 50 恶性 / 50 非恶性签名基因 |
| iterative_params.json | 主流程 | 本轮参数记录（种子/收敛阈值/上皮口径等） |
| run_summary.csv | summarize_runs.R | 9 Runs 每 Run 一行（n / 恶性% / A 侧% / B 侧% / 对照%） |
| per_patient_malignant.csv | summarize_runs.R | Run × 患者 × 侧别恶性长表（468 行，paired 患者 A/B 拆两行；对照患者仅有 A 行） |
| seed_comparison.csv | compare_seeds.R | 双种子 run1/run6 配对对比（label_agree + 签名 top10 Jaccard + A/B 侧恶性率回连注释重算） |

### data/ 各 Run 子目录（run1-6 + testseed_run1/6，每目录四件套）

| 文件 | 内容 |
|------|------|
| labels.csv | 该 Run 逐细胞标签 |
| iteration_log.csv | 该 Run 迭代轨迹 |
| final_signatures.csv | 该 Run 最终签名 |
| run_params.json | 该 Run 参数记录 |

### figures/（主流程图件）

| 文件 | 内容 |
|------|------|
| fig_scores_kmeans.png | 恶性/非恶性模块评分散点，K-means 标签着色 |
| fig_convergence.png | 收敛轨迹（agreement vs 轮次，红线 = 90% 阈值） |
