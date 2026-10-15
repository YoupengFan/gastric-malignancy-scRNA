#!/usr/bin/env Rscript
# ============================================================
# 迭代聚类结果汇总：9 份 labels.csv（主 05 + 主种子 6 Run
# + 测试种子 run1/run6）回连注释，统一拆出
#   总 / A、B 侧 / 对照（Pat25-29）/ 逐患者 的恶性率
#
# 输入：05_Iterative_Classify/data/ 各目录 labels.csv
#       03_Cluster_Annotate/data/cell_annotation.csv（回连 patient/tissue）
# 产出：05_Iterative_Classify/data/
#   run_summary.csv          每 Run 一行：n / 恶性% / A 侧% / B 侧% / 对照%
#   per_patient_malignant.csv Run × 患者 × 侧别 长表：n / 恶性数 / 恶性%
# ============================================================
suppressPackageStartupMessages(library(data.table))

REPO <- "."
DATA <- file.path(REPO, "05_Iterative_Classify/data")
ANNO <- file.path(REPO, "03_Cluster_Annotate/data/cell_annotation.csv")

anno <- fread(ANNO)[, .(barcode, patient, tissue)]

label_files <- c(
  main      = file.path(DATA, "labels.csv"),
  run1      = file.path(DATA, "run1_epi_chief_logFC", "labels.csv"),
  run2      = file.path(DATA, "run2_epi_4type_logFC", "labels.csv"),
  run3      = file.path(DATA, "run3_epi_chief_expr", "labels.csv"),
  run4      = file.path(DATA, "run4_epi_chief_logFC_blacklist", "labels.csv"),
  run5      = file.path(DATA, "run5_epi_only_logFC", "labels.csv"),
  run6      = file.path(DATA, "run6_epi_chief_logFC_lowconf_excl", "labels.csv"),
  test_run1 = file.path(DATA, "testseed_run1_epi_chief_logFC", "labels.csv"),
  test_run6 = file.path(DATA, "testseed_run6_epi_chief_logFC_lowconf_excl", "labels.csv")
)

CTRL_PAT <- paste0("Pat", 25:29)

# labels 回连注释（cell → patient/tissue）；漏连则警告
read_run <- function(path) {
  lab <- fread(path)
  m <- merge(lab, anno, by.x = "cell", by.y = "barcode")
  if (nrow(m) != nrow(lab))
    warning(sprintf("%s: %d/%d cells not found in annotation",
                    basename(path), nrow(lab) - nrow(m), nrow(lab)))
  m
}

# 子集恶性率（空集 → NA，避免 NaN 混入表格）
pct_of <- function(labels) {
  if (length(labels) == 0) NA_real_
  else round(mean(labels == "malignant") * 100, 1)
}

# ── 每 Run 一行汇总 ────────────────────────────────────────
run_summary <- rbindlist(lapply(names(label_files), function(rn) {
  m <- read_run(label_files[[rn]])
  data.table(
    run = rn,
    n_cells = nrow(m), n_malignant = sum(m$label == "malignant"),
    pct_malignant = pct_of(m$label),
    n_A = sum(m$tissue %in% "A"), pct_A = pct_of(m$label[m$tissue %in% "A"]),
    n_B = sum(m$tissue %in% "B"), pct_B = pct_of(m$label[m$tissue %in% "B"]),
    n_ctrl = sum(m$patient %in% CTRL_PAT),
    pct_ctrl = pct_of(m$label[m$patient %in% CTRL_PAT])
  )
}))
fwrite(run_summary, file.path(DATA, "run_summary.csv"))

# ── Run × 患者 × 侧别长表 ──────────────────────────────────
# A/B 拆开统计：paired 患者的 A（癌旁）与 B（肿瘤）是不同组织，
# 合并口径下的"恶性率"无法解读；对照 Pat25-29 与 Pat21 只有 A 侧，
# 不受影响（注：注释的 tissue 字段把对照标为 "A"，
# 此处 A 行对对照患者是样本侧而非癌旁）
per_patient <- rbindlist(lapply(names(label_files), function(rn) {
  m <- read_run(label_files[[rn]])
  m[, .(n_cells = .N, n_malignant = sum(label == "malignant")),
    by = .(patient, tissue)
    ][, `:=`(pct_malignant = round(n_malignant / n_cells * 100, 1),
             run = rn)]
}))
setcolorder(per_patient, c("run", "patient", "tissue", "n_cells", "n_malignant", "pct_malignant"))
per_patient <- per_patient[order(run, patient, tissue)]
fwrite(per_patient, file.path(DATA, "per_patient_malignant.csv"))

cat("\n=== run_summary.csv ===\n")
print(run_summary)
cat(sprintf("\nWrote: run_summary.csv (%d runs), per_patient_malignant.csv (%d rows)\n",
            nrow(run_summary), nrow(per_patient)))
