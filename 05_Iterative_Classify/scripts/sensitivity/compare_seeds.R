#!/usr/bin/env Rscript
# ============================================================
# 种子敏感性对比：测试种子（HiSeqV2+voom——早期误用数据源的种子）
# vs 主种子（Toil+limma+log2(TPM+1)）
#
# 测试种子只对比 run1 与 run6 两种配置：
#   run1 = 论文定义符合（Epithelial+Chief, logFC）
#   run6 = run1 + 剔除低置信注释细胞（低质量数据排除）
# 其余 run 与论文工作条件差异大（不同数据源 × 不同细胞类型
# 条件），不具说明力
#
# 输入（主种子与测试种子各自在 run1/run6 两种配置下的运行产物——共
# 4 个目录；跑 run 的是 run_iterative_runs.R（第二参数传种子目录），
# 本脚本只读产物做对比。种子基因集本身产自 04_Bulk_DEG/data/：
# 主种子 = top50_{tumor,normal}_limma.txt、测试种子 = testseed/
# 同名文件）：
#   05_Iterative_Classify/data/run{1,6}_*/（主种子 × run1/run6 产物）
#   05_Iterative_Classify/data/testseed_run{1,6}_*/（测试种子 × run1/run6 产物，同一数据版本）
# 产出：05_Iterative_Classify/data/seed_comparison.csv
# 对比维度：
#   - 每 Run：恶性率 / 收敛轮数与最终 agreement / 各组织恶性率
#   - 两种子间：逐细胞标签一致率、最终签名 top10 Jaccard
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

DATA <- "05_Iterative_Classify/data"
ANNO <- "03_Cluster_Annotate/data/cell_annotation.csv"
anno <- fread(ANNO)[, .(barcode, tissue)]

RUN_DIRS <- c(
  run1 = "run1_epi_chief_logFC",
  run6 = "run6_epi_chief_logFC_lowconf_excl"
)

read_run <- function(dir) {
  p <- fromJSON(file.path(dir, "run_params.json"))
  il <- fread(file.path(dir, "iteration_log.csv"))
  sig <- fread(file.path(dir, "final_signatures.csv"))
  list(
    n = p$n_epi_cells,
    pct = p$pct_malignant,
    n_iter = p$n_iterations,
    agree = if ("agreement" %in% names(il)) tail(il$agreement[!is.na(il$agreement)], 1) else NA,
    sig_mal = sig$malignant_genes[!is.na(sig$malignant_genes)],
    sig_non = sig$non_malignant_genes[!is.na(sig$non_malignant_genes)]
  )
}

# A/B 侧恶性率：labels 回连注释重算（tissue_summary.csv 已停产）
pct_tissue <- function(lab) {
  m <- merge(lab, anno, by.x = "cell", by.y = "barcode")
  round(c(A = mean(m$label[m$tissue == "A"] == "malignant") * 100,
          B = mean(m$label[m$tissue == "B"] == "malignant") * 100), 1)
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

rows <- list()
for (rn in names(RUN_DIRS)) {
  d_main <- file.path(DATA, RUN_DIRS[[rn]])
  d_test <- file.path(DATA, paste0("testseed_", RUN_DIRS[[rn]]))
  if (!dir.exists(d_main) || !dir.exists(d_test)) {
    cat(sprintf("SKIP %s: missing dir\n", rn))
    next
  }
  rn_main <- read_run(d_main)
  rn_test <- read_run(d_test)

  lab_main <- fread(file.path(d_main, "labels.csv"))
  lab_test <- fread(file.path(d_test, "labels.csv"))
  m <- merge(lab_main, lab_test, by = "cell", suffixes = c("_main", "_test"))
  label_agree <- mean(m$label_main == m$label_test)
  ts_main <- pct_tissue(lab_main)
  ts_test <- pct_tissue(lab_test)

  rows[[rn]] <- data.frame(
    run = rn,
    n_main = rn_main$n, pct_main = rn_main$pct, iter_main = rn_main$n_iter,
    agree_main = rn_main$agree,
    n_test = rn_test$n, pct_test = rn_test$pct, iter_test = rn_test$n_iter,
    agree_test = rn_test$agree,
    label_agree = label_agree,
    jac_mal_top10 = jaccard(head(rn_main$sig_mal, 10), head(rn_test$sig_mal, 10)),
    jac_non_top10 = jaccard(head(rn_main$sig_non, 10), head(rn_test$sig_non, 10)),
    adjA_main = unname(ts_main["A"]), adjA_test = unname(ts_test["A"]),
    tumorB_main = unname(ts_main["B"]), tumorB_test = unname(ts_test["B"])
  )
}

tab <- rbindlist(rows, use.names = TRUE)

out_csv <- file.path(DATA, "seed_comparison.csv")
fwrite(tab, out_csv)
cat("\n===== seed comparison (main vs test, run1 + run6) =====\n")
print(tab, row.names = FALSE)
cat(sprintf("\nwritten: %s\n", out_csv))
