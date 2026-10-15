#!/usr/bin/env Rscript
# ============================================================
# 检验敏感性扫描里的 run1 是否与主流程 05 的 labels 逐位一致
# 两者配置相同，理论应逐位一致（scores/iteration_log 为同一函数 +
# 同一随机流的中间产物，labels 一致即整体一致，无需单独比对）
# ============================================================
suppressPackageStartupMessages(library(data.table))

DATA <- "05_Iterative_Classify/data"

main_lab <- fread(file.path(DATA, "labels.csv"))
run1_lab <- fread(file.path(DATA, "run1_epi_chief_logFC", "labels.csv"))

cat(sprintf("main 05 labels: %d cells (malignant=%d, non=%d)\n",
            nrow(main_lab), sum(main_lab$label == "malignant"),
            sum(main_lab$label == "non_malignant")))
cat(sprintf("run1 labels:    %d cells (malignant=%d, non=%d)\n",
            nrow(run1_lab), sum(run1_lab$label == "malignant"),
            sum(run1_lab$label == "non_malignant")))

m <- merge(main_lab, run1_lab, by = "cell", suffixes = c("_main", "_run1"))
n_mis <- sum(m$label_main != m$label_run1)
cat(sprintf("\nrun1 anchor: %d cells matched, %d mismatches (bitwise)\n",
            nrow(m), n_mis))
if (nrow(m) != nrow(main_lab) || nrow(m) != nrow(run1_lab))
  cat("WARNING: cell sets differ between main 05 and run1!\n")
if (n_mis == 0) cat("ANCHOR PASS: run1 == main 05 (bit-identical)\n") else
  cat("ANCHOR FAIL\n")
