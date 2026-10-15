#!/usr/bin/env Rscript
# ============================================================
# 03_Cluster_Annotate 嫌疑污染枚举诊断 — contam_pairs_enumerate.R
#
# 用途：枚举嫌疑污染细胞的构成——vote × cluster-mode 配对、
#       谱系方向、患者分布；并 sanity-check 逐细胞导出计算的
#       列联表与导出的 annotation_qc_confusion.csv 一致
#
# 输入：03_Cluster_Annotate/data/annotation_qc_umap.csv +
#       annotation_qc_confusion.csv（annotation_qc_export.R 产出）
# 产出：无文件——统计打印到控制台；结论已固化进
#       本目录 NOTES_suspected_contamination.md
#
# 关键设计要点：
#   - 零重算：只读两个 QC CSV，data.table 即可，<1 min
#   - 谱系三分：epi / immune / stroma；dir = vote 谱系 →
#     cluster-mode 谱系（跨谱系 = 嫌疑，同谱系豁免）
#
# Usage (Windows native R 4.5.1)：
#   "E:\R-4.5.1\bin\Rscript.exe" 03_Cluster_Annotate\scripts\contamination_analysis\contam_pairs_enumerate.R
# ============================================================
suppressPackageStartupMessages(library(data.table))

BASE <- "03_Cluster_Annotate"
d <- fread(file.path(BASE, "data", "annotation_qc_umap.csv"),
           select = c("annotation", "cluster_annotation", "suspected_contam",
                      "low_conf", "patient"))
s <- d[suspected_contam == TRUE]
cat(sprintf("suspects: %d (%.2f%%)\n\n", nrow(s), 100 * nrow(s) / nrow(d)))

# ── vote × cluster-mode 配对（top 20）────────────────────
cat("== vote x cluster-mode pairs (top 20) ==\n")
print(s[, .N, by = .(annotation, cluster_annotation)][order(-N)][1:20])

cat("\n== by vote type ==\n")
print(s[, .N, by = annotation][order(-N)])

cat("\n== by cluster mode ==\n")
print(s[, .N, by = cluster_annotation][order(-N)])

# ── 谱系方向 ─────────────────────────────────────────────
cat("\n== by lineage direction ==\n")
LIN <- c(Epithelial = "epi", Chief = "epi", Parietal = "epi",
         Enteroendocrine = "epi", T_cells = "immune", B_cells = "immune",
         NK_cells = "immune", Myeloid = "immune", Mast_cells = "immune",
         Plasma_cells = "immune", Fibroblasts = "stroma",
         Endothelial = "stroma")
s[, dir := paste0(LIN[annotation], "->", LIN[cluster_annotation])]
print(s[, .N, by = dir][order(-N)])

cat("\n== by patient ==\n")
print(s[, .N, by = patient][order(-N)])

cat(sprintf("\nlow_conf rate: suspects %.1f%% | others %.1f%%\n",
            100 * mean(s$low_conf), 100 * mean(d$low_conf[!d$suspected_contam])))

# ── sanity：逐细胞导出计算的列联表 vs 导出的 qc_confusion.csv ──
conf_exp <- fread(file.path(BASE, "data", "annotation_qc_confusion.csv"))
conf_calc <- d[, .N, by = .(vote = annotation, cluster_mode = cluster_annotation)]
chk <- merge(conf_exp, conf_calc, by = c("vote", "cluster_mode"), all = TRUE)
chk[is.na(N.x), N.x := 0L]; chk[is.na(N.y), N.y := 0L]
cat(sprintf("\nsanity vs qc_confusion.csv: max|d| = %d (expect 0)\n",
            max(abs(chk$N.x - chk$N.y))))
