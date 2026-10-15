#!/usr/bin/env Rscript
# ============================================================
# 03_Cluster_Annotate 嫌疑污染患者分布诊断 — contam_patient_dist.R
#
# 用途：嫌疑污染细胞的 per-patient 分布——均匀散布还是集中于
#       少数患者？并给出 top 患者的谱系方向构成
#
# 输入：03_Cluster_Annotate/data/annotation_qc_umap.csv
#       （annotation_qc_export.R 产出）
# 产出：无文件——统计打印到控制台；结论已固化进
#       本目录 NOTES_suspected_contamination.md
#
# 关键设计要点：
#   - 零重算：只读一个 QC CSV，data.table 即可，<1 min
#   - 集中度度量：top-3/top-5/top-10 患者占全部嫌疑的比例 +
#     Gini（0 = 均匀，1 = 单患者全占）+ 嫌疑率与患者细胞数的
#     相关性（排除"患者大小驱动"的替代解释）
#
# Usage (Windows native R 4.5.1)：
#   "E:\R-4.5.1\bin\Rscript.exe" 03_Cluster_Annotate\scripts\contamination_analysis\contam_patient_dist.R
# ============================================================
suppressPackageStartupMessages(library(data.table))

BASE <- "03_Cluster_Annotate"
d <- fread(file.path(BASE, "data", "annotation_qc_umap.csv"),
           select = c("annotation", "cluster_annotation", "suspected_contam", "patient"))
setnames(d, "annotation", "vote")
s <- d[suspected_contam == TRUE]

cat("== per-patient suspects (ordered by count) ==\n")
p <- d[, .(total = .N, suspects = sum(suspected_contam == TRUE)), by = patient]
p[, rate_pct := 100 * suspects / total]
print(p[order(-suspects)])

cat("\n== spread stats ==\n")
cat(sprintf("n patients: %d\n", nrow(p)))
cat(sprintf("patients with >=1 suspect: %d / %d\n", sum(p$suspects > 0), nrow(p)))
cat(sprintf("rate range: %.2f - %.2f%% (median %.2f, mean %.2f)\n",
            min(p$rate_pct), max(p$rate_pct), median(p$rate_pct),
            mean(p$rate_pct)))
ord <- p[order(-suspects)]
cat(sprintf("top-3 patients hold %.1f%% of all suspects\n",
            100 * sum(ord$suspects[1:3]) / sum(p$suspects)))
cat(sprintf("top-5: %.1f%%; top-10: %.1f%%\n",
            100 * sum(ord$suspects[1:5]) / sum(p$suspects),
            100 * sum(ord$suspects[1:10]) / sum(p$suspects)))
# Gini（跨患者的嫌疑数）：0 = 完全均匀，1 = 单患者全占
x <- sort(p$suspects); n <- nrow(p)
cat(sprintf("Gini (suspects across patients): %.3f\n",
            (2 * sum((1:n) * x) - (n + 1) * sum(x)) / (n * sum(x))))
cat(sprintf("corr(rate_pct, patient cell count): %.3f (is rate size-driven?)\n",
            cor(p$rate_pct, p$total)))

cat("\n== per-patient lineage direction (top 6 patients) ==\n")
LIN <- c(Epithelial = "epi", Chief = "epi", Parietal = "epi",
         Enteroendocrine = "epi", T_cells = "immune", B_cells = "immune",
         NK_cells = "immune", Myeloid = "immune", Mast_cells = "immune",
         Plasma_cells = "immune", Fibroblasts = "stroma",
         Endothelial = "stroma")
s[, dir := paste0(LIN[vote], "->", LIN[cluster_annotation])]
top6 <- ord$patient[1:6]
print(s[patient %in% top6, .N, by = .(patient, dir)][order(patient, -N)])
