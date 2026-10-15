#!/usr/bin/env Rscript
# ============================================================
# 03_Cluster_Annotate 注释 QC — 导出部分（WSL）
#
# 从已注释 RDS 计算 QC 表（03 重跑之后），供 Windows 图件脚本
# annotation_qc_figures.R 使用。表回答注释质量审查的问题：
#   1. marker 是否恰当？          → 逐基因 violin 数据
#      （annotation_qc_violin.csv，分层抽样）
#   2. 投票 vs 簇 mode 一致度      → 混淆矩阵
#      （annotation_qc_confusion.csv）
#   3. 投票边际如何分布？          → 逐注释 z-delta 统计
#      （annotation_qc_stats.csv）
#   4. marker 是否具代表性？       → 逐基因 marker 审计：检出率
#      + 表达 + 特异性（annotation_qc_marker_stats.csv +
#      annotation_qc_marker_summary.csv）
#   5. UMAP feature-plot 数据：逐细胞UMAP 坐标 +
#      12 个代表 marker（log1p）+ 注释 meta（含suspected_contam flag）
#      ——供 Fig-E feature 图、注释总览、
#      簇混合诊断与下游 DEG 敏感性过滤使用
#      （annotation_qc_umap.csv）
#
# 输入：03_Cluster_Annotate/data/gse150290_annotated.rds
# 产出：03_Cluster_Annotate/data/annotation_qc_*.csv
#
# Usage (WSL):
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 03_Cluster_Annotate/scripts/annotation_qc_export.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

BASE     <- "."
RDS_PATH <- file.path(BASE, "03_Cluster_Annotate/data/gse150290_annotated.rds")
OUT_DIR  <- file.path(BASE, "03_Cluster_Annotate/data")
TMP_DIR  <- "/tmp/annot_qc_r"
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_STRAT     <- 600   # 每注释分层抽样数（violin 数据）
RANDOM_SEED <- 42

# 与 cluster_annotate.R 的 MARKER_DICT 保持同步。
MARKER_DICT <- list(
  T_cells         = c("CD3D", "CD3E", "CD2", "CD7", "LCK", "IL32", "CCL5"),
  B_cells         = c("CD79A", "CD79B", "MS4A1", "CD19", "PAX5"),
  NK_cells        = c("NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1"),
  Myeloid         = c("CD14", "CD68", "FCGR3A", "CSF1R", "ITGAM"),
  Mast_cells      = c("KIT", "TPSAB1", "CPA3", "HDC"),
  Plasma_cells    = c("MZB1", "SDC1"),
  Fibroblasts     = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "MMP2", "PDGFRA", "MYL9"),
  Endothelial     = c("CDH5", "VWF", "ENG", "CLDN5", "PLVAP", "KDR", "PTPRB"),
  Epithelial      = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "CLDN4", "CLDN7"),
  Enteroendocrine = c("CHGA", "GAST", "PROX1", "GHRL", "SST", "PCSK1N"),
  Parietal        = c("GIF", "ATP4B", "ATP4A"),
  Chief           = c("PGC", "LIPF", "PGA3", "PGA4")
)

log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
}

t0 <- Sys.time()
log(sprintf("Loading %s ...", RDS_PATH))
obj <- readRDS(RDS_PATH)
log(sprintf("  %d cells x %d genes [%.0fs]",
            ncol(obj), nrow(obj), difftime(Sys.time(), t0, units = "secs")))

data_m <- obj[["RNA"]]$data
cnts_m <- obj[["RNA"]]$counts

types <- sort(unique(obj$annotation))
log(sprintf("  %d annotations: %s", length(types), paste(types, collapse = ", ")))
for (p in names(MARKER_DICT)) {
  miss <- setdiff(MARKER_DICT[[p]], rownames(data_m))
  if (length(miss) > 0)
    log(sprintf("  %-16s: %d/%d genes found (missing: %s)",
                p, length(MARKER_DICT[[p]]) - length(miss),
                length(MARKER_DICT[[p]]), paste(miss, collapse = ",")))
}

# ── 1. Violin 数据：分层逐基因表达（长表）─────────────────────
set.seed(RANDOM_SEED)
keep <- unlist(lapply(types, function(ct) {
  bcs <- colnames(obj)[obj$annotation == ct]
  sample(bcs, min(N_STRAT, length(bcs)))
}))
panel_genes <- unique(unlist(MARKER_DICT))
panel_genes <- panel_genes[panel_genes %in% rownames(data_m)]
gene2panel <- unlist(lapply(names(MARKER_DICT), function(p)
  setNames(rep(p, length(MARKER_DICT[[p]])), MARKER_DICT[[p]])))

violin_rows <- vector("list", length(types))
for (i in seq_along(types)) {
  ct <- types[i]
  bcs <- keep[obj$annotation[keep] == ct]
  m <- as.matrix(data_m[panel_genes, bcs, drop = FALSE])
  violin_rows[[i]] <- data.table(
    cell       = rep(bcs, each = length(panel_genes)),
    annotation = ct,
    gene       = rep(panel_genes, times = length(bcs)),
    expr       = as.vector(m)
  )
}
violin <- rbindlist(violin_rows)
violin[, panel := gene2panel[gene]]
fwrite(violin, file.path(TMP_DIR, "annotation_qc_violin.csv"))
log(sprintf("  violin data: %d rows, %d cells [%.0fs]",
            nrow(violin), length(keep), difftime(Sys.time(), t0, units = "secs")))

# ── 2. 混淆：逐细胞投票 × 簇 mode ────────────────────────────
conf <- as.data.table(table(vote = obj$annotation,
                            cluster_mode = obj$cluster_annotation))
conf <- conf[N > 0]
fwrite(conf, file.path(TMP_DIR, "annotation_qc_confusion.csv"))
log(sprintf("  confusion: %d non-empty cells", nrow(conf)))

# ── 3. 逐注释投票边际统计 ─────────────────────────────────────
# z-delta（投票边际）= 本方法的一致性度量。标签转移式置信度
# （softmax 概率）不适用于 marker 投票注释，已退役（见 header）。
stats_dt <- data.table(
  annotation = obj$annotation,
  z_delta    = obj$cell_z_delta,
  low_conf   = obj$annotation_low_conf
)
stats <- stats_dt[, .(
  n             = .N,
  low_conf_pct  = round(100 * mean(low_conf), 2),
  z_delta_median = round(median(z_delta), 4),
  z_delta_q25   = round(quantile(z_delta, 0.25), 4),
  z_delta_q75   = round(quantile(z_delta, 0.75), 4)
), by = annotation][order(-n)]
fwrite(stats, file.path(TMP_DIR, "annotation_qc_stats.csv"))
log("  stats done")

# ── 4. 逐基因 marker 审计（逐类型检出率 + 表达）───────────────
# 回答 "marker 是否有代表性"：每个面板基因在 12 个注释中的
# 检出比例（raw counts，ambient 按设计计入）与 mean log1p 表达。
# 对 ambient 的操作性防御是 score_genes 的对照扣除，不是这张表
# ——广泛检出的基因（如 PGC）还必须显示大的表达水平差距
# （10-100×）才算有效。
marker_genes <- unique(unlist(MARKER_DICT))
marker_genes <- marker_genes[marker_genes %in% rownames(cnts_m)]
marker_rows <- lapply(types, function(ct) {
  bcs <- colnames(obj)[obj$annotation == ct]
  cm <- cnts_m[marker_genes, bcs, drop = FALSE]
  dm <- data_m[marker_genes, bcs, drop = FALSE]
  data.table(
    gene       = marker_genes,
    annotation = ct,
    det_pct    = round(100 * Matrix::rowMeans(cm > 0), 2),
    mean_log1p = round(Matrix::rowMeans(dm), 4)
  )
})
marker_stats <- rbindlist(marker_rows)
marker_stats[, panel := gene2panel[gene]]
fwrite(marker_stats, file.path(TMP_DIR, "annotation_qc_marker_stats.csv"))

# 逐基因汇总：自身类型检出/水平 vs 最佳其他类型
marker_sum <- marker_stats[, .(
  det_own        = det_pct[annotation == panel],
  mean_own       = mean_log1p[annotation == panel],
  det_max_other  = max(det_pct[annotation != panel]),
  mean_max_other = max(mean_log1p[annotation != panel])
), by = .(gene, panel)]
marker_sum[, spec_det := round(det_own / pmax(det_max_other, 0.1), 2)]
marker_sum[, expr_ratio := round(mean_own / pmax(mean_max_other, 0.01), 2)]
fwrite(marker_sum, file.path(TMP_DIR, "annotation_qc_marker_summary.csv"))
log(sprintf("  marker audit: %d genes x %d types", nrow(marker_sum), length(types)))

# ── 5. UMAP feature-plot 数据（12 个代表 marker）──────────────
# Zhang 2021 Fig E 风格：逐类型代表 marker 表达画在 UMAP 上。
# 宽表 CSV：cell, umap_1, umap_2 + 每个 marker 基因一列（log1p
# 归一化 data slot）。与 annotation_qc_figures.R 的 REP_MARKER
# （代表 marker violin 与 UMAP feature 图共用）保持同步。
REP_MARKER <- c(T_cells = "CD3D", B_cells = "MS4A1", NK_cells = "NKG7",
                Myeloid = "CD14", Mast_cells = "TPSAB1",
                Plasma_cells = "MZB1", Fibroblasts = "COL1A1",
                Endothelial = "VWF", Epithelial = "EPCAM",
                Enteroendocrine = "CHGA", Parietal = "ATP4A",
                Chief = "PGC")
rep_genes <- unname(REP_MARKER)
rep_genes <- rep_genes[rep_genes %in% rownames(data_m)]
um <- as.data.frame(Embeddings(obj, "umap"))
expr_rep <- as.matrix(data_m[rep_genes, , drop = FALSE])
# 跨谱系混合 flag：投票与簇 mode 属于不同谱系（epi / immune /
# stroma）。同谱系混合（如 Chief↔Epithelial）是分化连续谱，非
# 污染。跨谱系混合与转录组身份（簇/UMAP）矛盾；唯一自洽的解释
# 是外来 RNA 带入（裂解泄漏 / doublet）或投影错误——下游 DEG
# 应把两者都当疑似对象。Suspect = flag，永不裁决；无表达阈值。
LINEAGE <- c(Epithelial = "epi", Chief = "epi", Parietal = "epi",
             Enteroendocrine = "epi",
             T_cells = "immune", B_cells = "immune", NK_cells = "immune",
             Myeloid = "immune", Mast_cells = "immune", Plasma_cells = "immune",
             Fibroblasts = "stroma", Endothelial = "stroma")
lin_vote  <- LINEAGE[obj$annotation]
lin_clust <- LINEAGE[obj$cluster_annotation]
suspected_contam <- !is.na(lin_vote) & !is.na(lin_clust) & lin_vote != lin_clust
# meta 列与 colnames(obj) 同序 = rownames(um) 同序。
umap_dt <- data.table(cell = rownames(um), umap_1 = um[, 1], umap_2 = um[, 2],
                      annotation = obj$annotation,
                      cluster_annotation = obj$cluster_annotation,
                      low_conf = obj$annotation_low_conf,
                      z_delta = obj$cell_z_delta,
                      patient = obj$patient,
                      suspected_contam = suspected_contam)
for (g in rep_genes) set(umap_dt, j = g, value = expr_rep[g, ])
fwrite(umap_dt, file.path(TMP_DIR, "annotation_qc_umap.csv"))
log(sprintf("  umap feature data: %d cells x %d markers + meta",
            nrow(umap_dt), length(rep_genes)))

# ── 拷贝到仓库目录 ─────────────────────────────────────────────
for (f in list.files(TMP_DIR, pattern = "annotation_qc_.*\\.csv$")) {
  file.copy(file.path(TMP_DIR, f), file.path(OUT_DIR, f), overwrite = TRUE)
}
log(sprintf("Done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))
