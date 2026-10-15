#!/usr/bin/env Rscript
# ============================================================
# 03_Cluster_Annotate：Leiden 聚类 + 12 类 marker 注释
#
# 输入：02_Normalize_Batch/data/gse150290_norm_harmony.rds
# 产出：03_Cluster_Annotate/data/gse150290_annotated.rds
#       + cell_annotation.csv、cluster_annotation.csv、params json
#       （全部图件：annotation_qc_export.R → annotation_qc_figures.R）
#
# 关键设计要点：
#   - Leiden r=0.5 于 harmony 50 PC、k=15（igraph-Leiden 实现；
#     注释由 marker 驱动，簇号只作参照）
#   - marker 集（12 类）来源细分与 hg19 适配：详见
#     03_Cluster_Annotate/README.md
#   - 逐细胞 z-score（12 类 marker 均值）→ argmax 投票；低置信 =
#     z-delta（最优 − 次优）< 0.5
#   - 红细胞剔除：globin（HBA1/HBA2/HBB）计数占比 > 5%，在 Leiden
#     聚类之前完成（globin 仅凭 counts 即可计算，簇号对最终细胞集
#     自洽）
#   - 死细胞簇剔除：簇中位 MT% > 25 且 max 中位 marker z < 1.5
#     （cls7 高 MT 规则的数据驱动镜像；z 门槛 1→1.5 收紧：濒死上皮
#     保留身份 marker——如 Cls 7 Epithelial z≈1.42 旧门槛漏网，
#     Parietal 高 MT 系线粒体天然丰富（z=2.78）不受影响）
#   - 最终注释 = 逐细胞投票（非簇 mode）；簇 mode 仅参照列
#     （cluster_annotation + purity）
#   - CD3D 谱系救援：NK 面板投票 + CD3D+ → T_cells（CD3D 为确定性
#     T 谱系 marker；T 面板变薄把 CD3D+ T 细胞拉向 NK 面板——39%
#     的 NK 投票者 CD3D+）
#   - 图件分工：全部图由 annotation_qc_export.R（WSL 导出）+
#     annotation_qc_figures.R（Windows R + ragg, 300 dpi）产出，
#     本脚本不画图
#
# Usage:
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 03_Cluster_Annotate/scripts/cluster_annotate.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

# ── 路径 ─────────────────────────────────────────────────────
NORM_DIR <- "02_Normalize_Batch/data"
OUT_DIR  <- "03_Cluster_Annotate/data"
TMP_DIR  <- "/tmp/clust_r"
OUT_TMP  <- file.path(TMP_DIR, "out")

# ── 参数 ──────────────────────────────────────────────────────
N_NEIGHBORS <- 15
N_PCS       <- 50
LEIDEN_RES  <- 0.5
RANDOM_SEED <- 42
LOW_CONF_ZDELTA <- 0.5    # 逐细胞 z-delta < 0.5 = 低置信
GLOBIN_GENES     <- c("HBA1", "HBA2", "HBB")
GLOBIN_FRAC_MAX  <- 0.05  # globin 计数占比 > 5% = 红细胞
DEAD_MT_MEDIAN   <- 25    # 簇中位 MT% > 25 ...
DEAD_Z_MAX       <- 1.5   # ... 且 max 中位 marker z < 1.5 → 死细胞簇
# 收紧理由（2026-10-15）：旧门槛 z<1 放走全场最脏簇（Cls 7
# Epithelial：中位 MT 全场最高、纯度 51.8% 全场最低，z≈1.42 逃逸）——
# 濒死上皮仍表达上皮标记，z 判据区分不了"线粒体丰富的活类型"
# （Parietal z=2.78 系泌酸耗能，须保留）与"保留身份标记的濒死细胞"；
# 1.5 恰好只抓 Cls 7，Parietal 安全

# 12 类 marker 集——来源细分与 hg19 适配：见 03_Cluster_Annotate/
# README.md（同步于 params json 的 marker_source；与
# annotation_qc_export.R 的 MARKER_DICT 保持同步）
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

dir.create(OUT_TMP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
log_fh <- file(file.path(TMP_DIR, "cluster_annotate.log"), open = "a")
log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  writeLines(line, log_fh)
  flush(log_fh)
}

# ── 1. 读入 ───────────────────────────────────────────────────
t0 <- Sys.time()
log(sprintf("Loading %s ...", NORM_DIR))
obj <- readRDS(file.path(NORM_DIR, "gse150290_norm_harmony.rds"))
log(sprintf("  %d cells x %d genes [%.0fs]",
            ncol(obj), nrow(obj), difftime(Sys.time(), t0, units = "secs")))

# ── 2. 红细胞：globin 计数占比 > 5% ──────────────────────────
# （在聚类前：globin 仅凭 counts 即可计算，Leiden 跑在最终细胞集
# 上、簇号自洽）
t0 <- Sys.time()
counts_m <- obj[["RNA"]]$counts
globin_genes <- GLOBIN_GENES[GLOBIN_GENES %in% rownames(counts_m)]
globin_frac <- colSums(counts_m[globin_genes, , drop = FALSE]) / colSums(counts_m)
obj$globin_frac <- globin_frac
is_rbc <- globin_frac > GLOBIN_FRAC_MAX
log(sprintf("  RBC candidates (globin frac > %.2f): %d (%.2f%%)",
            GLOBIN_FRAC_MAX, sum(is_rbc), sum(is_rbc) / ncol(obj) * 100))
obj <- subset(obj, cells = colnames(obj)[!is_rbc])
log(sprintf("  Post-RBC: %d cells [%.0fs]", ncol(obj),
            difftime(Sys.time(), t0, units = "secs")))

# ── 3. harmony PC 上的 Leiden 聚类 ────────────────────────────
t0 <- Sys.time()
log(sprintf("FindNeighbors (harmony, dims=1:%d, k=%d) + Leiden (r=%.1f) ...",
            N_PCS, N_NEIGHBORS, LEIDEN_RES))
obj <- FindNeighbors(obj, reduction = "harmony", dims = seq_len(N_PCS),
                     k.param = N_NEIGHBORS, verbose = FALSE)
obj <- FindClusters(obj, algorithm = 4, resolution = LEIDEN_RES,
                    random.seed = RANDOM_SEED, verbose = FALSE)
n_cl <- length(unique(Idents(obj)))
log(sprintf("  %d clusters [%.0fs]", n_cl, difftime(Sys.time(), t0, units = "secs")))

# ── 3. UMAP（仅可视化）────────────────────────────────────────
t0 <- Sys.time()
log("RunUMAP (harmony) ...")
obj <- RunUMAP(obj, reduction = "harmony", dims = seq_len(N_PCS),
               n.neighbors = N_NEIGHBORS, verbose = FALSE)
log(sprintf("  UMAP done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 4. Marker 打分：25 bins × 100 同 bin 对照协议 ────────────
# score = 面板基因上 (expr_g − 100 个同 bin 对照基因均值) 的均值，
# 25 个表达 bin。对照扣除去掉 ambient 表达基线（例：PGC ambient
# RNA 在 57% 全部细胞检出）——裸 colMeans 会给每个细胞加常数
# Chief 偏移（已验证：裸均值投票把 Chief 从 9.9k 膨胀到 27.6k，
# 参考实现为 11.5k）。
t0 <- Sys.time()
log("Computing 12-type marker scores (score_genes protocol) ...")
data_m <- obj[["RNA"]]$data

# 打分协议：25 bins × 100 同 bin 对照（Puram 2017），
# 稀疏安全；移植自 05 step2_iterative_classify.R（Bug-2 修复：
# 逐细胞对照均值用 colMeans）
add_module_score <- function(expr, gene_list, ctrl_size = 100, n_bins = 25) {
  all_genes <- rownames(expr)
  valid <- intersect(gene_list, all_genes)
  if (length(valid) == 0) return(rep(0, ncol(expr)))

  gene_means <- rowMeans(expr)
  gene_bin <- setNames(as.integer(cut(gene_means, breaks = n_bins, labels = FALSE)),
                       all_genes)
  gene_bin[is.na(gene_bin) & gene_means <= min(gene_means)] <- 1
  gene_bin[is.na(gene_bin) & gene_means >= max(gene_means)] <- n_bins

  scores <- rep(0, ncol(expr))
  for (g in valid) {
    g_bin <- gene_bin[g]
    same_bin <- setdiff(which(gene_bin == g_bin), which(rownames(expr) == g))
    if (length(same_bin) == 0) {
      ctrl_idx <- sample.int(nrow(expr), min(ctrl_size, nrow(expr)))
    } else {
      ctrl_idx <- sample(same_bin, min(ctrl_size, length(same_bin)))
    }
    gene_expr <- expr[g, , drop = FALSE]
    # 同 bin 对照的逐细胞均值（colMeans：每细胞一个值；rowMeans 会
    # 给出每基因一个值并静默按细胞循环回收）
    ctrl_mean <- colMeans(expr[ctrl_idx, , drop = FALSE])
    scores <- scores + as.numeric(gene_expr) - ctrl_mean
  }
  scores / max(length(valid), 1)
}

set.seed(RANDOM_SEED)   # 对照基因抽样可复现性
score_mat <- matrix(0, nrow = ncol(obj), ncol = length(MARKER_DICT),
                    dimnames = list(colnames(obj), names(MARKER_DICT)))
for (ct in names(MARKER_DICT)) {
  genes <- MARKER_DICT[[ct]]
  found <- genes[genes %in% rownames(data_m)]
  missing <- genes[!genes %in% rownames(data_m)]
  if (length(missing) > 0)
    log(sprintf("  %-16s: %d/%d markers found (missing: %s)",
                ct, length(found), length(genes), paste(missing, collapse = ",")))
  score_mat[, ct] <- add_module_score(data_m, genes)
}
log(sprintf("  Scores done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 5. 逐细胞 z-score → argmax 投票 ──────────────────────────
t0 <- Sys.time()
log("Annotating: per-cell z-score -> argmax vote ...")
row_mean <- rowMeans(score_mat)
row_sd <- sqrt(rowSums((score_mat - row_mean)^2) / (ncol(score_mat) - 1))
z_mat <- (score_mat - row_mean) / row_sd
z_mat[!is.finite(z_mat)] <- 0    # 全零 marker 表达的细胞 → z=0
best_idx <- max.col(z_mat, ties.method = "first")
cell_vote <- colnames(z_mat)[best_idx]
# z-delta = 最优 − 次优
z_sorted <- apply(z_mat, 1, sort, decreasing = TRUE)
cell_z_delta <- z_sorted[1, ] - z_sorted[2, ]

# CD3D 谱系救援：NK 面板投票 + CD3D 检出 = 被误投 NK 的 T 谱系
# 细胞。T 面板失去 TRAC/TRBC1（本 hg19 参考缺失），把 CD3D+ T
# 细胞拉向 NK 面板（重投票实验：39% 的 NK 投票者 CD3D+）。CD3D
# 是确定性 T 谱系 marker，故救援误投、保留 CD3D- 的真 NK。
if ("CD3D" %in% rownames(obj[["RNA"]]$counts)) {
  cd3d_pos <- as.vector(obj[["RNA"]]$counts["CD3D", ] > 0)
  rescue <- cell_vote == "NK_cells" & cd3d_pos
  n_nk_votes <- sum(cell_vote == "NK_cells")
  if (n_nk_votes > 0) {
    log(sprintf("  CD3D rescue: %d/%d NK votes -> T_cells (%.1f%%)",
                sum(rescue), n_nk_votes, 100 * sum(rescue) / n_nk_votes))
  }
  cell_vote[rescue] <- "T_cells"
} else {
  log("  WARNING: CD3D not found in counts; CD3D rescue skipped")
}

obj$cell_vote <- cell_vote
obj$cell_z_delta <- cell_z_delta
obj$annotation <- cell_vote    # ★ 最终注释 = 逐细胞投票（CD3D 救援后）

clusters <- as.character(Idents(obj))
cl_anno <- tapply(cell_vote, clusters, function(v) {
  names(sort(table(v), decreasing = TRUE))[1]
})
cl_purity <- tapply(cell_vote, clusters, function(v) {
  max(table(v)) / length(v)
})
obj$cluster_annotation <- unname(cl_anno[clusters])   # 仅参照
obj$cluster_purity <- unname(cl_purity[clusters])
obj$annotation_low_conf <- cell_z_delta < LOW_CONF_ZDELTA

log("  Cluster annotation (mode vote, reference) + purity:")
for (cl in sort(unique(clusters))) {
  log(sprintf("    Cls %-3s -> %-16s purity %.3f", cl, cl_anno[[cl]], cl_purity[[cl]]))
}
n_low <- sum(obj$annotation_low_conf)
log(sprintf("  Per-cell low-confidence (z-delta < %.1f): %d (%.1f%%) [%.0fs]",
            LOW_CONF_ZDELTA, n_low, n_low / ncol(obj) * 100,
            difftime(Sys.time(), t0, units = "secs")))

# ── 6. 死细胞簇：中位 MT% > 25 且 max 中位 z < 1.5 ─────────────
cl_mt <- tapply(obj$percent.mt, clusters, median)
cl_zmax <- tapply(seq_len(ncol(obj)), clusters, function(idx) {
  max(colMeans(z_mat[idx, , drop = FALSE]))
})
dead_cl <- names(cl_mt)[cl_mt > DEAD_MT_MEDIAN & cl_zmax < DEAD_Z_MAX]
log(sprintf("  Dead-cell clusters (median MT%%>%.0f & max median z<%.1f): %s",
            DEAD_MT_MEDIAN, DEAD_Z_MAX,
            if (length(dead_cl) > 0) paste(dead_cl, collapse = ",") else "none"))

# ── 6b. 红细胞排查：逐簇 globin 分布 ─────────────────────────
# globin 占比 > 5% 已在聚类前剔除；此处核查剔除后是否仍有簇整体
# 偏高（部分红细胞污染漏网的信号）
cl_globin_med <- tapply(obj$globin_frac, clusters, median)
cl_globin_p90 <- tapply(obj$globin_frac, clusters, function(x) quantile(x, 0.9))
cl_globin_max <- tapply(obj$globin_frac, clusters, max)
n_globin_gt2  <- tapply(obj$globin_frac > 0.02, clusters, sum)
for (cl in sort(unique(clusters))) {
  log(sprintf("    Cls %-3s globin: med %.4f / p90 %.4f / max %.4f / >2%%: %d",
              cl, cl_globin_med[[cl]], cl_globin_p90[[cl]],
              cl_globin_max[[cl]], n_globin_gt2[[cl]]))
}

# ── 7. 应用剔除 ───────────────────────────────────────────────
exclude <- clusters %in% dead_cl
obj$exclude_reason <- ifelse(clusters %in% dead_cl, "dead_cluster", "")
log(sprintf("  Total excluded: %d (%.2f%%) [%.0fs]",
            sum(exclude), sum(exclude) / ncol(obj) * 100,
            difftime(Sys.time(), t0, units = "secs")))
obj <- subset(obj, cells = colnames(obj)[!exclude])
log(sprintf("  Remaining: %d cells in %d clusters", ncol(obj), length(unique(Idents(obj)))))
# clusters 更新为剔除后的簇号——cl_tab/cell_annotation 与剔除后
# 细胞集保持一致（cl_anno/cl_purity/cl_mt/cl_zmax 是按簇聚合的
# 向量，剔除前算好、按存活簇号索引即可；死细胞簇被剔除后不再
# 出现在表中，由 params json 的 dead_clusters 字段记录）
clusters <- as.character(Idents(obj))

# ── 9. 保存 ──────────────────────────────────────────────────
t0 <- Sys.time()
log("Saving Seurat object (RDS) ...")
saveRDS(obj, file.path(OUT_TMP, "gse150290_annotated.rds"), compress = FALSE)

cell_anno <- data.frame(
  barcode = colnames(obj),
  sample = obj$sample,
  patient = obj$patient,
  tissue = obj$tissue,
  cluster = as.character(Idents(obj)),
  annotation = obj$annotation,
  cluster_annotation = obj$cluster_annotation,
  cluster_purity = obj$cluster_purity,
  cell_vote = obj$cell_vote,
  cell_z_delta = round(obj$cell_z_delta, 4),
  annotation_low_conf = obj$annotation_low_conf,
  percent_mt = obj$percent.mt,
  globin_frac = round(obj$globin_frac, 4),
  stringsAsFactors = FALSE
)
fwrite(cell_anno, file.path(OUT_TMP, "cell_annotation.csv"))

cl_tab <- data.frame(
  cluster = sort(unique(cell_anno$cluster)),
  annotation = cl_anno[sort(unique(clusters))],
  purity = unname(cl_purity[sort(unique(clusters))]),
  n_cells = as.numeric(table(clusters)[sort(unique(clusters))]),
  median_mt = unname(cl_mt[sort(unique(clusters))]),
  max_median_z = unname(cl_zmax[sort(unique(clusters))]),
  stringsAsFactors = FALSE
)
fwrite(cl_tab, file.path(OUT_TMP, "cluster_annotation.csv"))

params <- list(
  step = "03_Cluster_Annotate",
  clustering = list(algorithm = "Leiden (igraph, R)", resolution = LEIDEN_RES,
                    n_neighbors = N_NEIGHBORS, n_pcs = N_PCS,
                    reduction = "harmony", random_seed = RANDOM_SEED,
                    note = "cluster IDs are Leiden-resolution-specific; annotation is marker-driven"),
  annotation = list(
    method = "per-cell z-score across 12 marker means -> argmax vote (final annotation)",
    marker_source = paste0(
      "Immune panels = Zhang 2021 (Kim paper lists immune cell type names only, no marker genes); ",
      "Stromal = Zhang 2021 base + Kim 2022 additions (Fibroblasts += MMP2/PDGFRA/MYL9; Endothelial += PLVAP/KDR/PTPRB); ",
      "Gastric epithelial = Kim 2022 genes (Enteroendocrine CHGA/GAST/PROX1; Chief PGC) + ",
      "Parietal GIF/ATP4B/ATP4A (no literature source, data-discovered: per-cluster DEG of Cls 15 ",
      "misvoted Fibroblasts at 6.3% purity under Zhang 9-type annotation; Kim has no parietal type) + ",
      "generic gastric mucosa markers (Epithelial panel, LIPF/PGA3/PGA4, GHRL/SST/PCSK1N). ",
      "Adapted to hg19 ref: T += CD7/LCK/IL32/CCL5 (TRAC/TRBC1 absent from hg19 ref); ",
      "Plasma = MZB1+SDC1 only (IGHG1/IGKC absent from hg19 ref; JCHAIN present as alias IGJ but ",
      "85% detected = ambient RNA; TNFRSF17/DERL3/FKBP11 too global); Fibroblasts -= FAP (dead); ",
      "Endothelial -= PECAM1 (absent from hg19 ref)"),
    cd3d_rescue = "NK-panel vote + CD3D detected -> T_cells (T-lineage rescue)",
    cluster_mode = "reference column only (cluster_annotation + purity)",
    low_confidence_z_delta = LOW_CONF_ZDELTA),
  exclusion = list(
    rbc = list(globin_genes = GLOBIN_GENES, globin_frac_max = GLOBIN_FRAC_MAX,
               note = "applied BEFORE clustering (cluster IDs self-consistent); post-RBC per-cluster globin diagnostics in log (section 6b)"),
    dead_clusters = list(median_mt_gt = DEAD_MT_MEDIAN, max_median_z_lt = DEAD_Z_MAX,
                         clusters = dead_cl),
    scrublet = "applied upstream in 02 (per-sample, threshold 0.25)"
  ),
  output = list(n_cells_after = ncol(obj), n_clusters = length(unique(Idents(obj))),
                n_excluded = sum(exclude))
)
writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(OUT_TMP, "cluster_annotate_params.json"))
log(sprintf("  Saved [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

log("Copying outputs to repo dir ...")
for (f in list.files(OUT_TMP)) {
  file.copy(file.path(OUT_TMP, f), file.path(OUT_DIR, f), overwrite = TRUE)
}

log(sprintf("\n%s", paste(rep("=", 60), collapse = "")))
log("Cluster + annotate complete.")
log(sprintf("  %d cells, %d clusters, %d types, %.1f%% low-confidence",
            ncol(obj), length(unique(Idents(obj))),
            length(unique(obj$annotation)),
            sum(obj$annotation_low_conf) / ncol(obj) * 100))
log(sprintf("  Type counts (per-cell vote):"))
for (ct in sort(names(table(obj$annotation)))) {
  log(sprintf("    %-16s %d", ct, sum(obj$annotation == ct)))
}
close(log_fh)
# 在收尾摘要写入之后再拷贝日志——存档副本才包含 "complete" +
# 类型计数段（此前拷贝先执行，存档日志永远停在 "Copying ..."）。
file.copy(file.path(TMP_DIR, "cluster_annotate.log"),
          file.path(OUT_DIR, "cluster_annotate.log"), overwrite = TRUE)
