#!/usr/bin/env Rscript
# ============================================================
# 05_Iterative_Classify: step2_iterative_classify — Zhang 2021 Gut
# 迭代 K-means 恶性分类的忠实复现（主流程 = run1 配置）
#
# 输入：03_Cluster_Annotate/data/gse150290_epi_counts.rds（上皮谱系极简
#       子集：RAW counts + annotation + patient；从全量注释对象导出）
#       04_Bulk_DEG/data/top50_{tumor,normal}_limma.txt（种子基因）
# 产出：05_Iterative_Classify/data/（labels/scores/iteration_log/
#       final_signatures + iterative_params.json）
#       + 05_Iterative_Classify/figures/
#
# Usage (WSL):
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 05_Iterative_Classify/scripts/step2_iterative_classify.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

# ══════════════════════════════════════════════════════════════
# 函数：step2_iterative_classify —— issue 要求的"第 2-3 步 = 一个函数"
# ══════════════════════════════════════════════════════════════
# counts       : 基因 × 细胞的 RAW UMI 稀疏计数矩阵——
#                必须是原始 counts（未归一化）。函数内部自行
#                log1p(CP10K) 归一化，不读 02 的归一化结果
# cell_meta    : data.frame(cell, patient) — cell = colnames(counts)
# tumor_genes  : 肿瘤上调种子基因（来自 step1_bulk_deg）
# normal_genes : 正常上调种子基因（来自 step1_bulk_deg）
# n_top        : 每轮 DEG 每方向基因数
# convergence  : 停步用的标签一致率阈值
# max_iter     : 迭代轮数安全上限
# 返回: list(labels, scores, iteration_log, final_signatures)
step2_iterative_classify <- function(counts, cell_meta, tumor_genes,
                                     normal_genes, n_top = 50,
                                     convergence = 0.90, max_iter = 10) {

  # ── 内部函数：AddModuleScore（Puram 2017：25 bins × 100 同 bin 对照）──
  # expr: 基因 × 细胞（稀疏或稠密）；得分为各基因
  # (expr_g − 100 个随机同 bin 对照均值) 的均值
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
      # 同 bin 对照的逐细胞均值（colMeans：每细胞一个值；
      # rowMeans 会给每基因一个值并按列静默循环复用）
      ctrl_mean <- colMeans(expr[ctrl_idx, , drop = FALSE])
      scores <- scores + as.numeric(gene_expr) - ctrl_mean
    }
    scores / max(length(valid), 1)
  }

  # ── 内部函数：pseudobulk 配对 Wilcoxon DEG ────────────────
  # 按 label 构建 per-patient 均值 profile，再逐基因做配对 Wilcoxon
  pseudobulk_wilcoxon <- function(counts, labels, patients, n_top, padj = 0.01) {
    mal_cells  <- names(labels)[labels == "malignant"]
    non_cells  <- names(labels)[labels == "non_malignant"]
    if (length(mal_cells) < 3 || length(non_cells) < 3)
      stop("groups too small for DEG")

    pt_mal  <- split(mal_cells, patients[mal_cells])
    pt_non  <- split(non_cells, patients[non_cells])
    shared  <- intersect(names(pt_mal), names(pt_non))
    if (length(shared) < 3) stop("fewer than 3 paired patients")

    mal_prof  <- sapply(shared, function(p) rowMeans(counts[, pt_mal[[p]], drop = FALSE]))
    non_prof  <- sapply(shared, function(p) rowMeans(counts[, pt_non[[p]], drop = FALSE]))

    # 基因过滤：至少一组中 ≥3 个患者表达（均值 > 0）
    ok <- rowSums(mal_prof > 0) >= 3 | rowSums(non_prof > 0) >= 3
    genes <- rownames(counts)[ok]

    pvals <- numeric(length(genes)); lfc <- numeric(length(genes))
    for (i in seq_along(genes)) {
      g <- genes[i]
      w <- tryCatch(suppressWarnings(
        wilcox.test(mal_prof[g, ], non_prof[g, ], paired = TRUE, exact = FALSE)),
        error = function(e) NULL)
      pvals[i] <- if (is.null(w)) 1 else w$p.value
      lfc[i] <- mean(mal_prof[g, ]) - mean(non_prof[g, ])
    }
    padjs <- p.adjust(pvals, method = "BH")

    # DEG 过滤：MT-/核糖体/MALAT1/globin
    hygiene <- grepl("^MT-|^MTRNR|^RPL|^RPS|^MRPL|^MRPS|^MALAT1$|^HBA[12]$|^HBB$", genes)

    sig <- data.frame(gene = genes, logFC = lfc, p_value = pvals,
                      p_adj = padjs, stringsAsFactors = FALSE)
    sig <- sig[!hygiene & sig$p_adj < padj, , drop = FALSE]
    if (nrow(sig) == 0) {
      sig <- data.frame(gene = genes, logFC = lfc, p_value = pvals,
                        p_adj = padjs, stringsAsFactors = FALSE)
      sig <- sig[!hygiene & sig$p_adj < 0.05, , drop = FALSE]
    }
    mal_genes  <- head(sig$gene[sig$logFC > 0][order(-sig$logFC[sig$logFC > 0])], n_top)
    non_genes  <- head(sig$gene[sig$logFC < 0][order(sig$logFC[sig$logFC < 0])], n_top)
    list(malig = mal_genes, nonmalig = non_genes, deg = sig)
  }

  # ── 准备 ───────────────────────────────────────────────────
  cells <- intersect(colnames(counts), cell_meta$cell)
  counts <- counts[, cells, drop = FALSE]
  patients <- setNames(cell_meta$patient[match(cells, cell_meta$cell)], cells)

  lib <- colSums(counts)
  # CP10K：log1p(x * 10K / lib)——右乘对角矩阵实现逐列缩放
  expr_norm <- counts %*% Diagonal(x = 10000 / lib)
  expr_norm@x <- log1p(expr_norm@x)
  # %*% 的 dimnames 规则：乘积行名取左操作数、列名取右操作数；
  # Diagonal(x=...) 本身无 dimnames → 乘积列名（细胞 barcode）丢失，
  # 行名（基因）保留。列名是下游 labels/scores 的细胞身份键，必须补回
  colnames(expr_norm) <- colnames(counts)

  mal_genes <- intersect(tumor_genes, rownames(counts))
  non_genes <- intersect(normal_genes, rownames(counts))
  log(sprintf("Seed genes: %d tumor / %d normal (data-matched)", length(mal_genes), length(non_genes)))

  prev_labels <- NULL
  iteration_log <- list()
  labels <- NULL; scores <- NULL

  for (iter in seq_len(max_iter)) {
    log(sprintf("===== Iteration %d ===== (mal %d, nonmal %d genes)",
                iter, length(mal_genes), length(non_genes)))

    # A. 打分
    mal_score <- add_module_score(expr_norm, mal_genes)
    non_score <- add_module_score(expr_norm, non_genes)

    # B. 缩放后的 2D 分数上 K-means k=2（nstart=25 组随机起点取最优）
    scores_scaled <- scale(cbind(mal_score, non_score))
    km <- kmeans(scores_scaled, centers = 2, nstart = 25)
    mal_mean <- tapply(mal_score, km$cluster, mean)
    mal_cl <- which.max(mal_mean)
    labels <- setNames(ifelse(km$cluster == mal_cl, "malignant", "non_malignant"),
                       colnames(expr_norm))   # kmeans() 的 cluster 无名字——显式命名
    n_mal <- sum(labels == "malignant")
    log(sprintf("  K-means: malignant cluster = %d -> %d malignant / %d non-malignant cells",
                mal_cl, n_mal, length(labels) - n_mal))

    # C. 收敛检查
    agreement <- if (is.null(prev_labels)) NA else mean(labels == prev_labels)
    converged <- !is.na(agreement) && agreement >= convergence
    log(sprintf("  Label agreement: %s (threshold %.2f)",
                if (is.na(agreement)) "NA" else sprintf("%.4f", agreement), convergence))

    iteration_log[[iter]] <- list(
      iteration = iter,
      n_malig_genes = length(mal_genes),
      n_nonmalig_genes = length(non_genes),
      n_malignant = n_mal,
      n_non_malignant = length(labels) - n_mal,
      agreement = agreement,
      converged = converged
    )

    if (converged) {
      log(sprintf("  >>> Converged (%.2f%% >= %.0f%%)", agreement * 100, convergence * 100))
      break
    }
    prev_labels <- labels

    # D. pseudobulk 配对 Wilcoxon DEG → 新签名
    deg <- tryCatch(pseudobulk_wilcoxon(counts, labels, patients, n_top = n_top),
                    error = function(e) { log(sprintf("  DEG failed: %s", e$message)); NULL })
    if (is.null(deg) || length(deg$malig) == 0 || length(deg$nonmalig) == 0) {
      log("  DEG empty -> stopping")
      break
    }
    mal_genes <- deg$malig
    non_genes <- deg$nonmalig
    log(sprintf("  New signatures: %d malignant / %d non-malignant genes", length(mal_genes), length(non_genes)))
  }

  list(
    labels = labels,
    scores = data.frame(cell = cells, mal_score = mal_score, nonmal_score = non_score,
                        stringsAsFactors = FALSE),
    iteration_log = iteration_log,
    final_signatures = list(malignant = mal_genes, non_malignant = non_genes)
  )
}

# ══════════════════════════════════════════════════════════════
# 主流程（R 管线的上皮谱系运行）
# ══════════════════════════════════════════════════════════════
ANNO_RDS  <- "03_Cluster_Annotate/data/gse150290_epi_counts.rds"
DEG_DIR   <- "04_Bulk_DEG/data"
OUT_DIR   <- "05_Iterative_Classify/data"
FIG_DIR   <- "05_Iterative_Classify/figures"
OUT_TMP   <- "/tmp/itercls_r"

# 上皮定义 = 论文口径：Epithelial + Chief，排除 Parietal/Enteroendocrine
EPI_TYPES <- c("Epithelial", "Chief")
N_TOP      <- 50
CONVERGENCE <- 0.90   # issue 原文
MAX_ITER   <- 10
SEED       <- 42      # 固定随机种子：K-means 起点 + AddModuleScore 对照抽样可复现
set.seed(SEED)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
if (dir.exists(OUT_TMP)) unlink(OUT_TMP, recursive = TRUE)
dir.create(OUT_TMP, recursive = TRUE, showWarnings = FALSE)
log_fh <- file(file.path(OUT_TMP, "iterative.log"), open = "a")
log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  writeLines(line, log_fh)
  flush(log_fh)
}

log(sprintf("Loading %s ...", ANNO_RDS))
epi_in <- readRDS(ANNO_RDS)   # list(counts, meta)：03 export_epi_counts 导出的上皮谱系极简子集
log(sprintf("  %d epithelial-lineage cells, annotations: %s", ncol(epi_in$counts),
            paste(sort(unique(epi_in$meta$annotation)), collapse = ", ")))

epi_cells <- epi_in$meta$cell[epi_in$meta$annotation %in% EPI_TYPES]
log(sprintf("Epithelial lineage (%s): %d cells",
            paste(EPI_TYPES, collapse = " + "), length(epi_cells)))

counts <- epi_in$counts[, epi_cells]
cell_meta <- data.frame(cell = epi_cells,
                        patient = epi_in$meta$patient[match(epi_cells, epi_in$meta$cell)],
                        stringsAsFactors = FALSE)

tumor_seed <- readLines(file.path(DEG_DIR, "top50_tumor_limma.txt"))
normal_seed <- readLines(file.path(DEG_DIR, "top50_normal_limma.txt"))
log(sprintf("Seed genes loaded: %d tumor / %d normal", length(tumor_seed), length(normal_seed)))

log(sprintf("Running step2_iterative_classify (n_top=%d, convergence=%.2f, max_iter=%d) ...",
            N_TOP, CONVERGENCE, MAX_ITER))
res <- step2_iterative_classify(counts, cell_meta, tumor_seed, normal_seed,
                                n_top = N_TOP, convergence = CONVERGENCE,
                                max_iter = MAX_ITER)

# ── 保存结果 ─────────────────────────────────────────────
fwrite(res$scores, file.path(OUT_TMP, "scores.csv"))
fwrite(data.frame(cell = names(res$labels), label = unname(res$labels),
                  stringsAsFactors = FALSE),
       file.path(OUT_TMP, "labels.csv"))

ilog <- rbindlist(lapply(res$iteration_log, as.data.table), fill = TRUE)
fwrite(ilog, file.path(OUT_TMP, "iteration_log.csv"))

final_df <- data.frame(
  malignant_genes = c(res$final_signatures$malignant,
                      rep(NA, max(0, length(res$final_signatures$non_malignant) -
                                    length(res$final_signatures$malignant)))),
  non_malignant_genes = c(res$final_signatures$non_malignant,
                          rep(NA, max(0, length(res$final_signatures$malignant) -
                                        length(res$final_signatures$non_malignant)))),
  stringsAsFactors = FALSE
)
fwrite(final_df, file.path(OUT_TMP, "final_signatures.csv"), na = "")

n_mal <- sum(res$labels == "malignant")
log(sprintf("\nFinal: %d malignant / %d non-malignant (%.1f%% malignant)",
            n_mal, length(res$labels) - n_mal, n_mal / length(res$labels) * 100))
log(sprintf("Final malignant signature (top 8): %s",
            paste(head(res$final_signatures$malignant, 8), collapse = ", ")))
log(sprintf("Final non-malignant signature (top 8): %s",
            paste(head(res$final_signatures$non_malignant, 8), collapse = ", ")))

# ── 图件 ──────────────────────────────────────────────────
log("Saving figures ...")
sc <- res$scores
sc$label <- res$labels[sc$cell]
p1 <- ggplot(sc, aes(x = mal_score, y = nonmal_score, color = label)) +
  geom_point(size = 0.4) +
  scale_color_manual(values = c(malignant = "#c51b7d", non_malignant = "#4d9221")) +
  labs(title = "Iterative classification: K-means k=2 on module scores",
       x = "Malignant module score", y = "Non-malignant module score") +
  theme_classic()
ggsave(file.path(OUT_TMP, "fig_scores_kmeans.png"), p1,
       device = ragg::agg_png, width = 7, height = 6, dpi = 300)

# 收敛轨迹：x = 迭代轮次（离散），首轮 agreement=NA 无对照
# 可比——图上从第 2 轮起，刻度取整数
if (nrow(ilog) > 1) {
  it_df <- data.frame(iter = ilog$iteration, agreement = ilog$agreement)
  it_df <- it_df[!is.na(it_df$agreement), ]
  p2 <- ggplot(it_df, aes(x = iter, y = agreement)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept = CONVERGENCE, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = it_df$iter) +
    labs(title = sprintf("Convergence trajectory (threshold %.0f%%)", CONVERGENCE * 100),
         x = "Iteration", y = "Label agreement") +
    ylim(0, 1) + theme_classic()
  ggsave(file.path(OUT_TMP, "fig_convergence.png"), p2,
         device = ragg::agg_png, width = 6, height = 4, dpi = 300)
}

params <- list(
  step = "05_Iterative_Classify",
  function_name = "step2_iterative_classify",
  method = "AddModuleScore (25 bins x 100 ctrl) -> K-means k=2 (nstart=25) -> pseudobulk paired Wilcoxon DEG -> iterate",
  seed = SEED,
  convergence = CONVERGENCE, max_iter = MAX_ITER, n_top = N_TOP,
  epithelial_definition = EPI_TYPES,
  seed_source = "04_Bulk_DEG top50_*_limma.txt (limma only)",
  output = list(n_epi_cells = length(epi_cells), n_malignant = n_mal,
                pct_malignant = round(n_mal / length(epi_cells) * 100, 1),
                n_iterations = nrow(ilog))
)
writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(OUT_TMP, "iterative_params.json"))

for (f in list.files(OUT_TMP, pattern = "\\.(csv|json)$")) {
  file.copy(file.path(OUT_TMP, f), file.path(OUT_DIR, f), overwrite = TRUE)
}
for (f in list.files(OUT_TMP, pattern = "\\.png$")) {
  file.copy(file.path(OUT_TMP, f), file.path(FIG_DIR, f), overwrite = TRUE)
}
log("Done.")
close(log_fh)
