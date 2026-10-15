#!/usr/bin/env Rscript
# ============================================================
# 六 Run 敏感性分析——星形设计执行器（种子目录由第二参数指定，
# 主/测试种子都经本脚本跑：主种子跑全 6 Run，测试种子只跑 run1 + run6）
#
# 输入：03_Cluster_Annotate/data/gse150290_epi_counts.rds（上皮谱系极简
#       子集 RAW counts — 未归一化；step2_iterative_classify 内部自行
#       log1p(CP10K) 归一化；由 03 的 export_epi_counts.R 导出）
#       03_Cluster_Annotate/data/cell_annotation.csv
#       04_Bulk_DEG/data/top50_{tumor,normal}_limma.txt（主种子）
#       或 testseed/ 下同名测试种子——由第二参数指定哪个）
# 产出：05_Iterative_Classify/data/
#   run{1..6}_*/（每 Run 四件套：labels/iteration_log/
#     final_signatures + run_params.json；
#     执行日志 runs.log 只留 /tmp，不拷入 data）
#
# Run 矩阵——星形设计：锚定 run1，run2–6 每个只变一个条件，
# 任何 Run 与 run1 的差异都可单独归因（对照不混淆）：
#
#   Run  上皮口径                       基因选择     管家黑名单  低置信排除  对照轴
#   ───  ────────────────────────────  ──────────  ──────────  ──────────  ──────────
#   1    Epithelial+Chief              logFC       无          否          论文迭代聚类流程
#   2    4 类型（Ep+Chief+EE+Parietal） logFC       无          否          上皮口径↑（广义）
#   3    Epithelial+Chief              expression  无          否          基因选择轴
#   4    Epithelial+Chief              logFC       有          否          管家黑名单轴
#   5    仅 Epithelial                 logFC       无          否          上皮口径↓（窄）
#   6    Epithelial+Chief              logFC       无          是          低置信排除轴
#
# Usage (WSL):
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 05_Iterative_Classify/scripts/sensitivity/run_iterative_runs.R
#   Rscript 05_Iterative_Classify/scripts/sensitivity/run_iterative_runs.R run6
#   （第一参数：单 Run 名或 "all"，缺省 = 全扫描）
#   Rscript 05_Iterative_Classify/scripts/sensitivity/run_iterative_runs.R all <种子目录>
#   （第二参数：种子目录，缺省 = 04_Bulk_DEG/data；测试种子传
#     04_Bulk_DEG/data/testseed，只跑 run1 + run6，输出加 testseed_ 前缀）
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

# ── 路径 ─────────────────────────────────────────────────────
REPO    <- "."
EPI_RDS <- file.path(REPO, "03_Cluster_Annotate/data/gse150290_epi_counts.rds")
ANNO_CSV<- file.path(REPO, "03_Cluster_Annotate/data/cell_annotation.csv")
DEG_DIR <- file.path(REPO, "04_Bulk_DEG/data")
OUT_DIR <- file.path(REPO, "05_Iterative_Classify/data")
TMP_DIR <- "/tmp/itercls_runs"
# 每次全扫描前清空临时目录——避免上次中断/旧配置残留的 Run 目录混入
if (dir.exists(TMP_DIR)) unlink(TMP_DIR, recursive = TRUE)
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 参数 ────────────────────────────────────────────────
N_TOP        <- 50
CONVERGENCE  <- 0.90    # issue 原文
MAX_ITER     <- 10
PADJ         <- 0.01
SEED         <- 42      # 固定随机种子：K-means 起点 + AddModuleScore 对照抽样可复现
set.seed(SEED)          # 全扫描共享一条随机流：同一调用模式重跑逐位一致；
                        # 单 Run 模式与全扫描模式的同一 Run 不保证相同（流位置不同）

# 12 个管家基因（use_blacklist=TRUE 时额外排除，即星形设计的 run4 黑名单轴；
# MT-/核糖体/MALAT1/globin 卫生过滤始终开启）
HK_GENES <- c("B2M", "EEF1A1", "TMSB4X", "TMSB10", "ACTB", "ACTG1",
              "GAPDH", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA8")

log_fh <- file(file.path(TMP_DIR, "runs.log"), open = "a")
log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  writeLines(line, log_fh)
  flush(log_fh)
}

# ══════════════════════════════════════════════════════════════
# 函数：step2_iterative_classify（05 的参数化副本；
# 新增 gene_selection + use_blacklist；卫生过滤始终开启）
# ══════════════════════════════════════════════════════════════
step2_iterative_classify <- function(counts, cell_meta, tumor_genes,
                                     normal_genes, n_top = 50,
                                     convergence = 0.90, max_iter = 10,
                                     gene_selection = "logFC",
                                     use_blacklist = FALSE) {

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

    # 卫生过滤（始终）：MT-/核糖体/MALAT1/globin
    hygiene <- grepl("^MT-|^MTRNR|^RPL|^RPS|^MRPL|^MRPS|^MALAT1$|^HBA[12]$|^HBB$",
                     genes)
    # 可选：12 个管家基因（早期黑名单构建，超出卫生过滤）
    if (use_blacklist) hygiene <- hygiene | genes %in% HK_GENES

    sig <- data.frame(gene = genes, logFC = lfc, p_value = pvals,
                      p_adj = padjs, stringsAsFactors = FALSE)
    # 分组均值 profile（基因 × 患者 → 各组基因均值）
    sig$mean_mal  <- rowMeans(mal_prof[genes, , drop = FALSE])
    sig$mean_non  <- rowMeans(non_prof[genes, , drop = FALSE])
    sig <- sig[!hygiene & sig$p_adj < padj, , drop = FALSE]
    if (nrow(sig) == 0) {
      sig <- data.frame(gene = genes, logFC = lfc, p_value = pvals,
                        p_adj = padjs, stringsAsFactors = FALSE)
      sig$mean_mal <- rowMeans(mal_prof[genes, , drop = FALSE])
      sig$mean_non <- rowMeans(non_prof[genes, , drop = FALSE])
      sig <- sig[!hygiene & sig$p_adj < 0.05, , drop = FALSE]
    }

    up_mal <- sig$gene[sig$logFC > 0]
    up_non <- sig$gene[sig$logFC < 0]
    # top-50 排序方式论文未精确规定，两种都跑：
    # logFC = 差异幅度；expression = 目标组平均表达
    if (gene_selection == "expression") {
      mal_genes <- head(up_mal[order(-sig$mean_mal[sig$logFC > 0])], n_top)
      non_genes <- head(up_non[order(-sig$mean_non[sig$logFC < 0])], n_top)
    } else {
      mal_genes <- head(up_mal[order(-sig$logFC[sig$logFC > 0])], n_top)
      non_genes <- head(up_non[order(sig$logFC[sig$logFC < 0])], n_top)
    }
    list(malig = mal_genes, nonmalig = non_genes, deg = sig)
  }

  cells <- intersect(colnames(counts), cell_meta$cell)
  counts <- counts[, cells, drop = FALSE]
  patients <- setNames(cell_meta$patient[match(cells, cell_meta$cell)], cells)

  lib <- colSums(counts)
  # CP10K：log1p(x * 10K / lib) — 社区默认 scale factor
  expr_norm <- counts %*% Diagonal(x = 10000 / lib)
  expr_norm@x <- log1p(expr_norm@x)
  colnames(expr_norm) <- colnames(counts)  # %*% Diagonal 会丢 colnames（右操作数无 dimnames）

  mal_genes <- intersect(tumor_genes, rownames(counts))
  non_genes <- intersect(normal_genes, rownames(counts))

  prev_labels <- NULL
  iteration_log <- list()
  labels <- NULL; scores <- NULL

  for (iter in seq_len(max_iter)) {
    log(sprintf("    Iter %d: %d mal / %d non genes", iter,
                length(mal_genes), length(non_genes)))

    mal_score <- add_module_score(expr_norm, mal_genes)
    non_score <- add_module_score(expr_norm, non_genes)

    scores_scaled <- scale(cbind(mal_score, non_score))
    km <- kmeans(scores_scaled, centers = 2, nstart = 25)
    mal_mean <- tapply(mal_score, km$cluster, mean)
    mal_cl <- which.max(mal_mean)
    labels <- setNames(ifelse(km$cluster == mal_cl, "malignant", "non_malignant"),
                       colnames(expr_norm))   # kmeans() 的 cluster 无名字——显式命名
    n_mal <- sum(labels == "malignant")

    agreement <- if (is.null(prev_labels)) NA else mean(labels == prev_labels)
    converged <- !is.na(agreement) && agreement >= convergence
    log(sprintf("    K-means: %d malignant | agreement %s",
                n_mal, if (is.na(agreement)) "NA" else sprintf("%.4f", agreement)))

    iteration_log[[iter]] <- list(
      iteration = iter,
      n_malig_genes = length(mal_genes),
      n_nonmalig_genes = length(non_genes),
      n_malignant = n_mal,
      n_non_malignant = length(labels) - n_mal,
      agreement = agreement,
      converged = converged
    )

    if (converged) break
    prev_labels <- labels

    deg <- tryCatch(pseudobulk_wilcoxon(counts, labels, patients, n_top = n_top),
                    error = function(e) { log(sprintf("    DEG failed: %s", e$message)); NULL })
    if (is.null(deg) || length(deg$malig) == 0 || length(deg$nonmalig) == 0) break
    mal_genes <- deg$malig
    non_genes <- deg$nonmalig
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
# 主流程：加载一次，横扫 6 个 Run
# ══════════════════════════════════════════════════════════════
# 命令行参数解析（种子目录必须在加载种子之前解析）：
#   [2] 种子目录（可选，缺省 = 04_Bulk_DEG/data；传其他种子目录
#       做种子敏感性对比时，输出目录/文件加 "testseed_" 前缀）
args <- commandArgs(trailingOnly = TRUE)
seed_dir <- DEG_DIR
run_tag <- ""
if (length(args) >= 2 && nzchar(args[2])) {
  # 相对路径按仓库根解析（与文档示例一致）；
  # 绝对路径原样使用
  seed_dir <- if (startsWith(args[2], "/")) args[2] else
    file.path(REPO, args[2])
  run_tag <- "testseed_"
}

log(sprintf("Loading counts from %s ...", EPI_RDS))
epi_in <- readRDS(EPI_RDS)          # list(counts, meta)：03 export_epi_counts 导出的上皮谱系极简子集
counts <- epi_in$counts
rownames(counts) <- make.unique(rownames(counts))   # 幂等；与 02 相同
anno <- fread(ANNO_CSV)
log(sprintf("  %d cells, annotations: %s", nrow(anno),
            paste(sort(unique(anno$annotation)), collapse = ", ")))

tumor_seed <- readLines(file.path(seed_dir, "top50_tumor_limma.txt"))
normal_seed <- readLines(file.path(seed_dir, "top50_normal_limma.txt"))
log(sprintf("Seeds: %d tumor / %d normal (limma)", length(tumor_seed), length(normal_seed)))

RUNS <- list(
  run1 = list(epi = c("Epithelial", "Chief"),         sel = "logFC",      bl = FALSE,
              dir = "run1_epi_chief_logFC"),           # 锚 = 主流程（逐位一致复现论文）
  run2 = list(epi = c("Epithelial", "Chief", "Enteroendocrine", "Parietal"),
              sel = "logFC", bl = FALSE,
              dir = "run2_epi_4type_logFC"),           # 上皮口径↑：广义（vs run1）
  run3 = list(epi = c("Epithelial", "Chief"),         sel = "expression", bl = FALSE,
              dir = "run3_epi_chief_expr"),            # 基因选择轴（vs run1）
  run4 = list(epi = c("Epithelial", "Chief"),         sel = "logFC",      bl = TRUE,
              dir = "run4_epi_chief_logFC_blacklist"), # 管家黑名单轴（vs run1；机制剖析见 docs/blacklist_ab_analysis.md）
  run5 = list(epi = c("Epithelial"),                  sel = "logFC",      bl = FALSE,
              dir = "run5_epi_only_logFC"),            # 上皮口径↓：窄（vs run1）
  run6 = list(epi = c("Epithelial", "Chief"),
              sel = "logFC", bl = FALSE, low_conf_excl = TRUE,
              dir = "run6_epi_chief_logFC_lowconf_excl")  # 低置信排除轴（vs run1）
)

# 命令行参数解析（续）：[1] 单 Run 名或 "all"（可选，缺省 = 全扫描）
if (length(args) >= 1 && nzchar(args[1]) && args[1] != "all") {
  if (!args[1] %in% names(RUNS))
    stop(sprintf("unknown run '%s' (available: %s)", args[1],
                 paste(names(RUNS), collapse = ", ")))
  RUNS <- RUNS[args[1]]
}

for (rn in names(RUNS)) {
  cfg <- RUNS[[rn]]
  out_run <- file.path(TMP_DIR, paste0(run_tag, cfg$dir))
  dir.create(out_run, recursive = TRUE, showWarnings = FALSE)
  log(sprintf("===== %s: epi=%s, selection=%s, HK blacklist=%s =====",
              rn, paste(cfg$epi, collapse = "+"), cfg$sel, cfg$bl))

  epi_cells <- anno$barcode[anno$annotation %in% cfg$epi]
  if (isTRUE(cfg$low_conf_excl)) {
    n_before <- length(epi_cells)
    epi_cells <- intersect(epi_cells, anno$barcode[!anno$annotation_low_conf])
    log(sprintf("  Epithelial cells: %d (excluded %d low-confidence)",
                length(epi_cells), n_before - length(epi_cells)))
  } else {
    log(sprintf("  Epithelial cells: %d", length(epi_cells)))
  }
  cnt <- counts[, epi_cells, drop = FALSE]
  meta <- data.frame(cell = epi_cells,
                     patient = anno$patient[match(epi_cells, anno$barcode)],
                     stringsAsFactors = FALSE)

  t0 <- Sys.time()
  res <- step2_iterative_classify(cnt, meta, tumor_seed, normal_seed,
                                  n_top = N_TOP, convergence = CONVERGENCE,
                                  max_iter = MAX_ITER,
                                  gene_selection = cfg$sel,
                                  use_blacklist = cfg$bl)
  log(sprintf("  Done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

  n_mal <- sum(res$labels == "malignant")
  n_cells <- length(res$labels)
  log(sprintf("  RESULT: %d malignant / %d non-malignant (%.1f%% malignant)",
              n_mal, n_cells - n_mal, n_mal / n_cells * 100))
  log(sprintf("  Final malignant signature (top 8): %s",
              paste(head(res$final_signatures$malignant, 8), collapse = ", ")))

  fwrite(data.frame(cell = names(res$labels), label = unname(res$labels),
                    stringsAsFactors = FALSE),
         file.path(out_run, "labels.csv"))

  ilog <- rbindlist(lapply(res$iteration_log, as.data.table), fill = TRUE)
  fwrite(ilog, file.path(out_run, "iteration_log.csv"))

  final_df <- data.frame(
    malignant_genes = c(res$final_signatures$malignant,
                        rep(NA, max(0, length(res$final_signatures$non_malignant) -
                                      length(res$final_signatures$malignant)))),
    non_malignant_genes = c(res$final_signatures$non_malignant,
                            rep(NA, max(0, length(res$final_signatures$malignant) -
                                          length(res$final_signatures$non_malignant)))),
    stringsAsFactors = FALSE
  )
  fwrite(final_df, file.path(out_run, "final_signatures.csv"), na = "")

  params <- list(
    run = rn,
    epithelial_definition = cfg$epi,
    gene_selection = cfg$sel,
    hk_blacklist = cfg$bl,
    hygiene = "MT-/MTRNR/RPL/RPS/MRPL/MRPS/MALAT1/globin (always on)",
    seeds = seed_dir,
    seed = SEED,
    convergence = CONVERGENCE, max_iter = MAX_ITER, n_top = N_TOP,
    n_epi_cells = n_cells, n_malignant = n_mal,
    low_conf_excluded = isTRUE(cfg$low_conf_excl),
    pct_malignant = round(n_mal / n_cells * 100, 1),
    n_iterations = nrow(ilog)
  )
  writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE),
             file.path(out_run, "run_params.json"))
  rm(cnt); gc()
}

# ── 全部拷回仓库目录 ────────────────────────────────────────
log("Copying runs to repo dir ...")
for (d in list.files(TMP_DIR, pattern = "run[0-9]", full.names = TRUE)) {
  out_d <- file.path(OUT_DIR, basename(d))
  dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
  for (f in list.files(d)) file.copy(file.path(d, f), file.path(out_d, f), overwrite = TRUE)
}
log("All runs complete.")
close(log_fh)
