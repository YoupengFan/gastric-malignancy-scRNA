#!/usr/bin/env Rscript
# ============================================================
# 01_QC：GSE150290 raw 10x 矩阵 → QC 过滤后的合并计数矩阵
#
# 输入：01_Data/raw/（52 个 GSM*_PatXX-{A,B} 原始 10x 矩阵；数据
#       不随仓库分发，从 GEO GSE150290 下载解压后放置于此）
#
# 输出：合并 dgCMatrix RDS + 逐细胞 metadata CSV + 逐样本统计 +
# 参数 JSON（Seurat 对象由 02_Normalize_Batch 从这些产物构建）。
#
# 下界逐字沿用 Zhang 2021（UMI>=1001 / genes>=501 / MT<=50%）。
# 与论文的出入共两处，均在下文阈值注释与 qc_params.json 留档：
#   - 上界放宽 6000 → 15000
#   - 看家基因过滤省略：论文 curated list 未随文公开，无法复现
#
# Usage (WSL)：
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 01_QC/scripts/qc_filter.R
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

# ── 路径 ─────────────────────────────────────────────────────
RAW       <- "01_Data/raw"
OUT_DIR   <- "01_QC/data"
# 中间产物与最终产物先写 /tmp 再复制到 OUT_DIR：/tmp 是 WSL 的
# ext4，写大文件远快于 NTFS 挂载，故收尾统一复制一次
TMP_DIR   <- "/tmp/qc_r"
SAMP_TMP  <- file.path(TMP_DIR, "samples")   # 逐样本过滤矩阵 (ext4)
OUT_TMP   <- file.path(TMP_DIR, "out")       # 最终产物，收尾时复制到 OUT_DIR

# ── QC 阈值（论文下界逐字沿用 + 上界放宽）──────────────
MIN_UMIS   <- 1001
MAX_UMIS   <- 15000   # 从论文 6000 放宽
MIN_GENES  <- 501
PCT_MT_MAX <- 50      # Zhang 2021 原文："over 50% UMIs ... mitochondrial genome"
# 看家基因过滤：省略（论文未公开 curated list）

# 目录三件套：逐样本中间产物 + 最终产物暂存（均在 ext4 上）+ 正式输出
dir.create(SAMP_TMP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TMP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# log()：一行日志同时写终端与文件。flush 保证崩溃时已写行不丢失；
# 文件 append 模式，多次运行的历史都保留在同一份 qc_filter.log
log_fh <- file(file.path(TMP_DIR, "qc_filter.log"), open = "a")
log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  writeLines(line, log_fh)
  flush(log_fh)
}

# ── 定位 hg19 矩阵目录（先试平铺路径，再递归兜底）────
find_hg19 <- function(sample_dir) {
  flat <- file.path(sample_dir, "raw_gene_bc_matrices", "hg19")
  req <- c("matrix.mtx", "genes.tsv", "barcodes.tsv")
  if (all(file.exists(file.path(flat, req)))) return(flat)
  dirs <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
  for (d in dirs) {
    if (basename(d) == "hg19" && all(file.exists(file.path(d, req)))) return(d)
  }
  return(NULL)
}

# ── 发现样本 ──────────────────────────────────────────────
sample_dirs <- sort(list.dirs(RAW, recursive = FALSE, full.names = TRUE))
sample_dirs <- sample_dirs[grepl("^GSM", basename(sample_dirs))]
sample_names <- basename(sample_dirs)
log(sprintf("Found %d sample directories in %s", length(sample_names), RAW))

parse_sample <- function(name) {
  # GSM4546300_Pat01-A -> patient=Pat01, tissue=A
  m <- regmatches(name, regexec("^GSM\\d+_(Pat\\d+)-([AB])$", name))[[1]]
  if (length(m) == 3) return(list(patient = m[2], tissue = m[3]))
  list(patient = name, tissue = "?")
}

# ── 单样本处理 ──────────────────────────────────────────────
# 读入一个 raw 10x hg19 矩阵、施加 QC 阈值，返回过滤后的 dgCMatrix +
# 逐细胞 QC 向量。样本缺矩阵目录 / 行列数不齐 / 无细胞通过时返回
# ok=FALSE（并记录原因）——调用方统计 raw barcode 的口径与旧版一致。
process_sample <- function(sname, sdir) {
  mtx_dir <- find_hg19(sdir)
  if (is.null(mtx_dir)) {
    log(sprintf("  SKIP %s: no hg19 matrix dir found", sname))
    return(list(ok = FALSE, raw = 0L, msg = "no hg19 matrix dir"))
  }

  # genes.tsv：10x v2，两列、无表头（ENSG id \t symbol）
  genes <- fread(file.path(mtx_dir, "genes.tsv"), header = FALSE,
                 col.names = c("gene_id", "symbol"), sep = "\t",
                 stringsAsFactors = FALSE, showProgress = FALSE)
  barcodes <- fread(file.path(mtx_dir, "barcodes.tsv"), header = FALSE,
                    col.names = "barcode", sep = "\t",
                    stringsAsFactors = FALSE, showProgress = FALSE)

  mtx <- as(readMM(file.path(mtx_dir, "matrix.mtx")), "dgCMatrix")
  raw <- ncol(mtx)
  if (nrow(mtx) != nrow(genes)) {
    log(sprintf("  FAIL %s: matrix rows (%d) != genes (%d)",
                sname, nrow(mtx), nrow(genes)))
    return(list(ok = FALSE, raw = raw, msg = "matrix rows != genes.tsv"))
  }

  # QC 指标（全部按列向量化，无逐细胞循环）
  total_umis <- colSums(mtx)
  n_genes    <- diff(mtx@p)          # dgCMatrix：列非零元数 = 表达基因数
  mt_rows    <- grepl("^MT-", genes$symbol)
  pct_mt     <- colSums(mtx[mt_rows, , drop = FALSE]) / total_umis * 100

  keep <- total_umis >= MIN_UMIS & total_umis <= MAX_UMIS &
          n_genes >= MIN_GENES & pct_mt <= PCT_MT_MAX
  if (sum(keep) == 0) {
    return(list(ok = FALSE, raw = raw, msg = "0 cells pass QC"))
  }

  mtx_f <- mtx[, keep, drop = FALSE]
  colnames(mtx_f) <- paste0(sname, "-", barcodes$barcode[keep])
  rownames(mtx_f) <- genes$symbol

  list(ok = TRUE, raw = raw, genes = genes$symbol, mtx = mtx_f,
       n_pass = sum(keep),
       total = total_umis[keep], n_genes = n_genes[keep], pct_mt = pct_mt[keep])
}

# ── 逐样本处理 ─────────────────────────────────────────────
# 四个累积器：逐样本统计行 / 过滤矩阵（末尾 cbind 合并用）/
# 逐细胞 meta / 第一样本的基因顺序（全样本一致性参照）
stats_rows <- list()
mats <- list()
metas <- list()
global_genes <- NULL
raw_barcodes_total <- 0

for (i in seq_along(sample_dirs)) {
  sname <- sample_names[i]
  t0 <- Sys.time()
  log(sprintf("[%d/%d] %s ...", i, length(sample_dirs), sname))

  res <- process_sample(sname, sample_dirs[i])
  raw_barcodes_total <- raw_barcodes_total + res$raw
  if (!res$ok) {
    log(sprintf("  -> %s, skipping sample [%.0fs]",
                res$msg, difftime(Sys.time(), t0, units = "secs")))
    next
  }

  # genes.tsv 一致性：52 个矩阵必须共享同一基因顺序才能安全 cbind
  # （全向量检查，而非只比前几行）
  if (is.null(global_genes)) {
    global_genes <- res$genes
  } else if (!identical(res$genes, global_genes)) {
    log(sprintf("  WARN %s: genes.tsv differs from sample 1 (rows=%d vs %d)",
                sname, length(res$genes), length(global_genes)))
  }

  # 逐样本 RDS 写 ext4 = 崩溃检查点（未实现断点续跑；全量约 15 分钟，
  # 重启代价不高，故仅作保险）
  saveRDS(res$mtx, file.path(SAMP_TMP, paste0(sname, ".rds")))
  mats[[sname]] <- res$mtx

  # 逐细胞 metadata：02 环节细胞注释的基础——QC 指标随细胞带走，
  # 下游任何过滤/诊断都能回溯到这里的原始 QC 值
  info <- parse_sample(sname)
  metas[[sname]] <- data.frame(
    barcode = colnames(res$mtx),
    sample  = sname,
    patient = info$patient,
    tissue  = info$tissue,
    n_umis  = res$total,
    n_genes = res$n_genes,
    pct_mt  = round(res$pct_mt, 3),
    stringsAsFactors = FALSE
  )

  # 逐样本统计：QC 报告表格的数据源（通过率 + 三个中位数）
  stats_rows[[sname]] <- data.frame(
    sample = sname, patient = info$patient, tissue = info$tissue,
    raw_barcodes = res$raw, n_pass = res$n_pass,
    pct_pass = round(res$n_pass / res$raw * 100, 2),
    median_genes = round(median(res$n_genes), 1),
    median_counts = round(median(res$total), 0),
    median_pct_mt = round(median(res$pct_mt), 2),
    stringsAsFactors = FALSE
  )

  log(sprintf("  -> %6d cells passed (%.1f%% of %d raw) [%.0fs]",
              res$n_pass, res$n_pass / res$raw * 100, res$raw,
              difftime(Sys.time(), t0, units = "secs")))
}

log(sprintf("\nProcessed %d samples, raw barcodes total: %d", length(mats), raw_barcodes_total))
if (length(mats) == 0) {
  log("ERROR: No samples passed QC. Aborting.")
  quit(status = 1)
}

# ── 跨样本合并 ─────────────────────────────────────────────
# cbind 的前提是 52 个矩阵基因顺序一致（上方循环已用 identical()
# 全向量校验，不一致会 WARN）；barcode 在 process_sample 里已加
# 样本前缀，合并后逐细胞身份不冲突
t0 <- Sys.time()
log("Merging sample matrices (cbind) ...")
merged <- do.call(cbind, mats)
log(sprintf("  merged: %d genes x %d cells [%.0fs]",
            nrow(merged), ncol(merged), difftime(Sys.time(), t0, units = "secs")))

meta <- rbindlist(metas)

stats_df <- rbindlist(stats_rows)

# ══ 导出：5 个产物 ═══════════════════════════════════════════
# 各产物服务不同下游（分工表见 01_QC/README.md）：
#   ① 合并矩阵 RDS    —— 02_Normalize_Batch 的输入（全部表达数据）
#   ② 逐细胞 meta CSV —— 02 环节细胞注释（QC 指标随细胞走）
#   ③ 逐样本 stats CSV — QC 报告表格数据源（通过率/中位数）
#   ④ 参数 JSON       —— 留档：阈值 + 与论文出入 + 全局汇总
#   ⑤ 运行日志        —— 完整运行记录（append，含末尾汇总）
# ①–④ 先写 /tmp (ext4) 再统一复制，原因见文件头路径注释。
t0 <- Sys.time()
log("Saving merged matrix (RDS) ...")
# compress=FALSE：02 环节要整体读入这个大矩阵，不压缩换取更快读入；
# 且无压缩的序列化结果确定，便于逐位复现校验
saveRDS(merged, file.path(OUT_TMP, "gse150290_qc_matrix.rds"), compress = FALSE)
log("Saving cell metadata (CSV) ...")
fwrite(meta, file.path(OUT_TMP, "gse150290_qc_meta.csv"))
log("Saving per-sample stats (CSV) ...")
fwrite(stats_df, file.path(OUT_TMP, "qc_stats.csv"))

# qc_params.json 四部分：
#   thresholds        → 本次实际阈值（机器可读，复现核对用）
#   paper_comparison  → 论文原值 vs 我们的改动
#   output            → 输出规模（细胞/基因/样本数）
#   summary_stats     → 质量摘要（末尾日志汇总的 JSON 版）
params <- list(
  step = "01_QC",
  pipeline = "R replication of Zhang 2021 Gut (issue #26)",
  thresholds = list(
    min_umis = MIN_UMIS, max_umis = MAX_UMIS,
    min_genes = MIN_GENES, pct_mt_max = PCT_MT_MAX,
    housekeeping_filter = "OMITTED: paper's curated list not published",
    paper_comparison = "Zhang 2021: 1001-6000 UMIs / >=501 genes / MT<50% / housekeeping mean<2"
  ),
  output = list(
    total_cells = ncol(merged), total_genes = nrow(merged),
    total_raw_barcodes = raw_barcodes_total, n_samples = length(mats)
  ),
  summary_stats = list(
    median_genes_per_cell = round(median(stats_df$median_genes), 1),
    median_counts_per_cell = round(median(stats_df$median_counts), 0),
    median_pct_mt = round(median(stats_df$median_pct_mt), 2),
    cells_per_sample_range = sprintf("%d-%d", min(stats_df$n_pass), max(stats_df$n_pass))
  )
)
writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(OUT_TMP, "qc_params.json"))
log(sprintf("Saved /tmp products [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# 收尾：ext4 → NTFS 统一复制一次。四个数据产物 + 运行日志
# （日志带回 Windows 侧供查阅，append 保留历史运行记录）
t0 <- Sys.time()
log("Copying outputs to repo dir ...")
for (f in c("gse150290_qc_matrix.rds", "gse150290_qc_meta.csv",
            "qc_stats.csv", "qc_params.json")) {
  file.copy(file.path(OUT_TMP, f), file.path(OUT_DIR, f), overwrite = TRUE)
}
file.copy(file.path(TMP_DIR, "qc_filter.log"), file.path(OUT_DIR, "qc_filter.log"), overwrite = TRUE)
log(sprintf("Copy done [%.0fs]", difftime(Sys.time(), t0, units = "secs")))

# ── 末尾汇总：日志的"封面页" ───────────────────────────────
# 任务跑完后日志最后一段就是本次 run 的摘要——阈值/输入/产出/得率，
# tail 一下即可汇报，不用翻 52 行逐样本记录
log(sprintf("\n%s", paste(rep("=", 60), collapse = "")))
log("QC Complete.")
log(sprintf("  Thresholds: UMI [%d, %d], genes>=%d, MT%%<=%d (no housekeeping filter)",
            MIN_UMIS, MAX_UMIS, MIN_GENES, PCT_MT_MAX))
log(sprintf("  Input:  %d samples, %d raw barcodes", length(sample_dirs), raw_barcodes_total))
log(sprintf("  Output: %d cells, %d genes", ncol(merged), nrow(merged)))
log(sprintf("  Yield:  %.2f%%", ncol(merged) / raw_barcodes_total * 100))
log(sprintf("  Files:  %s", file.path(OUT_DIR, "gse150290_qc_matrix.rds")))
close(log_fh)
