#!/usr/bin/env Rscript
# ============================================================
# 01_QC 图 — 计算部分：过滤前 QC 聚合量
#
# 输入：01_Data/raw/（52 个 GSM*_PatXX-{A,B} 原始 10x 矩阵；数据
#       不随仓库分发，从 GEO GSE150290 下载解压后放置于此）
#
# 设计要点：
#   - 计算/绘图拆分：raw 矩阵重活在本脚本（WSL）；图由
#     qc_plot_figure.R 在 Windows 原生 R 上渲染（ragg device）
#   - 自包含：过滤前指标直接从 raw 10x 矩阵计算，不依赖其他环节
#   - 阈值线 = 本管线 QC：min 1,001 / max 15,000 / genes 501 /
#     MT 50%；论文的 6,000 UMI 上限另画一条虚线
#   - Panel C knee 图按 rank 顺序抽稀到 ~200k 点以提速绘图；
#     完整 38M rank 向量不保存（~300 MB），需要时重跑本脚本
#
# 产出 (01_QC/data/)：qc_plot_hist_total.csv, qc_plot_hist_genes.csv,
#   qc_plot_knee.csv, qc_plot_scatter.csv, qc_plot_summary.json
#
# Usage (WSL)：
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 01_QC/scripts/qc_plot_compute.R   # ~12 min
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

# ── 路径 ─────────────────────────────────────────────────────
RAW      <- "01_Data/raw"
BASE     <- "01_QC"
DATA_DIR <- file.path(BASE, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# ── QC 阈值（本管线 + 论文对照线）─────────────────────────
MIN_UMIS   <- 1001
MAX_UMIS   <- 15000
PAPER_MAX  <- 6000    # Zhang 2021 上限
MIN_GENES  <- 501
PCT_MT_MAX <- 50

# 绘图参数：直方图 x 轴上限（截断长尾超界部分）、散点面板采样数、直方图分箱数
XMAX_COUNTS <- 20000
XMAX_GENES  <- 6000
N_SCATTER   <- 50000
N_BINS      <- 300

t0 <- Sys.time()
log <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# 与 qc_filter.R / qc_doublets.py 同款；各脚本自包含，故保留副本
find_hg19 <- function(sample_dir) {
  flat <- file.path(sample_dir, "raw_gene_bc_matrices", "hg19")
  req <- c("matrix.mtx", "genes.tsv", "barcodes.tsv")
  if (all(file.exists(file.path(flat, req)))) return(flat)
  dirs <- list.dirs(sample_dir, recursive = TRUE, full.names = TRUE)
  for (d in dirs) {
    if (basename(d) == "hg19" && all(file.exists(file.path(d, req)))) return(d)
  }
  NULL
}

# ══════════════════════════════════════════════════════════════
# PASS 1：从 raw 矩阵计算过滤前指标
# ══════════════════════════════════════════════════════════════
sample_dirs <- sort(list.dirs(RAW, recursive = FALSE, full.names = TRUE))
sample_dirs <- sample_dirs[grepl("^GSM", basename(sample_dirs))]
log(sprintf("Scanning %d samples for pre-filter metrics ...", length(sample_dirs)))

# 累积器：两个直方图用固定公共分箱（breaks 跨样本一致，52 个样本
# 才能直接相加）；rank_all 存每样本 total 向量（knee 图原料）；
# scatter_rows 存散点面板的候选行
brk_a <- seq(0, XMAX_COUNTS, length.out = N_BINS + 1)
brk_b <- seq(0, XMAX_GENES,  length.out = N_BINS + 1)
hist_a <- hist_b <- rep(0, N_BINS)
rank_all <- list()
scatter_rows <- list()
raw_total <- 0

for (i in seq_along(sample_dirs)) {
  sname <- basename(sample_dirs[i])
  mtx_dir <- find_hg19(sample_dirs[i])
  if (is.null(mtx_dir)) { log(sprintf("  SKIP %s: no hg19 dir", sname)); next }

  genes <- fread(file.path(mtx_dir, "genes.tsv"), header = FALSE,
                 col.names = c("gene_id", "symbol"), sep = "\t",
                 stringsAsFactors = FALSE, showProgress = FALSE)
  mtx <- as(readMM(file.path(mtx_dir, "matrix.mtx")), "dgCMatrix")
  raw_total <- raw_total + ncol(mtx)

  total <- colSums(mtx)
  ng    <- diff(mtx@p)
  mt    <- colSums(mtx[grepl("^MT-", genes$symbol), , drop = FALSE])
  pct   <- mt / total * 100          # total==0 的空液滴 -> NaN

  # 直方图只累加 <= x 轴上限的值（长尾截断，超界部分不进 bin；
  # 总量另计入 raw_total）
  hist_a <- hist_a + hist(total[total <= XMAX_COUNTS], breaks = brk_a,
                          plot = FALSE)$counts
  hist_b <- hist_b + hist(ng[ng <= XMAX_GENES], breaks = brk_b,
                          plot = FALSE)$counts

  # 散点面板候选：只收"有内容"的 barcode（ng>=50、depth<=20k、
  # MT<=50）——38M barcode 全量存散点装不下，且空液滴在图中只是
  # 原点噪音；NaN（空液滴）显式排除
  sc_mask <- ng >= 50 & total <= XMAX_COUNTS & pct <= PCT_MT_MAX
  sc_mask[is.na(sc_mask)] <- FALSE
  if (any(sc_mask)) {
    scatter_rows[[sname]] <- data.frame(
      total_counts = total[sc_mask], n_genes = ng[sc_mask], pct_mt = pct[sc_mask])
  }
  # knee 图原料：每样本 total 全量保留（合计 38M 向量，排序放聚合段）
  rank_all[[sname]] <- total
  rm(mtx); gc()   # 每样本用完即释放矩阵，控制峰值内存

  if (i %% 5 == 0 || i == length(sample_dirs)) {
    log(sprintf("  %d/%d done [elapsed %s]", i, length(sample_dirs),
                format(difftime(Sys.time(), t0, units = "mins"), digits = 3)))
  }
}

# ══════════════════════════════════════════════════════════════
# 聚合：rank 向量（knee 图）、scatter 子集
# ══════════════════════════════════════════════════════════════
log("Building rank vector (38M barcodes) ...")
# unlist + 降序排序一步完成，中间向量在函数返回后即被释放；
# 峰值内存 ≈ 列表 + unlist 结果 + 排序结果 ≈ 3 × 300 MB
rank_sorted <- sort(unlist(rank_all, use.names = FALSE), decreasing = TRUE)
rm(rank_all); gc()
n_barcodes <- length(rank_sorted)
step <- max(1, floor(n_barcodes / 200000))   # 抽稀步长：38M → ~200k 点
knee <- data.frame(rank = seq(1, n_barcodes, by = step),
                   depth = rank_sorted[seq(1, n_barcodes, by = step)])
# 步长不一定整除总长，rbind 补上末尾一点（最小 depth 处），
# 保证 log-log 曲线右端完整
knee <- rbind(knee, data.frame(rank = n_barcodes, depth = rank_sorted[n_barcodes]))

# Panel C 三条阈值线：本管线下限 1,001 / 论文上限 6,000 /
# 本管线上限 15,000。n_above = 排序向量上的精确计数（sum(>=th)），
# 渲染端用它校验自己插值出的 n_above 误差 < 抽稀步长
thr_knee <- data.frame(y = c(MIN_UMIS, PAPER_MAX, MAX_UMIS),
                       lab = c("1,001 UMI", "6,000 UMI", "15,000 UMI"),
                       col = c("red", "purple", "green"),
                       n_above = sapply(c(MIN_UMIS, PAPER_MAX, MAX_UMIS),
                                        function(th) sum(rank_sorted >= th)))
rm(rank_sorted); gc()

log("Sampling scatter panel ...")
scatter <- rbindlist(scatter_rows)
# 散点面板最多 5 万点；
# set.seed 固定抽样，保证每次运行抽到同一批点（逐位复现）
set.seed(42)
if (nrow(scatter) > N_SCATTER) scatter <- scatter[sample.int(nrow(scatter), N_SCATTER), ]

log("Saving aggregates ...")
fwrite(data.frame(mid = (brk_a[-1] + brk_a[-length(brk_a)]) / 2, count = hist_a),
       file.path(DATA_DIR, "qc_plot_hist_total.csv"))
fwrite(data.frame(mid = (brk_b[-1] + brk_b[-length(brk_b)]) / 2, count = hist_b),
       file.path(DATA_DIR, "qc_plot_hist_genes.csv"))
fwrite(knee, file.path(DATA_DIR, "qc_plot_knee.csv"))
fwrite(scatter, file.path(DATA_DIR, "qc_plot_scatter.csv"))
writeLines(toJSON(list(
  raw_barcodes_total = raw_total, n_scatter = nrow(scatter),
  thresholds = as.list(setNames(thr_knee$n_above, thr_knee$lab)),
  bins = N_BINS, decimation_step = step
), pretty = TRUE, auto_unbox = TRUE),
  file.path(DATA_DIR, "qc_plot_summary.json"))

log(sprintf("Done [elapsed %s]. Raw barcodes: %d",
            format(difftime(Sys.time(), t0, units = "mins"), digits = 3), raw_total))
