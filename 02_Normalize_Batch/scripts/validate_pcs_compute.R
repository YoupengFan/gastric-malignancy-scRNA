#!/usr/bin/env Rscript
# ============================================================
# 02_Normalize_Batch PC 验证 — 计算部分（WSL）
#
# 输入：03_Cluster_Annotate/data/gse150290_annotated.rds（含 pca +
#       harmony embeddings 与 CC/MT meta）
# 产出（02_Normalize_Batch/data/）：
#   pc_validation.csv               per-PC 指标（pre/post 两套系列）
#   pc_validation_summary.json      阈值检查 + Harmony 前后对比汇总
#   fnn_mixing_percell.csv          每细胞 k=20 同样本比例 pre/post
#   （图由 validate_pcs_figure.R 在 Windows 渲染，本脚本只算不画）
#
# 关键设计要点：
#   - 方差口径：stdev_i^2 / 2000 × 100（2000 HVG）
#   - LINEAGE_MARKERS = 谱系 marker 列表 + 胃上皮谱系扩充（PGC/LIPF/
#     ATP4A/ATP4B/GIF/CHGA/GHRL/SST/REG1A/LYZ/TFF1/TFF3/GAST/GKN2）
#   - Harmony 校验（sample 键）：PC1-30 平均 R²、PC1 R²、FNN k=20
#     混合率 pre/post
#   - 后-Harmony "方差解释" = 各维边际方差（与 pre 同分母同口径；
#     Harmony 移除的批次方差体现为总方差下降与分布变平）
#   - per-PC 指标在慢速 FNN 段之前提前落盘——timeout 无法清空结果
#
# Usage (WSL):
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 02_Normalize_Batch/scripts/validate_pcs_compute.R
#   # ~40-90 min：6.8 GB RDS 加载 (~2 min) + FNN get.knn 50 维 ×
#   # 160k 细胞（KD-tree 高维退化，慢但可完成，已验证）
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(FNN)
})

# ── 路径 ─────────────────────────────────────────────────────
BASE     <- "."
RDS_PATH <- file.path(BASE, "03_Cluster_Annotate/data/gse150290_annotated.rds")
DATA_DIR <- file.path(BASE, "02_Normalize_Batch/data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 阈值 ────────────────────────────────────────────────────
ETA2_BATCH_MAX <- 0.15    # PC1-3：eta2(patient~PC) 应 < 0.15
R2_CC_MT_MAX   <- 0.05    # PC1-10：R2(CC/MT~PC) 应 < 0.05
N_PCS          <- 30      # 详细评估的 PC 数
N_HVG          <- 2000    # 已缩放基因数（方差分母）

# ── loading 分类用基因类别 ───────────────────────────────────
LINEAGE_MARKERS <- c(
  'EPCAM','CDH1','KRT8','KRT18','KRT19',          # 上皮
  'CD3D','CD3E','CD3G','CD4','CD8A','CD8B',       # T 细胞
  'CD79A','CD79B','MS4A1','CD19','MZB1','JCHAIN','DERL3',  # B/浆细胞
  'CD14','CD68','FCGR3A','CSF1R','CST3',          # 髓系
  'NKG7','GNLY','KLRD1','KLRF1',                  # NK
  'MKI67','TOP2A','PCNA',                         # 增殖
  'TFF2','MUC5AC','MUC6','GKN1','PGA4',           # 胃上皮
  'COL1A1','COL1A2','ACTA2','RGS5',               # 间质
  'PECAM1','VWF','CDH5',                          # 内皮
  'KIT','TPSAB1','CPA3',                          # 肥大
  'PGC','LIPF','REG1A','LYZ',                     # Chief 主细胞
  'ATP4A','ATP4B','GIF',                          # Parietal 壁细胞
  'CHGA','GHRL','SST','GAST',                     # Enteroendocrine 肠内分泌
  'TFF1','TFF3','GKN2')                           # 表面黏液/小凹上皮
TECH_PREFIXES <- c('RPL','RPS','MRPL','MRPS','MT-')
TECH_EXACT    <- c('MALAT1')

t0 <- Sys.time()
log <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# ══════════════════════════════════════════════════════════════
# 读入已注释 Seurat 对象
# ══════════════════════════════════════════════════════════════
log("Loading annotated RDS (6.8 GB) ...")
obj <- readRDS(RDS_PATH)
log(sprintf("Loaded: %d cells x %d features", ncol(obj), nrow(obj)))

pca      <- Embeddings(obj, "pca")      # pre-Harmony PC（cells × 50）
harm     <- Embeddings(obj, "harmony")  # post-Harmony embedding
stdev    <- obj[["pca"]]@stdev
loadings <- obj[["pca"]]@feature.loadings   # genes × PCs
meta     <- obj@meta.data

patient <- as.factor(meta$patient)
sample  <- as.factor(meta$sample)
tissue  <- meta$tissue
if (is.null(tissue) || all(is.na(tissue))) {
  tissue <- sub(".*-(.)$", "\\1", meta$sample)   # 回退：从 sample 后缀解析
}
tissue <- as.factor(tissue)

n_cells <- ncol(obj)
rm(obj)            # 释放 Seurat 对象；下文只需 embedding/meta
invisible(gc())

# ══════════════════════════════════════════════════════════════
# 逐 PC 指标
# ══════════════════════════════════════════════════════════════
eta2 <- function(x, g) {   # eta² = 组间 SS / 总 SS
  g <- as.factor(g)
  gm <- mean(x)
  ss_t <- sum((x - gm)^2)
  ss_b <- sum(tapply(x, g, function(v) length(v) * (mean(v) - gm)^2))
  ss_b / ss_t
}

X_ccmt <- cbind(1, meta$S.Score, meta$G2M.Score, meta$percent.mt)
r2_ccmt <- function(y) {   # y ~ S.Score + G2M.Score + percent.mt 的 R²
  beta <- solve(crossprod(X_ccmt), crossprod(X_ccmt, y))
  1 - sum((y - X_ccmt %*% beta)^2) / sum((y - mean(y))^2)
}

log(sprintf("Evaluating PC1-%d ...", N_PCS))
res <- lapply(seq_len(N_PCS), function(i) {
  pcv <- pca[, i]
  ld  <- loadings[, i]
  ord <- order(abs(ld), decreasing = TRUE)
  top <- rownames(loadings)[ord[seq_len(20)]]
  is_lineage <- top %in% LINEAGE_MARKERS
  is_tech    <- vapply(top, function(g)
    any(startsWith(g, TECH_PREFIXES)) || g %in% TECH_EXACT, logical(1))
  list(
    pc            = i,
    variance_pct  = stdev[i]^2 / N_HVG * 100,
    eta2_patient  = eta2(pcv, patient),
    eta2_tissue   = eta2(pcv, tissue),
    r2_cc_mt      = r2_ccmt(pcv),
    n_lineage     = sum(is_lineage),
    n_tech        = sum(is_tech),
    n_other       = 20L - sum(is_lineage) - sum(is_tech),
    top_genes     = paste(top[1:5], collapse = ",")
  )
})
res <- rbindlist(res)
res[, cum_var_pct := cumsum(variance_pct)]

# ── post-Harmony per-PC 系列（before/after 图的矫正后半边）────
# 后-Harmony "方差解释" = 各维边际方差（与 pre 的 stdev²/2000 同
# 分母同口径；Harmony 移除的批次方差体现为总方差的下降与分布的
# 变平）；批次信号 = eta2(sample ~ dim)，pre/post 成对对比。
# eta2_sample_pre 的 PC1-30 均值即下方 harmony validation 的
# r2_sample_pre（同一计算，逐 PC 持久化）
harm_var <- apply(harm[, seq_len(N_PCS)], 2, var)
res[, variance_pct_post := harm_var / N_HVG * 100]
res[, cum_var_pct_post  := cumsum(variance_pct_post)]
res[, eta2_sample_pre  := vapply(seq_len(N_PCS),
                                 function(i) eta2(pca[, i],  sample), numeric(1))]
res[, eta2_sample_post := vapply(seq_len(N_PCS),
                                 function(i) eta2(harm[, i], sample), numeric(1))]

# 在慢速 FNN 段之前先落盘 per-PC 指标——timeout 无法清空结果
fwrite(res, file.path(DATA_DIR, "pc_validation.csv"))
log("Saved (partial): pc_validation.csv")

# ══════════════════════════════════════════════════════════════
# Harmony 校验（sample 键）
# ══════════════════════════════════════════════════════════════
log("Harmony validation (sample R^2 pre/post, FNN k=20 mixing) ...")
r2_sample_pre  <- mean(vapply(seq_len(N_PCS), function(i) eta2(pca[, i],  sample), numeric(1)))
r2_sample_post <- mean(vapply(seq_len(N_PCS), function(i) eta2(harm[, i], sample), numeric(1)))
pc1_r2_pre  <- eta2(pca[, 1],  sample)
pc1_r2_post <- eta2(harm[, 1], sample)

nn_metrics <- function(emb) {   # 一次 get.knn 调用 → 汇总 + 逐细胞两用
  nn  <- get.knn(emb, k = 20)$nn.index
  lab <- as.integer(sample)
  nnl <- matrix(lab[nn], nrow = nrow(nn))
  is_diff <- nnl != lab[row(nnl)]        # cells x 20: 邻居是否不同样本
  list(mean_diff = mean(is_diff),        # 不同样本邻居比例均值（JSON 用）
       same      = 1 - rowMeans(is_diff))  # 每细胞同样本比例（直方图用）
}
log("  FNN k=20 on pca (pre) — 160k x 50 dims, allow minutes ...")
pre_m  <- nn_metrics(pca)
log("  FNN k=20 on harmony (post) ...")
post_m <- nn_metrics(harm)
mix_pre  <- pre_m$mean_diff
mix_post <- post_m$mean_diff

log(sprintf("  mean R2(sample~PC1-3): %.4f -> %.4f", r2_sample_pre, r2_sample_post))
log(sprintf("  PC1 R2(sample):        %.4f -> %.4f", pc1_r2_pre, pc1_r2_post))
log(sprintf("  FNN k=20 mixing:       %.4f -> %.4f", mix_pre, mix_post))

# per-cell same-sample ratio（渲染端 batch mixing 直方图面板用）
fnn_pc <- data.frame(barcode         = rownames(pca),
                     same_ratio_pre  = pre_m$same,
                     same_ratio_post = post_m$same)
fwrite(fnn_pc, file.path(DATA_DIR, "fnn_mixing_percell.csv"))
log("Saved: fnn_mixing_percell.csv")

# ══════════════════════════════════════════════════════════════
# 产出（图由 validate_pcs_figure.R 在 Windows 渲染）
# ══════════════════════════════════════════════════════════════
writeLines(toJSON(list(
  n_cells            = n_cells,
  n_pcs_evaluated    = N_PCS,
  n_hvg              = N_HVG,
  variance_pct_pc1   = res$variance_pct[1],
  cum_var_pc30       = res$cum_var_pct[N_PCS],
  max_eta2_patient_pc1_3 = max(res$eta2_patient[1:3]),
  max_eta2_tissue_pc1_3  = max(res$eta2_tissue[1:3]),
  max_r2_cc_mt_pc1_10    = max(res$r2_cc_mt[1:10]),
  thresholds = list(eta2_batch_max = ETA2_BATCH_MAX, r2_cc_mt_max = R2_CC_MT_MAX),
  checks = list(
    eta2_patient_ok = max(res$eta2_patient[1:3]) < ETA2_BATCH_MAX,
    eta2_tissue_ok  = max(res$eta2_tissue[1:3])  < ETA2_BATCH_MAX,
    r2_cc_mt_ok     = max(res$r2_cc_mt[1:10])    < R2_CC_MT_MAX),
  harmony = list(
    mean_r2_sample_pc1_3_pre  = r2_sample_pre,
    mean_r2_sample_pc1_3_post = r2_sample_post,
    pc1_r2_pre  = pc1_r2_pre,
    pc1_r2_post = pc1_r2_post,
    fnn_k20_mixing_pre  = mix_pre,
    fnn_k20_mixing_post = mix_post)
), pretty = TRUE, auto_unbox = TRUE), file.path(DATA_DIR, "pc_validation_summary.json"))
log("Saved: pc_validation.csv, pc_validation_summary.json")
log(sprintf("Done [elapsed %s]",
            format(difftime(Sys.time(), t0, units = "mins"), digits = 3)))
