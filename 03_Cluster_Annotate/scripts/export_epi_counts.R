#!/usr/bin/env Rscript
# ============================================================
# 03_Cluster_Annotate: export_epi_counts — 导出 05 的上皮谱系极简输入
#
# 背景：gse150290_annotated.rds（~7 GB Seurat 全对象：15.9 万细胞 +
#       归一化矩阵 + PCA/Harmony/UMAP + 图）无法随仓库分发；而 05 的
#       step2_iterative_classify.R 与 sensitivity/run_iterative_runs.R
#       只用其中上皮谱系细胞的 RAW counts + annotation + patient。
# 产出：03_Cluster_Annotate/data/gse150290_epi_counts.rds（xz 压缩）
#       list(counts = dgCMatrix（基因 × 上皮谱系细胞，RAW UMI counts）,
#            meta   = data.frame(cell, annotation, patient,
#                                annotation_low_conf 若存在)）
#       上皮谱系 = Epithelial / Chief / Enteroendocrine / Parietal 四类型；
#       annotation 列完整保留——05 内部的 EPI_TYPES 筛选逻辑原样生效，
#       任意上皮口径（run1/2/5 等）均可在此文件上重放。
#
# 为仓库瘦身，把 05 的输入从 ~7 GB 全对象替换为本极简子集；
#
# Usage (WSL):
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 03_Cluster_Annotate/scripts/export_epi_counts.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# ── 路径 ─────────────────────────────────────────────────────
ANNO_RDS <- "03_Cluster_Annotate/data/gse150290_annotated.rds"
OUT_RDS  <- "03_Cluster_Annotate/data/gse150290_epi_counts.rds"

# 上皮谱系口径：广义上皮四类型（05 各 Run 的具体口径在其脚本内自行筛选）
EPI_TYPES <- c("Epithelial", "Chief", "Enteroendocrine", "Parietal")

log <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

log(sprintf("Loading %s ...", ANNO_RDS))
obj <- readRDS(ANNO_RDS)
log(sprintf("  %d cells, annotations: %s", ncol(obj),
            paste(sort(unique(obj$annotation)), collapse = ", ")))

md <- obj@meta.data
keep <- md$annotation %in% EPI_TYPES
cells <- rownames(md)[keep]
counts <- obj[["RNA"]]$counts[, cells]

meta_cols <- intersect(c("annotation", "patient", "annotation_low_conf"),
                       colnames(md))
meta <- md[cells, meta_cols, drop = FALSE]
meta$cell <- rownames(meta)

log(sprintf("Epithelial lineage: %d / %d cells (%.1f%%)",
            ncol(counts), nrow(md), 100 * ncol(counts) / nrow(md)))
log(sprintf("  nonzeros: %d, sparsity: %.1f%%", length(counts@x),
            100 * (1 - length(counts@x) / prod(dim(counts)))))
log(sprintf("  meta columns: %s", paste(meta_cols, collapse = ", ")))

log(sprintf("Saving %s (xz) ...", OUT_RDS))
saveRDS(list(counts = counts, meta = meta), OUT_RDS, compress = "xz")
log(sprintf("Done. size = %s", format(file.info(OUT_RDS)$size, big.mark = ",")))
