#!/usr/bin/env Rscript
# ============================================================
# 04_Bulk_DEG: step1_bulk_deg — TCGA-STAD 配对 limma DEG（主种子）
#
# 输入：01_Data/tcga_RSEM_gene_tpm.gz（论文原引数据源，值 =
#       log2(TPM+0.001)）
#       + 01_Data/TCGA.STAD.sampleMap_HiSeqV2.gz（仅取 STAD 患者
#         前缀集合）
# 产出：04_Bulk_DEG/data/top50_{tumor,normal}_limma.txt
#       + 全表 deg_limma.csv + bulk_deg_params.json
#       + figures/fig_bulk_deg.png（火山图，种子基因高亮）
#
# 数据源核实与论文 Table S1 逐位复现细节：见 04_Bulk_DEG/README.md
#
# 关键设计要点：
#   - 函数化：step1_bulk_deg() 接收表达矩阵 + 样本表，返回
#     list(tumor_genes, normal_genes) —— "第一步 = 一个函数"
#   - plain limma（论文原文 "using the limma package"）
#   - DEG 口径 = padj<0.01，无 FC 截断；|logFC|>1 仅火山图虚线
#     与着色边界，不参与筛选
#   - Ensembl → symbol 本地映射（01_Data/GSE150290_genes_hg19.tsv），
#     未映射基因丢弃、重名 symbol 取平均表达最大行
#
# Usage (WSL):
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 04_Bulk_DEG/scripts/step1_bulk_deg.R
# ============================================================

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(jsonlite)
  library(ggplot2) # 火山图
})

# ══════════════════════════════════════════════════════════════
# 函数：step1_bulk_deg（plain limma，配对设计，tumor vs normal）
# ══════════════════════════════════════════════════════════════
# bulk_expr   : 基因 × 样本的表达矩阵（已是 log2 尺度，如 log2(TPM+0.001)）
# sample_table: 数据框，含 `sample`（bulk_expr 列名）、
#               `patient`（配对 ID）、`condition`（"Tumor"/"Normal"）
# n_top       : 每方向种子基因数
# padj        : 显著性阈值（Benjamini-Hochberg adj.P.Val）——DEG 口径只有 padj，
#               无 FC 截断（论文method正文口径）
# 返回       : list(tumor_genes, normal_genes, deg_table, design_info)
step1_bulk_deg <- function(bulk_expr, sample_table, n_top = 50, padj = 0.01,
                           min_tpm = 1, min_frac = 0.2, ...) {

  # ── 对齐样本 ───────────────────────────────────────────
  keep_s <- intersect(sample_table$sample, colnames(bulk_expr))
  if (length(keep_s) < 4) stop("too few matched samples")
  st <- sample_table[match(keep_s, sample_table$sample), , drop = FALSE]
  expr <- bulk_expr[, keep_s, drop = FALSE]
  n_genes_in <- nrow(expr)

  # ── 输入变换：log2(TPM + 1) ──────────────────────────────
  # 文件值 = log2(TPM + 0.001)
  tpm <- 2^expr - 0.001
  expr <- log2(tpm + 1)

  # ── 低表达过滤 TPM>1 in ≥20%（论文方法未写明，Table S1 佐证，
  #    细节见 README）──────────────────────────────────────
  keep_g <- rowSums(tpm > min_tpm) >= min_frac * ncol(expr)
  expr <- expr[keep_g, , drop = FALSE]

  condition <- factor(st$condition, levels = c("Normal", "Tumor"))
  patient <- factor(st$patient)

  # ── limma：配对设计 ~ patient + condition（论文原文口径）──
  design <- model.matrix(~ patient + condition)
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)

  deg <- topTable(fit, coef = "conditionTumor", number = Inf, sort.by = "none")
  deg$gene <- rownames(deg)

  # DEG 口径 = 论文方法正文：padj < 0.01，无 FC 截断
  sig <- deg[!is.na(deg$adj.P.Val) & deg$adj.P.Val < padj, , drop = FALSE]
  up <- sig$gene[sig$logFC > 0]
  dn <- sig$gene[sig$logFC < 0]

  # ── 每方向按 |logFC| 取前 n_top ──────────────────────
  tumor_genes <- head(up[order(-sig$logFC[sig$logFC > 0])], n_top)
  normal_genes <- head(dn[order(sig$logFC[sig$logFC < 0])], n_top)

  list(
    tumor_genes = tumor_genes,
    normal_genes = normal_genes,
    deg_table = deg,
    design_info = list(
      n_patients = length(unique(patient)),
      n_tumor = sum(condition == "Tumor"),
      n_normal = sum(condition == "Normal"),
      n_sig = nrow(sig),
      n_up = length(up),
      n_down = length(dn),
      padj_cutoff = padj,
      n_genes_in = n_genes_in,
      n_genes_filtered = nrow(expr),
      transform = "log2(TPM+1)",
      low_expr_filter = sprintf("TPM>%g in >=%.0f%% paired samples",
                                min_tpm, min_frac * 100),
      method = "limma, paired ~ patient + condition (log2(TPM+1), no voom)"
    )
  )
}

# ══════════════════════════════════════════════════════════════
# 函数：plot_bulk_deg_volcano（火山图，top50 种子基因高亮）
# ══════════════════════════════════════════════════════════════
# deg_table   : limma topTable 全表（含 logFC/adj.P.Val/gene）
# tumor_genes : 肿瘤上调种子基因
# normal_genes: 正常上调种子基因
# padj        : 显著性阈值（与 DEG 同款）
# fc_cutoff   : 仅画图——垂直虚线位置 ±fc_cutoff，
#               且线内显著点不上色（着色以线为界）；不参与基因筛选
# out_png     : 输出 PNG 路径
plot_bulk_deg_volcano <- function(deg_table, tumor_genes, normal_genes,
                                  padj, fc_cutoff, out_png, ...) {
  deg_plot <- deg_table[!is.na(deg_table$adj.P.Val), ]
  up_lab <- sprintf("Tumor up (padj<%.3g, |logFC|>%.3g)", padj, fc_cutoff)
  dn_lab <- sprintf("Normal up (padj<%.3g, |logFC|>%.3g)", padj, fc_cutoff)
  # 上色 = padj 显著且超出 ±fc_cutoff 线（论文 Fig 3A 画法）；线内显著点
  # 不上色——着色仅为图注口径，基因筛选仍是 padj-only
  deg_plot$col <- ifelse(deg_plot$gene %in% tumor_genes, "Tumor top50 (seed)",
                  ifelse(deg_plot$gene %in% normal_genes, "Normal top50 (seed)",
                  ifelse(deg_plot$adj.P.Val < padj & abs(deg_plot$logFC) > fc_cutoff,
                         ifelse(deg_plot$logFC > 0, up_lab, dn_lab), "ns")))
  # 种子基因最后绘制，压在上层
  deg_plot <- deg_plot[order(deg_plot$col %in% c("Tumor top50 (seed)", "Normal top50 (seed)")), ]
  deg_pal <- setNames(c("#c51b7d", "#4d9221", "#f1b6da", "#b8e186", "grey85"),
                      c("Tumor top50 (seed)", "Normal top50 (seed)", up_lab, dn_lab, "ns"))
  # 阈值虚线：padj 水平线（显著性）+ |logFC| 垂直线 ×2（论文 supp007 Fig A 同款，
  # 着色同样以线为界）；图例移除，分组由象限内文字标注说明
  # （深色点 = top50 种子基因）
  x_max <- max(abs(deg_plot$logFC))
  y_max <- max(-log10(pmax(deg_plot$adj.P.Val, 1e-300)))
  p <- ggplot(deg_plot, aes(x = logFC, y = -log10(pmax(adj.P.Val, 1e-300)),
                            color = col)) +
    geom_point(size = 0.3) +
    geom_hline(yintercept = -log10(padj), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed", color = "grey40") +
    annotate("text", x = x_max * 0.97, y = y_max * 0.96, label = "Tumor up",
             hjust = 1, vjust = 1, size = 4, fontface = "bold", color = "#c51b7d") +
    annotate("text", x = -x_max * 0.97, y = y_max * 0.96, label = "Normal up",
             hjust = 0, vjust = 1, size = 4, fontface = "bold", color = "#4d9221") +
    annotate("text", x = 0, y = 0, label = "ns",
             hjust = 0.5, vjust = -0.4, size = 4, color = "grey55") +
    annotate("text", x = x_max * 0.97, y = -log10(padj),
             label = sprintf("padj = %.3g", padj),
             hjust = 1, vjust = -0.6, size = 3, color = "grey40") +
    scale_color_manual(values = deg_pal) +
    labs(title = "Bulk DEG: main seed genes",
         x = "logFC (Tumor vs Normal)", y = "-log10 adj.P.Val") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  ggsave(out_png, p, device = ragg::agg_png, width = 7, height = 6, dpi = 300)
  invisible(out_png)
}

# ══════════════════════════════════════════════════════════════
# 主流程（TCGA-STAD 演示运行）
# ══════════════════════════════════════════════════════════════
REPO     <- "."
EXPR_GZ  <- file.path(REPO, "01_Data", "tcga_RSEM_gene_tpm.gz")
STAD_REF <- file.path(REPO, "01_Data", "TCGA.STAD.sampleMap_HiSeqV2.gz")
GENES_TSV <- file.path(REPO, "01_Data", "GSE150290_genes_hg19.tsv")
OUT_DIR  <- file.path(REPO, "04_Bulk_DEG", "data")
OUT_TMP  <- "/tmp/bulkdeg_r"

N_TOP <- 50
PADJ  <- 0.01
FC_LINE <- 1   # 仅火山图用（论文 supp007 Fig A 画法）：±1 虚线 + 线内显著点不上色；不参与基因筛选

unlink(OUT_TMP, recursive = TRUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TMP, recursive = TRUE, showWarnings = FALSE)
log_fh <- file(file.path(OUT_TMP, "step1_bulk_deg.log"), open = "a")
log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  writeLines(line, log_fh)
  flush(log_fh)
}

# ── 读取 Toil tpm 的 STAD 子集（Toil 文件自身全部 STAD 前缀样本）──
log("Loading tcga_RSEM_gene_tpm.gz (Xena Toil hub, log2 TPM) ...")
hdr_toil <- strsplit(readLines(gzfile(EXPR_GZ), n = 1), "\t", fixed = TRUE)[[1]]
hdr_stad <- strsplit(readLines(gzfile(STAD_REF), n = 1), "\t", fixed = TRUE)[[1]][-1]
stad_prefix <- unique(vapply(strsplit(hdr_stad, "-", fixed = TRUE),
                             function(p) paste(p[1:3], collapse = "-"), ""))
is_stad <- vapply(strsplit(hdr_toil, "-", fixed = TRUE), function(p)
  length(p) >= 3 && paste(p[1:3], collapse = "-") %in% stad_prefix, logical(1))
keep_idx <- which(is_stad)
log(sprintf("  STAD subset: %d samples (of %d total columns)",
            length(keep_idx), length(hdr_toil) - 1))

toil <- fread(cmd = paste("gzip -dc", shQuote(EXPR_GZ)), select = c(1L, keep_idx))
expr <- as.matrix(toil[, -1, drop = FALSE])
rownames(expr) <- toil[[1]]
colnames(expr) <- hdr_toil[keep_idx]
log(sprintf("  Matrix: %d genes x %d samples", nrow(expr), ncol(expr)))

# ── Ensembl → symbol 本地映射（01_Data/GSE150290_genes_hg19.tsv）──
log(sprintf("Mapping Ensembl -> symbol via %s", GENES_TSV))
gmap <- fread(GENES_TSV, header = FALSE, select = 1:2)
setnames(gmap, c("ensg", "symbol"))
gmap <- gmap[!duplicated(ensg)]
sym <- gmap$symbol[match(sub("\\..*$", "", rownames(expr)), gmap$ensg)]

n_unmapped <- sum(is.na(sym))
expr <- expr[!is.na(sym), , drop = FALSE]
sym <- sym[!is.na(sym)]
rownames(expr) <- sym
# 重名 symbol：取平均表达最大的一行（limma 常规去重口径）
ord <- order(rowMeans(expr), decreasing = TRUE)
expr <- expr[ord, , drop = FALSE]
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
log(sprintf("  %d genes mapped, %d unmapped dropped, %d unique symbols",
            nrow(expr), n_unmapped, nrow(expr)))

# ── 解析 barcode → 配对 Tumor/Normal ────────────────────
parse_tcga_barcode <- function(barcode) {
  parts <- strsplit(barcode, "-", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(c(patient = NA, code = NA))
  c(patient = paste(parts[1:3], collapse = "-"), code = substr(parts[4], 1, 2))
}
info <- t(sapply(colnames(expr), parse_tcga_barcode))
sample_type <- ifelse(info[, "code"] == "01", "Tumor",
                      ifelse(info[, "code"] == "11", "Normal", NA))

patient_list <- split(data.frame(sample = colnames(expr),
                                 condition = sample_type,
                                 stringsAsFactors = FALSE),
                      info[, "patient"])
paired <- Filter(function(x) all(c("Tumor", "Normal") %in% x$condition), patient_list)
log(sprintf("  Paired patients: %d", length(paired)))

sample_table <- do.call(rbind, lapply(names(paired), function(pid) {
  d <- paired[[pid]]
  data.frame(sample = c(d$sample[d$condition == "Tumor"][1],
                        d$sample[d$condition == "Normal"][1]),
             patient = pid,
             condition = c("Tumor", "Normal"),
             stringsAsFactors = FALSE)
}))

# ── 运行函数 ─────────────────────────────────────────────
log(sprintf("Running step1_bulk_deg (plain limma, log2(TPM+1), n_top=%d, padj=%g) ...",
            N_TOP, PADJ))
res <- step1_bulk_deg(expr, sample_table, n_top = N_TOP, padj = PADJ)

log(sprintf("  Genes: %d -> %d after low-expr filter (%s)",
            res$design_info$n_genes_in, res$design_info$n_genes_filtered,
            res$design_info$low_expr_filter))
log(sprintf("  Design: %d patients (%d tumor + %d normal)",
            res$design_info$n_patients, res$design_info$n_tumor, res$design_info$n_normal))
log(sprintf("  Significant (padj < %g): %d (T-up %d, N-up %d)",
            PADJ, res$design_info$n_sig,
            res$design_info$n_up, res$design_info$n_down))
log(sprintf("  Top %d tumor genes:  %s", N_TOP, paste(head(res$tumor_genes, 8), collapse = ", ")))
log(sprintf("  Top %d normal genes: %s", N_TOP, paste(head(res$normal_genes, 8), collapse = ", ")))

# ── 保存 ─────────────────────────────────────────────────────
writeLines(res$tumor_genes, file.path(OUT_TMP, "top50_tumor_limma.txt"))
writeLines(res$normal_genes, file.path(OUT_TMP, "top50_normal_limma.txt"))
fwrite(res$deg_table[order(res$deg_table$adj.P.Val), ],
       file.path(OUT_TMP, "deg_limma.csv"))

params <- list(
  step = "04_Bulk_DEG",
  function_name = "step1_bulk_deg",
  method = "limma, paired ~ patient + condition, input log2(TPM+1) (paper: 'DEGs ... using the limma package')",
  data = "tcga_RSEM_gene_tpm.gz (Xena Toil hub, log2(TPM+0.001)), STAD subset = Toil columns with STAD patient prefix",
  design = "~ patient + condition (paired), ref = Normal",
  thresholds = list(n_top = N_TOP, padj = PADJ),
  low_expr_filter = res$design_info$low_expr_filter,
  output = list(n_sig = res$design_info$n_sig, n_up = res$design_info$n_up,
                n_down = res$design_info$n_down, n_patients = res$design_info$n_patients,
                n_genes_in = res$design_info$n_genes_in,
                n_genes_filtered = res$design_info$n_genes_filtered)
)
writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(OUT_TMP, "bulk_deg_params.json"))

for (f in list.files(OUT_TMP, pattern = "\\.(txt|csv|json)$")) {
  file.copy(file.path(OUT_TMP, f), file.path(OUT_DIR, f), overwrite = TRUE)
}

# ── 图件 ─────────────────────────────────────────────────────
FIG_DIR <- file.path(REPO, "04_Bulk_DEG", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
plot_bulk_deg_volcano(res$deg_table, res$tumor_genes, res$normal_genes,
                      PADJ, FC_LINE, file.path(OUT_TMP, "fig_bulk_deg.png"))
for (f in list.files(OUT_TMP, pattern = "\\.png$")) {
  file.copy(file.path(OUT_TMP, f), file.path(FIG_DIR, f), overwrite = TRUE)
}

log(sprintf("Done. Output: %s", OUT_DIR))
close(log_fh)
