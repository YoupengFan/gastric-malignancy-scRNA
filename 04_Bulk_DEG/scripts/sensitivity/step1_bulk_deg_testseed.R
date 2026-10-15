#!/usr/bin/env Rscript
# ============================================================
# 04_Bulk_DEG (sensitivity): step1_bulk_deg_testseed — 测试种子生成器
#
# 用途：按 limma-voom 配方用 HiSeqV2 矩阵生成测试种子，与主种子
#       （../step1_bulk_deg.R，Toil + plain limma）构成双种子敏感性
#       对照。
#
# 输入：01_Data/TCGA.STAD.sampleMap_HiSeqV2.gz（值 = log2(RSEM+1)）
# 产出：04_Bulk_DEG/data/testseed/top50_{tumor,normal}_limma.txt
#       + bulk_deg_testseed_params.json
#       + figures/fig_bulk_deg_testseed.png（火山图，与主种子同款画法）
#
# 关键设计要点：
#   - 函数化：step1_testseed_deg() —— "第一步 = 一个函数"；
#     方法（anti-log → counts>1 in ≥20% → TMM → voom →
#     eBayes trend/robust）逐步骤注释在函数体内
#
# Usage (WSL):
#   export PATH="$HOME/softwares/bin:$PATH"
#   Rscript 04_Bulk_DEG/scripts/sensitivity/step1_bulk_deg_testseed.R
# ============================================================

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)   # 仅 calcNormFactors(TMM) 归一化；DEG 检验仅 limma
  library(ggplot2) # 火山图
  library(data.table)
  library(jsonlite)
})

# ══════════════════════════════════════════════════════════════
# 函数：step1_testseed_deg（limma-voom，配对设计，tumor vs normal）
# ══════════════════════════════════════════════════════════════
# bulk_expr   : 基因 × 样本表达矩阵（已是 log2 尺度 = log2(RSEM+1)）
# sample_table: 数据框，含 sample（bulk_expr 列名）、patient（配对
#               ID）、condition（"Tumor"/"Normal"）
# n_top/padj  : 每方向种子基因数 / BH 显著性阈值
# 返回        : list(tumor_genes, normal_genes, deg_table, design_info)
step1_testseed_deg <- function(bulk_expr, sample_table, n_top = 50, padj = 0.01) {
  keep_s <- intersect(sample_table$sample, colnames(bulk_expr))
  if (length(keep_s) < 4) stop("too few matched samples")
  st <- sample_table[match(keep_s, sample_table$sample), , drop = FALSE]
  expr <- bulk_expr[, keep_s, drop = FALSE]

  condition <- factor(st$condition, levels = c("Normal", "Tumor"))
  patient <- factor(st$patient)
  design <- model.matrix(~ patient + condition)

  # ── anti-log → pseudo-counts（取整 + 负值截零）──
  pseudo_counts <- round(2^expr - 1)
  pseudo_counts[pseudo_counts < 0] <- 0

  # ── 低表达过滤：counts > 1 在 ≥20% 配对样本──
  min_samples <- ceiling(ncol(pseudo_counts) * 0.2)
  keep <- rowSums(pseudo_counts > 1) >= min_samples
  counts <- pseudo_counts[keep, , drop = FALSE]

  # ── limma-voom：TMM + voom + eBayes(trend, robust)──
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design)
  fit <- eBayes(lmFit(v, design), trend = TRUE, robust = TRUE)

  deg <- topTable(fit, coef = "conditionTumor", number = Inf, sort.by = "none")
  deg$gene <- rownames(deg)

  sig <- deg[!is.na(deg$adj.P.Val) & deg$adj.P.Val < padj, , drop = FALSE]
  up <- sig$gene[sig$logFC > 0]
  dn <- sig$gene[sig$logFC < 0]
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
      n_genes_in = nrow(bulk_expr),
      n_genes_filtered = nrow(counts),
      filter = sprintf("pseudo_counts>1 in >=%.0f%% paired samples", 20),
      method = "limma-voom: round(2^x-1) -> filter -> TMM -> voom -> eBayes(trend, robust), paired"
    )
  )
}

# ══════════════════════════════════════════════════════════════
# 函数：plot_bulk_deg_volcano（火山图，top50 种子基因高亮）
# 与主种子 ../step1_bulk_deg.R 同款画法，仅标题注明数据源/检验链差异
# ══════════════════════════════════════════════════════════════
plot_bulk_deg_volcano <- function(deg_table, tumor_genes, normal_genes,
                                  padj, fc_cutoff, out_png, ...) {
  deg_plot <- deg_table[!is.na(deg_table$adj.P.Val), ]
  up_lab <- sprintf("Tumor up (padj<%.3g, |logFC|>%.3g)", padj, fc_cutoff)
  dn_lab <- sprintf("Normal up (padj<%.3g, |logFC|>%.3g)", padj, fc_cutoff)
  # 上色 = padj 显著且超出 ±fc_cutoff 线；线内显著点不上色——
  # 着色仅为图注口径，基因筛选仍是 padj-only
  deg_plot$col <- ifelse(deg_plot$gene %in% tumor_genes, "Tumor top50 (seed)",
                  ifelse(deg_plot$gene %in% normal_genes, "Normal top50 (seed)",
                  ifelse(deg_plot$adj.P.Val < padj & abs(deg_plot$logFC) > fc_cutoff,
                         ifelse(deg_plot$logFC > 0, up_lab, dn_lab), "ns")))
  # 种子基因最后绘制，压在上层
  deg_plot <- deg_plot[order(deg_plot$col %in% c("Tumor top50 (seed)", "Normal top50 (seed)")), ]
  deg_pal <- setNames(c("#c51b7d", "#4d9221", "#f1b6da", "#b8e186", "grey85"),
                      c("Tumor top50 (seed)", "Normal top50 (seed)", up_lab, dn_lab, "ns"))
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
    labs(title = "Bulk DEG: test seed genes",
         x = "logFC (Tumor vs Normal)", y = "-log10 adj.P.Val") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  ggsave(out_png, p, device = ragg::agg_png, width = 7, height = 6, dpi = 300)
  invisible(out_png)
}

# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════
REPO    <- "."
HIV2_GZ <- file.path(REPO, "01_Data", "TCGA.STAD.sampleMap_HiSeqV2.gz")
OUT_DIR <- file.path(REPO, "04_Bulk_DEG", "data", "testseed")
OUT_TMP <- "/tmp/testseed_deg_r"
FIG_DIR <- file.path(REPO, "04_Bulk_DEG", "figures")
N_TOP <- 50
PADJ  <- 0.01
FC_LINE <- 1   # 仅火山图：±1 虚线 + 线内显著点不上色；不参与基因筛选

unlink(OUT_TMP, recursive = TRUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TMP, recursive = TRUE, showWarnings = FALSE)
log_fh <- file(file.path(OUT_TMP, "step1_bulk_deg_testseed.log"), open = "a")
log <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n")
  writeLines(line, log_fh)
  flush(log_fh)
}

# ── 读取 HiSeqV2（Xena sampleMap：行 = 基因符号，列 = STAD barcode）──
log(sprintf("Loading %s ...", HIV2_GZ))
hdr <- strsplit(readLines(gzfile(HIV2_GZ), n = 1), "\t", fixed = TRUE)[[1]]
mat <- fread(cmd = paste("gzip -dc", shQuote(HIV2_GZ)))
expr <- as.matrix(mat[, -1, drop = FALSE])
rownames(expr) <- mat[[1]]
colnames(expr) <- hdr[-1]
log(sprintf("  Matrix: %d genes x %d samples", nrow(expr), ncol(expr)))

# 重名 symbol 去重（sampleMap 矩阵通常已唯一；保险起见按平均表达取最大行）
ord <- order(rowMeans(expr), decreasing = TRUE)
expr <- expr[ord, , drop = FALSE]
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
log(sprintf("  %d unique symbols after dedup", nrow(expr)))

# ── 解析 barcode → 配对 Tumor/Normal（HiSeqV2 内 01/11 齐全 = 32 对）──
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

# ── 运行函数 ─────────────────────────────────────────────────
log(sprintf("Running step1_testseed_deg (n_top=%d, padj=%g) ...", N_TOP, PADJ))
res <- step1_testseed_deg(expr, sample_table, n_top = N_TOP, padj = PADJ)
log(sprintf("  Design: %d patients (%d tumor + %d normal)",
            res$design_info$n_patients, res$design_info$n_tumor, res$design_info$n_normal))
log(sprintf("  Genes: %d -> %d after filter (%s)",
            res$design_info$n_genes_in, res$design_info$n_genes_filtered,
            res$design_info$filter))
log(sprintf("  Significant (padj < %g): %d (T-up %d, N-up %d)",
            PADJ, res$design_info$n_sig, res$design_info$n_up, res$design_info$n_down))
log(sprintf("  Top %d tumor genes:  %s", N_TOP, paste(head(res$tumor_genes, 8), collapse = ", ")))
log(sprintf("  Top %d normal genes: %s", N_TOP, paste(head(res$normal_genes, 8), collapse = ", ")))

# ── 保存 ─────────────────────────────────────────────────────
writeLines(res$tumor_genes, file.path(OUT_TMP, "top50_tumor_limma.txt"))
writeLines(res$normal_genes, file.path(OUT_TMP, "top50_normal_limma.txt"))

params <- list(
  step = "04_Bulk_DEG (sensitivity)",
  function_name = "step1_testseed_deg",
  role = "testseed — 早期误用数据源的种子，保留作稳健性/敏感性对照（历史见 04_Bulk_DEG/README.md）",
  method = res$design_info$method,
  data = "TCGA.STAD.sampleMap_HiSeqV2.gz (Xena sampleMap, log2(RSEM+1))",
  design = "~ patient + condition (paired), ref = Normal",
  thresholds = list(n_top = N_TOP, padj = PADJ),
  filter = res$design_info$filter,
  output = list(n_sig = res$design_info$n_sig, n_up = res$design_info$n_up,
                n_down = res$design_info$n_down, n_patients = res$design_info$n_patients,
                n_genes_in = res$design_info$n_genes_in,
                n_genes_filtered = res$design_info$n_genes_filtered)
)
writeLines(toJSON(params, pretty = TRUE, auto_unbox = TRUE),
           file.path(OUT_TMP, "bulk_deg_testseed_params.json"))

for (f in list.files(OUT_TMP, pattern = "\\.(txt|csv|json)$")) {
  file.copy(file.path(OUT_TMP, f), file.path(OUT_DIR, f), overwrite = TRUE)
}

# ── 火山图（与主种子同款画法，标题注明数据源/检验链差异；补充材料）──
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
plot_bulk_deg_volcano(res$deg_table, res$tumor_genes, res$normal_genes,
                      PADJ, FC_LINE, file.path(OUT_TMP, "fig_bulk_deg_testseed.png"))
file.copy(file.path(OUT_TMP, "fig_bulk_deg_testseed.png"),
          file.path(FIG_DIR, "fig_bulk_deg_testseed.png"), overwrite = TRUE)

log(sprintf("Done. Output: %s", OUT_DIR))
close(log_fh)
