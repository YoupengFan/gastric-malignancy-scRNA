#!/usr/bin/env Rscript
# ============================================================
# 恶性/非恶性两池构成：4 份 labels.csv（主 05 + 主 run6
# + 测试种子 run1/run6）回连注释，按 对照（Pat25-29）/
# A 癌旁 / B 肿瘤 拆出每池的细胞数与池内占比
#
# 用途：量化"参照池自污"论点——非恶性参照池中来自癌症患者
# （A+B）与直接来自肿瘤组织（B）的细胞占比；反向看"恶性"池中
# 对照样本细胞的占比。该构成对种子/低置信剔除稳定 = 结构性。
#
# 输入：05_Iterative_Classify/data/ 各目录 labels.csv
#       03_Cluster_Annotate/data/cell_annotation.csv（回连 patient/tissue）
# 产出：05_Iterative_Classify/data/pool_composition.csv
#       （run × pool × group 长表：n_cells + pct_of_pool）
#
# Usage（仓库根目录）：
#   "E:\R-4.5.1\bin\Rscript.exe" 05_Iterative_Classify\scripts\sensitivity\pool_composition.R
# ============================================================
suppressPackageStartupMessages(library(data.table))

REPO <- "."
DATA <- file.path(REPO, "05_Iterative_Classify/data")
ANNO <- file.path(REPO, "03_Cluster_Annotate/data/cell_annotation.csv")

anno <- fread(ANNO)[, .(barcode, patient, tissue)]

label_files <- c(
  main      = file.path(DATA, "labels.csv"),
  run6      = file.path(DATA, "run6_epi_chief_logFC_lowconf_excl", "labels.csv"),
  test_run1 = file.path(DATA, "testseed_run1_epi_chief_logFC", "labels.csv"),
  test_run6 = file.path(DATA, "testseed_run6_epi_chief_logFC_lowconf_excl", "labels.csv")
)

CTRL_PAT <- paste0("Pat", 25:29)

# labels 回连注释；group = 对照（优先于 tissue A）/ A / B
read_run <- function(path) {
  lab <- fread(path)
  m <- merge(lab, anno, by.x = "cell", by.y = "barcode")
  if (nrow(m) != nrow(lab))
    warning(sprintf("%s: %d/%d cells not found in annotation",
                    basename(path), nrow(lab) - nrow(m), nrow(lab)))
  m[, group := fifelse(patient %in% CTRL_PAT, "ctrl",
                       fifelse(tissue %in% "A", "A", "B"))]
  m
}

# ── run × pool × group 长表 ────────────────────────────────
pool_comp <- rbindlist(lapply(names(label_files), function(rn) {
  m <- read_run(label_files[[rn]])
  m[, .(n_cells = .N), by = .(pool = label, group)
    ][, `:=`(pct_of_pool = round(n_cells / sum(n_cells) * 100, 1),
             run = rn), by = pool]
}))
setcolorder(pool_comp, c("run", "pool", "group", "n_cells", "pct_of_pool"))
pool_comp <- pool_comp[order(run, pool, group)]
fwrite(pool_comp, file.path(DATA, "pool_composition.csv"))

cat("\n=== pool_composition.csv ===\n")
print(pool_comp)
cat(sprintf("\nWrote: pool_composition.csv (%d rows)\n", nrow(pool_comp)))
