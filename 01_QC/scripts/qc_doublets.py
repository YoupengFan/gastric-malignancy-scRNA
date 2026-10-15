#!/usr/bin/env python3
# ============================================================
# 01_QC (stage 2)：逐样本 Scrublet doublet 检测
#
# 输入：01_Data/raw/（52 个 GSM*_PatXX-{A,B} 原始 10x 矩阵；数据
#       不随仓库分发，从 GEO GSE150290 下载解压后放置于此）
#
# 仅对 QC 通过的细胞运行（barcode 来自 qc_meta.csv，阈值来自
# qc_filter.R）。逐样本：读 raw mtx → 取 QC 通过的 barcode 子集
# → 跑 Scrublet → 收集分数。
#
# 设计理由：Scrublet（k-NN 图歧义）用于剔除异型 doublet，相对
# Zhang 2021 QC 做与时俱进的改进（原论文无提及 doublet 处理）。
#
# Checkpoint/resume：逐样本结果以 append 方式写 /tmp/scrublet_ckpt/
# （ext4 写速优于 NTFS 挂载），每样本后 fsync；重启时已在
# checkpoint 中的样本跳过。最终 CSV/JSON 在收尾时从 checkpoint
# 组装。
#
# Usage (WSL)：
#   source ~/softwares/miniforge3/etc/profile.d/conda.sh
#   conda activate bioinfo
#   python 10_R_pipeline/01_QC/scripts/qc_doublets.py
# ============================================================

import importlib.metadata
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

import scrublet as scr

RAW_DIR = "/mnt/e/gastric-msda-upload/01_Data/raw"
QC_DIR = "/mnt/e/gastric-msda-upload/10_R_pipeline/01_QC/data"
OUT_CSV = os.path.join(QC_DIR, "gse150290_doublet_scores.csv")
OUT_JSON = os.path.join(QC_DIR, "qc_doublets_summary.json")
CKPT_DIR = "/tmp/scrublet_ckpt"
CKPT_CSV = os.path.join(CKPT_DIR, "gse150290_doublet_ckpt.csv")

THRESHOLD = 0.25  # Scrublet 默认分数阈值
DOUBLET_RATE_PER_1K = 0.008  # 10x 经验法则：每 1,000 个回收细胞约 0.8%

t0 = time.time()


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)


# ── QC 通过的 barcode ─────────────────────────────────────────
meta = pd.read_csv(os.path.join(QC_DIR, "gse150290_qc_meta.csv"))
meta["sample"] = meta["sample"].astype(str)
log(f"qc_meta: {len(meta)} cells across {meta['sample'].nunique()} samples")

# ── Checkpoint（断点续跑）────────────────────────────────────
# 完整性校验：写入中途被杀可能截断最后一个样本的行；其名字仍会
# 出现在 checkpoint 里，所以逐样本核对行数与 qc_meta，不完整的重做。
done_samples = set()
if os.path.exists(CKPT_CSV):
    try:
        ckpt = pd.read_csv(CKPT_CSV, on_bad_lines="skip")
        counts_ckpt = ckpt["sample"].astype(str).value_counts()
        counts_meta = meta["sample"].value_counts()
        complete = [
            s for s in counts_ckpt.index if counts_ckpt[s] == counts_meta.get(s, -1)
        ]
        for s in counts_ckpt.index:
            if s not in complete:
                log(
                    f"  Checkpoint incomplete for {s} "
                    f"({counts_ckpt[s]} vs {counts_meta.get(s, -1)} rows) — redo"
                )
        if len(complete) < len(counts_ckpt):
            ckpt = ckpt[ckpt["sample"].astype(str).isin(complete)]
            ckpt.to_csv(CKPT_CSV, index=False)  # rewrite trimmed checkpoint
        done_samples = set(complete)
        log(
            f"Checkpoint: {len(ckpt)} rows, {len(done_samples)} samples done — resuming"
        )
    except Exception as e:
        log(f"Checkpoint unreadable ({e}) — starting fresh")
else:
    os.makedirs(CKPT_DIR, exist_ok=True)

sample_dirs = sorted(
    p for p in Path(RAW_DIR).iterdir() if p.is_dir() and p.name.startswith("GSM")
)


def find_hg19(sdir):
    flat = sdir / "raw_gene_bc_matrices" / "hg19"
    req = ["matrix.mtx", "genes.tsv", "barcodes.tsv"]
    if all((flat / f).exists() for f in req):
        return flat
    for d in sdir.rglob("hg19"):
        if all((d / f).exists() for f in req):
            return d
    return None


for i, sdir in enumerate(sample_dirs, 1):
    sname = sdir.name
    if sname in done_samples:
        log(f"[{i}/{len(sample_dirs)}] {sname}: in checkpoint — skip")
        continue
    mtx_dir = find_hg19(sdir)
    if mtx_dir is None:
        log(f"  SKIP {sname}: no hg19 dir")
        continue
    sub = meta.loc[meta["sample"] == sname]
    if len(sub) == 0:
        log(f"  SKIP {sname}: no QC-passed cells")
        continue

    t1 = time.time()
    # 读 raw 矩阵（genes × cells，CSC 利于按列取子集）
    mtx = scipy.io.mmread(str(mtx_dir / "matrix.mtx")).tocsc()
    barcodes = pd.read_csv(mtx_dir / "barcodes.tsv", sep="\t", header=None)[0].astype(
        str
    )
    raw_name = barcodes.to_numpy()
    # 拼"样本-条码"全名（与 qc_filter.R 的命名规则一致），isin 匹配
    # 本样本 QC 通过的细胞位置
    full_name = np.array([f"{sname}-{b}" for b in raw_name])
    keep_mask = np.isin(full_name, sub["barcode"].to_numpy())
    # Scrublet 要求 cells × genes 输入；mtx 是 genes × cells
    # → 先按列（细胞）取子集，再转置
    mtx_pass = mtx[:, keep_mask].tocsc().T.tocsr()

    n_cells = mtx_pass.shape[0]
    # 期望 doublet 率随回收细胞数缩放（10x 经验法则 ~0.8%/千细胞）；
    # Scrublet 按此比例合成模拟 doublet 来标定 k-NN 分数
    expected_rate = DOUBLET_RATE_PER_1K * n_cells / 1000.0
    scrub = scr.Scrublet(
        mtx_pass, expected_doublet_rate=expected_rate, sim_doublet_ratio=2
    )
    # scrublet 0.2.3：scrub_doublets() 内部已用默认阈值调用过
    # call_doublets()；再调一次以应用我们的阈值，然后读 *_obs_ 属性
    scrub.scrub_doublets(verbose=False)
    scrub.call_doublets(threshold=THRESHOLD, verbose=False)
    scores = scrub.doublet_scores_obs_
    predicted = scrub.predicted_doublets_

    # 逐细胞结果表：score 为连续分数（下游可自定义阈值），
    # predicted 为 0.25 阈值下的二值判定
    df = pd.DataFrame(
        {
            "barcode": full_name[keep_mask],
            "sample": sname,
            "doublet_score": scores,
            "predicted_doublet": predicted,
        }
    )
    log(
        f"[{i}/{len(sample_dirs)}] {sname}: {n_cells} cells, "
        f"{int(predicted.sum())} doublets ({100*predicted.mean():.2f}%) "
        f"[{time.time()-t1:.0f}s]"
    )
    # checkpoint 以 append 方式写 (ext4) + fsync，被杀也不丢样本
    df.to_csv(CKPT_CSV, mode="a", header=not os.path.exists(CKPT_CSV), index=False)
    fd = os.open(CKPT_CSV, os.O_RDONLY)
    os.fsync(fd)
    os.close(fd)

# ── 从 checkpoint 组装最终产物 ────────────────────────────────
# 最终 CSV 直接由 checkpoint 重排写出（同源，保证断点续跑与全量
# 直跑结果一致）；JSON 为参数留档 + 逐样本汇总
out = pd.read_csv(CKPT_CSV, on_bad_lines="skip")
out.to_csv(OUT_CSV, index=False)

per_sample = {}
for sname, g in out.groupby("sample"):
    per_sample[sname] = {
        "n_cells": int(len(g)),
        "n_doublets": int(g["predicted_doublet"].sum()),
        "pct_doublets": round(100 * g["predicted_doublet"].mean(), 2),
        "median_score": round(float(np.median(g["doublet_score"])), 4),
    }

summary = {
    "step": "01_QC_doublets",
    "tool": f"scrublet {importlib.metadata.version('scrublet')}",
    "input": "QC-passed cells per qc_meta.csv (thresholds from qc_filter.R)",
    "parameters": {
        "threshold": THRESHOLD,
        "expected_doublet_rate": f"{DOUBLET_RATE_PER_1K} * n_cells / 1000",
        "sim_doublet_ratio": 2,
    },
    "totals": {
        "cells_scored": int(len(out)),
        "doublets": int(out["predicted_doublet"].sum()),
        "pct_doublets": round(100 * out["predicted_doublet"].mean(), 2),
    },
    "per_sample": per_sample,
    "note": (
        "Python original Tier 2 Scrublet removed 1036/163747 (0.63%) "
        "on a slightly different QC cell set; not used as a tuning target."
    ),
}
with open(OUT_JSON, "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

log(
    f"TOTAL: {len(out)} scored, {int(out['predicted_doublet'].sum())} doublets "
    f"({100*out['predicted_doublet'].mean():.2f}%) [{time.time()-t0:.0f}s]"
)
log(f"Output: {OUT_CSV}")
