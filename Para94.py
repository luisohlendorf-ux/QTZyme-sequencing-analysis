#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Inline-barcode demultiplex
from pathlib import Path
import gzip
import csv
from collections import Counter

# ---------- paths ----------
BASE = Path()
R1_IN = BASE / "Undetermined_S0_L001_R1_001.fastq.gz"
R2_IN = BASE / "Undetermined_S0_L001_R2_001.fastq.gz"

OUTDIR = BASE / "HDV_inline_demux_pos128_seedAGATCG_HDV2fixed_jupyter"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ---------- demux settings ----------
BARCODE_TO_LABEL = {
    "TGACCA": "HDV4",
    "CGATGT": "HDV2",   
}
ADAPTER_SEED = "AGATCG"  
WINDOW_START = 120         
WINDOW_END = 145           

labels = ["HDV4", "HDV2", "unassigned"]

# ---------- open writers ----------
writers = {}
for lab in labels:
    writers[(lab, "R1")] = gzip.open(OUTDIR / f"{lab}_inline_R1.fastq.gz", "wt")
    writers[(lab, "R2")] = gzip.open(OUTDIR / f"{lab}_inline_R2.fastq.gz", "wt")

# ---------- demultiplex ----------
counts = Counter()

with gzip.open(R1_IN, "rt") as f1, gzip.open(R2_IN, "rt") as f2:
    while True:
        h1 = f1.readline()
        h2 = f2.readline()
        if not h1 or not h2:
            break

        s1 = f1.readline().strip()
        s2 = f2.readline().strip()
        p1 = f1.readline()
        p2 = f2.readline()
        q1 = f1.readline()
        q2 = f2.readline()

        cls = "unassigned"

        max_i = min(len(s1) - 6 - 3 - len(ADAPTER_SEED) + 1, WINDOW_END)
        for i in range(max(0, WINDOW_START), max_i):
            bc = s1[i:i+6]
            if bc in BARCODE_TO_LABEL:
                j = i + 6 + 3  # barcode + NNN
                if s1[j:j+len(ADAPTER_SEED)] == ADAPTER_SEED:  # exact
                    cls = BARCODE_TO_LABEL[bc]
                    break

        counts[cls] += 1
        counts["total"] += 1

        writers[(cls, "R1")].write(h1); writers[(cls, "R1")].write(s1 + "\n"); writers[(cls, "R1")].write(p1); writers[(cls, "R1")].write(q1)
        writers[(cls, "R2")].write(h2); writers[(cls, "R2")].write(s2 + "\n"); writers[(cls, "R2")].write(p2); writers[(cls, "R2")].write(q2)

for w in writers.values():
    w.close()

# ---------- write summary ----------
summary_tsv = OUTDIR / "inline_demux_summary.tsv"
with open(summary_tsv, "w", newline="") as fh:
    wr = csv.writer(fh, delimiter="\t")
    wr.writerow(["source", "HDV4", "HDV2", "unassigned", "total"])
    wr.writerow([
        "Undetermined_S0_L001",
        counts["HDV4"],
        counts["HDV2"],
        counts["unassigned"],
        counts["total"]
    ])

print("Done.")
print("Output folder:", OUTDIR)
print("Summary file:", summary_tsv)
print(dict(counts))


# In[7]:


# trim + merge
import subprocess
import shutil

TRIM_DIR = BASE / "HDV_trim_merge_cutadaptAGATCGG_vsearch_jupyter"
for d in [TRIM_DIR, TRIM_DIR/"cutadapt", TRIM_DIR/"merged", TRIM_DIR/"logs", TRIM_DIR/"tables"]:
    d.mkdir(parents=True, exist_ok=True)

if shutil.which("cutadapt") is None:
    raise RuntimeError("cutadapt not found in PATH")
if shutil.which("vsearch") is None:
    raise RuntimeError("vsearch not found in PATH")

def run(cmd, log_file):
    with open(log_file, "w") as lf:
        subprocess.run(cmd, check=True, stdout=lf, stderr=lf)

def count_fastq(path):
    op = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with op(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            f.readline(); f.readline(); f.readline()
            n += 1
    return n

samples = ["HDV2", "HDV4"]

for s in samples:
    r1 = OUTDIR / f"{s}_inline_R1.fastq.gz"
    r2 = OUTDIR / f"{s}_inline_R2.fastq.gz"

    c1 = TRIM_DIR / "cutadapt" / f"{s}_cut_R1.fastq.gz"
    c2 = TRIM_DIR / "cutadapt" / f"{s}_cut_R2.fastq.gz"

    merged = TRIM_DIR / "merged" / f"{s}_merged.fastq"
    notm1 = TRIM_DIR / "merged" / f"{s}_notmerged_R1.fastq"
    notm2 = TRIM_DIR / "merged" / f"{s}_notmerged_R2.fastq"

    cut_log = TRIM_DIR / "logs" / f"{s}_cutadapt.log"
    merge_log = TRIM_DIR / "logs" / f"{s}_vsearch_merge.log"

    # trimming
    cmd_cut = [
        "cutadapt",
        "-a", "AGATCGG",
        "-A", "AGATCGG",
        "-e", "0.1",
        "-O", "10",
        "-m", "30",
        "-o", str(c1),
        "-p", str(c2),
        str(r1), str(r2),
    ]
    run(cmd_cut, cut_log)

    # merging
    cmd_merge = [
        "vsearch",
        "--fastq_mergepairs", str(c1),
        "--reverse", str(c2),
        "--fastqout", str(merged),
        "--fastqout_notmerged_fwd", str(notm1),
        "--fastqout_notmerged_rev", str(notm2),
        "--fastq_minovlen", "25",
        "--fastq_maxdiffs", "10",
        "--fastq_allowmergestagger",
        "--threads", "8",
    ]
    run(cmd_merge, merge_log)

# summary table
rows = []
for s in samples:
    raw = count_fastq(OUTDIR / f"{s}_inline_R1.fastq.gz")
    cut = count_fastq(TRIM_DIR / "cutadapt" / f"{s}_cut_R1.fastq.gz")
    mrg = count_fastq(TRIM_DIR / "merged" / f"{s}_merged.fastq")
    rows.append([s, raw, cut, mrg, round(mrg/cut, 4) if cut else 0, round(mrg/raw, 4) if raw else 0])

summary_tsv = TRIM_DIR / "tables" / "trim_merge_summary.tsv"
with open(summary_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample", "raw_pairs", "pairs_after_cutadapt", "merged_reads", "merge_ratio_vs_cutadapt", "merge_ratio_vs_raw"])
    w.writerows(rows)

print("Done:", TRIM_DIR)
print("Summary:", summary_tsv)
for r in rows:
    print(r)


# In[6]:


import os
os.environ["PATH"] += ":/opt/homebrew/bin"


# In[8]:


# Cell 3: Core95 extraction
# Uses BASE and TRIM_DIR from your Cell 1/2

from collections import Counter
import pandas as pd
import json

EXACT_DIR = TRIM_DIR / "exact95_profile_shifted"
for d in [EXACT_DIR, EXACT_DIR/"extracted_fastq", EXACT_DIR/"tables", EXACT_DIR/"figures"]:
    d.mkdir(parents=True, exist_ok=True)

rules = {
    "HDV2": {"black_len": 111, "trim5": 16},
    "HDV4": {"black_len": 118, "trim5": 23},
}
samples = ["HDV2", "HDV4"]
anchor = "GGGTCGGCATGGCATC"

def parse_fastq(path):
    with open(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().strip().upper()
            p = f.readline()
            q = f.readline().strip()
            yield h, s, p, q

summary_rows = []

for s in samples:
    in_fq = TRIM_DIR / "merged" / f"{s}_merged.fastq"
    out_fq = EXACT_DIR / "extracted_fastq" / f"{s}_core95.fastq"

    need_black = rules[s]["black_len"]
    trim5 = rules[s]["trim5"]

    total = 0
    with_anchor = 0
    without_anchor = 0
    req_black = 0
    wrong_black = 0
    core95 = 0
    black_len_counts = Counter()

    with open(out_fq, "wt") as fo:
        for h, seq, p, q in parse_fastq(in_fq):
            total += 1
            i = seq.find(anchor)
            if i == -1:
                without_anchor += 1
                continue

            with_anchor += 1
            black = seq[:i]
            black_len_counts[len(black)] += 1

            if len(black) != need_black:
                wrong_black += 1
                continue

            req_black += 1
            core = black[trim5:trim5+95]
            if len(core) != 95:
                continue

            core95 += 1
            qq = q[trim5:trim5+95] if len(q) >= trim5+95 else ("I" * 95)

            fo.write(h)
            fo.write(core + "\n")
            fo.write("+\n")
            fo.write(qq + "\n")

    summary_rows.append({
        "sample": s,
        "merged_reads_total": total,
        "reads_with_anchor": with_anchor,
        "reads_without_anchor": without_anchor,
        "required_black_len": need_black,
        "reads_with_required_black_len": req_black,
        "reads_wrong_black_len": wrong_black,
        "core95_reads": core95,
        "core95_ratio_vs_total": round(core95 / total if total else 0, 4),
        "core95_ratio_vs_anchor": round(core95 / with_anchor if with_anchor else 0, 4),
        "top_black_lengths_from_anchor": ";".join([f"{k}:{v}" for k, v in black_len_counts.most_common(10)])
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(EXACT_DIR / "tables" / "core95_extraction_summary.tsv", sep="\t", index=False)

cfg = {"samples": samples, "rules": rules, "anchor": anchor}
with open(EXACT_DIR / "tables" / "core95_extraction_config.json", "w") as fh:
    json.dump(cfg, fh, indent=2)

print("Extraction done:", EXACT_DIR / "extracted_fastq")
display(summary_df)


# ## Cell 4: Build mutation tables 
# 
# import numpy as np
# from collections import Counter
# 
# REF_111 = "CTTCGGAGTCTAGGAGAGGTGCAAACAAACAACAAAACAAACAAACAGGGCAGCGCAGTCGGACATTGATAACGGAATCGAGCCCTATTGTTGTTGATATGCACCTCTCTT"
# REF_95 = REF_111[16:]
# L = 95
# samples = ["HDV2", "HDV4"]
# 
# def parse_fastq(path):
#     with open(path, "rt") as f:
#         while True:
#             h = f.readline()
#             if not h:
#                 break
#             s = f.readline().strip().upper()
#             p = f.readline()
#             q = f.readline().strip()
#             yield h, s, p, q
# 
# pos_rows = []
# 
# for s in samples:
#     fq = EXACT_DIR / "extracted_fastq" / f"{s}_core95.fastq"
#     base_counts = [Counter() for _ in range(L)]
# 
#     for h, seq, p, q in parse_fastq(fq):
#         if len(seq) != L:
#             continue
#         for i, b in enumerate(seq):
#             base_counts[i][b] += 1
# 
#     for i in range(L):
#         t = sum(base_counts[i].values())
#         rb = REF_95[i]
#         mut_reads = t - base_counts[i][rb]
# 
#         pos_rows.append({
#             "sample": s,
#             "position_1based": i + 1,
#             "ref_base": rb,
#             "total_reads": t,
#             "mut_reads": mut_reads,
#             "mut_ratio": (mut_reads / t) if t else 0.0,
#             "A_count": base_counts[i]["A"],
#             "C_count": base_counts[i]["C"],
#             "G_count": base_counts[i]["G"],
#             "T_count": base_counts[i]["T"],
#             "N_count": base_counts[i]["N"],
#         })
# 
# pos_df = pd.DataFrame(pos_rows)
# 
# # base tables
# pos_df[["sample","position_1based","ref_base","total_reads","mut_reads","mut_ratio"]].to_csv(
#     EXACT_DIR / "tables" / "mutation_rate_by_position_core95.tsv", sep="\t", index=False
# )
# pos_df.to_csv(EXACT_DIR / "tables" / "mutation_details_by_position_core95.tsv", sep="\t", index=False)
# 
# wide = pos_df.pivot(index="position_1based", columns="sample", values="mut_ratio").reset_index()
# wide.insert(1, "ref_base", [REF_95[i-1] for i in wide["position_1based"]])
# wide.to_csv(EXACT_DIR / "tables" / "mutation_rate_by_position_core95_wide.tsv", sep="\t", index=False)
# 
# mean_df = (
#     pos_df.groupby("sample", as_index=False)["mut_ratio"]
#     .mean()
#     .rename(columns={"mut_ratio": "mean_mut_ratio_across_95_positions"})
# )
# mean_df.to_csv(EXACT_DIR / "tables" / "mean_mutation_ratio_core95.tsv", sep="\t", index=False)
# 
# # detailed mutation-target table
# detail_df = pos_df.copy()
# for b in ["A", "C", "G", "T"]:
#     detail_df[f"mut_to_{b}"] = detail_df[f"{b}_count"]
# for b in ["A", "C", "G", "T"]:
#     detail_df.loc[detail_df["ref_base"] == b, f"mut_to_{b}"] = 0
# 
# target_cols = ["mut_to_A", "mut_to_C", "mut_to_G", "mut_to_T"]
# detail_df["top_mut_to"] = detail_df[target_cols].idxmax(axis=1).str.replace("mut_to_", "", regex=False)
# detail_df["top_mut_count"] = detail_df[target_cols].max(axis=1)
# detail_df["top_mut_ratio_of_mut"] = np.where(
#     detail_df["mut_reads"] > 0,
#     detail_df["top_mut_count"] / detail_df["mut_reads"],
#     0.0
# )
# 
# detail_cols = [
#     "sample","position_1based","ref_base","total_reads","mut_reads","mut_ratio",
#     "A_count","C_count","G_count","T_count","N_count",
#     "mut_to_A","mut_to_C","mut_to_G","mut_to_T",
#     "top_mut_to","top_mut_count","top_mut_ratio_of_mut"
# ]
# detail_df = detail_df[detail_cols]
# detail_df.to_csv(EXACT_DIR / "tables" / "mutation_details_by_position_core95_with_targets.tsv", sep="\t", index=False)
# 
# print("Tables written:", EXACT_DIR / "tables")
# display(mean_df)
# display(detail_df.head(20))

# In[13]:


# Cell 5: Heatmap + line plot for last 90 nt only (positions 6..95)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

pos_df = pd.read_csv(EXACT_DIR / "tables" / "mutation_rate_by_position_core95.tsv", sep="\t")

start_pos = 6   # last 90 nt of 95
samples = ["HDV2", "HDV4"]
pos_sub = pos_df[pos_df["position_1based"] >= start_pos].copy()

# Heatmap
mat = np.vstack([pos_sub[pos_sub["sample"] == s]["mut_ratio"].to_numpy() for s in samples])
x_labels = pos_sub[pos_sub["sample"] == samples[0]]["position_1based"].to_numpy()

plt.figure(figsize=(12, 2.8))
im = plt.imshow(mat, aspect="auto", cmap="magma", interpolation="nearest")
plt.yticks([0, 1], samples)
tick_idx = np.arange(0, len(x_labels), 5)
plt.xticks(tick_idx, x_labels[tick_idx], fontsize=8)
plt.xlabel(f"Position in normalized core95 region (1-based, {start_pos}-95)")
plt.title("Mutation-ratio heatmap (last 90 nt)")
cbar = plt.colorbar(im)
cbar.set_label("Mutation ratio")
plt.tight_layout()
plt.savefig(EXACT_DIR / "figures" / "mutation_ratio_heatmap_core95_last90nt.png", dpi=240)
plt.show()

# Line plot
plt.figure(figsize=(12, 4.5))
for s, col in [("Selected", "#ff7f0e"), ("Unselected", "#1f77b4")]:
    d = pos_sub[pos_sub["sample"] == s]
    plt.plot(d["position_1based"], d["mut_ratio"], label=s, color=col, lw=1.8)

plt.xlabel(f"Position in normalized core95 region (1-based, {start_pos}-95)")
plt.ylabel("Mutation ratio")
plt.title("Per-position mutation ratio (last 90 nt)")
plt.ylim(bottom=0)
plt.grid(alpha=0.25)
plt.legend()
plt.tight_layout()
plt.savefig(EXACT_DIR / "figures" / "mutation_ratio_lineplot_core95_last90nt.png", dpi=240)
plt.show()

print("Figures saved in:", EXACT_DIR / "figures")


# In[11]:


# Cell 7: Build tables (including mutation-target detail table)

pos_df = pd.read_csv(INT_DIR / "tables" / "pos_df_integrated_core95.tsv", sep="\t")

# Basic rate table
pos_df[["sample","position_1based","ref_base","total_reads","mut_reads","mut_ratio"]].to_csv(
    INT_DIR / "tables" / "mutation_rate_by_position_core95_integrated.tsv", sep="\t", index=False
)

# Wide table
wide = pos_df.pivot(index="position_1based", columns="sample", values="mut_ratio").reset_index()
wide.insert(1, "ref_base", [REF_95[i-1] for i in wide["position_1based"]])
wide.to_csv(INT_DIR / "tables" / "mutation_rate_by_position_core95_integrated_wide.tsv", sep="\t", index=False)

# Mutation-target table
detail_df = pos_df.copy()
for b in ["A", "C", "G", "T"]:
    detail_df[f"mut_to_{b}"] = detail_df[f"{b}_count"]
for b in ["A", "C", "G", "T"]:
    detail_df.loc[detail_df["ref_base"] == b, f"mut_to_{b}"] = 0

target_cols = ["mut_to_A","mut_to_C","mut_to_G","mut_to_T"]
detail_df["top_mut_to"] = detail_df[target_cols].idxmax(axis=1).str.replace("mut_to_", "", regex=False)
detail_df["top_mut_count"] = detail_df[target_cols].max(axis=1)
detail_df["top_mut_ratio_of_mut"] = np.where(
    detail_df["mut_reads"] > 0,
    detail_df["top_mut_count"] / detail_df["mut_reads"],
    0.0
)

detail_cols = [
    "sample","position_1based","ref_base","total_reads","mut_reads","mut_ratio",
    "A_count","C_count","G_count","T_count","N_count",
    "mut_to_A","mut_to_C","mut_to_G","mut_to_T",
    "top_mut_to","top_mut_count","top_mut_ratio_of_mut"
]
detail_df = detail_df[detail_cols]
detail_df.to_csv(INT_DIR / "tables" / "mutation_details_by_position_core95_integrated_with_targets.tsv", sep="\t", index=False)

# Last-90 subset + mean table
start_pos = 6
last90 = pos_df[pos_df["position_1based"] >= start_pos].copy()
last90.to_csv(INT_DIR / "tables" / "mutation_rate_by_position_last90_integrated.tsv", sep="\t", index=False)

means = (
    last90.groupby("sample", as_index=False)["mut_ratio"]
    .mean()
    .rename(columns={"mut_ratio": "mean_mut_ratio_last90"})
)
means.to_csv(INT_DIR / "tables" / "mean_mutation_ratio_last90_integrated.tsv", sep="\t", index=False)

print("Tables written:", INT_DIR / "tables")
display(means)


# In[ ]:




