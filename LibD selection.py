#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pathlib import Path
import gzip
import csv

BASE = Path()
OUT = BASE / "N50_extracted"
OUT.mkdir(parents=True, exist_ok=True)

FLANK = "GGCAGGGCAGCGCAGTC"
N50_LENGTH = 50

r1_files = sorted(BASE.glob("*_R1_001.fastq.gz"))

if not r1_files:
    raise FileNotFoundError(f"No R1 FASTQ files found in {BASE}")

for fastq in r1_files:
    output_tsv = OUT / fastq.name.replace(".fastq.gz", "_N50.tsv")

    total_reads = 0
    flank_found = 0
    n50_extracted = 0

    with gzip.open(fastq, "rt") as fin, open(output_tsv, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["read_id", "N50"])

        while True:
            header = fin.readline().strip()
            if not header:
                break

            sequence = fin.readline().strip().upper()
            fin.readline()  # +
            fin.readline()  # quality

            total_reads += 1
            flank_start = sequence.find(FLANK)

            if flank_start == -1:
                continue

            flank_found += 1

            if flank_start < N50_LENGTH:
                continue

            n50 = sequence[flank_start - N50_LENGTH:flank_start]

            if len(n50) == N50_LENGTH:
                writer.writerow([header.removeprefix("@"), n50])
                n50_extracted += 1

    print(
        f"{fastq.name}: total={total_reads:,}, "
        f"flank_found={flank_found:,}, "
        f"N50_extracted={n50_extracted:,}"
    )

print("Output folder:", OUT)


# In[2]:


from pathlib import Path
import pandas as pd

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
IN_DIR = BASE / "N50_extracted"
OUT_DIR = BASE / "N50_collapsed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

for infile in sorted(IN_DIR.glob("*_N50.tsv")):
    df = pd.read_csv(infile, sep="\t")

    counts = (
        df["N50"]
        .value_counts()
        .rename_axis("N50")
        .reset_index(name="count")
    )

    total = counts["count"].sum()
    counts["fraction"] = counts["count"] / total if total else 0

    outfile = OUT_DIR / infile.name.replace("_N50.tsv", "_collapsed.tsv")
    counts.to_csv(outfile, sep="\t", index=False)

    print(
        f"{infile.name}: reads={total:,}, "
        f"unique_N50={len(counts):,}"
    )

print("Output folder:", OUT_DIR)


# In[3]:


from pathlib import Path
import pandas as pd
import subprocess
import shutil

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
IN_DIR = BASE / "N50_collapsed"
OUT_DIR = BASE / "N50_clusters_97"
OUT_DIR.mkdir(parents=True, exist_ok=True)

if shutil.which("vsearch") is None:
    raise RuntimeError("vsearch not found in PATH")

# Combine identical N50s across all rounds and sum their read counts
global_counts = {}

for infile in sorted(IN_DIR.glob("*_collapsed.tsv")):
    df = pd.read_csv(infile, sep="\t")

    for sequence, count in zip(df["N50"], df["count"]):
        global_counts[sequence] = global_counts.get(sequence, 0) + int(count)

# VSEARCH cluster_size expects sequences sorted by decreasing abundance
sorted_sequences = sorted(
    global_counts.items(),
    key=lambda x: x[1],
    reverse=True
)

input_fasta = OUT_DIR / "all_N50_global.fasta"
mapping_tsv = OUT_DIR / "sequence_ids.tsv"

with open(input_fasta, "w") as fasta, open(mapping_tsv, "w") as mapping:
    mapping.write("sequence_id\tN50\ttotal_count\n")

    for i, (sequence, count) in enumerate(sorted_sequences, start=1):
        sequence_id = f"seq_{i}"
        fasta.write(f">{sequence_id};size={count}\n{sequence}\n")
        mapping.write(f"{sequence_id}\t{sequence}\t{count}\n")

cmd = [
    "vsearch",
    "--cluster_size", str(input_fasta),
    "--id", "0.97",
    "--sizein",
    "--sizeout",
    "--strand", "plus",
    "--centroids", str(OUT_DIR / "family_centroids.fasta"),
    "--uc", str(OUT_DIR / "cluster_assignments.uc"),
]

log_file = OUT_DIR / "vsearch_clustering.log"

with open(log_file, "w") as log:
    subprocess.run(
        cmd,
        check=True,
        stdout=log,
        stderr=subprocess.STDOUT
    )

print(f"Globally unique N50 sequences: {len(global_counts):,}")
print("Cluster assignments:", OUT_DIR / "cluster_assignments.uc")
print("Family centroids:", OUT_DIR / "family_centroids.fasta")


# In[4]:


from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_97"

# Load sequence IDs and global counts
seq_info = pd.read_csv(
    CLUSTER_DIR / "sequence_ids.tsv",
    sep="\t"
)

seq_count = dict(zip(seq_info["sequence_id"], seq_info["total_count"]))
seq_sequence = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

# Parse VSEARCH cluster assignments
cluster_counts = {}
cluster_centroids = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] not in {"S", "H"}:
            continue

        cluster_id = int(fields[1])
        query_id = fields[8].split(";")[0]

        cluster_counts[cluster_id] = (
            cluster_counts.get(cluster_id, 0)
            + seq_count.get(query_id, 0)
        )

        if fields[0] == "S":
            cluster_centroids[cluster_id] = seq_sequence.get(query_id, "")

# Create ranked table
clusters = pd.DataFrame([
    {
        "cluster": f"Cluster {cluster_id + 1}",
        "centroid": cluster_centroids.get(cluster_id, ""),
        "read_count": count
    }
    for cluster_id, count in cluster_counts.items()
]).sort_values("read_count", ascending=False)

clusters["fraction"] = clusters["read_count"] / clusters["read_count"].sum()

clusters.to_csv(
    CLUSTER_DIR / "cluster_abundance.tsv",
    sep="\t",
    index=False
)

# Plot top clusters
top = clusters.head(20).sort_values("fraction")

plt.figure(figsize=(9, 7))
plt.barh(top["cluster"], top["fraction"] * 100)

plt.xlabel("Share of all extracted reads (%)")
plt.ylabel("Sequence family")
plt.title("Top N50 sequence clusters at 97% identity")

plt.tight_layout()
plt.savefig(
    CLUSTER_DIR / "top_20_clusters.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

display(clusters.head(20))


# In[5]:


from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
COLLAPSED_DIR = BASE / "N50_collapsed"
CLUSTER_DIR = BASE / "N50_clusters_97"

# Load sequence IDs and sequences
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

# Map each N50 sequence to its cluster
sequence_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] not in {"S", "H"}:
            continue

        cluster_id = int(fields[1])
        sequence_id = fields[8].split(";")[0]

        if sequence_id in id_to_seq:
            sequence_to_cluster[id_to_seq[sequence_id]] = cluster_id

# Identify the globally top 10 clusters
cluster_abundance = pd.read_csv(
    CLUSTER_DIR / "cluster_abundance.tsv",
    sep="\t"
)

top10_clusters = [
    int(name.replace("Cluster ", "")) - 1
    for name in cluster_abundance.head(10)["cluster"]
]

# Calculate cluster shares in each round
round_cluster_counts = {}
round_totals = {}

for infile in sorted(COLLAPSED_DIR.glob("*_collapsed.tsv")):
    match = re.search(r"(?:LUR|^R)(1|4|5|6|7|8)", infile.name)

    if not match:
        print(f"Could not identify round: {infile.name}")
        continue

    round_name = f"R{match.group(1)}"
    df = pd.read_csv(infile, sep="\t")

    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"])
    df["cluster_id"] = df["cluster_id"].astype(int)

    round_totals[round_name] = df["count"].sum()
    round_cluster_counts[round_name] = df.groupby("cluster_id")["count"].sum()

rounds = ["R1", "R4", "R5", "R6", "R7", "R8"]

shares = pd.DataFrame(
    index=[f"Cluster {cluster_id + 1}" for cluster_id in top10_clusters],
    columns=rounds,
    dtype=float
)

for cluster_id in top10_clusters:
    for round_name in rounds:
        count = round_cluster_counts[round_name].get(cluster_id, 0)
        shares.loc[f"Cluster {cluster_id + 1}", round_name] = (
            100 * count / round_totals[round_name]
        )

shares.to_csv(
    CLUSTER_DIR / "top10_cluster_share_by_round.tsv",
    sep="\t"
)

plt.figure(figsize=(10, 6))

for cluster_name in shares.index:
    plt.plot(
        rounds,
        shares.loc[cluster_name],
        marker="o",
        label=cluster_name
    )

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Abundance trajectories of the top 10 sequence clusters")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top10_cluster_trajectories.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
display(shares.round(4))


# In[6]:


# Combine Clusters 1 and 4, rank by R8 abundance, and replot

combined_shares = shares.copy()

combined_shares.loc["Cluster 1 + 4"] = (
    combined_shares.loc["Cluster 1"] +
    combined_shares.loc["Cluster 4"]
)

combined_shares = combined_shares.drop(
    index=["Cluster 1", "Cluster 4"]
)

# Select and rank the top 10 by round 8 abundance
top10_combined = (
    combined_shares
    .sort_values("R8", ascending=False)
    .head(10)
)

top10_combined.index = [
    f"Clone {i}" for i in range(1, len(top10_combined) + 1)
]

# Save the renamed abundance table
top10_combined.to_csv(
    CLUSTER_DIR / "top10_clones_cluster1_and4_combined.tsv",
    sep="\t"
)

# Plot
plt.figure(figsize=(10, 6))

for clone in top10_combined.index:
    plt.plot(
        rounds,
        top10_combined.loc[clone],
        marker="o",
        label=clone
    )

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Abundance trajectories of the top 10 clones")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top10_clone_trajectories_combined.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
display(top10_combined.round(4))


# In[7]:


import numpy as np
import matplotlib.pyplot as plt

all_rounds = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]

plot_data = top10_combined.reindex(columns=all_rounds)

plt.figure(figsize=(10, 6))

for clone in plot_data.index:
    plt.plot(
        all_rounds,
        plot_data.loc[clone],
        marker="o",
        label=clone
    )

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Abundance trajectories of the top 10 clones")
plt.xticks(all_rounds)
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top10_clone_trajectories_all_rounds.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[8]:


import matplotlib.pyplot as plt

round_positions = [1, 4, 5, 6, 7, 8]
round_labels = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]

plt.figure(figsize=(10, 6))

for clone in top10_combined.index:
    plt.plot(
        round_positions,
        top10_combined.loc[clone, ["R1", "R4", "R5", "R6", "R7", "R8"]],
        marker="o",
        label=clone
    )

plt.xticks(range(1, 9), round_labels)
plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Abundance trajectories of the top 10 clones")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top10_clone_trajectories_all_rounds_connected.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[9]:


from pathlib import Path
import pandas as pd
import subprocess
import shutil

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
IN_DIR = BASE / "N50_collapsed"
OUT_DIR = BASE / "N50_clusters_94"
OUT_DIR.mkdir(parents=True, exist_ok=True)

if shutil.which("vsearch") is None:
    raise RuntimeError("vsearch not found in PATH")

# Pool identical N50 sequences across all rounds
global_counts = {}

for infile in sorted(IN_DIR.glob("*_collapsed.tsv")):
    df = pd.read_csv(infile, sep="\t")

    for sequence, count in zip(df["N50"], df["count"]):
        global_counts[sequence] = global_counts.get(sequence, 0) + int(count)

sorted_sequences = sorted(
    global_counts.items(),
    key=lambda item: item[1],
    reverse=True
)

input_fasta = OUT_DIR / "all_N50_global.fasta"
mapping_tsv = OUT_DIR / "sequence_ids.tsv"

with open(input_fasta, "w") as fasta, open(mapping_tsv, "w") as mapping:
    mapping.write("sequence_id\tN50\ttotal_count\n")

    for i, (sequence, count) in enumerate(sorted_sequences, start=1):
        sequence_id = f"seq_{i}"
        fasta.write(f">{sequence_id};size={count}\n{sequence}\n")
        mapping.write(f"{sequence_id}\t{sequence}\t{count}\n")

with open(OUT_DIR / "vsearch_clustering.log", "w") as log:
    subprocess.run(
        [
            "vsearch",
            "--cluster_size", str(input_fasta),
            "--id", "0.94",
            "--sizein",
            "--sizeout",
            "--strand", "plus",
            "--centroids", str(OUT_DIR / "family_centroids.fasta"),
            "--uc", str(OUT_DIR / "cluster_assignments.uc"),
        ],
        check=True,
        stdout=log,
        stderr=subprocess.STDOUT
    )

print(f"Unique N50 sequences clustered: {len(global_counts):,}")
print("94% clustering complete:", OUT_DIR)


# In[10]:


from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
COLLAPSED_DIR = BASE / "N50_collapsed"
CLUSTER_DIR = BASE / "N50_clusters_94"

# Map each N50 sequence to its 94% cluster
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

sequence_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] not in {"S", "H"}:
            continue

        sequence_id = fields[8].split(";")[0]
        cluster_id = int(fields[1])

        if sequence_id in id_to_seq:
            sequence_to_cluster[id_to_seq[sequence_id]] = cluster_id

# Calculate each cluster's percentage of each round
data = {}
rounds_with_data = ["R1", "R4", "R5", "R6", "R7", "R8"]

for infile in sorted(COLLAPSED_DIR.glob("*_collapsed.tsv")):
    match = re.search(r"(?:LUR|^R)(1|4|5|6|7|8)", infile.name)

    if not match:
        continue

    round_name = f"R{match.group(1)}"
    df = pd.read_csv(infile, sep="\t")

    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"])
    df["cluster_id"] = df["cluster_id"].astype(int)

    total_reads = df["count"].sum()
    data[round_name] = (
        df.groupby("cluster_id")["count"].sum() / total_reads * 100
    )

shares_94 = pd.DataFrame(data).fillna(0)
shares_94 = shares_94.reindex(columns=rounds_with_data, fill_value=0)

# Rank by round 8 and rename
top10_94 = shares_94.sort_values("R8", ascending=False).head(10).copy()
top10_94.index = [f"Clone {i}" for i in range(1, 11)]

top10_94.to_csv(
    CLUSTER_DIR / "top10_clones_share_by_round.tsv",
    sep="\t"
)

# Plot with rounds 2 and 3 shown, while connecting R1 directly to R4
x_positions = [1, 4, 5, 6, 7, 8]

plt.figure(figsize=(10, 6))

for clone in top10_94.index:
    plt.plot(
        x_positions,
        top10_94.loc[clone, rounds_with_data],
        marker="o",
        label=clone
    )

plt.xticks(range(1, 9), [f"R{i}" for i in range(1, 9)])
plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Top 10 clone trajectories at 94% identity")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top10_clone_trajectories_94.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
display(top10_94.round(4))


# In[2]:


from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
COLLAPSED_DIR = BASE / "N50_collapsed"
CLUSTER_DIR = BASE / "N50_clusters_94"

# Map each N50 sequence to its 94% cluster
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

sequence_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] in {"S", "H"}:
            seq_id = fields[8].split(";")[0]
            if seq_id in id_to_seq:
                sequence_to_cluster[id_to_seq[seq_id]] = int(fields[1])

# Load round 8 and identify the second-most abundant cluster
r8_file = next(COLLAPSED_DIR.glob("*R8*_collapsed.tsv"))
r8 = pd.read_csv(r8_file, sep="\t")
r8["cluster_id"] = r8["N50"].map(sequence_to_cluster)
r8 = r8.dropna(subset=["cluster_id"])
r8["cluster_id"] = r8["cluster_id"].astype(int)

clone2_cluster = (
    r8.groupby("cluster_id")["count"]
    .sum()
    .sort_values(ascending=False)
    .index[1]
)

clone2 = r8[r8["cluster_id"] == clone2_cluster].copy()

sequences = clone2["N50"].tolist()
weights = clone2["count"].to_numpy(dtype=float)
weights /= weights.sum()

# Read-count-weighted mutual information
length = 50
bases = "ACGT"
mi = np.zeros((length, length))

for i in range(length):
    for j in range(i + 1, length):
        joint = {
            (a, b): sum(
                w for seq, w in zip(sequences, weights)
                if seq[i] == a and seq[j] == b
            )
            for a in bases for b in bases
        }

        pi = {a: sum(joint[a, b] for b in bases) for a in bases}
        pj = {b: sum(joint[a, b] for a in bases) for b in bases}

        for a in bases:
            for b in bases:
                p = joint[a, b]
                if p > 0 and pi[a] > 0 and pj[b] > 0:
                    mi[i, j] += p * np.log2(p / (pi[a] * pj[b]))

        mi[j, i] = mi[i, j]

plt.figure(figsize=(8, 7))
plt.imshow(mi, origin="lower", aspect="equal")
plt.colorbar(label="Mutual information (bits)")
plt.xlabel("N50 position")
plt.ylabel("N50 position")
plt.title("Clone 2 covariation in round 8")
plt.xticks(range(0, 50, 5), range(1, 51, 5))
plt.yticks(range(0, 50, 5), range(1, 51, 5))
plt.tight_layout()
plt.show()

print(f"Clone 2 cluster ID: {clone2_cluster}")
print(f"Unique sequences: {len(clone2):,}")
print(f"Total reads: {clone2['count'].sum():,}")


# In[3]:


# Create the read-count-weighted consensus sequence
bases = "ACGT"
consensus = ""

for position in range(50):
    base_counts = {
        base: sum(
            count
            for sequence, count in zip(clone2["N50"], clone2["count"])
            if sequence[position] == base
        )
        for base in bases
    }
    consensus += max(base_counts, key=base_counts.get)

axis_labels = [
    f"{position + 1}\n{base}"
    for position, base in enumerate(consensus)
]

plt.figure(figsize=(12, 10))
plt.imshow(mi, origin="lower", aspect="equal")
plt.colorbar(label="Mutual information (bits)")

plt.xticks(
    range(50),
    axis_labels,
    fontsize=7
)
plt.yticks(
    range(50),
    axis_labels,
    fontsize=7
)

plt.xlabel("Clone 2 consensus position and nucleotide")
plt.ylabel("Clone 2 consensus position and nucleotide")
plt.title("Clone 2 covariation in round 8")

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "clone2_R8_covariation_MI_sequence_axes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Clone 2 consensus sequence:")
print(consensus)


# In[4]:


mask = np.triu(np.ones_like(mi, dtype=bool), k=1)
mi_lower = np.ma.array(mi, mask=mask)

plt.figure(figsize=(12, 10))
plt.imshow(mi_lower, origin="lower", aspect="equal")
plt.colorbar(label="Mutual information (bits)")

plt.xticks(
    range(50),
    axis_labels,
    fontsize=7
)
plt.yticks(
    range(50),
    axis_labels,
    fontsize=7
)

plt.xlabel("Clone 2 consensus position and nucleotide")
plt.ylabel("Clone 2 consensus position and nucleotide")
plt.title("Clone 2 covariation in round 8")

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "clone2_R8_covariation_MI_sequence_axes_lower_triangle.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[5]:


from pathlib import Path
import gzip
import csv

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
OUT = BASE / "N50_extracted_repeats"
OUT.mkdir(parents=True, exist_ok=True)

FLANK = "GGCAGGGCAGCGCAGTC"
N50_LENGTH = 50

samples = {
    "R5_repeat": BASE / "QT5_S4_L001_R1_001.fastq.gz",
    "R6_repeat": BASE / "QT6_S5_L001_R1_001.fastq.gz",
}

for sample_name, fastq in samples.items():
    output_tsv = OUT / f"{sample_name}_N50.tsv"

    total_reads = 0
    flank_found = 0
    extracted = 0

    with gzip.open(fastq, "rt") as fin, open(output_tsv, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["read_id", "N50"])

        while True:
            header = fin.readline().strip()
            if not header:
                break

            sequence = fin.readline().strip().upper()
            fin.readline()
            fin.readline()

            total_reads += 1
            flank_start = sequence.find(FLANK)

            if flank_start == -1:
                continue

            flank_found += 1

            if flank_start >= N50_LENGTH:
                n50 = sequence[flank_start - N50_LENGTH:flank_start]
                writer.writerow([header.removeprefix("@"), n50])
                extracted += 1

    print(
        f"{sample_name}: total={total_reads:,}, "
        f"flank_found={flank_found:,}, "
        f"N50_extracted={extracted:,}"
    )


# In[6]:


from pathlib import Path
import pandas as pd

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
IN_DIR = BASE / "N50_extracted_repeats"
OUT_DIR = BASE / "N50_collapsed_repeats"
OUT_DIR.mkdir(parents=True, exist_ok=True)

for infile in sorted(IN_DIR.glob("*_N50.tsv")):
    df = pd.read_csv(infile, sep="\t")

    collapsed = (
        df["N50"]
        .value_counts()
        .rename_axis("N50")
        .reset_index(name="count")
    )

    total = collapsed["count"].sum()
    collapsed["fraction"] = collapsed["count"] / total if total else 0

    outfile = OUT_DIR / infile.name.replace("_N50.tsv", "_collapsed.tsv")
    collapsed.to_csv(outfile, sep="\t", index=False)

    print(
        f"{infile.name}: reads={total:,}, "
        f"unique_N50={len(collapsed):,}"
    )


# In[7]:


from pathlib import Path
import pandas as pd

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
REPEAT_DIR = BASE / "N50_collapsed_repeats"
CLUSTER_DIR = BASE / "N50_clusters_94"

# Rebuild N50 → cluster mapping
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

sequence_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] in {"S", "H"}:
            sequence_id = fields[8].split(";")[0]
            cluster_id = int(fields[1])

            if sequence_id in id_to_seq:
                sequence_to_cluster[id_to_seq[sequence_id]] = cluster_id

# Assign each repeat to existing clusters
for infile in sorted(REPEAT_DIR.glob("*_collapsed.tsv")):
    df = pd.read_csv(infile, sep="\t")
    df["cluster_id"] = df["N50"].map(sequence_to_cluster)

    assigned_reads = df.loc[df["cluster_id"].notna(), "count"].sum()
    total_reads = df["count"].sum()

    outfile = infile.with_name(
        infile.name.replace("_collapsed.tsv", "_clustered94.tsv")
    )
    df.to_csv(outfile, sep="\t", index=False)

    print(
        f"{infile.name}: assigned={assigned_reads:,}/{total_reads:,} "
        f"({100 * assigned_reads / total_reads:.2f}%)"
    )


# In[8]:


from pathlib import Path
import pandas as pd
import subprocess

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
REPEAT_DIR = BASE / "N50_collapsed_repeats"
CLUSTER_DIR = BASE / "N50_clusters_94"
OUT_DIR = BASE / "N50_repeat_clusters_94"
OUT_DIR.mkdir(exist_ok=True)

for infile in sorted(REPEAT_DIR.glob("*_collapsed.tsv")):
    sample = infile.name.replace("_collapsed.tsv", "")
    df = pd.read_csv(infile, sep="\t")

    fasta = OUT_DIR / f"{sample}.fasta"
    uc_file = OUT_DIR / f"{sample}_assignments.uc"

    with open(fasta, "w") as fh:
        for i, row in df.iterrows():
            fh.write(f">seq_{i};size={int(row['count'])}\n{row['N50']}\n")

    subprocess.run([
        "vsearch",
        "--usearch_global", str(fasta),
        "--db", str(CLUSTER_DIR / "family_centroids.fasta"),
        "--id", "0.94",
        "--strand", "plus",
        "--sizein",
        "--uc", str(uc_file),
        "--notmatched", str(OUT_DIR / f"{sample}_unmatched.fasta")
    ], check=True)

    print(f"{sample}: similarity-based assignment complete")


# In[10]:


from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
REPEAT_DIR = BASE / "N50_repeat_clusters_94"

# Original top-10 clone trajectories
top10 = pd.read_csv(
    CLUSTER_DIR / "top10_clones_share_by_round.tsv",
    sep="\t",
    index_col=0
)

# Map original centroid FASTA labels such as seq_2 -> VSEARCH cluster ID
centroid_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] == "S":
            cluster_id = int(fields[1])
            centroid_name = fields[8].split(";")[0]
            centroid_to_cluster[centroid_name] = cluster_id

# Recover which cluster IDs correspond to Clone 1–10
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

sequence_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] in {"S", "H"}:
            seq_id = fields[8].split(";")[0]

            if seq_id in id_to_seq:
                sequence_to_cluster[id_to_seq[seq_id]] = int(fields[1])

r8_file = next((BASE / "N50_collapsed").glob("*R8*_collapsed.tsv"))
r8 = pd.read_csv(r8_file, sep="\t")
r8["cluster_id"] = r8["N50"].map(sequence_to_cluster)
r8 = r8.dropna(subset=["cluster_id"])
r8["cluster_id"] = r8["cluster_id"].astype(int)

top_cluster_ids = (
    r8.groupby("cluster_id")["count"]
    .sum()
    .sort_values(ascending=False)
    .head(10)
    .index
    .tolist()
)

clone_by_cluster = {
    cluster_id: f"Clone {i}"
    for i, cluster_id in enumerate(top_cluster_ids, start=1)
}

# Calculate repeat-round shares
repeat_files = {
    "R5 repeat": REPEAT_DIR / "R5_repeat_assignments.uc",
    "R6 repeat": REPEAT_DIR / "R6_repeat_assignments.uc",
}

repeat_shares = pd.DataFrame(
    0.0,
    index=[f"Clone {i}" for i in range(1, 11)],
    columns=repeat_files
)

for round_name, uc_file in repeat_files.items():
    counts = {}
    total_assigned = 0

    with open(uc_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")

            if fields[0] != "H":
                continue

            query_label = fields[8]
            centroid_label = fields[9].split(";")[0]

            size = 1
            for part in query_label.split(";"):
                if part.startswith("size="):
                    size = int(part.split("=")[1])
                    break

            cluster_id = centroid_to_cluster.get(centroid_label)

            if cluster_id is not None:
                counts[cluster_id] = counts.get(cluster_id, 0) + size
                total_assigned += size

    for cluster_id, clone_name in clone_by_cluster.items():
        repeat_shares.loc[clone_name, round_name] = (
            100 * counts.get(cluster_id, 0) / total_assigned
            if total_assigned else 0
        )

# Plot repeat rounds together
plt.figure(figsize=(8, 5))

for clone in repeat_shares.index:
    plt.plot(
        [5, 6],
        repeat_shares.loc[clone],
        marker="o",
        label=clone
    )

plt.xticks([5, 6], ["R5 repeat", "R6 repeat"])
plt.xlabel("Selection round")
plt.ylabel("Share of assigned pool (%)")
plt.title("Top 10 clones in repeat rounds")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.show()

# Integrate repeats with original rounds
full = top10.copy()
full["R5 repeat"] = repeat_shares["R5 repeat"]
full["R6 repeat"] = repeat_shares["R6 repeat"]

plot_columns = ["R1", "R4", "R5", "R5 repeat", "R6", "R6 repeat", "R7", "R8"]
x_positions = [1, 4, 5, 5.2, 6, 6.2, 7, 8]

plt.figure(figsize=(11, 6))

for clone in full.index:
    plt.plot(
        x_positions,
        full.loc[clone, plot_columns],
        marker="o",
        label=clone
    )

plt.xticks(
    [1, 2, 3, 4, 5, 5.2, 6, 6.2, 7, 8],
    ["R1", "R2", "R3", "R4", "R5", "R5\nrepeat",
     "R6", "R6\nrepeat", "R7", "R8"]
)

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Top 10 clone trajectories with repeat rounds")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.show()

display(repeat_shares.round(4))


# In[11]:


from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"

FULL_SEQUENCE = (
    "TGATACGAGTTGCAGATATCTCTGTAGAACGTCTTCGCCTTGATTGACTA"
    "GGCAGGGCAGCGCAGTCGG"
)
CENTROID_N50 = FULL_SEQUENCE[:50]

# Rebuild N50 → cluster mapping
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

sequence_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] in {"S", "H"}:
            seq_id = fields[8].split(";")[0]
            if seq_id in id_to_seq:
                sequence_to_cluster[id_to_seq[seq_id]] = int(fields[1])

if CENTROID_N50 not in sequence_to_cluster:
    raise ValueError("The supplied centroid N50 was not found in the 94% clustering.")

clone8_cluster = sequence_to_cluster[CENTROID_N50]

# Extract all round-8 members of this cluster
r8_file = next(COLLAPSED_DIR.glob("*R8*_collapsed.tsv"))
r8 = pd.read_csv(r8_file, sep="\t")
r8["cluster_id"] = r8["N50"].map(sequence_to_cluster)

clone8 = r8[r8["cluster_id"] == clone8_cluster].copy()

sequences = clone8["N50"].tolist()
weights = clone8["count"].to_numpy(dtype=float)
weights /= weights.sum()

# Read-count-weighted mutual information
bases = "ACGT"
mi = np.zeros((50, 50))

for i in range(50):
    for j in range(i + 1, 50):
        joint = {
            (a, b): sum(
                w for seq, w in zip(sequences, weights)
                if seq[i] == a and seq[j] == b
            )
            for a in bases for b in bases
        }

        pi = {a: sum(joint[a, b] for b in bases) for a in bases}
        pj = {b: sum(joint[a, b] for a in bases) for b in bases}

        for a in bases:
            for b in bases:
                p = joint[a, b]
                if p > 0 and pi[a] > 0 and pj[b] > 0:
                    mi[i, j] += p * np.log2(p / (pi[a] * pj[b]))

        mi[j, i] = mi[i, j]

# Keep only one nonredundant triangle
mi_triangle = np.ma.array(
    mi,
    mask=np.triu(np.ones_like(mi, dtype=bool), k=0)
)

axis_labels = [
    f"{position + 1}\n{base}"
    for position, base in enumerate(CENTROID_N50)
]

plt.figure(figsize=(12, 10))
plt.imshow(mi_triangle, origin="lower", aspect="equal")
plt.colorbar(label="Mutual information (bits)")
plt.xticks(range(50), axis_labels, fontsize=7)
plt.yticks(range(50), axis_labels, fontsize=7)
plt.xlabel("Clone 8 centroid position and nucleotide")
plt.ylabel("Clone 8 centroid position and nucleotide")
plt.title("Clone 8 nucleotide covariation in round 8")
plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "clone8_R8_covariation_MI.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Clone 8 cluster ID:", clone8_cluster)
print("Clone 8 centroid N50:", CENTROID_N50)
print("Unique sequences analysed:", len(clone8))
print("Total reads analysed:", int(clone8["count"].sum()))


# In[12]:


from matplotlib.colors import LinearSegmentedColormap

# Display the opposite triangular half
mi_triangle = np.ma.array(
    mi,
    mask=np.tril(np.ones_like(mi, dtype=bool), k=0)
)

blue_orange = LinearSegmentedColormap.from_list(
    "blue_orange",
    ["#2166AC", "#F7F7F7", "#E66101"]
)

plt.figure(figsize=(12, 10))

plt.imshow(
    mi_triangle,
    origin="lower",
    aspect="equal",
    cmap=blue_orange,
    vmin=0,
    vmax=mi.max()
)

plt.colorbar(label="Mutual information (bits)")
plt.xticks(range(50), axis_labels, fontsize=7)
plt.yticks(range(50), axis_labels, fontsize=7)

plt.xlabel("Clone 8 centroid position and nucleotide")
plt.ylabel("Clone 8 centroid position and nucleotide")
plt.title("Clone 8 nucleotide covariation in round 8")

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "clone8_R8_covariation_MI_upper_blue_orange.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[13]:


from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"

# Rebuild sequence-to-cluster mapping
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

sequence_to_cluster = {}
cluster_centroid = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] not in {"S", "H"}:
            continue

        seq_id = fields[8].split(";")[0]
        cluster_id = int(fields[1])

        if seq_id in id_to_seq:
            sequence_to_cluster[id_to_seq[seq_id]] = cluster_id

            if fields[0] == "S":
                cluster_centroid[cluster_id] = id_to_seq[seq_id]

# Load round 8
r8 = pd.read_csv(
    next(COLLAPSED_DIR.glob("*R8*_collapsed.tsv")),
    sep="\t"
)
r8["cluster_id"] = r8["N50"].map(sequence_to_cluster)
r8 = r8.dropna(subset=["cluster_id"])
r8["cluster_id"] = r8["cluster_id"].astype(int)

# Clone 2 = second-most abundant cluster in R8
ranked_clusters = (
    r8.groupby("cluster_id")["count"]
    .sum()
    .sort_values(ascending=False)
)

clone2_cluster = ranked_clusters.index[1]
clone2_centroid = cluster_centroid[clone2_cluster]

# Clone 8 defined by the supplied centroid sequence
clone8_centroid = "TGATACGAGTTGCAGATATCTCTGTAGAACGTCTTCGCCTTGATTGACTA"
clone8_cluster = sequence_to_cluster[clone8_centroid]

bases = "ACGT"

def single_mutation_matrix(cluster_id, centroid):
    family = r8[r8["cluster_id"] == cluster_id]
    total_reads = family["count"].sum()
    matrix = np.zeros((4, 50))

    for sequence, count in zip(family["N50"], family["count"]):
        differences = [
            i for i, (observed, reference) in enumerate(zip(sequence, centroid))
            if observed != reference
        ]

        if len(differences) == 1:
            position = differences[0]
            mutant_base = sequence[position]

            if mutant_base in bases:
                matrix[bases.index(mutant_base), position] += (
                    100 * count / total_reads
                )

    # Hide the centroid nucleotide at each position
    for position, reference_base in enumerate(centroid):
        matrix[bases.index(reference_base), position] = np.nan

    return matrix

for clone_name, cluster_id, centroid in [
    ("Clone 2", clone2_cluster, clone2_centroid),
    ("Clone 8", clone8_cluster, clone8_centroid),
]:
    matrix = single_mutation_matrix(cluster_id, centroid)

    plt.figure(figsize=(14, 3.5))
    image = plt.imshow(matrix, aspect="auto", cmap="Oranges")

    plt.colorbar(image, label="Share of clone-family reads (%)")
    plt.yticks(range(4), list(bases))
    plt.xticks(
        range(50),
        [f"{i + 1}\n{base}" for i, base in enumerate(centroid)],
        fontsize=7
    )

    plt.xlabel("Centroid position and nucleotide")
    plt.ylabel("Mutant nucleotide")
    plt.title(f"{clone_name}: single-nucleotide mutations in round 8")
    plt.tight_layout()

    plt.savefig(
        CLUSTER_DIR / f"{clone_name.lower().replace(' ', '')}_R8_single_mutations.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


# In[14]:


from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
IN_DIR = BASE / "N50_collapsed"
OUT_DIR = BASE / "diversity_analysis"
OUT_DIR.mkdir(exist_ok=True)

round_order = ["R1", "R4", "R5", "R6", "R7", "R8"]
results = []

for infile in IN_DIR.glob("*_collapsed.tsv"):
    match = re.search(r"(?:LUR|^R)(1|4|5|6|7|8)", infile.name)
    if not match:
        continue

    round_name = f"R{match.group(1)}"
    df = pd.read_csv(infile, sep="\t")

    counts = df["count"].to_numpy(dtype=float)
    frequencies = counts / counts.sum()
    nonzero = frequencies[frequencies > 0]

    shannon = -np.sum(nonzero * np.log(nonzero))
    simpson = np.sum(nonzero ** 2)

    results.append({
        "Round": round_name,
        "Total reads": int(counts.sum()),
        "Unique N50s": len(counts),
        "Shannon diversity": np.exp(shannon),
        "Simpson diversity": 1 / simpson
    })

diversity = (
    pd.DataFrame(results)
    .set_index("Round")
    .reindex(round_order)
)

diversity.to_csv(
    OUT_DIR / "sequence_diversity_by_round.tsv",
    sep="\t"
)

for metric, filename in [
    ("Unique N50s", "unique_N50_diversity.png"),
    ("Shannon diversity", "shannon_effective_diversity.png"),
    ("Simpson diversity", "simpson_effective_diversity.png")
]:
    plt.figure(figsize=(7, 4.5))
    plt.plot(
        [1, 4, 5, 6, 7, 8],
        diversity[metric],
        marker="o"
    )
    plt.xticks(range(1, 9), [f"R{i}" for i in range(1, 9)])
    plt.xlabel("Selection round")
    plt.ylabel(metric)
    plt.title(f"{metric} throughout selection")
    plt.tight_layout()
    plt.savefig(OUT_DIR / filename, dpi=300, bbox_inches="tight")
    plt.show()

display(diversity.round(2))


# In[15]:


for metric, filename in [
    ("Unique N50s", "unique_N50_diversity_log.png"),
    ("Shannon diversity", "shannon_effective_diversity_log.png"),
    ("Simpson diversity", "simpson_effective_diversity_log.png")
]:
    plt.figure(figsize=(7, 4.5))

    plt.plot(
        [1, 4, 5, 6, 7, 8],
        diversity[metric],
        marker="o"
    )

    plt.yscale("log")
    plt.xticks(range(1, 9), [f"R{i}" for i in range(1, 9)])
    plt.xlabel("Selection round")
    plt.ylabel(metric)
    plt.title(f"{metric} throughout selection")

    plt.tight_layout()
    plt.savefig(
        OUT_DIR / filename,
        dpi=300,
        bbox_inches="tight"
    )
    plt.show()


# In[17]:


def find_original_file(round_name):
    matches = [
        path for path in ORIGINAL_DIR.glob("*_collapsed.tsv")
        if round_name.lower() in path.stem.lower()
        and "repeat" not in path.stem.lower()
    ]

    if len(matches) != 1:
        print("Available files:")
        for path in sorted(ORIGINAL_DIR.glob("*_collapsed.tsv")):
            print(" ", path.name)

        raise ValueError(
            f"Expected one original {round_name} file, found: "
            f"{[p.name for p in matches]}"
        )

    return matches[0]


# In[19]:


original = {}


# In[20]:


def find_original_file(round_name):
    matches = [
        path for path in ORIGINAL_DIR.glob("*_collapsed.tsv")
        if round_name.lower() in path.stem.lower()
        and "repeat" not in path.stem.lower()
    ]

    if len(matches) != 1:
        print("Available files:")
        for path in sorted(ORIGINAL_DIR.glob("*_collapsed.tsv")):
            print(" ", path.name)

        raise ValueError(
            f"Expected one original {round_name} file, found: "
            f"{[p.name for p in matches]}"
        )

    return matches[0]


# In[22]:


from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
ORIGINAL_DIR = BASE / "N50_collapsed"
REPEAT_DIR = BASE / "N50_repeat_clusters_94"
OUT_DIR = BASE / "enrichment_top100"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Build mappings for the existing 94% sequence-family clusters
# ------------------------------------------------------------

seq_info = pd.read_csv(
    CLUSTER_DIR / "sequence_ids.tsv",
    sep="\t"
)

id_to_seq = dict(
    zip(seq_info["sequence_id"], seq_info["N50"])
)

sequence_to_cluster = {}
centroid_name_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] not in {"S", "H"}:
            continue

        cluster_id = int(fields[1])
        sequence_id = fields[8].split(";")[0]

        if sequence_id in id_to_seq:
            sequence_to_cluster[id_to_seq[sequence_id]] = cluster_id

        if fields[0] == "S":
            centroid_name_to_cluster[sequence_id] = cluster_id

# ------------------------------------------------------------
# Find and load original selection-round files
# ------------------------------------------------------------

def find_original_file(round_name):
    round_number = round_name.removeprefix("R")

    patterns = [
        rf"(^|[_-])R{round_number}([_-]|$)",
        rf"(^|[_-])LUR{round_number}([_-]|$)"
    ]

    matches = []

    for path in ORIGINAL_DIR.glob("*_collapsed.tsv"):
        filename = path.name

        if "repeat" in filename.lower():
            continue

        if any(
            re.search(pattern, filename, flags=re.IGNORECASE)
            for pattern in patterns
        ):
            matches.append(path)

    if len(matches) != 1:
        print("\nAvailable collapsed files:")
        for path in sorted(ORIGINAL_DIR.glob("*_collapsed.tsv")):
            print(" ", path.name)

        raise ValueError(
            f"\nExpected exactly one original {round_name} file, "
            f"but found {len(matches)}: {[p.name for p in matches]}"
        )

    return matches[0]


def load_original_counts(round_name):
    infile = find_original_file(round_name)
    df = pd.read_csv(infile, sep="\t")

    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"]).copy()
    df["cluster_id"] = df["cluster_id"].astype(int)

    counts = df.groupby("cluster_id")["count"].sum()
    total = int(df["count"].sum())

    print(
        f"{round_name}: {infile.name} | "
        f"assigned reads={total:,} | "
        f"clusters={len(counts):,}"
    )

    return counts, total


original = {
    round_name: load_original_counts(round_name)
    for round_name in ["R4", "R5", "R6"]
}

# ------------------------------------------------------------
# Load repeat assignments from VSEARCH UC files
# ------------------------------------------------------------

def parse_size(label):
    match = re.search(
        r"(?:^|;)size=(\d+)(?:;|$)",
        label
    )
    return int(match.group(1)) if match else 1


def load_repeat_counts(uc_file):
    if not uc_file.exists():
        raise FileNotFoundError(f"Missing repeat UC file: {uc_file}")

    counts = {}
    assigned_total = 0

    with open(uc_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")

            # usearch_global matches are H records
            if fields[0] != "H":
                continue

            query_label = fields[8]
            centroid_label = fields[9].split(";")[0]

            cluster_id = centroid_name_to_cluster.get(centroid_label)

            if cluster_id is None:
                continue

            read_count = parse_size(query_label)

            counts[cluster_id] = (
                counts.get(cluster_id, 0) + read_count
            )
            assigned_total += read_count

    return pd.Series(counts, dtype=float), assigned_total


repeat = {
    "R5_repeat": load_repeat_counts(
        REPEAT_DIR / "R5_repeat_assignments.uc"
    ),
    "R6_repeat": load_repeat_counts(
        REPEAT_DIR / "R6_repeat_assignments.uc"
    )
}

print(
    f"R5 repeat assigned reads: {repeat['R5_repeat'][1]:,}"
)
print(
    f"R6 repeat assigned reads: {repeat['R6_repeat'][1]:,}"
)

# ------------------------------------------------------------
# Select the 100 most abundant clusters in original R6
# ------------------------------------------------------------

top100 = (
    original["R6"][0]
    .sort_values(ascending=False)
    .head(100)
    .index
)

results = pd.DataFrame(index=top100)
results.index.name = "cluster_id"
results["R6_rank"] = np.arange(1, len(results) + 1)

# Half-read pseudocount prevents infinite enrichment values
PSEUDOCOUNT = 0.5


def calculate_frequencies(counts, total, cluster_ids):
    adjusted_counts = (
        counts.reindex(cluster_ids, fill_value=0)
        + PSEUDOCOUNT
    )

    adjusted_total = (
        total + PSEUDOCOUNT * len(cluster_ids)
    )

    return adjusted_counts / adjusted_total

# Original rounds
for sample in ["R4", "R5", "R6"]:
    counts, total = original[sample]

    results[f"{sample}_count"] = counts.reindex(
        top100,
        fill_value=0
    )

    results[f"{sample}_fraction"] = calculate_frequencies(
        counts,
        total,
        top100
    )

# Repeat rounds
for sample in ["R5_repeat", "R6_repeat"]:
    counts, total = repeat[sample]

    results[f"{sample}_count"] = counts.reindex(
        top100,
        fill_value=0
    )

    results[f"{sample}_fraction"] = calculate_frequencies(
        counts,
        total,
        top100
    )

# ------------------------------------------------------------
# Calculate log2 enrichment
# ------------------------------------------------------------

results["log2_enrichment_R4_to_R5"] = np.log2(
    results["R5_fraction"]
    / results["R4_fraction"]
)

results["log2_enrichment_R5_to_R6"] = np.log2(
    results["R6_fraction"]
    / results["R5_fraction"]
)

# R5 repeat is compared with the same original R4 input
results["log2_enrichment_R4_to_R5_repeat"] = np.log2(
    results["R5_repeat_fraction"]
    / results["R4_fraction"]
)

results["log2_enrichment_R5repeat_to_R6repeat"] = np.log2(
    results["R6_repeat_fraction"]
    / results["R5_repeat_fraction"]
)

results.to_csv(
    OUT_DIR / "top100_cluster_enrichment_original_vs_repeats.tsv",
    sep="\t"
)

# ------------------------------------------------------------
# Scatter plots: original enrichment versus repeat enrichment
# ------------------------------------------------------------

comparisons = [
    (
        "log2_enrichment_R4_to_R5",
        "log2_enrichment_R4_to_R5_repeat",
        "R4_to_R5",
        "R4 → R5"
    ),
    (
        "log2_enrichment_R5_to_R6",
        "log2_enrichment_R5repeat_to_R6repeat",
        "R5_to_R6",
        "R5 → R6"
    )
]

for original_col, repeat_col, filename_label, title_label in comparisons:
    plot_data = results[
        [original_col, repeat_col]
    ].replace([np.inf, -np.inf], np.nan).dropna()

    x = plot_data[original_col]
    y = plot_data[repeat_col]

    pearson_r = x.corr(y, method="pearson")
    spearman_r = x.corr(y, method="spearman")

    axis_limit = max(
        abs(x).max(),
        abs(y).max()
    ) * 1.05

    if not np.isfinite(axis_limit) or axis_limit == 0:
        axis_limit = 1

    plt.figure(figsize=(6, 6))

    plt.scatter(
        x,
        y,
        alpha=0.75
    )

    plt.axhline(
        0,
        linewidth=0.8
    )
    plt.axvline(
        0,
        linewidth=0.8
    )

    plt.plot(
        [-axis_limit, axis_limit],
        [-axis_limit, axis_limit],
        linestyle="--",
        linewidth=1
    )

    plt.xlim(-axis_limit, axis_limit)
    plt.ylim(-axis_limit, axis_limit)

    plt.xlabel(
        f"Original {title_label} log₂ enrichment"
    )
    plt.ylabel(
        f"Repeat {title_label} log₂ enrichment"
    )

    plt.title(
        f"Top 100 cluster enrichment: {title_label}\n"
        f"Pearson r = {pearson_r:.3f}; "
        f"Spearman ρ = {spearman_r:.3f}"
    )

    plt.tight_layout()

    plt.savefig(
        OUT_DIR
        / f"top100_{filename_label}_original_vs_repeat.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

display(
    results[
        [
            "R6_rank",
            "R4_count",
            "R5_count",
            "R5_repeat_count",
            "R6_count",
            "R6_repeat_count",
            "log2_enrichment_R4_to_R5",
            "log2_enrichment_R4_to_R5_repeat",
            "log2_enrichment_R5_to_R6",
            "log2_enrichment_R5repeat_to_R6repeat"
        ]
    ].round(3)
)


# In[23]:


# Calculate how much more strongly each cluster enriched in repeats

results["repeat_advantage_R4_to_R5"] = (
    results["log2_enrichment_R4_to_R5_repeat"]
    - results["log2_enrichment_R4_to_R5"]
)

results["repeat_advantage_R5_to_R6"] = (
    results["log2_enrichment_R5repeat_to_R6repeat"]
    - results["log2_enrichment_R5_to_R6"]
)

# Combined score across both transitions
results["combined_repeat_advantage"] = (
    results["repeat_advantage_R4_to_R5"]
    + results["repeat_advantage_R5_to_R6"]
)

# Keep clusters that were more enriched in repeats in at least one transition
repeat_enriched = results[
    (results["repeat_advantage_R4_to_R5"] > 0)
    | (results["repeat_advantage_R5_to_R6"] > 0)
].copy()

repeat_enriched = repeat_enriched.sort_values(
    "combined_repeat_advantage",
    ascending=False
)

repeat_enriched.to_csv(
    OUT_DIR / "clusters_more_enriched_in_repeats.tsv",
    sep="\t"
)

display(
    repeat_enriched[
        [
            "R6_rank",
            "log2_enrichment_R4_to_R5",
            "log2_enrichment_R4_to_R5_repeat",
            "repeat_advantage_R4_to_R5",
            "log2_enrichment_R5_to_R6",
            "log2_enrichment_R5repeat_to_R6repeat",
            "repeat_advantage_R5_to_R6",
            "combined_repeat_advantage"
        ]
    ].head(30).round(3)
)


# In[24]:


import numpy as np
import matplotlib.pyplot as plt

# Combine rounds 5 and 6 into one R4 → R6 enrichment measurement
results["original_R4_to_R6_enrichment"] = (
    results["log2_enrichment_R4_to_R5"]
    + results["log2_enrichment_R5_to_R6"]
)

results["repeat_R4_to_R6_enrichment"] = (
    results["log2_enrichment_R4_to_R5_repeat"]
    + results["log2_enrichment_R5repeat_to_R6repeat"]
)

# Positive score = enriched more strongly in repeat selection
results["unified_repeat_advantage"] = (
    results["repeat_R4_to_R6_enrichment"]
    - results["original_R4_to_R6_enrichment"]
)

# Select the 50 clusters with the strongest repeat advantage
top50_advantage = (
    results
    .replace([np.inf, -np.inf], np.nan)
    .dropna(subset=["unified_repeat_advantage"])
    .sort_values("unified_repeat_advantage", ascending=False)
    .head(50)
    .copy()
)

top50_advantage["cluster_label"] = [
    f"Cluster {cluster_id}"
    for cluster_id in top50_advantage.index
]

# Sort ascending so the strongest cluster appears at the top
plot_data = top50_advantage.sort_values(
    "unified_repeat_advantage",
    ascending=True
)

plt.figure(figsize=(9, 12))

plt.barh(
    plot_data["cluster_label"],
    plot_data["unified_repeat_advantage"]
)

plt.axvline(0, linewidth=1)

plt.xlabel(
    "Repeat advantage score\n"
    "(repeat R4→R6 log₂ enrichment − original R4→R6 log₂ enrichment)"
)
plt.ylabel("Sequence cluster")
plt.title("Top 50 clusters enriched more strongly in repeat selection")

plt.tight_layout()

plt.savefig(
    OUT_DIR / "top50_unified_repeat_advantage_R4_to_R6.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

top50_advantage.to_csv(
    OUT_DIR / "top50_unified_repeat_advantage_R4_to_R6.tsv",
    sep="\t"
)

display(
    top50_advantage[
        [
            "R6_rank",
            "original_R4_to_R6_enrichment",
            "repeat_R4_to_R6_enrichment",
            "unified_repeat_advantage"
        ]
    ].round(3)
)


# In[25]:


import pandas as pd

# Rebuild cluster ID → centroid N50 mapping
cluster_to_centroid = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] != "S":
            continue

        cluster_id = int(fields[1])
        centroid_id = fields[8].split(";")[0]

        if centroid_id in id_to_seq:
            cluster_to_centroid[cluster_id] = id_to_seq[centroid_id]

# Keep all clusters with positive repeat advantage
positive_advantage = (
    results[results["unified_repeat_advantage"] > 0]
    .copy()
    .sort_values("unified_repeat_advantage", ascending=False)
)

positive_advantage["centroid_N50"] = (
    positive_advantage.index.map(cluster_to_centroid)
)

output_columns = [
    "centroid_N50",
    "R6_rank",
    "original_R4_to_R6_enrichment",
    "repeat_R4_to_R6_enrichment",
    "unified_repeat_advantage"
]

positive_sequences = positive_advantage[output_columns]

positive_sequences.to_csv(
    OUT_DIR / "positive_repeat_advantage_cluster_N50_sequences.tsv",
    sep="\t",
    index=True
)

display(positive_sequences.round(3))

print(
    f"Positive-repeat-advantage clusters: "
    f"{len(positive_sequences):,}"
)


# In[26]:


from pathlib import Path
import pandas as pd

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")

files_to_check = {
    "Extracted R8": next((BASE / "N50_extracted").glob("*R8*_N50.tsv")),
    "Collapsed R8": next((BASE / "N50_collapsed").glob("*R8*_collapsed.tsv")),
    "Cluster sequence IDs": BASE / "N50_clusters_94" / "sequence_ids.tsv",
}

for label, path in files_to_check.items():
    df = pd.read_csv(path, sep="\t")

    sequence_column = "N50"
    lengths = df[sequence_column].astype(str).str.strip().str.len()

    print(f"\n{label}")
    print("File:", path.name)
    print("Length counts:")
    print(lengths.value_counts().sort_index())
    print("Example:", df[sequence_column].iloc[0])


# In[27]:


import pandas as pd

pd.set_option("display.max_colwidth", None)

positive_sequences = pd.read_csv(
    OUT_DIR / "positive_repeat_advantage_cluster_N50_sequences.tsv",
    sep="\t"
)

positive_sequences["N50_length"] = (
    positive_sequences["centroid_N50"]
    .astype(str)
    .str.strip()
    .str.len()
)

print("Length distribution:")
print(positive_sequences["N50_length"].value_counts().sort_index())

display(
    positive_sequences[
        [
            "cluster_id",
            "centroid_N50",
            "N50_length",
            "unified_repeat_advantage"
        ]
    ].head(20)
)


# In[28]:


# Count unique N50 sequence members in every 94% cluster
cluster_member_counts = (
    pd.Series(sequence_to_cluster, name="cluster_id")
    .rename_axis("N50")
    .reset_index()
    .groupby("cluster_id")
    .size()
    .rename("unique_N50_members")
)

positive_sequences["unique_N50_members"] = (
    positive_sequences["cluster_id"]
    .map(cluster_member_counts)
)

positive_sequences = positive_sequences.sort_values(
    "unified_repeat_advantage",
    ascending=False
)

positive_sequences.to_csv(
    OUT_DIR / "positive_repeat_advantage_clusters_with_member_counts.tsv",
    sep="\t",
    index=False
)

display(
    positive_sequences[
        [
            "cluster_id",
            "centroid_N50",
            "unique_N50_members",
            "unified_repeat_advantage"
        ]
    ]
)


# In[29]:


from pathlib import Path
import re
import pandas as pd

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
COLLAPSED_DIR = BASE / "N50_collapsed"

round_order = ["R1", "R4", "R5", "R6", "R7", "R8"]
read_counts_by_cluster = pd.DataFrame()

for infile in COLLAPSED_DIR.glob("*_collapsed.tsv"):
    match = re.search(r"(?:LUR|^R)(1|4|5|6|7|8)", infile.name)
    if not match:
        continue

    round_name = f"R{match.group(1)}"
    df = pd.read_csv(infile, sep="\t")

    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"]).copy()
    df["cluster_id"] = df["cluster_id"].astype(int)

    read_counts_by_cluster[round_name] = (
        df.groupby("cluster_id")["count"]
        .sum()
    )

read_counts_by_cluster = (
    read_counts_by_cluster
    .fillna(0)
    .astype(int)
    .reindex(columns=round_order)
)

read_counts_by_cluster["total_original_reads"] = (
    read_counts_by_cluster.sum(axis=1)
)

positive_sequences = positive_sequences.drop(
    columns=round_order + ["total_original_reads"],
    errors="ignore"
)

positive_sequences = positive_sequences.merge(
    read_counts_by_cluster.reset_index(),
    on="cluster_id",
    how="left"
)

positive_sequences.to_csv(
    OUT_DIR / "positive_repeat_advantage_clusters_with_read_counts.tsv",
    sep="\t",
    index=False
)

display(
    positive_sequences[
        [
            "cluster_id",
            "centroid_N50",
            "unique_N50_members",
            *round_order,
            "total_original_reads",
            "unified_repeat_advantage"
        ]
    ].sort_values(
        "unified_repeat_advantage",
        ascending=False
    )
)


# In[30]:


import matplotlib.pyplot as plt
import numpy as np

# Assumes:
# top10 contains original rounds with rows = Clone 1...Clone 10
# repeat_shares contains columns "R5 repeat" and "R6 repeat"

plt.figure(figsize=(11, 7))

colors = plt.cm.tab10(np.linspace(0, 1, 10))

for i, clone in enumerate(top10.index):
    color = colors[i]

    # Shared trajectory up to R4
    plt.plot(
        [1, 4],
        [top10.loc[clone, "R1"], top10.loc[clone, "R4"]],
        marker="o",
        color=color,
        linewidth=2,
        label=clone
    )

    # Original branch: R4 -> R5 -> R6 -> R7 -> R8
    plt.plot(
        [4, 5, 6, 7, 8],
        [
            top10.loc[clone, "R4"],
            top10.loc[clone, "R5"],
            top10.loc[clone, "R6"],
            top10.loc[clone, "R7"],
            top10.loc[clone, "R8"]
        ],
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5
    )

    # Repeat branch: R4 -> R5 repeat -> R6 repeat
    plt.plot(
        [4, 5, 6],
        [
            top10.loc[clone, "R4"],
            repeat_shares.loc[clone, "R5 repeat"],
            repeat_shares.loc[clone, "R6 repeat"]
        ],
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2
    )

plt.xticks(
    range(1, 9),
    [f"R{i}" for i in range(1, 9)]
)

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Top 10 clone trajectories: original vs repeat branches")

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top10_clone_trajectories_split_original_repeat.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[31]:


import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(8, 6))

colors = plt.cm.tab10(np.linspace(0, 1, 10))

for i, clone in enumerate(top10.index):
    color = colors[i]

    # Original branch: R4 -> R5 -> R6
    plt.plot(
        [4, 5, 6],
        [
            top10.loc[clone, "R4"],
            top10.loc[clone, "R5"],
            top10.loc[clone, "R6"]
        ],
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5,
        label=clone
    )

    # Repeat branch: R4 -> R5 repeat -> R6 repeat
    plt.plot(
        [4, 5, 6],
        [
            top10.loc[clone, "R4"],
            repeat_shares.loc[clone, "R5 repeat"],
            repeat_shares.loc[clone, "R6 repeat"]
        ],
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2
    )

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Top 10 clone trajectories: R4–R6 original vs repeat")

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top10_clone_trajectories_R4_R6_split_original_repeat.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[32]:


import matplotlib.pyplot as plt
import numpy as np

exclude = {"Clone 1", "Clone 2", "Clone 5"}
clones_to_plot = [c for c in top10.index if c not in exclude]

plt.figure(figsize=(8, 6))

colors = plt.cm.tab10(np.linspace(0, 1, 10))
color_map = dict(zip(top10.index, colors))

all_values = []

for clone in clones_to_plot:
    color = color_map[clone]

    original_values = [
        top10.loc[clone, "R4"],
        top10.loc[clone, "R5"],
        top10.loc[clone, "R6"]
    ]

    repeat_values = [
        top10.loc[clone, "R4"],
        repeat_shares.loc[clone, "R5 repeat"],
        repeat_shares.loc[clone, "R6 repeat"]
    ]

    all_values.extend(original_values)
    all_values.extend(repeat_values)

    # Original branch
    plt.plot(
        [4, 5, 6],
        original_values,
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5,
        label=clone
    )

    # Repeat branch
    plt.plot(
        [4, 5, 6],
        repeat_values,
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2
    )

# Automatically tighten y-axis around remaining clones
ymax = max(all_values)
plt.ylim(0, ymax * 1.08)

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Clone trajectories: R4–R6 original vs repeat")

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "clone_trajectories_R4_R6_excluding_1_2_5.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[33]:


import matplotlib.pyplot as plt
import numpy as np

# Choose the full original-round share table if available.
# Falls back to top10 only if that is all that exists.
if "all_cluster_shares" in globals():
    original_shares = all_cluster_shares
elif "cluster_shares" in globals():
    original_shares = cluster_shares
else:
    original_shares = top10

# Top 20 clusters by R6 repeat abundance
top20_repeat = (
    repeat_shares["R6 repeat"]
    .sort_values(ascending=False)
    .head(20)
    .index
)

# Keep only clusters present in both datasets
clusters_to_plot = [
    c for c in top20_repeat
    if c in original_shares.index and c in repeat_shares.index
]

if len(clusters_to_plot) < 20:
    print(
        f"Warning: only {len(clusters_to_plot)} of the top 20 R6-repeat "
        "clusters are present in the original-round share table."
    )

plt.figure(figsize=(10, 7))

colors = plt.cm.tab20(np.linspace(0, 1, len(clusters_to_plot)))

all_values = []

for rank, (clone, color) in enumerate(
    zip(clusters_to_plot, colors),
    start=1
):
    original_values = [
        original_shares.loc[clone, "R4"],
        original_shares.loc[clone, "R5"],
        original_shares.loc[clone, "R6"]
    ]

    repeat_values = [
        original_shares.loc[clone, "R4"],
        repeat_shares.loc[clone, "R5 repeat"],
        repeat_shares.loc[clone, "R6 repeat"]
    ]

    all_values.extend(original_values)
    all_values.extend(repeat_values)

    # Original branch — 50% transparent
    plt.plot(
        [4, 5, 6],
        original_values,
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5
    )

    # Repeat branch
    plt.plot(
        [4, 5, 6],
        repeat_values,
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2,
        label=f"{rank}. {clone}"
    )

# Adjust y-axis to displayed clusters
plt.ylim(0, max(all_values) * 1.08)

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Top 20 clusters in R6 repeat: original vs repeat trajectories")

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=8,
    title="R6 repeat rank"
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top20_R6_repeat_cluster_trajectories_R4_R6.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Top 20 clusters ranked by R6 repeat abundance:")
for rank, clone in enumerate(top20_repeat, start=1):
    print(
        f"{rank:2d}. {clone}: "
        f"{repeat_shares.loc[clone, 'R6 repeat']:.4f}%"
    )


# In[34]:


import matplotlib.pyplot as plt
import numpy as np

exclude = {"Clone 1", "Clone 2", "Clone 5"}

# Choose the full original-round share table if available
if "all_cluster_shares" in globals():
    original_shares = all_cluster_shares
elif "cluster_shares" in globals():
    original_shares = cluster_shares
else:
    original_shares = top10

# Rank by R6 repeat abundance, then exclude clones 1, 2, and 5
top20_repeat = (
    repeat_shares["R6 repeat"]
    .sort_values(ascending=False)
    .index
)

clusters_to_plot = [
    c for c in top20_repeat
    if c not in exclude
    and c in original_shares.index
    and c in repeat_shares.index
][:20]

plt.figure(figsize=(10, 7))

colors = plt.cm.tab20(np.linspace(0, 1, len(clusters_to_plot)))
all_values = []

for rank, (clone, color) in enumerate(zip(clusters_to_plot, colors), start=1):

    original_values = [
        original_shares.loc[clone, "R4"],
        original_shares.loc[clone, "R5"],
        original_shares.loc[clone, "R6"]
    ]

    repeat_values = [
        original_shares.loc[clone, "R4"],
        repeat_shares.loc[clone, "R5 repeat"],
        repeat_shares.loc[clone, "R6 repeat"]
    ]

    all_values.extend(original_values)
    all_values.extend(repeat_values)

    # Original branch — 50% transparent
    plt.plot(
        [4, 5, 6],
        original_values,
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5
    )

    # Repeat branch
    plt.plot(
        [4, 5, 6],
        repeat_values,
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2,
        label=f"{clone}"
    )

plt.ylim(0, max(all_values) * 1.08)

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Top R6-repeat clusters excluding Clones 1, 2, and 5")

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=8
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "top20_R6_repeat_excluding_clones_1_2_5.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Clusters plotted:")
for i, clone in enumerate(clusters_to_plot, start=1):
    print(
        f"{i:2d}. {clone}: "
        f"{repeat_shares.loc[clone, 'R6 repeat']:.4f}% in R6 repeat"
    )


# In[35]:


import matplotlib.pyplot as plt
import numpy as np

exclude = {"Clone 1", "Clone 2", "Clone 5"}

# Full original-round share table
if "all_cluster_shares" in globals():
    original_shares = all_cluster_shares
elif "cluster_shares" in globals():
    original_shares = cluster_shares
else:
    original_shares = top10


# ------------------------------------------------------------
# 1. Top 20 R6-repeat clusters, excluding clones 1, 2, and 5
# ------------------------------------------------------------

top20_repeat = [
    c for c in repeat_shares["R6 repeat"]
    .sort_values(ascending=False)
    .index
    if c not in exclude
    and c in original_shares.index
][:20]


# ------------------------------------------------------------
# 2. All clusters with positive unified repeat advantage
# ------------------------------------------------------------

positive_clusters = results.loc[
    results["unified_repeat_advantage"] > 0
].index.tolist()

positive_clusters = [
    c for c in positive_clusters
    if c not in exclude
    and c in original_shares.index
    and c in repeat_shares.index
]


# ------------------------------------------------------------
# 3. Union: top 20 + every positive-advantage cluster
# ------------------------------------------------------------

clusters_to_plot = list(top20_repeat)

for cluster in positive_clusters:
    if cluster not in clusters_to_plot:
        clusters_to_plot.append(cluster)


print(f"Top-20 R6-repeat clusters:             {len(top20_repeat)}")
print(f"Positive repeat-advantage clusters:    {len(positive_clusters)}")
print(f"Total unique clusters shown:           {len(clusters_to_plot)}")


# ------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------

plt.figure(figsize=(12, 8))

colors = plt.cm.nipy_spectral(
    np.linspace(0, 1, len(clusters_to_plot))
)

all_values = []

for cluster, color in zip(clusters_to_plot, colors):

    original_values = [
        original_shares.loc[cluster, "R4"],
        original_shares.loc[cluster, "R5"],
        original_shares.loc[cluster, "R6"]
    ]

    repeat_values = [
        original_shares.loc[cluster, "R4"],
        repeat_shares.loc[cluster, "R5 repeat"],
        repeat_shares.loc[cluster, "R6 repeat"]
    ]

    all_values.extend(original_values)
    all_values.extend(repeat_values)

    is_positive = cluster in positive_clusters

    # Original trajectory — 50% transparent
    plt.plot(
        [4, 5, 6],
        original_values,
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5
    )

    # Repeat trajectory
    plt.plot(
        [4, 5, 6],
        repeat_values,
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2.5 if is_positive else 2,
        label=f"{cluster}*" if is_positive else str(cluster)
    )


# Automatically rescale y-axis
plt.ylim(0, max(all_values) * 1.08)

plt.xticks(
    [4, 5, 6],
    ["R4", "R5", "R6"]
)

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")

plt.title(
    "R4–R6 trajectories: top R6-repeat clusters + positive repeat advantage"
)

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=7,
    title="* Positive repeat advantage"
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR /
    "R4_R6_top20_plus_all_positive_repeat_advantage.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[36]:


print("RESULTS")
print("Index example:", list(results.index[:10]))
print("Index type:", type(results.index[0]))
print("Columns:")
print(list(results.columns))

print("\nREPEAT_SHARES")
print("Index example:", list(repeat_shares.index[:10]))
print("Index type:", type(repeat_shares.index[0]))
print("Columns:")
print(list(repeat_shares.columns))

if "all_cluster_shares" in globals():
    print("\nALL_CLUSTER_SHARES")
    print("Index example:", list(all_cluster_shares.index[:10]))
    print("Index type:", type(all_cluster_shares.index[0]))
    print("Columns:")
    print(list(all_cluster_shares.columns))

elif "cluster_shares" in globals():
    print("\nCLUSTER_SHARES")
    print("Index example:", list(cluster_shares.index[:10]))
    print("Index type:", type(cluster_shares.index[0]))
    print("Columns:")
    print(list(cluster_shares.columns))

print("\nPOSITIVE REPEAT-ADVANTAGE CLUSTERS")
positive = results.loc[results["unified_repeat_advantage"] > 0]

print("Number:", len(positive))
print("IDs:")
print(list(positive.index))


# In[37]:


from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"
REPEAT_DIR = BASE / "N50_repeat_clusters_94"

# ---------- helpers ----------
def parse_size(label):
    m = re.search(r"(?:^|;)size=(\d+)(?:;|$)", label)
    return int(m.group(1)) if m else 1

def find_original_file(round_name):
    round_number = round_name.removeprefix("R")
    patterns = [
        rf"(^|[_-])R{round_number}([_-]|$)",
        rf"(^|[_-])LUR{round_number}([_-]|$)"
    ]
    matches = []
    for path in COLLAPSED_DIR.glob("*_collapsed.tsv"):
        if "repeat" in path.name.lower():
            continue
        if any(re.search(p, path.name, flags=re.IGNORECASE) for p in patterns):
            matches.append(path)
    if len(matches) != 1:
        raise ValueError(f"{round_name}: expected 1 file, found {[p.name for p in matches]}")
    return matches[0]

# ---------- rebuild cluster mappings ----------
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

sequence_to_cluster = {}
centroid_name_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] not in {"S", "H"}:
            continue
        cluster_id = int(fields[1])
        seq_id = fields[8].split(";")[0]

        if seq_id in id_to_seq:
            sequence_to_cluster[id_to_seq[seq_id]] = cluster_id

        if fields[0] == "S":
            centroid_name_to_cluster[seq_id] = cluster_id

# ---------- map cluster IDs to Clone 1..10 from original R8 ----------
r8 = pd.read_csv(find_original_file("R8"), sep="\t")
r8["cluster_id"] = r8["N50"].map(sequence_to_cluster)
r8 = r8.dropna(subset=["cluster_id"]).copy()
r8["cluster_id"] = r8["cluster_id"].astype(int)

top10_r8_cluster_ids = (
    r8.groupby("cluster_id")["count"]
    .sum()
    .sort_values(ascending=False)
    .head(10)
    .index
    .tolist()
)

cluster_to_clone = {
    cluster_id: f"Clone {i}"
    for i, cluster_id in enumerate(top10_r8_cluster_ids, start=1)
}

exclude_ids = {
    cid for cid, name in cluster_to_clone.items()
    if name in {"Clone 1", "Clone 2", "Clone 5"}
}

# ---------- full original shares for R4, R5, R6 ----------
original_counts = {}
original_totals = {}

for round_name in ["R4", "R5", "R6"]:
    df = pd.read_csv(find_original_file(round_name), sep="\t")
    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"]).copy()
    df["cluster_id"] = df["cluster_id"].astype(int)

    original_counts[round_name] = df.groupby("cluster_id")["count"].sum()
    original_totals[round_name] = df["count"].sum()

original_shares_num = pd.DataFrame({
    round_name: 100 * original_counts[round_name] / original_totals[round_name]
    for round_name in ["R4", "R5", "R6"]
}).fillna(0)

# ---------- full repeat shares for R5 repeat, R6 repeat ----------
def load_repeat_share(uc_file):
    counts = {}
    total = 0
    with open(uc_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if fields[0] != "H":
                continue
            query_label = fields[8]
            centroid_label = fields[9].split(";")[0]
            cluster_id = centroid_name_to_cluster.get(centroid_label)
            if cluster_id is None:
                continue
            n = parse_size(query_label)
            counts[cluster_id] = counts.get(cluster_id, 0) + n
            total += n
    return pd.Series(counts, dtype=float) * 100 / total

repeat_shares_num = pd.DataFrame({
    "R5 repeat": load_repeat_share(REPEAT_DIR / "R5_repeat_assignments.uc"),
    "R6 repeat": load_repeat_share(REPEAT_DIR / "R6_repeat_assignments.uc"),
}).fillna(0)

# ---------- choose clusters to plot ----------
top20_repeat_ids = [
    cid for cid in repeat_shares_num["R6 repeat"].sort_values(ascending=False).index
    if cid not in exclude_ids
][:20]

positive_ids = [
    int(cid) for cid in results.index
    if results.loc[cid, "unified_repeat_advantage"] > 0 and cid not in exclude_ids
]

cluster_ids_to_plot = []
for cid in top20_repeat_ids + positive_ids:
    if cid not in cluster_ids_to_plot:
        cluster_ids_to_plot.append(cid)

# keep only clusters with all needed values available
cluster_ids_to_plot = [
    cid for cid in cluster_ids_to_plot
    if cid in original_shares_num.index and cid in repeat_shares_num.index
]

print(f"Top 20 by R6 repeat (after exclusions): {len(top20_repeat_ids)}")
print(f"Positive repeat-advantage clusters:     {len(positive_ids)}")
print(f"Total clusters plotted:                 {len(cluster_ids_to_plot)}")

# ---------- labels ----------
def cluster_label(cid):
    base = cluster_to_clone.get(cid, f"Cluster {cid}")
    if cid in positive_ids:
        base += " *"
    return base

# ---------- plot ----------
plt.figure(figsize=(12, 9))
colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(cluster_ids_to_plot)))
all_values = []

for cid, color in zip(cluster_ids_to_plot, colors):
    original_values = [
        original_shares_num.loc[cid, "R4"],
        original_shares_num.loc[cid, "R5"],
        original_shares_num.loc[cid, "R6"],
    ]
    repeat_values = [
        original_shares_num.loc[cid, "R4"],
        repeat_shares_num.loc[cid, "R5 repeat"],
        repeat_shares_num.loc[cid, "R6 repeat"],
    ]

    all_values.extend(original_values + repeat_values)

    # original branch
    plt.plot(
        [4, 5, 6],
        original_values,
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5
    )

    # repeat branch
    plt.plot(
        [4, 5, 6],
        repeat_values,
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2.5 if cid in positive_ids else 2,
        label=cluster_label(cid)
    )

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])
plt.ylim(0, max(all_values) * 1.08)
plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("R4–R6 trajectories: top R6-repeat clusters + positive repeat advantage")
plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=7,
    title="* positive repeat advantage"
)
plt.tight_layout()
plt.savefig(
    CLUSTER_DIR / "R4_R6_top20_repeat_plus_positive_advantage_FIXED.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()


# In[39]:


from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"
REPEAT_DIR = BASE / "N50_repeat_clusters_94"

# ---------- helpers ----------
def parse_size(label):
    m = re.search(r"(?:^|;)size=(\d+)(?:;|$)", label)
    return int(m.group(1)) if m else 1

def find_original_file(round_name):
    round_number = round_name.removeprefix("R")
    patterns = [
        rf"(^|[_-])R{round_number}([_-]|$)",
        rf"(^|[_-])LUR{round_number}([_-]|$)"
    ]
    matches = []
    for path in COLLAPSED_DIR.glob("*_collapsed.tsv"):
        if "repeat" in path.name.lower():
            continue
        if any(re.search(p, path.name, flags=re.IGNORECASE) for p in patterns):
            matches.append(path)
    if len(matches) != 1:
        raise ValueError(f"{round_name}: expected 1 file, found {[p.name for p in matches]}")
    return matches[0]

# ---------- rebuild cluster mappings ----------
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

sequence_to_cluster = {}
centroid_name_to_cluster = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] not in {"S", "H"}:
            continue
        cluster_id = int(fields[1])
        seq_id = fields[8].split(";")[0]

        if seq_id in id_to_seq:
            sequence_to_cluster[id_to_seq[seq_id]] = cluster_id

        if fields[0] == "S":
            centroid_name_to_cluster[seq_id] = cluster_id

# ---------- map cluster IDs to Clone 1..10 from original R8 ----------
r8 = pd.read_csv(find_original_file("R8"), sep="\t")
r8["cluster_id"] = r8["N50"].map(sequence_to_cluster)
r8 = r8.dropna(subset=["cluster_id"]).copy()
r8["cluster_id"] = r8["cluster_id"].astype(int)

top10_r8_cluster_ids = (
    r8.groupby("cluster_id")["count"]
    .sum()
    .sort_values(ascending=False)
    .head(10)
    .index
    .tolist()
)

cluster_to_clone = {
    cluster_id: f"Clone {i}"
    for i, cluster_id in enumerate(top10_r8_cluster_ids, start=1)
}

exclude_ids = {
    cid for cid, name in cluster_to_clone.items()
    if name in {"Clone 1", "Clone 2", "Clone 5"}
}

# ---------- full original shares for R4, R5, R6 ----------
original_counts = {}
original_totals = {}

for round_name in ["R4", "R5", "R6"]:
    df = pd.read_csv(find_original_file(round_name), sep="\t")
    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"]).copy()
    df["cluster_id"] = df["cluster_id"].astype(int)

    original_counts[round_name] = df.groupby("cluster_id")["count"].sum()
    original_totals[round_name] = df["count"].sum()

original_shares_num = pd.DataFrame({
    round_name: 100 * original_counts[round_name] / original_totals[round_name]
    for round_name in ["R4", "R5", "R6"]
}).fillna(0)

# ---------- full repeat shares for R5 repeat, R6 repeat ----------
def load_repeat_share(uc_file):
    counts = {}
    total = 0
    with open(uc_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if fields[0] != "H":
                continue
            query_label = fields[8]
            centroid_label = fields[9].split(";")[0]
            cluster_id = centroid_name_to_cluster.get(centroid_label)
            if cluster_id is None:
                continue
            n = parse_size(query_label)
            counts[cluster_id] = counts.get(cluster_id, 0) + n
            total += n
    return pd.Series(counts, dtype=float) * 100 / total

repeat_shares_num = pd.DataFrame({
    "R5 repeat": load_repeat_share(REPEAT_DIR / "R5_repeat_assignments.uc"),
    "R6 repeat": load_repeat_share(REPEAT_DIR / "R6_repeat_assignments.uc"),
}).fillna(0)

# ---------- choose clusters to plot ----------
top20_repeat_ids = [
    cid for cid in repeat_shares_num["R6 repeat"].sort_values(ascending=False).index
    if cid not in exclude_ids
][:20]

positive_ids = [
    int(cid) for cid in results.index
    if results.loc[cid, "unified_repeat_advantage"] > 0 and cid not in exclude_ids
]

cluster_ids_to_plot = []
for cid in top20_repeat_ids + positive_ids:
    if cid not in cluster_ids_to_plot:
        cluster_ids_to_plot.append(cid)

# keep only clusters with all needed values available
cluster_ids_to_plot = [
    cid for cid in cluster_ids_to_plot
    if cid in original_shares_num.index and cid in repeat_shares_num.index
]

print(f"Top 20 by R6 repeat (after exclusions): {len(top20_repeat_ids)}")
print(f"Positive repeat-advantage clusters:     {len(positive_ids)}")
print(f"Total clusters plotted:                 {len(cluster_ids_to_plot)}")

# ---------- labels ----------
def cluster_label(cid):
    base = cluster_to_clone.get(cid, f"Cluster {cid}")
    if cid in positive_ids:
        base += " *"
    return base

# ---------- plot ----------
plt.figure(figsize=(12, 9))
colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(cluster_ids_to_plot)))
all_values = []

for cid, color in zip(cluster_ids_to_plot, colors):
    original_values = [
        original_shares_num.loc[cid, "R4"],
        original_shares_num.loc[cid, "R5"],
        original_shares_num.loc[cid, "R6"],
    ]
    repeat_values = [
        original_shares_num.loc[cid, "R4"],
        repeat_shares_num.loc[cid, "R5 repeat"],
        repeat_shares_num.loc[cid, "R6 repeat"],
    ]

    all_values.extend(original_values + repeat_values)

    # original branch
    plt.plot(
        [4, 5, 6],
        original_values,
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5
    )

    # repeat branch
    plt.plot(
        [4, 5, 6],
        repeat_values,
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2.5 if cid in positive_ids else 2,
        label=cluster_label(cid)
    )

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])
plt.ylim(0, max(all_values) * 1.08)
plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("R4–R6 trajectories: top R6-repeat clusters + positive repeat advantage")
plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=7,
    title="* positive repeat advantage"
)
plt.tight_layout()
plt.savefig(
    CLUSTER_DIR / "R4_R6_top20_repeat_plus_positive_advantage_FIXED.png",
    dpi=300,
    bbox_inches="tight"
)
plt.yscale("log")
plt.show()


# In[40]:


from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"
REPEAT_DIR = BASE / "N50_repeat_clusters_94"

clusters_to_plot = [4, 14, 6, 7, 23, 72, 93, 18]

def parse_size(label):
    m = re.search(r"(?:^|;)size=(\d+)(?:;|$)", label)
    return int(m.group(1)) if m else 1

def find_original_file(round_name):
    round_number = round_name.removeprefix("R")
    patterns = [
        rf"(^|[_-])R{round_number}([_-]|$)",
        rf"(^|[_-])LUR{round_number}([_-]|$)"
    ]
    matches = []
    for path in COLLAPSED_DIR.glob("*_collapsed.tsv"):
        if "repeat" in path.name.lower():
            continue
        if any(re.search(p, path.name, flags=re.IGNORECASE) for p in patterns):
            matches.append(path)
    if len(matches) != 1:
        raise ValueError(f"{round_name}: expected 1 file, found {[p.name for p in matches]}")
    return matches[0]

# sequence ID -> N50
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

# maps
sequence_to_cluster = {}
centroid_name_to_cluster = {}
cluster_to_centroid = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] not in {"S", "H"}:
            continue

        cluster_id = int(fields[1])
        seq_id = fields[8].split(";")[0]

        if seq_id in id_to_seq:
            sequence_to_cluster[id_to_seq[seq_id]] = cluster_id

        if fields[0] == "S":
            centroid_name_to_cluster[seq_id] = cluster_id
            if seq_id in id_to_seq:
                cluster_to_centroid[cluster_id] = id_to_seq[seq_id]

# original shares
original_counts = {}
original_totals = {}

for round_name in ["R4", "R5", "R6"]:
    df = pd.read_csv(find_original_file(round_name), sep="\t")
    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"]).copy()
    df["cluster_id"] = df["cluster_id"].astype(int)

    original_counts[round_name] = df.groupby("cluster_id")["count"].sum()
    original_totals[round_name] = df["count"].sum()

original_shares = pd.DataFrame({
    round_name: 100 * original_counts[round_name] / original_totals[round_name]
    for round_name in ["R4", "R5", "R6"]
}).fillna(0)

# repeat shares
def load_repeat_share(uc_file):
    counts = {}
    total = 0
    with open(uc_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if fields[0] != "H":
                continue

            query_label = fields[8]
            centroid_label = fields[9].split(";")[0]
            cluster_id = centroid_name_to_cluster.get(centroid_label)
            if cluster_id is None:
                continue

            n = parse_size(query_label)
            counts[cluster_id] = counts.get(cluster_id, 0) + n
            total += n

    return pd.Series(counts, dtype=float) * 100 / total

repeat_shares = pd.DataFrame({
    "R5 repeat": load_repeat_share(REPEAT_DIR / "R5_repeat_assignments.uc"),
    "R6 repeat": load_repeat_share(REPEAT_DIR / "R6_repeat_assignments.uc"),
}).fillna(0)

# table
table = pd.DataFrame(index=clusters_to_plot)
table.index.name = "cluster_id"
table["centroid_N50"] = [cluster_to_centroid.get(c, "") for c in clusters_to_plot]

for col in ["R4", "R5", "R6"]:
    table[col] = [original_shares.loc[c, col] if c in original_shares.index else 0 for c in clusters_to_plot]

for col in ["R5 repeat", "R6 repeat"]:
    table[col] = [repeat_shares.loc[c, col] if c in repeat_shares.index else 0 for c in clusters_to_plot]

table = table.reset_index()
table["label"] = table["cluster_id"].apply(lambda x: f"Cluster {x}")

# plot
plt.figure(figsize=(9, 6))
colors = plt.cm.tab10(np.linspace(0, 1, len(clusters_to_plot)))
eps = 1e-4  # for log plotting of zeros

for color, cluster_id in zip(colors, clusters_to_plot):
    r4 = table.loc[table["cluster_id"] == cluster_id, "R4"].iloc[0]
    r5 = table.loc[table["cluster_id"] == cluster_id, "R5"].iloc[0]
    r6 = table.loc[table["cluster_id"] == cluster_id, "R6"].iloc[0]
    r5r = table.loc[table["cluster_id"] == cluster_id, "R5 repeat"].iloc[0]
    r6r = table.loc[table["cluster_id"] == cluster_id, "R6 repeat"].iloc[0]

    orig = [max(r4, eps), max(r5, eps), max(r6, eps)]
    rep = [max(r4, eps), max(r5r, eps), max(r6r, eps)]

    plt.plot(
        [4, 5, 6], orig,
        marker="o", color=color, linewidth=2, alpha=0.5
    )
    plt.plot(
        [4, 5, 6], rep,
        marker="s", linestyle="--", color=color, linewidth=2.2,
        label=f"Cluster {cluster_id}"
    )

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])
plt.yscale("log")
plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Selected clusters: original vs repeat trajectories")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.savefig(
    CLUSTER_DIR / "selected_clusters_R4_R6_log.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

display(
    table[
        ["label", "centroid_N50", "R4", "R5", "R5 repeat", "R6", "R6 repeat"]
    ].round(4)
)


# In[42]:


from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"
REPEAT_DIR = BASE / "N50_repeat_clusters_94"

clusters_to_plot = [4, 14, 6, 7, 23, 72, 18]

def parse_size(label):
    m = re.search(r"(?:^|;)size=(\d+)(?:;|$)", label)
    return int(m.group(1)) if m else 1

def find_original_file(round_name):
    round_number = round_name.removeprefix("R")
    patterns = [
        rf"(^|[_-])R{round_number}([_-]|$)",
        rf"(^|[_-])LUR{round_number}([_-]|$)"
    ]
    matches = []
    for path in COLLAPSED_DIR.glob("*_collapsed.tsv"):
        if "repeat" in path.name.lower():
            continue
        if any(re.search(p, path.name, flags=re.IGNORECASE) for p in patterns):
            matches.append(path)
    if len(matches) != 1:
        raise ValueError(f"{round_name}: expected 1 file, found {[p.name for p in matches]}")
    return matches[0]

# sequence ID -> N50
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

# maps
sequence_to_cluster = {}
centroid_name_to_cluster = {}
cluster_to_centroid = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] not in {"S", "H"}:
            continue

        cluster_id = int(fields[1])
        seq_id = fields[8].split(";")[0]

        if seq_id in id_to_seq:
            sequence_to_cluster[id_to_seq[seq_id]] = cluster_id

        if fields[0] == "S":
            centroid_name_to_cluster[seq_id] = cluster_id
            if seq_id in id_to_seq:
                cluster_to_centroid[cluster_id] = id_to_seq[seq_id]

# original shares
original_counts = {}
original_totals = {}

for round_name in ["R4", "R5", "R6"]:
    df = pd.read_csv(find_original_file(round_name), sep="\t")
    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"]).copy()
    df["cluster_id"] = df["cluster_id"].astype(int)

    original_counts[round_name] = df.groupby("cluster_id")["count"].sum()
    original_totals[round_name] = df["count"].sum()

original_shares = pd.DataFrame({
    round_name: 100 * original_counts[round_name] / original_totals[round_name]
    for round_name in ["R4", "R5", "R6"]
}).fillna(0)

# repeat shares
def load_repeat_share(uc_file):
    counts = {}
    total = 0
    with open(uc_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if fields[0] != "H":
                continue

            query_label = fields[8]
            centroid_label = fields[9].split(";")[0]
            cluster_id = centroid_name_to_cluster.get(centroid_label)
            if cluster_id is None:
                continue

            n = parse_size(query_label)
            counts[cluster_id] = counts.get(cluster_id, 0) + n
            total += n

    return pd.Series(counts, dtype=float) * 100 / total

repeat_shares = pd.DataFrame({
    "R5 repeat": load_repeat_share(REPEAT_DIR / "R5_repeat_assignments.uc"),
    "R6 repeat": load_repeat_share(REPEAT_DIR / "R6_repeat_assignments.uc"),
}).fillna(0)

# table
table = pd.DataFrame(index=clusters_to_plot)
table.index.name = "cluster_id"
table["centroid_N50"] = [cluster_to_centroid.get(c, "") for c in clusters_to_plot]

for col in ["R4", "R5", "R6"]:
    table[col] = [original_shares.loc[c, col] if c in original_shares.index else 0 for c in clusters_to_plot]

for col in ["R5 repeat", "R6 repeat"]:
    table[col] = [repeat_shares.loc[c, col] if c in repeat_shares.index else 0 for c in clusters_to_plot]

table = table.reset_index()
table["label"] = table["cluster_id"].apply(lambda x: f"Cluster {x}")

# plot
plt.figure(figsize=(9, 6))
colors = plt.cm.tab10(np.linspace(0, 1, len(clusters_to_plot)))
eps = 1e-4  # for log plotting of zeros

for color, cluster_id in zip(colors, clusters_to_plot):
    r4 = table.loc[table["cluster_id"] == cluster_id, "R4"].iloc[0]
    r5 = table.loc[table["cluster_id"] == cluster_id, "R5"].iloc[0]
    r6 = table.loc[table["cluster_id"] == cluster_id, "R6"].iloc[0]
    r5r = table.loc[table["cluster_id"] == cluster_id, "R5 repeat"].iloc[0]
    r6r = table.loc[table["cluster_id"] == cluster_id, "R6 repeat"].iloc[0]

    orig = [max(r4, eps), max(r5, eps), max(r6, eps)]
    rep = [max(r4, eps), max(r5r, eps), max(r6r, eps)]

    plt.plot(
        [4, 5, 6], orig,
        marker="o", color=color, linewidth=2, alpha=0.5
    )
    plt.plot(
        [4, 5, 6], rep,
        marker="s", linestyle="--", color=color, linewidth=2.2,
        label=f"Cluster {cluster_id}"
    )

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])
plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Selected clusters: original vs repeat trajectories")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.savefig(
    CLUSTER_DIR / "selected_clusters_R4_R6_log.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

display(
    table[
        ["label", "centroid_N50", "R4", "R5", "R5 repeat", "R6", "R6 repeat"]
    ].round(4)
)


# In[43]:


clone_name = "Clone 7"

cluster_id = next(
    cid for cid, name in cluster_to_clone.items()
    if name == clone_name
)

print(f"{clone_name} = cluster ID {cluster_id}")


# In[44]:


from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"
REPEAT_DIR = BASE / "N50_repeat_clusters_94"

clusters_to_plot = [4, 14, 6, 7, 9, 23, 72, 18]

def parse_size(label):
    m = re.search(r"(?:^|;)size=(\d+)(?:;|$)", label)
    return int(m.group(1)) if m else 1

def find_original_file(round_name):
    round_number = round_name.removeprefix("R")
    patterns = [
        rf"(^|[_-])R{round_number}([_-]|$)",
        rf"(^|[_-])LUR{round_number}([_-]|$)"
    ]
    matches = []
    for path in COLLAPSED_DIR.glob("*_collapsed.tsv"):
        if "repeat" in path.name.lower():
            continue
        if any(re.search(p, path.name, flags=re.IGNORECASE) for p in patterns):
            matches.append(path)
    if len(matches) != 1:
        raise ValueError(f"{round_name}: expected 1 file, found {[p.name for p in matches]}")
    return matches[0]

# sequence ID -> N50
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

# maps
sequence_to_cluster = {}
centroid_name_to_cluster = {}
cluster_to_centroid = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] not in {"S", "H"}:
            continue

        cluster_id = int(fields[1])
        seq_id = fields[8].split(";")[0]

        if seq_id in id_to_seq:
            sequence_to_cluster[id_to_seq[seq_id]] = cluster_id

        if fields[0] == "S":
            centroid_name_to_cluster[seq_id] = cluster_id
            if seq_id in id_to_seq:
                cluster_to_centroid[cluster_id] = id_to_seq[seq_id]

# original shares
original_counts = {}
original_totals = {}

for round_name in ["R4", "R5", "R6"]:
    df = pd.read_csv(find_original_file(round_name), sep="\t")
    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"]).copy()
    df["cluster_id"] = df["cluster_id"].astype(int)

    original_counts[round_name] = df.groupby("cluster_id")["count"].sum()
    original_totals[round_name] = df["count"].sum()

original_shares = pd.DataFrame({
    round_name: 100 * original_counts[round_name] / original_totals[round_name]
    for round_name in ["R4", "R5", "R6"]
}).fillna(0)

# repeat shares
def load_repeat_share(uc_file):
    counts = {}
    total = 0
    with open(uc_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if fields[0] != "H":
                continue

            query_label = fields[8]
            centroid_label = fields[9].split(";")[0]
            cluster_id = centroid_name_to_cluster.get(centroid_label)
            if cluster_id is None:
                continue

            n = parse_size(query_label)
            counts[cluster_id] = counts.get(cluster_id, 0) + n
            total += n

    return pd.Series(counts, dtype=float) * 100 / total

repeat_shares = pd.DataFrame({
    "R5 repeat": load_repeat_share(REPEAT_DIR / "R5_repeat_assignments.uc"),
    "R6 repeat": load_repeat_share(REPEAT_DIR / "R6_repeat_assignments.uc"),
}).fillna(0)

# table
table = pd.DataFrame(index=clusters_to_plot)
table.index.name = "cluster_id"
table["centroid_N50"] = [cluster_to_centroid.get(c, "") for c in clusters_to_plot]

for col in ["R4", "R5", "R6"]:
    table[col] = [original_shares.loc[c, col] if c in original_shares.index else 0 for c in clusters_to_plot]

for col in ["R5 repeat", "R6 repeat"]:
    table[col] = [repeat_shares.loc[c, col] if c in repeat_shares.index else 0 for c in clusters_to_plot]

table = table.reset_index()
table["label"] = table["cluster_id"].apply(lambda x: f"Cluster {x}")

# plot
plt.figure(figsize=(9, 6))
colors = plt.cm.tab10(np.linspace(0, 1, len(clusters_to_plot)))
eps = 1e-4  # for log plotting of zeros

for color, cluster_id in zip(colors, clusters_to_plot):
    r4 = table.loc[table["cluster_id"] == cluster_id, "R4"].iloc[0]
    r5 = table.loc[table["cluster_id"] == cluster_id, "R5"].iloc[0]
    r6 = table.loc[table["cluster_id"] == cluster_id, "R6"].iloc[0]
    r5r = table.loc[table["cluster_id"] == cluster_id, "R5 repeat"].iloc[0]
    r6r = table.loc[table["cluster_id"] == cluster_id, "R6 repeat"].iloc[0]

    orig = [max(r4, eps), max(r5, eps), max(r6, eps)]
    rep = [max(r4, eps), max(r5r, eps), max(r6r, eps)]

    plt.plot(
        [4, 5, 6], orig,
        marker="o", color=color, linewidth=2, alpha=0.5
    )
    plt.plot(
        [4, 5, 6], rep,
        marker="s", linestyle="--", color=color, linewidth=2.2,
        label=f"Cluster {cluster_id}"
    )

plt.xticks([4, 5, 6], ["R4", "R5", "R6"])
plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Selected clusters: original vs repeat trajectories")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.savefig(
    CLUSTER_DIR / "selected_clusters_R4_R6_log.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()

display(
    table[
        ["label", "centroid_N50", "R4", "R5", "R5 repeat", "R6", "R6 repeat"]
    ].round(4)
)


# In[46]:


def find_original_file(round_name):
    round_number = round_name.removeprefix("R")

    matches = [
        path for path in COLLAPSED_DIR.glob("*_collapsed.tsv")
        if path.name.startswith(f"LUR{round_number}_")
        and "repeat" not in path.name.lower()
    ]

    if len(matches) != 1:
        raise ValueError(
            f"{round_name}: expected 1 file, found {[p.name for p in matches]}"
        )

    return matches[0]


# Check that all rounds resolve correctly
for rnd in ["R1", "R4", "R5", "R6", "R7", "R8"]:
    print(rnd, "->", find_original_file(rnd).name)


# In[47]:


# Add R1, R7 and R8 to the selected-cluster trajectory plot

# Build shares for the additional original rounds
for round_name in ["R1", "R7", "R8"]:
    df = pd.read_csv(find_original_file(round_name), sep="\t")
    df["cluster_id"] = df["N50"].map(sequence_to_cluster)
    df = df.dropna(subset=["cluster_id"]).copy()
    df["cluster_id"] = df["cluster_id"].astype(int)

    counts = df.groupby("cluster_id")["count"].sum()
    total = df["count"].sum()

    original_shares[round_name] = 100 * counts / total

original_shares = original_shares.fillna(0)

# Plot
plt.figure(figsize=(11, 7))

colors = plt.cm.tab10(np.linspace(0, 1, len(clusters_to_plot)))

eps = 1e-4  # prevents zeros breaking the log scale

for color, cluster_id in zip(colors, clusters_to_plot):

    # Full ORIGINAL trajectory
    original_rounds = [1, 4, 5, 6, 7, 8]
    original_values = [
        original_shares.loc[cluster_id, "R1"],
        original_shares.loc[cluster_id, "R4"],
        original_shares.loc[cluster_id, "R5"],
        original_shares.loc[cluster_id, "R6"],
        original_shares.loc[cluster_id, "R7"],
        original_shares.loc[cluster_id, "R8"]
    ]

    original_values = [max(x, eps) for x in original_values]

    # Original trajectory — 50% transparent
    plt.plot(
        original_rounds,
        original_values,
        marker="o",
        color=color,
        linewidth=2,
        alpha=0.5
    )

    # REPEAT branch starts from R4
    repeat_values = [
        original_shares.loc[cluster_id, "R4"],
        repeat_shares.loc[cluster_id, "R5 repeat"],
        repeat_shares.loc[cluster_id, "R6 repeat"]
    ]

    repeat_values = [max(x, eps) for x in repeat_values]

    plt.plot(
        [4, 5, 6],
        repeat_values,
        marker="s",
        linestyle="--",
        color=color,
        linewidth=2.5,
        label=f"Cluster {cluster_id}"
    )

plt.xticks(
    [1, 4, 5, 6, 7, 8],
    ["R1", "R4", "R5", "R6", "R7", "R8"]
)

plt.yscale("log")

plt.xlabel("Selection round")
plt.ylabel("Share of pool (%)")
plt.title("Selected cluster trajectories: original vs repeat")

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "selected_clusters_R1_R4_R8_original_repeat_log.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# In[48]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# All clusters observed in R6 repeat
# ------------------------------------------------------------

all_r6_repeat_clusters = repeat_shares_num.index[
    repeat_shares_num["R6 repeat"] > 0
].astype(int)

print(f"Clusters observed in R6 repeat: {len(all_r6_repeat_clusters):,}")


# ------------------------------------------------------------
# Get READ COUNTS for original R4/R5/R6
# ------------------------------------------------------------

orig_counts = pd.DataFrame({
    "R4_count": original_counts["R4"],
    "R5_count": original_counts["R5"],
    "R6_count": original_counts["R6"]
}).fillna(0)


# ------------------------------------------------------------
# Rebuild repeat READ COUNTS from UC files
# ------------------------------------------------------------

def load_repeat_counts(uc_file):
    counts = {}

    with open(uc_file) as fh:
        for line in fh:
            fields = line.rstrip().split("\t")

            if fields[0] != "H":
                continue

            query_label = fields[8]
            centroid_label = fields[9].split(";")[0]

            cluster_id = centroid_name_to_cluster.get(centroid_label)

            if cluster_id is None:
                continue

            n = parse_size(query_label)

            counts[cluster_id] = counts.get(cluster_id, 0) + n

    return pd.Series(counts, dtype=float)


r5_repeat_counts = load_repeat_counts(
    REPEAT_DIR / "R5_repeat_assignments.uc"
)

r6_repeat_counts = load_repeat_counts(
    REPEAT_DIR / "R6_repeat_assignments.uc"
)


# ------------------------------------------------------------
# Total reads in each population
# ------------------------------------------------------------

R4_total = original_totals["R4"]
R5_total = original_totals["R5"]
R6_total = original_totals["R6"]

R5_repeat_total = r5_repeat_counts.sum()
R6_repeat_total = r6_repeat_counts.sum()


# ------------------------------------------------------------
# Calculate enrichment for EVERY R6-repeat cluster
# ------------------------------------------------------------

PSEUDOCOUNT = 0.5

rows = []

for cid in all_r6_repeat_clusters:

    r4 = orig_counts["R4_count"].get(cid, 0)
    r5 = orig_counts["R5_count"].get(cid, 0)
    r6 = orig_counts["R6_count"].get(cid, 0)

    r5r = r5_repeat_counts.get(cid, 0)
    r6r = r6_repeat_counts.get(cid, 0)

    # Frequencies with half-read pseudocount
    f4 = (r4 + PSEUDOCOUNT) / (R4_total + PSEUDOCOUNT)
    f5 = (r5 + PSEUDOCOUNT) / (R5_total + PSEUDOCOUNT)
    f6 = (r6 + PSEUDOCOUNT) / (R6_total + PSEUDOCOUNT)

    f5r = (r5r + PSEUDOCOUNT) / (R5_repeat_total + PSEUDOCOUNT)
    f6r = (r6r + PSEUDOCOUNT) / (R6_repeat_total + PSEUDOCOUNT)

    # Original overall enrichment R4 -> R6
    original_enrichment = np.log2(f6 / f4)

    # Repeat overall enrichment:
    # R4 -> R5 repeat -> R6 repeat
    repeat_enrichment = np.log2(f6r / f4)

    # Difference between repeat and original
    repeat_advantage = repeat_enrichment - original_enrichment

    rows.append({
        "cluster_id": cid,
        "R4_count": r4,
        "R5_count": r5,
        "R6_count": r6,
        "R5_repeat_count": r5r,
        "R6_repeat_count": r6r,
        "original_R4_to_R6_enrichment": original_enrichment,
        "repeat_R4_to_R6_enrichment": repeat_enrichment,
        "unified_repeat_advantage": repeat_advantage
    })


all_repeat_advantage = (
    pd.DataFrame(rows)
    .set_index("cluster_id")
    .sort_values("unified_repeat_advantage", ascending=False)
)


# ------------------------------------------------------------
# Save full table
# ------------------------------------------------------------

OUT = CLUSTER_DIR / "all_R6_repeat_clusters_repeat_advantage.tsv"

all_repeat_advantage.to_csv(
    OUT,
    sep="\t"
)

print(f"Calculated: {len(all_repeat_advantage):,} clusters")
print(
    f"Positive advantage: "
    f"{(all_repeat_advantage['unified_repeat_advantage'] > 0).sum():,}"
)
print(
    f"Negative advantage: "
    f"{(all_repeat_advantage['unified_repeat_advantage'] < 0).sum():,}"
)

print(f"Saved: {OUT}")


# ------------------------------------------------------------
# Histogram
# ------------------------------------------------------------

plt.figure(figsize=(9, 6))

plt.hist(
    all_repeat_advantage["unified_repeat_advantage"],
    bins=60,
    edgecolor="black",
    linewidth=0.5
)

plt.axvline(
    0,
    linestyle="--",
    linewidth=1.5,
    color="black"
)

plt.xlabel("Repeat advantage enrichment (log2)")
plt.ylabel("Number of clusters")

plt.title(
    "Distribution of repeat enrichment advantage\n"
    "All clusters observed in R6 repeat"
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR / "all_R6_repeat_clusters_repeat_advantage_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

display(all_repeat_advantage.head(30))


# In[49]:


# Normalise repeat advantage by each cluster's R4 pool share

PSEUDOCOUNT = 0.5

norm = all_repeat_advantage.copy()

# R4 frequency using the same pseudocount logic as the enrichment calculation
norm["R4_fraction"] = (
    (norm["R4_count"] + PSEUDOCOUNT)
    / (R4_total + PSEUDOCOUNT)
)

# Normalised repeat advantage
norm["repeat_advantage_normalised_R4"] = (
    norm["unified_repeat_advantage"]
    / norm["R4_fraction"]
)

# Save
OUT = CLUSTER_DIR / "all_R6_repeat_clusters_repeat_advantage_normalised_R4.tsv"

norm.to_csv(
    OUT,
    sep="\t"
)

print(f"Clusters analysed: {len(norm):,}")
print(f"Saved: {OUT}")

# ------------------------------------------------------------
# Plot distribution
# ------------------------------------------------------------

plt.figure(figsize=(9, 6))

plt.hist(
    norm["repeat_advantage_normalised_R4"],
    bins=60,
    edgecolor="black",
    linewidth=0.5
)

plt.axvline(
    0,
    linestyle="--",
    linewidth=1.5,
    color="black"
)

plt.xlabel(
    "Repeat advantage enrichment / R4 pool fraction"
)
plt.ylabel("Number of clusters")

plt.title(
    "Repeat enrichment advantage normalised by R4 abundance"
)

plt.tight_layout()

plt.savefig(
    CLUSTER_DIR /
    "all_R6_repeat_clusters_repeat_advantage_normalised_R4_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

display(
    norm.sort_values(
        "repeat_advantage_normalised_R4",
        ascending=False
    ).head(30)
)


# In[2]:


from pathlib import Path
import itertools
import pandas as pd

# Paths
BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"

# ------------------------------------------------------------
# 1. Load sequence ID -> N50 sequence mapping
# ------------------------------------------------------------

seq_info = pd.read_csv(
    CLUSTER_DIR / "sequence_ids.tsv",
    sep="\t"
)

id_to_seq = dict(
    zip(seq_info["sequence_id"], seq_info["N50"])
)

# ------------------------------------------------------------
# 2. Rebuild sequence -> cluster and cluster -> centroid maps
# ------------------------------------------------------------

sequence_to_cluster = {}
cluster_to_centroid = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")

        if fields[0] not in {"S", "H"}:
            continue

        cluster_id = int(fields[1])
        seq_id = fields[8].split(";")[0]

        if seq_id in id_to_seq:
            sequence_to_cluster[id_to_seq[seq_id]] = cluster_id

        # S records define the cluster centroid
        if fields[0] == "S" and seq_id in id_to_seq:
            cluster_to_centroid[cluster_id] = id_to_seq[seq_id]

# ------------------------------------------------------------
# 3. Determine Clone 1–10 from R8 abundance
# ------------------------------------------------------------

r8 = pd.read_csv(
    COLLAPSED_DIR / "LUR8_S6_L001_R1_001_collapsed.tsv",
    sep="\t"
)

r8["cluster_id"] = r8["N50"].map(sequence_to_cluster)

r8 = r8.dropna(subset=["cluster_id"]).copy()
r8["cluster_id"] = r8["cluster_id"].astype(int)

top10_cluster_ids = (
    r8.groupby("cluster_id")["count"]
    .sum()
    .sort_values(ascending=False)
    .head(10)
    .index
    .tolist()
)

cluster_to_clone = {
    cid: f"Clone {i}"
    for i, cid in enumerate(top10_cluster_ids, start=1)
}

# ------------------------------------------------------------
# 4. Get centroid sequences
# ------------------------------------------------------------

clone_to_sequence = {
    clone: cluster_to_centroid[cid]
    for cid, clone in cluster_to_clone.items()
}

# ------------------------------------------------------------
# 5. Calculate all 45 pairwise identities
# ------------------------------------------------------------

pairwise = []

for (clone_a, seq_a), (clone_b, seq_b) in itertools.combinations(
    clone_to_sequence.items(), 2
):
    identity = 100 * sum(
        a == b for a, b in zip(seq_a, seq_b)
    ) / len(seq_a)

    pairwise.append({
        "clone_1": clone_a,
        "clone_2": clone_b,
        "percent_identity": identity
    })

pairwise_df = pd.DataFrame(pairwise)

# ------------------------------------------------------------
# Results
# ------------------------------------------------------------

print("Clone / cluster mapping:")
for cid, clone in cluster_to_clone.items():
    print(f"{clone}: Cluster {cid}")

print()
print(f"Pairwise comparisons:     {len(pairwise_df)}")
print(f"Mean pairwise identity:   {pairwise_df['percent_identity'].mean():.2f}%")
print(f"Median pairwise identity: {pairwise_df['percent_identity'].median():.2f}%")
print(f"Minimum identity:         {pairwise_df['percent_identity'].min():.2f}%")
print(f"Maximum identity:         {pairwise_df['percent_identity'].max():.2f}%")

display(
    pairwise_df.sort_values(
        "percent_identity",
        ascending=False
    )
)


# In[3]:


from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Paths
BASE = Path("/Users/luisohlendorf/Documents/PhD/QT_reanalysis")
CLUSTER_DIR = BASE / "N50_clusters_94"
COLLAPSED_DIR = BASE / "N50_collapsed"

# Load sequence ID -> N50
seq_info = pd.read_csv(CLUSTER_DIR / "sequence_ids.tsv", sep="\t")
id_to_seq = dict(zip(seq_info["sequence_id"], seq_info["N50"]))

# Rebuild maps
sequence_to_cluster = {}
cluster_to_centroid = {}

with open(CLUSTER_DIR / "cluster_assignments.uc") as fh:
    for line in fh:
        fields = line.rstrip().split("\t")
        if fields[0] not in {"S", "H"}:
            continue

        cluster_id = int(fields[1])
        seq_id = fields[8].split(";")[0]

        if seq_id in id_to_seq:
            sequence_to_cluster[id_to_seq[seq_id]] = cluster_id

        if fields[0] == "S" and seq_id in id_to_seq:
            cluster_to_centroid[cluster_id] = id_to_seq[seq_id]

# Determine Clone 1 and Clone 2 from R8 abundance
r8 = pd.read_csv(COLLAPSED_DIR / "LUR8_S6_L001_R1_001_collapsed.tsv", sep="\t")
r8["cluster_id"] = r8["N50"].map(sequence_to_cluster)
r8 = r8.dropna(subset=["cluster_id"]).copy()
r8["cluster_id"] = r8["cluster_id"].astype(int)

top10_cluster_ids = (
    r8.groupby("cluster_id")["count"]
    .sum()
    .sort_values(ascending=False)
    .head(10)
    .index
    .tolist()
)

clone1_cluster = top10_cluster_ids[0]
clone2_cluster = top10_cluster_ids[1]

seq1 = cluster_to_centroid[clone1_cluster]
seq2 = cluster_to_centroid[clone2_cluster]

matches = [a == b for a, b in zip(seq1, seq2)]
identity = 100 * sum(matches) / len(seq1)

# Plot
fig, ax = plt.subplots(figsize=(16, 3.8))
ax.set_xlim(0.5, len(seq1) + 0.5)
ax.set_ylim(0, 3)
ax.axis("off")

for i, (a, b, m) in enumerate(zip(seq1, seq2, matches), start=1):
    color = "#c7f0c2" if m else "#f7c6c6"  # green for match, red for mismatch

    # background boxes
    ax.add_patch(plt.Rectangle((i - 0.5, 1.6), 1, 0.7, facecolor=color, edgecolor="black", linewidth=0.5))
    ax.add_patch(plt.Rectangle((i - 0.5, 0.7), 1, 0.7, facecolor=color, edgecolor="black", linewidth=0.5))

    # letters
    ax.text(i, 1.95, a, ha="center", va="center", fontsize=9, family="monospace")
    ax.text(i, 1.05, b, ha="center", va="center", fontsize=9, family="monospace")

    # position label
    if i % 5 == 0 or i == 1:
        ax.text(i, 2.45, str(i), ha="center", va="bottom", fontsize=8)

# row labels
ax.text(-1.5, 1.95, f"Clone 1\n(Cluster {clone1_cluster})", ha="right", va="center", fontsize=10)
ax.text(-1.5, 1.05, f"Clone 2\n(Cluster {clone2_cluster})", ha="right", va="center", fontsize=10)

# legend / summary
ax.text(
    1, 2.8,
    f"Pairwise comparison of centroid N50 sequences\n"
    f"Green = match, Red = mismatch   |   Identity = {identity:.1f}% "
    f"({sum(matches)}/{len(seq1)} positions match)",
    ha="left", va="top", fontsize=11
)

plt.tight_layout()
plt.show()

print("Clone 1 sequence:", seq1)
print("Clone 2 sequence:", seq2)
print(f"Percent identity: {identity:.2f}%")


# In[4]:


# Align and display the top 10 clone centroid sequences

from collections import Counter

# Get Clone 1–10 sequences
aligned_clones = {
    f"Clone {i}": cluster_to_centroid[cluster_id]
    for i, cluster_id in enumerate(top10_cluster_ids, start=1)
}

# Build majority-rule consensus
L = len(next(iter(aligned_clones.values())))

consensus = ""

for pos in range(L):
    column = [seq[pos] for seq in aligned_clones.values()]
    consensus += Counter(column).most_common(1)[0][0]

# Mark fully conserved positions
conservation_line = ""

for pos in range(L):
    column = [seq[pos] for seq in aligned_clones.values()]
    conservation_line += "*" if len(set(column)) == 1 else " "

# Display alignment
print("Top 10 clone alignment\n")

for clone, seq in aligned_clones.items():
    cluster_id = top10_cluster_ids[int(clone.split()[-1]) - 1]
    print(f"{clone:8s} C{cluster_id:<4d}  {seq}")

print(" " * 15 + conservation_line)
print(f"{'Consensus':13s}  {consensus}")

# Save as FASTA
OUT = CLUSTER_DIR / "top10_clone_centroid_alignment.fasta"

with open(OUT, "w") as fh:
    for clone, seq in aligned_clones.items():
        cluster_id = top10_cluster_ids[int(clone.split()[-1]) - 1]
        fh.write(f">{clone.replace(' ', '_')}_cluster_{cluster_id}\n{seq}\n")

print(f"\nSaved to: {OUT}")


# In[5]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# aligned_clones is assumed to come from the previous cell
clone_names = list(aligned_clones.keys())
seqs = [aligned_clones[name] for name in clone_names]
L = len(seqs[0])

base_to_num = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}
matrix = np.array([[base_to_num.get(base, 4) for base in seq] for seq in seqs])

cmap = ListedColormap([
    "#4daf4a",  # A
    "#377eb8",  # C
    "#ff7f00",  # G
    "#e41a1c",  # T
    "#dddddd"   # gap
])

plt.figure(figsize=(14, 4))

plt.imshow(
    matrix,
    aspect="auto",
    interpolation="nearest",
    cmap=cmap,
    vmin=0,
    vmax=4
)

plt.yticks(range(len(clone_names)), clone_names)
plt.xticks(range(L), [f"{i+1}\n{consensus[i]}" for i in range(L)], fontsize=8)

plt.xlabel("Alignment position and consensus nucleotide")
plt.ylabel("Clone")
plt.title("Alignment of top 10 clone centroid sequences")

plt.tight_layout()
plt.show()


# In[6]:


# Align Clones 2–10 to Clone 1, allowing insertions/deletions

from pathlib import Path

# Uses aligned_clones from the previous top-10 clone cell
# e.g. {"Clone 1": "...", ..., "Clone 10": "..."}

reference = aligned_clones["Clone 1"]

# ------------------------------------------------------------
# Needleman-Wunsch global pairwise alignment
# ------------------------------------------------------------

def global_align(ref, query, match=2, mismatch=-1, gap=-2):
    n, m = len(ref), len(query)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = i * gap
        trace[i][0] = "U"

    for j in range(1, m + 1):
        score[0][j] = j * gap
        trace[0][j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            diag = score[i-1][j-1] + (
                match if ref[i-1] == query[j-1] else mismatch
            )
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap

            best = max(diag, up, left)
            score[i][j] = best

            # Prefer substitutions/matches over gaps in ties
            if best == diag:
                trace[i][j] = "D"
            elif best == up:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

    aligned_ref = []
    aligned_query = []

    i, j = n, m

    while i > 0 or j > 0:
        move = trace[i][j]

        if move == "D":
            aligned_ref.append(ref[i-1])
            aligned_query.append(query[j-1])
            i -= 1
            j -= 1

        elif move == "U":
            aligned_ref.append(ref[i-1])
            aligned_query.append("-")
            i -= 1

        else:
            aligned_ref.append("-")
            aligned_query.append(query[j-1])
            j -= 1

    return (
        "".join(reversed(aligned_ref)),
        "".join(reversed(aligned_query))
    )


# ------------------------------------------------------------
# Convert each pairwise alignment into reference coordinates
# ------------------------------------------------------------

def decompose_reference_alignment(aligned_ref, aligned_query, ref_length):
    # insertions[b] = bases inserted before reference position b
    # b = 0 means before first base; b = ref_length means after last
    insertions = [""] * (ref_length + 1)
    ref_bases = ["-"] * ref_length

    ref_pos = 0

    for r, q in zip(aligned_ref, aligned_query):

        if r == "-":
            # insertion relative to Clone 1
            insertions[ref_pos] += q

        else:
            # nucleotide or deletion aligned to this Clone 1 position
            ref_bases[ref_pos] = q
            ref_pos += 1

    return insertions, ref_bases


decomposed = {}

# Clone 1 itself
decomposed["Clone 1"] = (
    [""] * (len(reference) + 1),
    list(reference)
)

# Align all other clones to Clone 1
for clone in [f"Clone {i}" for i in range(2, 11)]:

    aligned_ref, aligned_query = global_align(
        reference,
        aligned_clones[clone]
    )

    decomposed[clone] = decompose_reference_alignment(
        aligned_ref,
        aligned_query,
        len(reference)
    )


# ------------------------------------------------------------
# Merge all pairwise alignments into one reference-anchored MSA
# ------------------------------------------------------------

max_insertions = []

for boundary in range(len(reference) + 1):
    max_insertions.append(
        max(
            len(decomposed[clone][0][boundary])
            for clone in decomposed
        )
    )


msa = {}

for clone in [f"Clone {i}" for i in range(1, 11)]:

    insertions, ref_bases = decomposed[clone]
    pieces = []

    for i in range(len(reference)):

        # insertion before this reference base
        ins = insertions[i]
        pieces.append(
            ins + "-" * (max_insertions[i] - len(ins))
        )

        # base aligned to reference position
        pieces.append(ref_bases[i])

    # insertion after final reference base
    ins = insertions[len(reference)]
    pieces.append(
        ins + "-" * (
            max_insertions[len(reference)] - len(ins)
        )
    )

    msa[clone] = "".join(pieces)


# ------------------------------------------------------------
# Display alignment
# ------------------------------------------------------------

print("Reference-anchored alignment\n")

for clone in [f"Clone {i}" for i in range(1, 11)]:
    print(f"{clone:8s}  {msa[clone]}")

print()
print(f"Clone 1 length:    {len(reference)} nt")
print(f"Alignment length:  {len(msa['Clone 1'])} columns")


# ------------------------------------------------------------
# Save as FASTA
# ------------------------------------------------------------

OUT = CLUSTER_DIR / "top10_clones_aligned_to_clone1_with_indels.fasta"

with open(OUT, "w") as fh:
    for clone in [f"Clone {i}" for i in range(1, 11)]:
        fh.write(f">{clone.replace(' ', '_')}\n{msa[clone]}\n")

print(f"Saved to: {OUT}")


# In[7]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

# Assumes `msa` was created in the previous cell
clone_order = [f"Clone {i}" for i in range(1, 11)]
aligned = [msa[c] for c in clone_order]
ref = msa["Clone 1"]
L = len(ref)

# Encode:
# 0 = match to Clone 1
# 1 = A difference
# 2 = C difference
# 3 = G difference
# 4 = T difference
# 5 = gap ("-")
code = {"A": 1, "C": 2, "G": 3, "T": 4, "-": 5}

matrix = np.zeros((len(clone_order), L), dtype=int)

for i, seq in enumerate(aligned):
    for j, base in enumerate(seq):
        if i == 0:
            matrix[i, j] = code.get(base, 5)
        else:
            if base == ref[j]:
                matrix[i, j] = 0
            else:
                matrix[i, j] = code.get(base, 5)

# Colours: match grey, then A/C/G/T/gap
cmap = ListedColormap([
    "#e6e6e6",  # match
    "#4daf4a",  # A
    "#377eb8",  # C
    "#ff7f00",  # G
    "#e41a1c",  # T
    "#000000",  # gap
])

norm = BoundaryNorm(np.arange(-0.5, 6.5, 1), cmap.N)

plt.figure(figsize=(16, 5))
plt.imshow(matrix, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

plt.yticks(range(len(clone_order)), clone_order)

# Show every 5th position on the x-axis, with the Clone 1 base
xticks = list(range(0, L, 5))
xticklabels = [f"{i+1}\n{ref[i]}" for i in xticks]
plt.xticks(xticks, xticklabels, fontsize=8)

plt.xlabel("Alignment column / Clone 1 reference base")
plt.ylabel("Clone")
plt.title("Top 10 clones aligned to Clone 1\nGrey = match to Clone 1, colours = differences, black = gap")

# Simple legend
from matplotlib.patches import Patch
legend_handles = [
    Patch(facecolor="#e6e6e6", edgecolor="black", label="Match to Clone 1"),
    Patch(facecolor="#4daf4a", edgecolor="black", label="A"),
    Patch(facecolor="#377eb8", edgecolor="black", label="C"),
    Patch(facecolor="#ff7f00", edgecolor="black", label="G"),
    Patch(facecolor="#e41a1c", edgecolor="black", label="T"),
    Patch(facecolor="#000000", edgecolor="black", label="Gap"),
]
plt.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left")

plt.tight_layout()
plt.show()


# In[8]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Uses clone_to_sequence from the previous analysis
# Reorder explicitly
clone_names = [f"Clone {i}" for i in range(1, 11)]
sequences = {c: clone_to_sequence[c] for c in clone_names}


# ------------------------------------------------------------
# Global alignment allowing substitutions + indels
# ------------------------------------------------------------

def global_align(ref, query, match=2, mismatch=-1, gap=-2):
    n, m = len(ref), len(query)

    score = np.zeros((n + 1, m + 1), dtype=float)
    trace = np.empty((n + 1, m + 1), dtype=object)

    for i in range(1, n + 1):
        score[i, 0] = i * gap
        trace[i, 0] = "U"

    for j in range(1, m + 1):
        score[0, j] = j * gap
        trace[0, j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            diag = score[i-1, j-1] + (
                match if ref[i-1] == query[j-1] else mismatch
            )
            up = score[i-1, j] + gap
            left = score[i, j-1] + gap

            best = max(diag, up, left)
            score[i, j] = best

            if best == diag:
                trace[i, j] = "D"
            elif best == up:
                trace[i, j] = "U"
            else:
                trace[i, j] = "L"

    aligned_ref = []
    aligned_query = []

    i, j = n, m

    while i > 0 or j > 0:

        if i > 0 and j > 0 and trace[i, j] == "D":
            aligned_ref.append(ref[i-1])
            aligned_query.append(query[j-1])
            i -= 1
            j -= 1

        elif i > 0 and trace[i, j] == "U":
            aligned_ref.append(ref[i-1])
            aligned_query.append("-")
            i -= 1

        else:
            aligned_ref.append("-")
            aligned_query.append(query[j-1])
            j -= 1

    return (
        "".join(reversed(aligned_ref)),
        "".join(reversed(aligned_query))
    )


# ------------------------------------------------------------
# Pairwise identity matrix
# identity = matching nucleotides / total alignment columns
# ------------------------------------------------------------

N = len(clone_names)

identity_matrix = np.eye(N) * 100
distance_matrix = np.zeros((N, N))

for i in range(N):
    for j in range(i + 1, N):

        a, b = global_align(
            sequences[clone_names[i]],
            sequences[clone_names[j]]
        )

        matches = sum(
            x == y and x != "-"
            for x, y in zip(a, b)
        )

        identity = 100 * matches / len(a)

        identity_matrix[i, j] = identity
        identity_matrix[j, i] = identity

        # distance from 0 to 1
        distance = 1 - identity / 100

        distance_matrix[i, j] = distance
        distance_matrix[j, i] = distance


identity_df = pd.DataFrame(
    identity_matrix,
    index=clone_names,
    columns=clone_names
)


# ------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------

upper = identity_matrix[np.triu_indices(N, k=1)]

print(f"Mean pairwise identity:   {upper.mean():.2f}%")
print(f"Median pairwise identity: {np.median(upper):.2f}%")
print(f"Minimum identity:         {upper.min():.2f}%")
print(f"Maximum identity:         {upper.max():.2f}%")

display(identity_df.round(1))


# ------------------------------------------------------------
# Identity heatmap
# ------------------------------------------------------------

plt.figure(figsize=(8, 7))

im = plt.imshow(
    identity_matrix,
    cmap="viridis",
    vmin=0,
    vmax=100
)

plt.colorbar(
    im,
    label="Global alignment identity (%)"
)

plt.xticks(
    range(N),
    clone_names,
    rotation=45,
    ha="right"
)

plt.yticks(
    range(N),
    clone_names
)

# Put identity values in cells
for i in range(N):
    for j in range(N):
        plt.text(
            j,
            i,
            f"{identity_matrix[i, j]:.0f}",
            ha="center",
            va="center",
            fontsize=8
        )

plt.title("Pairwise sequence identity of top 10 clones")
plt.tight_layout()
plt.show()


# ------------------------------------------------------------
# Hierarchical clustering dendrogram
# ------------------------------------------------------------

condensed = squareform(
    distance_matrix,
    checks=False
)

Z = linkage(
    condensed,
    method="average"
)

plt.figure(figsize=(9, 5))

dendrogram(
    Z,
    labels=clone_names,
    leaf_rotation=45,
    leaf_font_size=10
)

plt.ylabel("Sequence distance (1 − identity)")
plt.title("Sequence-similarity clustering of top 10 clones")

plt.tight_layout()
plt.show()


# In[9]:


# Recalculate pairwise identity excluding Clone 1

identity_no_clone1 = identity_matrix[1:, 1:]

upper_no_clone1 = identity_no_clone1[
    np.triu_indices(identity_no_clone1.shape[0], k=1)
]

print(f"Pairwise comparisons:     {len(upper_no_clone1)}")
print(f"Mean pairwise identity:   {upper_no_clone1.mean():.2f}%")
print(f"Median pairwise identity: {np.median(upper_no_clone1):.2f}%")
print(f"Minimum identity:         {upper_no_clone1.min():.2f}%")
print(f"Maximum identity:         {upper_no_clone1.max():.2f}%")


# In[10]:


for i, cluster_id in enumerate(top10_cluster_ids, start=1):
    print(f"Clone {i}: {cluster_to_centroid[cluster_id]}")


# In[ ]:




