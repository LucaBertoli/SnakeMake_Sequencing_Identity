#!/usr/bin/env python3
"""
Crea due subplot affiancati:
- Sinistra: qualità assegnata dal sequenziatore (median Qscore)
- Destra: error rate (mismatches + indels)
Stratificati per intervalli di 50 bp di insert size (150–800 bp).
Può prendere in input 4 TSV e plotta le 4 serie nello stesso grafico.
"""

import gzip
import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argparse

def load_tsv(tsv_path, max_reads_per_bin=10_000_000):
    bin_data = {}
    bin_full_flags = set()
    with gzip.open(tsv_path, "rt") if tsv_path.endswith(".gz") else open(tsv_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            insert_size = int(row["insert_size"])
            insert_bin = (insert_size // 50) * 50
            if insert_bin < 150 or insert_bin > 800:
                continue
            if insert_bin not in bin_data:
                bin_data[insert_bin] = []

            if len(bin_data[insert_bin]) < max_reads_per_bin:
                bases = int(row["aligned_read_length"])
                errors = int(row["mismatch_filtered"]) + int(row["indel_filtered"])
                bin_data[insert_bin].append({
                    "insert_bin": insert_bin,
                    "mean_qual_log": float(row["mean_qual_log"]),
                    "errors": errors,
                    "bases": bases
                })
                if len(bin_data[insert_bin]) == max_reads_per_bin and insert_bin not in bin_full_flags:
                    bin_full_flags.add(insert_bin)

            if all(len(v) >= max_reads_per_bin for v in bin_data.values()):
                break
    all_rows = [row for rows in bin_data.values() for row in rows]
    return pd.DataFrame(all_rows)

def compute_median_qscore(df):
    return df.groupby("insert_bin")["mean_qual_log"].median().reset_index()

def compute_error_rate(df):
    agg = df.groupby("insert_bin")[["errors", "bases"]].sum().reset_index()
    agg["Error_Rate"] = agg["errors"] / agg["bases"]
    return agg[["insert_bin", "Error_Rate"]]

def main(tsv_files, output_path, y_margin=0.05):
    colors = ["b", "g", "r", "cyan"]
    labels = ["Element", "GeneMind", "Illumina", "MGI"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    sns.set(style="whitegrid", font_scale=1.1)

    # ====== Qscore median plot ======
    ax1 = axes[0]
    for tsv, color, label in zip(tsv_files, colors, labels):
        df = load_tsv(tsv)
        medians = compute_median_qscore(df)
        ax1.plot(
            medians["insert_bin"], medians["mean_qual_log"],
            marker="o", markersize=4, linewidth=1.5,
            color=color, label=label
        )
    ax1.set_xlabel("Insert size (bp)")
    ax1.set_ylabel("Median read Qscore")
    ax1.set_title("Average Sequencer Qscore")
    ax1.tick_params(axis="x", rotation=45)
    ax1.set_ylim(0, 45)
    ax1.grid(True, linestyle="--", alpha=0.7)

    # ====== Error Rate plot ======
    ax2 = axes[1]
    for tsv, color, label in zip(tsv_files, colors, labels):
        df = load_tsv(tsv)
        error_df = compute_error_rate(df)
        ax2.plot(
            error_df["insert_bin"], error_df["Error_Rate"],
            marker="o", markersize=4, linewidth=1.5,
            color=color, label=label
        )
    ax2.set_xlabel("Insert size (bp)")
    ax2.set_ylabel("Error Rate (%)")
    ax2.set_title("Read Error Rate")
    ax2.tick_params(axis="x", rotation=45)
    ax2.grid(True, linestyle="--", alpha=0.7)
    ax2.set_ylim(0, 0.02)
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y*100:.2f}"))

    # ====== Legenda globale ======
    
    handles = [Line2D([0], [0], color=color, lw=2, marker='o') for color in colors]
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.05),
               ncol=4, frameon=False)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✅ Plot salvato in: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Median Qscore e Error Rate per 4 TSV")
    parser.add_argument("-i1", "--input1", required=True)
    parser.add_argument("-i2", "--input2", required=True)
    parser.add_argument("-i3", "--input3", required=True)
    parser.add_argument("-i4", "--input4", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--y-margin", type=float, default=0.05)
    args = parser.parse_args()

    main([args.input1, args.input2, args.input3, args.input4], args.output, args.y_margin)
