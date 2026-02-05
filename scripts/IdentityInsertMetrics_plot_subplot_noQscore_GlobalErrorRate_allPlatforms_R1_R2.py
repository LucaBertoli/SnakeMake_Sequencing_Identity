#!/usr/bin/env python3
"""
Crea un plot 2x2:
- Riga 1: Read 1 (R1)
- Riga 2: Read 2 (R2)
- Colonna sinistra: Median Qscore (mean_qual_log)
- Colonna destra: Error rate (mismatches + indels)

Stratificato per intervalli di insert size da 50 bp (150–800 bp).
Supporta 4 TSV (Element, GeneMind, Illumina, MGI).
"""

import gzip
import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argparse


# ============================================================
# Lettura TSV e separazione R1 / R2
# ============================================================

def load_tsv(tsv_path):
    bin_data = {
        "R1": {},
        "R2": {}
    }

    opener = gzip.open if tsv_path.endswith(".gz") else open
    with opener(tsv_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            mode = row["read_type"].strip().upper()
            if mode not in ("R1", "R2"):
                continue

            insert_size = int(row["insert_size"])
            insert_bin = (insert_size // 50) * 50
            if insert_bin < 150 or insert_bin > 800:
                continue

            if insert_bin not in bin_data[mode]:
                bin_data[mode][insert_bin] = []

            bases = int(row["aligned_read_length"])
            errors = int(row["mismatch_filtered"]) + int(row["indel_filtered"])

            bin_data[mode][insert_bin].append({
                "insert_bin": insert_bin,
                "mean_qual_log": float(row["mean_qual_log"]),
                "errors": errors,
                "bases": bases
            })

    def make_df(d):
        rows = [r for rows in d.values() for r in rows]
        return pd.DataFrame(rows) if rows else None

    return make_df(bin_data["R1"]), make_df(bin_data["R2"])


# ============================================================
# Metriche
# ============================================================

def compute_median_qscore(df):
    return (
        df.groupby("insert_bin")["mean_qual_log"]
        .median()
        .reset_index()
    )


def compute_error_rate(df):
    agg = (
        df.groupby("insert_bin")[["errors", "bases"]]
        .sum()
        .reset_index()
    )
    agg["Error_Rate"] = agg["errors"] / agg["bases"]
    return agg[["insert_bin", "Error_Rate"]]


# ============================================================
# Main
# ============================================================

def main(tsv_files, output_path, y_margin=0.05):
    labels = ["Element", "GeneMind", "Illumina", "MGI"]
    colors = ["#E69F00", "#009E73", "#004488", "#CC79A7"]

    sns.set(style="whitegrid", font_scale=1.1)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=False)

    for row, read_type in enumerate(["R1", "R2"]):
        for col, metric in enumerate(["Qscore", "ErrorRate"]):
            ax = axes[row, col]

            for tsv, label, color in zip(tsv_files, labels, colors):
                df_R1, df_R2 = load_tsv(tsv)
                df = df_R1 if read_type == "R1" else df_R2

                if df is None or df.empty:
                    continue

                if metric == "Qscore":
                    data = compute_median_qscore(df)
                    ax.plot(
                        data["insert_bin"],
                        data["mean_qual_log"],
                        marker="o",
                        markersize=4,
                        linewidth=1.5,
                        color=color,
                        label=label
                    )
                    ax.set_ylabel("Median read Qscore")
                    ax.set_ylim(0, 45)

                else:
                    data = compute_error_rate(df)
                    ax.plot(
                        data["insert_bin"],
                        data["Error_Rate"],
                        marker="o",
                        markersize=4,
                        linewidth=1.5,
                        color=color,
                        label=label
                    )
                    ax.set_ylabel("Error rate (%)")
                    ax.set_ylim(0, 0.02)
                    ax.yaxis.set_major_formatter(
                        plt.FuncFormatter(lambda y, _: f"{y*100:.2f}")
                    )

            ax.set_xlabel("Insert size (bp)")
            ax.tick_params(axis="x", rotation=45)
            ax.grid(True, linestyle="--", alpha=0.7)

            title = (
                f"{read_type} – Sequencer Qscore"
                if metric == "Qscore"
                else f"{read_type} – Read Error Rate"
            )
            ax.set_title(title)

    # ========================================================
    # Legenda globale
    # ========================================================

    handles = [
        Line2D([0], [0], color=color, lw=2, marker="o")
        for color in colors
    ]
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.03),
        ncol=4,
        frameon=False
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✅ Plot salvato in: {output_path}")


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot 2x2 Qscore ed Error Rate per R1/R2 da 4 TSV"
    )
    parser.add_argument("-i1", "--input1", required=True)
    parser.add_argument("-i2", "--input2", required=True)
    parser.add_argument("-i3", "--input3", required=True)
    parser.add_argument("-i4", "--input4", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--y-margin", type=float, default=0.05)

    args = parser.parse_args()

    main(
        [args.input1, args.input2, args.input3, args.input4],
        args.output,
        args.y_margin
    )
