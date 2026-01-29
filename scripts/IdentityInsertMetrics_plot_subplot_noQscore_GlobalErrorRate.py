#!/usr/bin/env python3
"""
Crea due subplot affiancati:
- Sinistra: qualità assegnata dal sequenziatore (mean_qual_log)
- Destra: identità delle read (identity_filtered_with_indels)
Entrambi stratificati per intervalli di 50 bp di insert size (150–800 bp).

Input:
    - TSV.gz generato da IdentityInsertMetrics.py
Output:
    - File PNG o PDF
"""

import gzip
import csv
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argparse
import math

def main(tsv_path, output_path, y_margin=0.05, max_reads_per_bin=10_000_000):
    bin_data = {}
    bin_full_flags = set()

    # --- Lettura TSV.gz ---
    with gzip.open(tsv_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            insert_size = int(row["insert_size"])
            insert_bin = (insert_size // 50) * 50
            if insert_bin < 150 or insert_bin > 800:
                continue
            if insert_bin not in bin_data:
                bin_data[insert_bin] = []

            if len(bin_data[insert_bin]) < max_reads_per_bin:
                bases=int(row["aligned_read_length"])
                errors = (int(row["mismatch_filtered"]) + int(row["indel_filtered"]))
                bin_data[insert_bin].append({
                    "insert_bin": insert_bin,
                    "mean_qual_log": float(row["mean_qual_log"]),
                    "errors": errors,
                    "bases": bases
                })
                if len(bin_data[insert_bin]) == max_reads_per_bin and insert_bin not in bin_full_flags:
                    print(f"⚡ Bin {insert_bin} ha raggiunto il limite di {max_reads_per_bin} read.")
                    bin_full_flags.add(insert_bin)

            # Stop se tutti i bin sono pieni
            if all(len(v) >= max_reads_per_bin for v in bin_data.values()):
                break

    # --- Crea DataFrame ---
    all_rows = [row for rows in bin_data.values() for row in rows]
    if not all_rows:
        raise ValueError("Nessuna read trovata nei bin richiesti.")
    df = pd.DataFrame(all_rows)

    # --- Impostazioni di stile ---
    sns.set(style="whitegrid", font_scale=1.1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    colors = {"mean_qual_log": "#36E2EE", "Error_Rate": "#43F05F"}

    # ================================
    # SUBPLOT 1 – Qualità del sequenziatore
    # ================================
    ax1 = axes[0]
    sns.violinplot(
        data=df, x="insert_bin", y="mean_qual_log",
        color=colors["mean_qual_log"],
        inner="quart", scale="width", cut=0, ax=ax1
    )
    ax1.set_ylabel("Average read Qscore")
    ax1.set_xlabel("Insert size (bp)")
    ax1.set_ylim(0, df["mean_qual_log"].max() * (1 + y_margin))
    ax1.tick_params(axis="x", rotation=45)
    ax1.set_title("Sequencer Qscore")

    medians1 = df.groupby("insert_bin")["mean_qual_log"].median()
    ax1.plot(
        range(len(medians1)),
        medians1.values,
        color="black",
        marker="o",
        markersize=3,
        linewidth=1,
        label="Mediana"
    )    

    # ================================
    # SUBPLOT 2 – Identità (non Phred-scaled)
    # ================================
    df_error_rate = df.groupby("insert_bin").agg({
        "errors": "sum",
        "bases": "sum"
        }).reset_index()
    df_error_rate["Error_Rate"] = df_error_rate["errors"] / df_error_rate["bases"]

    ax2 = axes[1]
    ax2.plot(
        df_error_rate["insert_bin"], df_error_rate["Error_Rate"],
        color=colors["Error_Rate"],
        marker="o", 
        markersize=5,
        linewidth=1
    )
    ax2.set_ylabel("Error Rate (%)")
    ax2.set_xlabel("Insert size (bp)")
    ax2.tick_params(axis="x", rotation=45)
    ax2.set_title("Read's Error Rate (mismatches + indels)")
    ax2.set_ylim(0, 0.02)
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y*100:.1f}"))

    # --- Griglia e ticks ---
    ax1.set_yticks(range(0, int(df["mean_qual_log"].max()) + 5, 5))
    ax1.grid(True, axis="both", linestyle="--", alpha=0.7)
    ax2.grid(True, axis="both", linestyle="--", alpha=0.7)

    # ================================
    # LEGENDA GLOBALE
    # ================================
    handles = [
        Line2D([0], [0], color=colors["mean_qual_log"], lw=6),
        Line2D([0], [0], color=colors["Error_Rate"], lw=6),
        Line2D([0], [0], color="black", lw=2, marker="o")
    ]
    labels = ["Qualità assegnata", "Error Rate", "Mediana"]
    fig.legend(handles, labels,
               loc="upper center", bbox_to_anchor=(0.5, -0.05),
               ncol=3, frameon=False)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✅ Plot salvato in: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Violin plot separati per qualità assegnata e identità")
    parser.add_argument("-i", "--input", required=True, help="TSV.gz da IdentityInsertMetrics.py")
    parser.add_argument("-o", "--output", required=True, help="File immagine di output (PNG, PDF, ecc.)")
    parser.add_argument("--y-margin", type=float, default=0.05, help="Percentuale margine extra sull'asse Y per Qscore")
    parser.add_argument("--max-reads-per-bin", type=int, default=10_000_000, help="Numero massimo di read per bin")
    args = parser.parse_args()

    main(args.input, args.output, args.y_margin, args.max_reads_per_bin)
