#!/usr/bin/env python3
"""
Crea due subplot affiancati:
- Sinistra: qualità assegnata dal sequenziatore (mean_qual_log)
- Destra: qualità calcolata dall'identità (Phred-scaled)
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

def main(tsv_path, output_path, max_q_identity=45, y_margin=0.05, max_reads_per_bin=100_000):
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
                identity = float(row["identity_filtered_with_indels"])
                q_identity = -10 * math.log10(max(1 - identity, 1e-6))
                q_identity = min(q_identity, max_q_identity)
                bin_data[insert_bin].append({
                    "insert_bin": insert_bin,
                    "mean_qual_log": float(row["mean_qual_log"]),
                    "Q_identity": q_identity
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
    colors = {"mean_qual_log": "#36E2EE", "Q_identity": "#43F05F"}

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
    # SUBPLOT 2 – Qualità calcolata dall'identità
    # ================================
    ax2 = axes[1]
    sns.violinplot(
        data=df, x="insert_bin", y="Q_identity",
        color=colors["Q_identity"],
        inner="quart", scale="width", cut=0, ax=ax2
    )
    ax2.set_ylabel("Qscore")
    ax2.set_xlabel("Insert size (bp)")
    ax2.set_ylim(0, df["Q_identity"].max() * (1 + y_margin))
    ax2.tick_params(axis="x", rotation=45)
    ax2.set_title("Identity-based Qscore")
    ax1.set_yticks(range(0, int(df["mean_qual_log"].max()) + 5, 5))
    ax2.set_yticks(range(0, int(df["Q_identity"].max()) + 5, 5))
    ax1.grid(True, axis="both", linestyle="--", alpha=0.7)
    ax2.grid(True, axis="both", linestyle="--", alpha=0.7)

    medians2 = df.groupby("insert_bin")["Q_identity"].median()
    ax2.plot(
        range(len(medians2)),
        medians2.values,
        color="black",
        marker="o",
        markersize=3,
        linewidth=1,
        label="Mediana"
    )

    # ================================
    # LEGENDA GLOBALE
    # ================================
    handles = [
        Line2D([0], [0], color=colors["mean_qual_log"], lw=6),
        Line2D([0], [0], color=colors["Q_identity"], lw=6),
        Line2D([0], [0], color="black", lw=2, marker="o")
    ]
    labels = ["Qualità assegnata", "Qualità calcolata", "Mediana"]
    fig.legend(handles, labels,
               loc="upper center", bbox_to_anchor=(0.5, -0.05),
               ncol=3, frameon=False)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✅ Plot salvato in: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Violin plot separati per qualità assegnata e calcolata")
    parser.add_argument("-i", "--input", required=True, help="TSV.gz da IdentityInsertMetrics.py")
    parser.add_argument("-o", "--output", required=True, help="File immagine di output (PNG, PDF, ecc.)")
    parser.add_argument("--max-q-identity", type=float, default=45, help="Limite superiore per Q_identity")
    parser.add_argument("--y-margin", type=float, default=0.05, help="Percentuale margine extra sull'asse Y")
    parser.add_argument("--max-reads-per-bin", type=int, default=100_000, help="Numero massimo di read per bin")
    args = parser.parse_args()

    main(args.input, args.output, args.max_q_identity, args.y_margin, args.max_reads_per_bin)
