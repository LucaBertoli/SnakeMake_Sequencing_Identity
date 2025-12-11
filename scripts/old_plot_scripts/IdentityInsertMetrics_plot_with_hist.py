#!/usr/bin/env python3
"""
Crea violin plot bipartiti mostrando:
- Qualità assegnata dal sequenziatore (mean_qual_log)
- Qualità calcolata dall'identità (Phred-scaled)
per intervalli di 50 bp di insert size (150–800 bp).
Mantiene solo le prime N read per bin e segnala quando un bin raggiunge il limite.

Input:
    - TSV.gz generato da IdentityInsertMetrics.py
Output:
    - File PNG
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
    bin_full_flags = set()  # traccia i bin già pieni

    # Lettura gzip riga per riga
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

            # Stop se tutti i bin hanno raggiunto max_reads_per_bin
            if all(len(v) >= max_reads_per_bin for v in bin_data.values()):
                break

    # Unisci tutte le righe
    all_rows = [row for rows in bin_data.values() for row in rows]
    if not all_rows:
        raise ValueError("Nessuna read trovata nei bin richiesti.")
    df = pd.DataFrame(all_rows)

    # Dati per violin plot
    df_long = df.melt(id_vars="insert_bin", value_vars=["mean_qual_log","Q_identity"],
                      var_name="QualityType", value_name="PhredScore")

    # Conteggio read per bin
    count_df = df.groupby("insert_bin").size().reset_index(name="count")

    # --- Plot ---
    sns.set(style="whitegrid", font_scale=1.1)
    fig, ax1 = plt.subplots(figsize=(12,6))
    ax2 = ax1.twinx()

    # Barplot sullo sfondo
    sns.barplot(data=count_df, x="insert_bin", y="count", color="lightgray", alpha=0.8, ax=ax1, zorder=1)
    ax1.set_ylabel("Numero di read per bin", color="gray")
    ax1.tick_params(axis="y", colors="gray")
    ax1.grid(False)
    ax1.set_ylim(0, count_df["count"].max() * (1 + y_margin))

    # Violin plot in primo piano
    sns.violinplot(data=df_long, x="insert_bin", y="PhredScore", hue="QualityType",
                   split=True, inner="quart",
                   palette={"mean_qual_log": "#4C72B0", "Q_identity": "#55A868"},
                   scale="width", cut=0, ax=ax2, zorder=2)
    ax2.set_xlabel("Insert size (binned every 50 bp, 150–800 bp)")
    ax2.set_ylabel("Phred-scaled quality score")
    ax2.set_ylim(0, df_long["PhredScore"].max() * (1 + y_margin))
    total_reads = len(df)
    ax2.set_title(f"Qualità assegnata vs osservata per bin di insert size (150–800 bp)\n(n = {max_reads_per_bin:,} read per bin)")
    ax2.tick_params(axis="x", rotation=45)

    # Legenda esterna sotto il grafico
    bar_handle = [Line2D([0],[0], color="lightgray", lw=6, alpha=0.8)]
    handles_violin, labels_violin = ax2.get_legend_handles_labels()
    ax2.legend(
        handles=bar_handle + handles_violin,
        labels=["Numero di read", "Qualità assegnata", "Qualità calcolata"],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=3,
        frameon=False
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✅ Plot salvato in: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Violin + barplot insert size quality con prime N read per bin")
    parser.add_argument("-i","--input", required=True, help="TSV.gz da IdentityInsertMetrics.py")
    parser.add_argument("-o","--output", required=True, help="File immagine di output (PNG, PDF, ecc.)")
    parser.add_argument("--max-q-identity", type=float, default=45, help="Limite superiore per Q_identity")
    parser.add_argument("--y-margin", type=float, default=0.05, help="Percentuale margine extra sulle assi Y")
    parser.add_argument("--max-reads-per-bin", type=int, default=100_000, help="Numero massimo read per bin")
    args = parser.parse_args()

    main(args.input, args.output, args.max_q_identity, args.y_margin, args.max_reads_per_bin)
