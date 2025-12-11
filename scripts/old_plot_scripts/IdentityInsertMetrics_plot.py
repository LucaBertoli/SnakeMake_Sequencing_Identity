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
from matplotlib.colors import to_rgba
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
    df_long["insert_bin"] = pd.Categorical(df_long["insert_bin"], ordered=True)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(9, 4))

    # Violin plot in primo piano
    sns.violinplot(data=df_long, 
                   x="insert_bin", 
                   y="PhredScore", 
                   hue="QualityType",
                   split=True, 
                   inner="quart", 
                   width=0.8,
                   palette={
                       "mean_qual_log": to_rgba("#36E2EE", 0.6),
                       "Q_identity": to_rgba("#43F05F", 0.6)
                   },
                   scale="width", 
                   linewidth=1,
                   cut=0, 
                   ax=ax
                   )

    sns.set_style("whitegrid", {'axes.grid': True})
    ax.yaxis.grid(True, linestyle="--", alpha=0.4, linewidth=0.6)
    ax.xaxis.grid(True, linestyle="--", alpha=0.4, linewidth=0.6)
    ax.set_xlabel("Insert size (bp)")
    ax.set_ylabel("Phred-scaled Qscore")
    ax.set_ylim(0, df_long["PhredScore"].max() * (1 + y_margin))
    ax.set_title(f"Insert size on Sequencer and Identity Qscore ({max_reads_per_bin:,} reads per bin)")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend_.remove()  # rimuove la legenda automatica di Seaborn
    label_map = {
        "mean_qual_log": "Average Sequencer Qscore",
        "Q_identity": "Identity-based Qscore"
    }
    labels = [label_map.get(l, l) for l in labels]
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),  # sotto il grafico
        ncol=2,                       # due colonne
        frameon=True
    )
    ax.tick_params(axis="x", rotation=45)
    ax.set_yticks(range(0, int(df_long["PhredScore"].max()) + 5, 5))

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
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
