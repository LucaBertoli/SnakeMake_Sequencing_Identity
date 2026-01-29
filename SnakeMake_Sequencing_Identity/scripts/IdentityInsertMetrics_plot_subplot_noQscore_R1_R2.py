#!/usr/bin/env python3
"""
Crea un plot 2x2:
- Riga 1: read1 (R1)
- Riga 2: read2 (R2)
- Colonna sinistra: qualità assegnata dal sequenziatore (mean_qual_log)
- Colonna destra: identità delle read (identity_filtered_with_indels)
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


def main(tsv_path, output_path, y_margin=0.05, max_reads_per_bin=100_000):
    # Dizionari separati per R1 e R2
    bin_data_R1 = {}
    bin_data_R2 = {}
    bin_full_flags = {"R1": set(), "R2": set()}

    # --- Lettura TSV.gz riga per riga ---
    with gzip.open(tsv_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            mode = row["read_type"].strip().upper()
            if mode not in ("R1", "R2"):  # ignora 'S' o altro
                continue

            insert_size = int(row["insert_size"])
            insert_bin = (insert_size // 50) * 50
            if insert_bin < 150 or insert_bin > 800:
                continue

            target_dict = bin_data_R1 if mode == "R1" else bin_data_R2

            if insert_bin not in target_dict:
                target_dict[insert_bin] = []

            if len(target_dict[insert_bin]) < max_reads_per_bin:
                identity = float(row["identity_filtered_with_indels"])
                target_dict[insert_bin].append({
                    "insert_bin": insert_bin,
                    "mean_qual_log": float(row["mean_qual_log"]),
                    "Identity": identity
                })
                if len(target_dict[insert_bin]) == max_reads_per_bin and insert_bin not in bin_full_flags[mode]:
                    print(f"⚡ Bin {insert_bin} ({mode}) ha raggiunto il limite di {max_reads_per_bin} read.")
                    bin_full_flags[mode].add(insert_bin)

            # Stop se tutti i bin R1 e R2 hanno raggiunto il limite
            if all(
                len(v) >= max_reads_per_bin
                for d in [bin_data_R1, bin_data_R2]
                for v in d.values()
            ):
                break

    # --- Crea DataFrame per R1 e R2 ---
    def make_df(data_dict):
        all_rows = [row for rows in data_dict.values() for row in rows]
        return pd.DataFrame(all_rows) if all_rows else None

    df_R1 = make_df(bin_data_R1)
    df_R2 = make_df(bin_data_R2)

    if df_R1 is None and df_R2 is None:
        raise ValueError("Nessuna read trovata per R1 o R2 nei bin richiesti.")

    # --- Stile plot ---
    sns.set(style="whitegrid", font_scale=1.1)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=False)
    colors = {"mean_qual_log": "#36E2EE", "Identity": "#43F05F"}

    def plot_pair(df, ax_q, ax_i, title_prefix):
        if df is None or df.empty:
            ax_q.set_visible(False)
            ax_i.set_visible(False)
            return

        # --- Subplot Qscore ---
        sns.violinplot(
            data=df, x="insert_bin", y="mean_qual_log",
            color=colors["mean_qual_log"],
            inner="quart", scale="width", cut=0, ax=ax_q
        )
        ax_q.set_ylabel("Average read Qscore")
        ax_q.set_xlabel("Insert size (bp)")
        ax_q.set_ylim(0, df["mean_qual_log"].max() * (1 + y_margin))
        ax_q.tick_params(axis="x", rotation=45)
        ax_q.set_title(f"{title_prefix} – Sequencer Qscore")

        medians_q = df.groupby("insert_bin")["mean_qual_log"].median()
        ax_q.plot(
            range(len(medians_q)),
            medians_q.values,
            color="black",
            marker="o",
            markersize=3,
            linewidth=1,
            label="Mediana"
        )

        # --- Subplot Identity ---
        sns.violinplot(
            data=df, x="insert_bin", y="Identity",
            color=colors["Identity"],
            inner="quart", scale="width", cut=0, ax=ax_i
        )
        ax_i.set_ylabel("Read identity (%)")
        ax_i.set_xlabel("Insert size (bp)")
        ax_i.tick_params(axis="x", rotation=45)
        ax_i.set_title(f"{title_prefix} – Alignment Identity")
        ax_i.set_ylim(0.98, 1.001)
        ax_i.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y*100:.1f}"))

        medians_i = df.groupby("insert_bin")["Identity"].median()
        ax_i.plot(
            range(len(medians_i)),
            medians_i.values,
            color="black",
            marker="o",
            markersize=3,
            linewidth=1,
            label="Mediana"
        )

        # --- Griglia ---
        ax_q.grid(True, axis="both", linestyle="--", alpha=0.7)
        ax_i.grid(True, axis="both", linestyle="--", alpha=0.7)

    # --- Plot R1 (riga 1) ---
    plot_pair(df_R1, axes[0, 0], axes[0, 1], "Read 1")

    # --- Plot R2 (riga 2) ---
    plot_pair(df_R2, axes[1, 0], axes[1, 1], "Read 2")

    # --- Legenda globale ---
    handles = [
        Line2D([0], [0], color=colors["mean_qual_log"], lw=6),
        Line2D([0], [0], color=colors["Identity"], lw=6),
        Line2D([0], [0], color="black", lw=2, marker="o")
    ]
    labels = ["Qualità assegnata", "Identità", "Mediana"]
    fig.legend(
        handles, labels,
        loc="upper center", bbox_to_anchor=(0.5, -0.01),
        ncol=3, frameon=False
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✅ Plot salvato in: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Violin plot 2x2 per qualità assegnata e identità (R1/R2)")
    parser.add_argument("-i", "--input", required=True, help="TSV.gz da IdentityInsertMetrics.py")
    parser.add_argument("-o", "--output", required=True, help="File immagine di output (PNG, PDF, ecc.)")
    parser.add_argument("--y-margin", type=float, default=0.05, help="Percentuale margine extra sull'asse Y per Qscore")
    parser.add_argument("--max-reads-per-bin", type=int, default=100_000, help="Numero massimo di read per bin per ciascun tipo")
    args = parser.parse_args()

    main(args.input, args.output, args.y_margin, args.max_reads_per_bin)
