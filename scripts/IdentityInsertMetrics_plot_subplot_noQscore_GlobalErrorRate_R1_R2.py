#!/usr/bin/env python3
"""
Crea un plot 2x2:
- Riga 1: Read 1 (R1)
- Riga 2: Read 2 (R2)
- Colonna sinistra: qualità assegnata dal sequenziatore (mean_qual_log)
- Colonna destra: error rate (mismatches + indels)

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


# ============================================================
# Main
# ============================================================

def main(tsv_path, output_path, y_margin=0.05, max_reads_per_bin=5_000_000):
    # Dizionari separati per R1 e R2
    bin_data = {"R1": {}, "R2": {}}
    bin_full_flags = {"R1": set(), "R2": set()}

    # --------------------------------------------------------
    # Lettura TSV.gz
    # --------------------------------------------------------
    with gzip.open(tsv_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            read_type = row["read_type"].strip().upper()
            if read_type not in ("R1", "R2"):
                continue

            insert_size = int(row["insert_size"])
            insert_bin = (insert_size // 50) * 50
            if insert_bin < 150 or insert_bin > 800:
                continue

            if insert_bin not in bin_data[read_type]:
                bin_data[read_type][insert_bin] = []

            if len(bin_data[read_type][insert_bin]) < max_reads_per_bin:
                bases = int(row["aligned_read_length"])
                errors = int(row["mismatch_filtered"]) + int(row["indel_filtered"])

                bin_data[read_type][insert_bin].append({
                    "insert_bin": insert_bin,
                    "mean_qual_log": float(row["mean_qual_log"]),
                    "errors": errors,
                    "bases": bases
                })

                if len(bin_data[read_type][insert_bin]) == max_reads_per_bin \
                        and insert_bin not in bin_full_flags[read_type]:
                    print(f"⚡ Bin {insert_bin} ({read_type}) ha raggiunto {max_reads_per_bin} read.")
                    bin_full_flags[read_type].add(insert_bin)

            # Stop se tutti i bin R1 e R2 sono pieni
            if all(
                len(v) >= max_reads_per_bin
                for d in bin_data.values()
                for v in d.values()
            ):
                break

    # --------------------------------------------------------
    # DataFrame
    # --------------------------------------------------------
    def make_df(d):
        rows = [r for rows in d.values() for r in rows]
        return pd.DataFrame(rows) if rows else None

    df_R1 = make_df(bin_data["R1"])
    df_R2 = make_df(bin_data["R2"])

    if df_R1 is None and df_R2 is None:
        raise ValueError("Nessuna read trovata per R1 o R2 nei bin richiesti.")

    # --------------------------------------------------------
    # Stile plot
    # --------------------------------------------------------
    sns.set(style="whitegrid", font_scale=1.1)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey=False)

    colors = {
        "mean_qual_log": "#36E2EE",
        "Error_Rate": "#43F05F"
    }

    # --------------------------------------------------------
    # Funzione di plotting per R1 / R2
    # --------------------------------------------------------
    def plot_row(df, ax_q, ax_e, title_prefix):
        if df is None or df.empty:
            ax_q.set_visible(False)
            ax_e.set_visible(False)
            return

        # ===== Qscore (violin) =====
        sns.violinplot(
            data=df,
            x="insert_bin",
            y="mean_qual_log",
            color=colors["mean_qual_log"],
            inner="quart",
            scale="width",
            cut=0,
            ax=ax_q
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
            linewidth=1
        )

        # ===== Error rate =====
        df_error = (
            df.groupby("insert_bin")[["errors", "bases"]]
            .sum()
            .reset_index()
        )
        df_error["Error_Rate"] = df_error["errors"] / df_error["bases"]

        ax_e.plot(
            df_error["insert_bin"],
            df_error["Error_Rate"],
            color=colors["Error_Rate"],
            marker="o",
            markersize=5,
            linewidth=1
        )
        ax_e.set_ylabel("Error Rate (%)")
        ax_e.set_xlabel("Insert size (bp)")
        ax_e.tick_params(axis="x", rotation=45)
        ax_e.set_title(f"{title_prefix} – Read Error Rate")
        ax_e.set_ylim(0, 0.02)
        ax_e.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda y, _: f"{y*100:.1f}")
        )

        ax_q.grid(True, axis="both", linestyle="--", alpha=0.7)
        ax_e.grid(True, axis="both", linestyle="--", alpha=0.7)

    # --------------------------------------------------------
    # Plot R1 / R2
    # --------------------------------------------------------
    plot_row(df_R1, axes[0, 0], axes[0, 1], "Read 1")
    plot_row(df_R2, axes[1, 0], axes[1, 1], "Read 2")

    # --------------------------------------------------------
    # Legenda globale
    # --------------------------------------------------------
    handles = [
        Line2D([0], [0], color=colors["mean_qual_log"], lw=6),
        Line2D([0], [0], color=colors["Error_Rate"], lw=6),
        Line2D([0], [0], color="black", lw=2, marker="o")
    ]
    labels = ["Qualità assegnata", "Error Rate", "Mediana"]

    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.03),
        ncol=3,
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
        description="Violin plot Qscore ed Error Rate separati per R1/R2"
    )
    parser.add_argument("-i", "--input", required=True, help="TSV.gz da IdentityInsertMetrics.py")
    parser.add_argument("-o", "--output", required=True, help="File immagine di output (PNG, PDF, ecc.)")
    parser.add_argument("--y-margin", type=float, default=0.05)
    parser.add_argument("--max-reads-per-bin", type=int, default=5_000_000)

    args = parser.parse_args()
    main(args.input, args.output, args.y_margin, args.max_reads_per_bin)
