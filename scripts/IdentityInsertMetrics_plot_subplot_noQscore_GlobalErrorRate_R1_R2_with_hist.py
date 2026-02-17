#!/usr/bin/env python3
"""
Creates a  2x2 insert-size stratified plot:
- Line 1: Read 1 (R1)
- Line 2: Read 2 (R2)
- Left Column: sequencer-assigned quality (mean_qual_log)
- Right Column: error rate (mismatches + indels) + Count Histogram (background)
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

def main(tsv_path, output_path, y_margin=0.05):
    # Dizionari separati per R1 e R2
    bin_data = {"R1": {}, "R2": {}}

    # --------------------------------------------------------
    # Lettura TSV.gz (TUTTI i frammenti, senza limiti)
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

            bases = int(row["aligned_read_length"])
            errors = int(row["mismatch_filtered"]) + int(row["indel_filtered"])

            bin_data[read_type][insert_bin].append({
                "insert_bin": insert_bin,
                "mean_qual_log": float(row["mean_qual_log"]),
                "errors": errors,
                "bases": bases
            })

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
    # Conteggi per istogrammi
    # --------------------------------------------------------
    counts_R1 = df_R1.groupby("insert_bin").size().reset_index(name="Count") if df_R1 is not None else pd.DataFrame()
    counts_R2 = df_R2.groupby("insert_bin").size().reset_index(name="Count") if df_R2 is not None else pd.DataFrame()

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
    def plot_row(df, counts_df, ax_q, ax_e, title_prefix):
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

        # ===== Error rate + ISTOGRAMMA (SOLO SU AX_E) =====
        df_error = (
            df.groupby("insert_bin")[["errors", "bases"]]
            .sum()
            .reset_index()
        )
        df_error["Error_Rate"] = df_error["errors"] / df_error["bases"]

        # Asse principale: Error Rate
        # Secondo asse: Conteggi (background)
        ax_eb = ax_e.twinx()
        ax_eb.grid(False)  # ← NO GRIDLINES sul secondo asse
        bin_width = 40

        # ISTOGRAMMA BACKGROUND (solo su ax_eb)
        ax_eb.bar(
            counts_df["insert_bin"],
            counts_df["Count"],
            width=bin_width,
            color="lightgray",
            alpha=0.3,
            zorder=1,
            label="Counts"
        )

        # LINEA Error Rate (primo piano)
        ax_e.plot(
            df_error["insert_bin"],
            df_error["Error_Rate"],
            color=colors["Error_Rate"],
            marker="o",
            markersize=5,
            linewidth=2,
            zorder=4
        )
        
        ax_e.set_ylabel("Error Rate (%)")
        ax_e.set_xlabel("Insert size (bp)")
        ax_e.tick_params(axis="x", rotation=45)
        ax_e.set_title(f"{title_prefix} – Error Rate + Counts")
        ax_e.set_ylim(0, 0.02)
        ax_e.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda y, _: f"{y*100:.1f}")
        )
        ax_eb.set_ylabel("Fragment count")

        ax_q.grid(True, axis="both", linestyle="--", alpha=0.7)
        ax_e.grid(True, axis="both", linestyle="--", alpha=0.7)

    # --------------------------------------------------------
    # Plot R1 / R2
    # --------------------------------------------------------
    plot_row(df_R1, counts_R1, axes[0, 0], axes[0, 1], "Read 1")
    plot_row(df_R2, counts_R2, axes[1, 0], axes[1, 1], "Read 2")

    # --------------------------------------------------------
    # Legenda globale
    # --------------------------------------------------------
    handles = [
        Line2D([0], [0], color=colors["mean_qual_log"], lw=6),
        Line2D([0], [0], color=colors["Error_Rate"], lw=6, marker="o"),
        Line2D([0], [0], color="black", lw=2, marker="o"),
        plt.Rectangle((0,0),1,1, color="lightgray", alpha=0.3),
    ]
    labels = ["Qscore", "Error Rate", "Mediana", "R1 Counts"]

    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.03),
        ncol=5,
        frameon=False
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✅ Plot salvato in: {output_path}")
    print(f"📊 R1: {len(df_R1):,} | R2: {len(df_R2):,} frammenti totali")
    print(f"📊 Bin: R1={len(bin_data['R1'])}, R2={len(bin_data['R2'])}")

# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Violin plot Qscore ed Error Rate separati per R1/R2 + Counts"
    )
    parser.add_argument("-i", "--input", required=True, help="TSV.gz da IdentityInsertMetrics.py")
    parser.add_argument("-o", "--output", required=True, help="File immagine di output")
    parser.add_argument("--y-margin", type=float, default=0.05)
    args = parser.parse_args()
    main(args.input, args.output, args.y_margin)
