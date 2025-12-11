#grafico per plottare la relazione tra i q score del sequenziatore e i q score corrispondenti all'identità filtrata per snv/indel e binnata.
# python ../scripts/plot_identity_BQ_stratifiction_binned_Q_seq_to_real_Q.py.py ../SALUS_WGS/IDENTITY_run_workflow_20250617_WGS_2_1_GIAB_trimmed_basequal_strat.tsv.gz.binned.stats ../ILLUMINA_NovaSeqX_WGS/NA12878_KAPA/IDENTITY_NA12878_NovaSeqX_KAPA_GIAB_basequal_strat.tsv.gz.binned.stats ../Element_WGS/RUN_4/IDENTITY_NA12878_Element_run4_GIAB_basequal_strat.tsv.gz.binned.stats ../GeneMind_WGS/IDENTITY_NA12878_GeneMind_GIAB_basequal_strat.tsv.gz.binned.stats ../MGI/SnakeMake_Sequencing_Identity_results/NA12878_KAPA/IDENTITY_NA12878_KAPA_base_quality_stratification.tsv.gz.binned.stats 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

colori_matplotlib = [
    "r", "b", "g", "c", "m", "y", "k",
    "royalblue", "limegreen", "tomato", "deepskyblue", "darkorchid", "goldenrod", "dimgray",
    "dodgerblue", "springgreen", "salmon", "skyblue", "mediumorchid", "khaki", "lightgray",
    "cornflowerblue", "palegreen", "lightcoral", "lightcyan", "lightpink", "lightyellow", "silver",
    "steelblue", "springgreen", "lightsalmon", "powderblue", "mediumvioletred", "gold", "gainsboro",
    "cornflowerblue", "chartreuse", "darkorange", "cadetblue", "mediumslateblue", "darkkhaki", "beige",
]

label = ["Illumina NovaSeqX", "Element AVITI", "GeneMind", "MGI G400"]

fig, ax = plt.subplots(figsize=(10, 8))

# ordine dei bin e loro range numerico
ordered_bins = ["3-17", "18-29", "30+"]
bin_ranges = {"3-17": (3, 17), "18-29": (18, 29), "30+": (30, 45)}  # range attesi per ciascun bin

# colori pastello chiari per i range (rettangoli)
range_colors = {
    "3-17": "salmon",
    "18-29": "lightblue",
    "30+": "lightgreen"
}

for i in range(1, len(sys.argv)):
    dir = sys.argv[i]
    df = pd.read_csv(dir, sep="\t")

    # converto in Phred
    df["identity_with_ins_filtered"] = df["identity_with_ins_filtered"].apply(lambda x: -10 * np.log10(1 - x))
    print(sys.argv[i])
    print(df["identity_with_ins_filtered"])
    df["BQ_bin"] = df["BQ_bin"].astype(str)

    # filtro e ordino
    df = df[df["BQ_bin"].isin(ordered_bins)]
    if df.empty:
        print(f"⚠️ Nessun bin valido in {dir}, skippo.")
        continue

    df["BQ_bin"] = pd.Categorical(df["BQ_bin"], categories=ordered_bins, ordered=True)
    color = colori_matplotlib[(i - 1) % len(colori_matplotlib)]

    ax.plot(
        df['BQ_bin'],
        df['identity_with_ins_filtered'],
        marker='o',
        color=color,
        label=label[i - 1],
        linewidth=0.8,
        markersize=5
    )

# --- Aggiunta dei quadrati colorati pastello ---
for i, b in enumerate(ordered_bins):
    y1, y2 = bin_ranges[b]
    ax.fill_between(
        [i - 0.4, i + 0.4],
        y1, y2,
        color=range_colors[b],
        alpha=0.25,
        zorder=0,
        edgecolor='gray',
        linewidth=0.8
    )

# sistemazione assi
ax.set_title('Sequencer-based Qscore vs identity-based Qscore binned (mismatch + indel)')
ax.set_xlabel('Sequencer Qscore (Phred)')
ax.set_ylabel('Identity-based Qscore (Phred)')

# asse X categoriale
ax.set_xticks(range(len(ordered_bins)))
ax.set_xticklabels(ordered_bins)

ax.set_ylim(0, 45)
ax.set_xlim(-0.5, len(ordered_bins) - 0.5)

# griglia e legenda
ax.grid(alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig("identity_plot_binned_identity_with_ins_filtered_Q_seq_to_real_Q.png", dpi=300)
plt.show()
