import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
import os

colori_matplotlib = [
    "r", "b", "g", "c", "m", "y", "k",
    "royalblue", "limegreen", "tomato", "deepskyblue", "darkorchid", "goldenrod", "dimgray",
    "dodgerblue", "springgreen", "salmon", "skyblue", "mediumorchid", "khaki", "lightgray",
    "cornflowerblue", "palegreen", "lightcoral", "lightcyan", "lightpink", "lightyellow", "silver",
    "steelblue", "springgreen", "lightsalmon", "powderblue", "mediumvioletred", "gold", "gainsboro",
    "cornflowerblue", "chartreuse", "darkorange", "cadetblue", "mediumslateblue", "darkkhaki", "beige",
    "mediumblue", "mediumseagreen", "crimson", "darkturquoise", "darkslategray", "olive", "floralwhite",
    "mediumslateblue", "chartreuse", "darkorange", "cadetblue", "mediumslateblue", "darkkhaki", "beige",
]


label = ["SALUS", "Illumina NovaSeqX", "Element AVITI", "GeneMind", "MGI G400"]

fig, ax = plt.subplots(figsize=(10, 8))

for i in range(1, len(sys.argv)):
    dir = sys.argv[i]

    df = pd.read_csv(dir, sep="\t")

    # converto in percentuale
    df["identity_with_ins_filtered"] = df["identity_with_ins_filtered"] * 100

    # imposto ordine dei bin (senza 0-2)
    ordered_bins = ["3-17", "18-29", "30+"]

    # forza BQ_bin a stringa
    df["BQ_bin"] = df["BQ_bin"].astype(str)

    # tieni solo i valori attesi
    df = df[df["BQ_bin"].isin(ordered_bins)]

    # assegna come categoria ordinata
    df["BQ_bin"] = pd.Categorical(df["BQ_bin"], categories=ordered_bins, ordered=True)

    # se non ci sono righe valide, passa oltre
    if df.empty:
        print(f"⚠️ Nessun bin valido trovato in {dir}, skippo.")
        continue

    color = colori_matplotlib[(i - 1) % len(colori_matplotlib)]

    ax.plot(
        df['BQ_bin'],
        df['identity_with_ins_filtered'],
        marker='o',
        color=color,
        label=label[i-1],
        linewidth=0.5,
        markersize=5
    )

ax.set_title('Stratified Identity by Base Quality Bin')
ax.set_xlabel('Base Quality Bin')
ax.set_ylabel('Identity (%)')

ax.set_yticks(np.arange(0, 105, 1))
ax.set_ylim(83, 103)
ax.grid(alpha=0.3)
ax.legend()
plt.tight_layout()
plt.savefig("identity_plot_binned_identity_with_ins_filtered.png")
