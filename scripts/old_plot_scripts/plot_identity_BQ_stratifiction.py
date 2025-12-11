import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
import sys
import numpy as np
import os

colori_matplotlib = [
    "r", "b", "g", "c", "m", "y", "k",  # Colori di base
    "royalblue", "limegreen", "tomato", "deepskyblue", "darkorchid", "goldenrod", "dimgray",  # Colori brillanti
    "dodgerblue", "springgreen", "salmon", "skyblue", "mediumorchid", "khaki", "lightgray",  # Altri colori brillanti
    "cornflowerblue", "palegreen", "lightcoral", "lightcyan", "lightpink", "lightyellow", "silver",  # Colori chiari
    "steelblue", "springgreen", "lightsalmon", "powderblue", "mediumvioletred", "gold", "gainsboro",  # Altri colori chiari
    "cornflowerblue", "chartreuse", "darkorange", "cadetblue", "mediumslateblue", "darkkhaki", "beige",  # Altri colori
    "mediumblue", "mediumseagreen", "crimson", "darkturquoise", "darkslategray", "olive", "floralwhite",  # Altri colori
    "mediumslateblue", "chartreuse", "darkorange", "cadetblue", "mediumslateblue", "darkkhaki", "beige",  # Altri colori
]

label=["Illumina NovaSeq X", "Element AVITI", "GeneMind", "MGI G400"]

fig, ax = plt.subplots(figsize=(10, 8))
for i in range(1, len(sys.argv)):
    dir = sys.argv[i]

    df = pd.read_csv(dir, sep="\t")
    df = df.dropna()
    df["identity_filtered"] = df["identity_filtered"] * 100

    color = colori_matplotlib[(i-1) % len(colori_matplotlib)]
    # label = os.path.basename(os.path.dirname(os.path.abspath(dir)))
    ax.plot(df['BaseQuality'], df['identity_filtered'], marker='o', color=color, label=label[i-1], linewidth=0.5, markersize=5)

ax.set_title('Stratified Identity by Base Quality (mismatch + indel)')
ax.set_xlabel('Base Quality (Phred)')
ax.set_ylabel('Identity (%)')

ax.set_xticks(np.arange(0, max(df['BaseQuality']) + 1, 5))
# ax.set_yticks(np.arange(0, 105, 5))

#y and x lim
ax.set_yticks(np.arange(0, 101, 5))
ax.set_xticks(np.arange(0, 45, 5))
# ax.set_ylim(98, 101)
# ax.set_xlim(20, 46)

ax.grid()
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig("identity_plot_only_mismatch.png")
