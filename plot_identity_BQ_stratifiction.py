import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
import sys
import numpy as np
import os

colori_matplotlib = [
    "r", "b", "g", "c", "m", "y", "k",  # Colori di base
]

legend_config = [
    Patch(color='r', label='seqWell IlluminaPrep_A02'),
    # Patch(color='g', label='Euroclone'),
    # Patch(color='r', label='19/02/2024'),
]

fig, ax = plt.subplots(figsize=(10, 8))

for i in range(1, len(sys.argv)):
    dir = sys.argv[i]
    print(dir)

    df = pd.read_csv(dir, sep="\t")
    df = df.dropna()
    df["identity_filtered"] = df["identity_filtered"] * 100

    color = colori_matplotlib[(i-1) % len(colori_matplotlib)]
    label = os.path.basename(os.path.dirname(dir))
    ax.plot(df['BaseQuality'], df['identity_filtered'], marker='o', color=color, label=label, linewidth=0.5, markersize=5)

ax.set_title('Stratified Identity by Base Quality')
ax.set_xlabel('Base Quality (Phred)')
ax.set_ylabel('Identity (%)')

ax.set_xticks(np.arange(0, max(df['BaseQuality']) + 1, 5))
ax.set_yticks(np.arange(0, 105, 5))

ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("identity_plot.png")
