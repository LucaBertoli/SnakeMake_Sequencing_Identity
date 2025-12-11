#python scripts/plot_identity_BQ_stratifiction_Q_seq_to_real_Q.py ILLUMINA_NovaSeqX_WGS/NA12878_KAPA/IDENTITY_NA12878_NovaSeqX_KAPA_GIAB_basequal_strat.tsv.gz.stats Element_WGS/RUN_4/IDENTITY_NA12878_Element_run4_GIAB_basequal_strat.tsv.gz.stats GeneMind_WGS/IDENTITY_NA12878_GeneMind_GIAB_basequal_strat.tsv.gz.stats MGI/SnakeMake_Sequencing_Identity_results_old/NA12878_KAPA/IDENTITY_NA12878_KAPA_base_quality_stratification.tsv.gz.stats

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

fig, ax = plt.subplots(figsize=(8, 8))
for i in range(1, len(sys.argv)):
    dir = sys.argv[i]

    df = pd.read_csv(dir, sep="\t")
    df = df.dropna()
    df["identity_with_ins_filtered"] = df["identity_with_ins_filtered"].apply(lambda x: -10*np.log10(1 - x))
    
    color = colori_matplotlib[(i-1) % len(colori_matplotlib)]
    # label = os.path.basename(os.path.dirname(os.path.abspath(dir)))
    ax.plot(df['BaseQuality'], df['identity_with_ins_filtered'], marker='o', color=color, label=label[i-1], linewidth=0.5, markersize=5)

ax.set_title('Sequencer-based Qscore vs identity-based Qscore (mismatch + indel)')
ax.set_xlabel('Sequencer Qscore (Phred)')
ax.set_ylabel('Identity-based Qscore (Phred)')

ax.set_xticks(np.arange(0, max(df['BaseQuality']) + 1, 5))
ax.set_yticks(np.arange(0, max(df['identity_with_ins_filtered']) + 1, 5))

ax.plot([-5, max(df['BaseQuality']) + 5], [-5, max(df['BaseQuality']) + 5], 'k--', linewidth=1)

ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("identity_plot_Q_seq_to_real_Q.png")
