import sys
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

#funzione per convertire percentuale di identità in punteggio Phred
# La formula è: Phred = -10 * log10(1 - percentuale_identità)
#La identità a 1 viene considerata come 0.999 per evitare log(0)
def percent_identity_to_phred(identity_percent):
    identity = min(float(identity_percent), 0.999)
    error_prob = 1.0 - identity
    phred = -10 * np.log10(error_prob)
    return round(phred, 2)

# Funzione per caricare il file di identità e convertire le identità in punteggi Phred con la funzione precedente (percent_identity_to_phred)
def load_identity_file(identita_file):
    with gzip.open(identita_file, 'rt') as f:
        df_id = pd.read_csv(f, sep='\t', skiprows=1, header=None, usecols=[4], names=['percent_identity'])
    phred_estimated = df_id['percent_identity'].apply(percent_identity_to_phred).values
    return phred_estimated

# Funzione per caricare il file dei valori Phred misurati dal sequenziatore
# Si assume che il file sia in formato tab-separated e contenga i valori Phred nella seconda colonna
# Il file deve avere un'intestazione, quindi si salta la prima riga
def load_phred_file(qvalori_file):
    with gzip.open(qvalori_file, 'rt') as f:
        df_q = pd.read_csv(f, sep='\t', skiprows=1, header=None, usecols=[1], names=['phred_measured'])
    return df_q['phred_measured'].values

# Funzione per convertire la lista di valori Phred in istogramma
def compute_histogram(values, bins):
    counts, _ = np.histogram(values, bins=bins)
    return counts

# Inizio main
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python Identity_plot.py identita.tsv.gz qvalori.tsv.gz [output_prefix]")
        sys.exit(1)

    identita_file = sys.argv[1]
    qvalori_file = sys.argv[2]
    output_prefix = sys.argv[3] if len(sys.argv) >= 4 else "phred_hist"

    

    #Se non esiste il fine CSV contenente l'istogramma delle Phred quality, lo calcolo e lo salvo.
    if not os.path.exists(f"{output_prefix}.csv"):
        #richiamo le funzioni per caricare i file di identità e Phred
        phred_estimated = load_identity_file(identita_file)
        phred_measured = load_phred_file(qvalori_file)

        # Definisco i bin dell'istogramma
        bins = np.arange(0, 60.01, 0.01)  # 0.01 step → 6000 bin

        # Calcolo gli istogrammi con le due liste di valori Phred
        counts_estimated = compute_histogram(phred_estimated, bins)
        counts_measured = compute_histogram(phred_measured, bins)

        # Salvo i valori degli istogrammi in un DataFrame e lo scrivo nel CSV in output
        # Il DataFrame conterrà le colonne: bin_start, bin_end, count_estimated (qualità calcolate sulla base dell'identità), count_measured (qualità assegnate dal sequenziatore)
        df_hist = pd.DataFrame({
            'bin_start': bins[:-1],
            'bin_end': bins[1:],
            'count_estimated': counts_estimated,
            'count_measured': counts_measured
        })
        df_hist.to_csv(f'{output_prefix}.csv', index=False)

    # Carico CSV per plotting
    df_hist = pd.read_csv(f'{output_prefix}.csv')

    # Raggruppo i bin 0.01 → 1
    df_hist['bin_center'] = (df_hist['bin_start'] + df_hist['bin_end']) / 2
    df_hist['bin_group'] = np.round(df_hist['bin_center']).astype(int)

    df_plot = df_hist.groupby('bin_group').agg({
        'count_estimated': 'sum',
        'count_measured': 'sum'
    }).reset_index().sort_values('bin_group')

    # # Plot istogramma
    # plt.figure(figsize=(12, 6))
    # width = 0.4

    # plt.bar(df_plot['bin_group'] - width/2, df_plot['count_estimated'], width=width, label='Estimated by Alignment', color='orange')
    # plt.bar(df_plot['bin_group'] + width/2, df_plot['count_measured'], width=width, label='Measured by Sequencer', color='steelblue')

    # plt.xlabel('Phred Score')
    # plt.ylabel('Read Count')
    # plt.title(f'Distribution of Estimated vs Measured Phred Scores ({output_prefix})')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f'{output_prefix}.png', dpi=300)
    # plt.close()

    
    # Subplot: istogrammi affiancati
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

    # Estimated by Alignment
    axes[0].bar(df_plot['bin_group'], df_plot['count_estimated'], color='orange', width=0.8)
    axes[0].set_title('Estimated by Alignment')
    axes[0].set_xlabel('Phred Score')
    axes[0].set_ylabel('Read Count')

    # Measured by Sequencer
    axes[1].bar(df_plot['bin_group'], df_plot['count_measured'], color='steelblue', width=0.8)
    axes[1].set_title('Measured by Sequencer')
    axes[1].set_xlabel('Phred Score')

    plt.suptitle(f'Distribution of Phred Scores ({output_prefix})')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'{output_prefix}_subplot.png', dpi=300)
    plt.close()



    # # Preparo dati per violin plot
    # df_long = pd.DataFrame({
    #     'phred_score': np.concatenate([phred_estimated, phred_measured]),
    #     'Score_Type': (['Estimated by Alignment'] * len(phred_estimated)) + (['Measured by Sequencer'] * len(phred_measured))
    # })

    # # Plot violin
    # plt.figure(figsize=(6, 6))
    # sns.violinplot(
    #     data=df_long,
    #     x='Score_Type',
    #     y='phred_score',
    #     palette={'Measured by Sequencer': 'steelblue', 'Estimated by Alignment': 'orange'}
    # )
    # plt.xlabel('Score Type')
    # plt.ylabel('Phred Score')
    # plt.title(f'Distribution of Estimated vs Measured Phred Scores ({output_prefix})')
    # plt.tight_layout()
    # plt.savefig(f"{output_prefix}_violin_plot.png", dpi=300)
    # plt.close()