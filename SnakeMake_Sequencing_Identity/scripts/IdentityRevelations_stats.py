# """
# IdentityRevelations_stats.py

# Questo script analizza un file TSV contenente statistiche di identità tra sequenze (calcolato con lo script IdentityRevelations.py), calcolando e stampando:
# - Il numero totale di sequenze presenti nel file di input.
# - Media  della percentuale di identità filtrata (identity_filtered_with_indels, che contiene l'identità escluse SNV ed INDEL presenti in un VCF).
# - Il numero e la percentuale di sequenze con 0, 1, 2, 3 e più di 3 errori (dove gli errori sono la somma delle basi corrispondenti a mismatches, inserzioni e delezioni filtrate).

# Utilizzo:
#     python IdentityRevelations_stats.py <input_identity_tsv_file>

# Argomenti:
#     input_identity_tsv_file: File TSV di input calcolato con IdentityRevelations.py con colonne tra cui 'identity_filtered_with_indels', 'mismatches_filtered', 'ins_filtered', 'del_filtered'.

# Output:
#     Statistiche riassuntive stampate a schermo.
# """


#!/usr/bin/env python3

import pandas as pd
import sys

if len(sys.argv) != 2:
    print("Usage: python IdentityRevelations_stats_optimized.py <input_identity_tsv_file>")
    sys.exit(1)

input_file = sys.argv[1]

# Inizializzazione dei conteggi delle read totali e della somma delle identità
total_reads = 0
sum_identity = 0.0

# Inizializzazione dei bucket per il conteggio degli errori
# I bucket sono: 0 errori, 1 errore, 2 errori, 3 errori, e più di 3 errori
error_buckets = {
    0: 0,
    1: 0,
    2: 0,
    3: 0,
    ">3": 0
}

# inizializzo una lista con le colonne da leggere dal file TSV in input
usecols = [
    'identity_filtered_with_indels',
    'mismatches_filtered',
    'ins_filtered_bases',
    'del_filtered_bases'
]

# Lettura del file TSV in input in chunks per ottimizzare l'uso della memoria
# Si assume che il file sia in formato tab-separated (TSV) e che le colonne desiderate siano presenti
# Si utilizza il parametro chunksize per leggere il file in blocchi di 100.000 righe alla volta
chunksize = 100_000
for chunk in pd.read_csv(input_file, sep='\t', usecols=usecols, chunksize=chunksize):
    #calcolo la colonna total_errors come somma delle colonne 'mismatches_filtered', 'ins_filtered_bases' e 'del_filtered_bases'.
    chunk['total_errors'] = (
        chunk['mismatches_filtered'] +
        chunk['ins_filtered_bases'] +
        chunk['del_filtered_bases']
    )

    # Sommo i valodi di identità di ogni read in ogni chunk
    sum_identity += chunk['identity_filtered_with_indels'].sum()

    # Conto il numero totale di read in questo chunk
    total_reads += len(chunk)

    # Sommo gli errori nel chunk
    error_counts = chunk['total_errors'].value_counts()

    # Aggiorno i conteggi degli errori nei bucket
    # Se il numero di errori nella read è 0, lo aggiungo al bucket 0, se è 1 lo aggiungo al bucket 1, e così via
    for err, count in error_counts.items():
        if err == 0:
            error_buckets[0] += count
        elif err == 1:
            error_buckets[1] += count
        elif err == 2:
            error_buckets[2] += count
        elif err == 3:
            error_buckets[3] += count
        elif err > 3:
            error_buckets[">3"] += count

#stampo i risultati in stdout
print(f"\nNumero totale di sequenze: {total_reads}")
mean_identity = sum_identity / total_reads if total_reads > 0 else float('nan')
print(f"Media della identità filtrata: {mean_identity:.6f}\n")

# Output errori
for key in [0, 1, 2, 3, ">3"]:
    count = error_buckets[key]
    perc = (count / total_reads) * 100 if total_reads > 0 else 0
    print(f"Sequenze con {key} errori: {count} ({perc:.4f}%)")