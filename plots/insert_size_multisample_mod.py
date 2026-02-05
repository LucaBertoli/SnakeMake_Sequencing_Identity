# by lucab

##plot istogramma dell'inserto, va lanciato con python dando in input un numero arbitrario di file ".output" prodotti da Picard CollectInsertSizeMetrics.

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys



## MODIFICARE ls SEGUENTe LISTa IN BASE AL NUMERO DI FILE IN INPUT, UTILIZZARE COLORI DI MATPLOTLIB
colori_matplotlib = [
    "#E69F00", 
    "#009E73",
    "#004488",
    "#CC79A7"
]
###### Open tsv files
##DA MODIFICARE PER IMPOSTARE LA DIMENSIONE DELLA FIGURA
fig, ax1 = plt.subplots(figsize=(8, 4))

hist_data_list = []
metrics_list = []


# Leggi i file e calcola i limiti degli assi
for i in range(1, len(sys.argv)):
    file_path = sys.argv[i]

    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Estrai i valori MEDIAN_INSERT_SIZE e MEAN_INSERT_SIZE
    for line in lines:
        if line.startswith('## METRICS CLASS'):
            metrics_line = line
        elif line.startswith('## HISTOGRAM'):
            histogram_start_idx = lines.index(line) + 1
            break

    metrics_data = lines[lines.index(metrics_line) + 2]
    metrics_values = metrics_data.split('\t')
    median_insert_size = int(metrics_values[0])
    mean_insert_size = float(metrics_values[5])
    mode_insert_size = float(metrics_values[1])
    metrics_list.append((median_insert_size, mean_insert_size, mode_insert_size))

    # Estrai i dati per l'istogramma
    histogram_data = lines[histogram_start_idx + 1:]  # Salta l'intestazione dell'istogramma
    insert_sizes = []
    counts = []

    for data in histogram_data:
        if not data.strip():
            continue
        try:
            size, count = map(int, data.split())
            insert_sizes.append(size)
            counts.append(count)
        except ValueError:
            raise ValueError(f"Errore di conversione dei dati in valori numerici nella riga: {data.strip()}")

    # for data in histogram_data:
    #     if not data.strip():
    #         continue
    #     try:
    #         parts = data.split()
    #         size = int(parts[0])
    #         if len(parts) == 2:
    #             count = int(parts[1])
    #         elif len(parts) == 3:
    #             count = int(parts[1]) + int(parts[2])
    #         else:
    #             continue  # Skip if the data doesn't match expected patterns
    #         insert_sizes.append(size)
    #         counts.append(count)
    #     except ValueError:
    #         # Handle or log the error if needed
    #         continue

    hist_data_list.append((insert_sizes, counts))

    # Determina i limiti degli assi X e Y
    x_min = min(min(insert_sizes) for insert_sizes, counts in hist_data_list)
    x_max = max(max(insert_sizes) for insert_sizes, counts in hist_data_list)
    y_max = max(max(counts) for insert_sizes, counts in hist_data_list)

    # Crea gli istogrammi con bin di 10 basi
    ##DA MODIFICARE PER IMPOSTARE LA DIMENSIONE DEI BIN
    bin_width = 1
    bins = range(x_min, x_max + bin_width, bin_width)
    line_legends = []
    line_labels = []

    ##DA COMMENTARE O MODIFICARE PER IMPOSTARE IL NOME DELLE CONDIZIONI PLOTTATE, USARE COLORI DELLA LEGGENDA IMPOSTATA IN PRECEDENZA 
    line_legends.append(Patch(color='#004488'))
    line_labels.append("Illumina NovaSeq X")
    line_legends.append(Patch(color='#E69F00'))
    line_labels.append("Element AVITI")
    line_legends.append(Patch(color='#009E73'))
    line_labels.append("GeneMind SurfSeq 5000")
    line_legends.append(Patch(color='#CC79A7'))
    line_labels.append("MGI T1")


    for j, (insert_sizes, counts) in enumerate(hist_data_list):

        #DECOMMENTARE SE SI VOGLIONO LE LINEE VERTICALI DI MEDIA, MEDIANA E MODA
        # Plot delle linee tratteggiate per le medie
        # mean_line = ax1.axvline(metrics_list[j][1], color=colori_matplotlib[j], linestyle='solid', linewidth=1.5, label=f'Mean Insert Size: {metrics_list[j][1]:.2f}')
        # line_legends.append(mean_line)
        line_labels.append(f'Mean: {metrics_list[j][1]:.2f}')

        # Plot delle linee punteggiate per le mediane
        # median_line = ax1.axvline(metrics_list[j][0], color=colori_matplotlib[j], linestyle='dashed', linewidth=1.5, label=f'Median Insert Size: {metrics_list[j][0]}')
        # line_legends.append(median_line)
        line_labels.append(f'Median: {metrics_list[j][0]}')

        # Plot delle linee tratteggiate per le medie
        # mode_line = ax1.axvline(metrics_list[j][2], color=colori_matplotlib[j], linestyle='dashdot', linewidth=1.5, label=f'Mean Insert Size: {metrics_list[j][1]:.2f}')
        # line_legends.append(mode_line)
        line_labels.append(f'Mode: {metrics_list[j][2]:.2f}')


        ##PER AVERA IL GRAFICO CON UNA LINEA SENZA RIEMPIMENTO, PER ESEMPIO SE SI PLOTTANO MOLTI CAMPIONI
        ax1.hist(insert_sizes, bins=bins, weights=counts, edgecolor=colori_matplotlib[j], histtype='step', linewidth=1.5, fill=False)

        ##PER AVERE IL GRAFICO CON L'AREA RIEMPITA, PER ESEMPIO SE PLOTTO UNO O DUE CAMPIONI, MODIFICARE IL PARAMETRO DI TRASPARENZA ALPHA A PIACIMENTO
        # ax1.hist(insert_sizes, bins=bins, weights=counts, color=colori_matplotlib[j], histtype='step', fill=True, alpha=0.4)


    ##DA MODIFICARE PER IMPOSTARE I LIMITI DEGLI ASSI IN BASE AI DATI
    ax1.set_xlim([0, 1000])
    ax1.set_ylim([0, 2100000])

    ax1.set_xlabel('Insert Size (bp)')
    ax1.set_ylabel('Number of Fragments')
    ax1.ticklabel_format(style='plain')

    ##DA MODIFICARE PER IMPOSTARE IL TITOLO DEL GRAFICO
    # ax1.set_title('Insert Size Histogram')

    plt.legend(line_legends, line_labels, loc='upper right', bbox_to_anchor=(1, 1), fontsize='small', prop={'size': 10})


##DA MODIFICARE PER IMPOSTARE IL NOME DEL GRAFICO IN OUTPUT
plt.savefig('insert_size_hist.png')
plt.show()
