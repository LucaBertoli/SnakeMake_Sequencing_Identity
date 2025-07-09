Questa repository contiene la collezione di script per il calcolo dell'identità su read sequenziate. 
Lo script principale à:
- Run_Identity_Workflow.py
    - questo script prende in input un BAM, un VCF ed un nome di file in output. Il bam deve essere generato allineando le read grezze senza eseguire trimming e senza eseguire deduplicazione/ricalibrazione. Il VCF deve essere normalizzato con bcftools norm, da utilizzare il GoldSeq GIAB se si usa NA12878 o il VCF chiamato sullo specifico campione. l'output_name verrà utilizzato come basename dei file in output.

Lo script Run_Identity_Workflow.py richiamerà i seguenti step:

- intersezione con bcftools intersect del BAM con un file BED di regioni ad alta confidenza (es. HG001_GRCh38_1_22_v4.2.1_benchmark.bed), da fare soltanto una volta per generare il BAM nelle high confidence regions (HCR)

- IdentityRevelations.py
    - questo script calcola a partire dal BAM e dal VCF un file tsv.gz contenente i mismatch, indel e valori di identità di ogni read (escludendo secondary, supplementay e unmapped fragments). Per le specifiche vedi sotto.

- IdentityRevelations_stats.py
    - questo step calcola alcune statistiche partendo dal file tsv.gz generato nello step precedente, in particolare numero di read, identità media (media delle identità di ciascuna read), numero di read con 0,1,2,3,>3 mismatch. Lo script stampa tutto in un file .stats.

- Scatter_BAM.py
    - lo script calcolara le qualità medie delle read in un file BAM. Le qualità considerate sono le qualità in phread score delle basi allineate (esclude le softclipped). L'output risulterà un file "_mean_quality.tsv.gz". Lo script può anche separare le read in due file BAM distinti in base alla loro qualità media.

- Identity_plot.py
    - lo script prende in input il tsv.gz con le identità (generato da IdentityRevelations.py) ed il tsv.gz con i phread score del BAM (generato da Scatter_BAM.py), converte le identità in Q-score con la formula: Q-score=-10*log10(1-identity). Per necessità la identità a 1 viene approssimata a 0.999. Infine crea un plot con i Q-score assegnati dal sequenziatore ed quelli calcolati sulla base dell'allineamento.

- identity_calculalculation_all_types.sh
    - lo script prende in input il BAM, VCF, fasta e le regioni target e un output_name. Genera un secondo report calcolando, tramite samtools mpileup, l'identità dei dati allineati. Per fare questo calcola il numero totale di basi allineate che non sono in posizioni varianti (nel vcf) e conta il numero di basi mismatch/indel. L'output è un file .identity.all.types.pileup.report.

# Esempio per lanciare il workflow:
    bash Run_Identity_Workflow.sh $PATH/start_sorted.bam $PATH/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz run_workflow_NOME_CAMPIONE >> IDENTITY_run_workflow_NOME_CAMPIONE.log



# 🧬 IdentityRevelations.py

Questo script è il cuore del calcolo dell'identità sulla base dell'allineamento. 
Calcola l'identità di ogni singola read allineata in un file BAM (escludendo unmapped, secondary e supplementary alignments) con il supporto di un file VCF contenente varianti (SNV e INDEL). Lo script calcola quattro metriche di identità per ciascuna read, tutte espresse come valori compresi tra 0 e 1, dove 1 indica identità perfetta con il riferimento:

- **identity**:  
  Identità basata solo sul numero totale di mismatch (ottenuti dal tag `MD`).  
  Formula:  
  identity = 1 - {{mismatches_total}/{aligned_length_total}}

- **identity_filtered**:  
  Come sopra, ma escludendo mismatch che coincidono con varianti SNV note nel VCF.  

  identity_filtered = 1 - {{mismatches_filtered}/{aligned_length_total}}
  

- **identity_with_indels**:  
  Identità considerando sia i mismatch che la lunghezza totale di inserzioni e delezioni.  
  
  identity_with_indels = 1 - {{mismatches_total} + {inserted_bases} + {deleted_bases}}/{aligned_length_total}
  

- **identity_filtered_with_indels**:  
  Come sopra, ma escludendo mismatch e indel sovrapposti a varianti note (SNV e INDEL).  

  identity_filtered_with_indels = 1 - {{mismatches_filtered} + {filtered_inserted_bases} + {filtered_deleted_bases}}/{aligned_length_total}


**Note**:
- `aligned_length_total` rappresenta la lunghezza totale dell’allineamento della read.


Quest'ultima è l'identità che viene utilizzata negli step successivi della pipeline.

Utilizzando la libreria pysam, script estrae mismatch tramite il tag MD e indel dalla CIGAR. Questi mismatch vengono confrontati con due dizionari contenenti rispettivamente le SNV e le INDEL. Se i mismatch o le indel identificati nella CIGAR e nel tag MD corrispondono a una variante nei due dizionari, quest'ultimi non vengono considerati come mismatch nel calcolo dell'identità. Il risultato è un file .tsv.gz contenente metriche per ciascuna read:

    Conteggio di mismatch totali e filtrati

    Conteggio e lunghezza di inserzioni e delezioni (filtrate e non)

    Lunghezza della read allineata

    Identità osservata nelle modalità sopracitate