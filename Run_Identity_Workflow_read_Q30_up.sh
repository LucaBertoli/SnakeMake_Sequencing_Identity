#script per lanciare in modo sequenziale il workflow per l'analisi dell'identità.
#in input richiede un bam allineato, un vcf ed un nume dei file in output.

bam_temporaneo=$1
vcf=$2
output_name=$3

bam="$(dirname "$bam_temporaneo")/$(basename "$bam_temporaneo" .bam)_HCR.bam"

echo "Directory di lavoro: $(pwd)"
echo "BAM:" $bam_temporaneo
echo "BAM finale:"  $(basename $bam)
echo "VCF:" $vcf
echo "Nome file di output:" $output_name

echo "Inizio esecuzione: $(date)" 
start_time=$(date +%s)

scripts_folder="/home/comparazione_kit/calcolo_identità_sequenze_Sequenziatori"
FASTA="/home/db/hg38_bwa2/hg38_chr1-22-XYM/Homo_sapiens_assembly38_noalt.fasta"
BED="/home/db/hg38_bwa2/hg38_chr1-22-XYM/Homo_sapiens_assembly38_noalt.bed"


# bedtools intersect -abam $bam_temporaneo -b /home/comparazione_kit/calcolo_identità_sequenze_Sequenziatori/HG001_GRCh38_1_22_v4.2.1_benchmark.bed -f 1 -u > $bam
# samtools index $bam

echo "Step3: INIZIO estrazione delle read con una qualità media mappata superiore a 30..."
python $scripts_folder/Scatter_BAM.py $bam BAM
samtools index $(basename $bam .bam)_reads_30_up.bam
echo "Step3: FINE estrazione delle read con una qualità media mappata superiore a 30..."

echo "Step1: INIZIO calcolo dell'identità sul bam..."
python $scripts_folder/IdentityRevelations.py $(basename $bam .bam)_reads_30_up.bam $vcf IDENTITY_${output_name}.tsv.gz
echo "Step1: FINE calcolo dell'identità sul bam..."

echo "Step2: INIZIO calcolo delle statistiche sull'identità del bam..."
python $scripts_folder/IdentityRevelations_stats.py IDENTITY_${output_name}.tsv.gz >> IDENTITY_${output_name}.tsv.gz.stats
echo "Step2: FINE calcolo delle statistiche sull'identità del bam..."

echo "Step3: INIZIO calcolo delle qualità medie delle read sequenziate..."
python $scripts_folder/Scatter_BAM.py $(basename $bam .bam)_reads_30_up.bam QUAL
echo "Step3: FINE calcolo delle qualità medie delle read sequenziate..."

echo "Step4: INIZIO plot delle qualità stimate dal sequenziatore e calcolate in base all'identità..."
python $scripts_folder/Identity_plot.py IDENTITY_${output_name}.tsv.gz $(dirname $bam)/$(basename $bam .bam)_reads_30_up.bam_mean_quality.tsv.gz IDENTITY_$output_name
echo "Step4: FINE plot delle qualità stimate dal sequenziatore e calcolate in base all'identità..."

echo "Step5: INIZIO calcolo dell'identità utilizzando tutti i dati assieme..."
bash $scripts_folder/identity_calculalculation_all_types.sh $FASTA $(basename $bam .bam)_reads_30_up.bam $vcf $BED IDENTITY_${output_name}
echo "Step5: FINE calcolo dell'identità utilizzando tutti i dati assieme..."


#### I seguenti comandi non vanno utilizzati a meno che non si voglia calcolare l'identità separando read 1 e read 2.
# echo "Step1.1: INIZIO calcolo dell'identità sul bam read1..."
# python $scripts_folder/IdentityRevelations.py $bam $vcf IDENTITY_${output_name}_read1.tsv.gz read1
# echo "Step1.1: FINE calcolo dell'identità sul bam read1..."

# echo "Step2.1: INIZIO calcolo delle statistiche sull'identità del bam read1..."
# python $scripts_folder/IdentityRevelations_stats.py IDENTITY_${output_name}_read1.tsv.gz >> IDENTITY_${output_name}_read1.tsv.gz.stats
# echo "Step2.1: FINE calcolo delle statistiche sull'identità del bam read1..."

# echo "Step1.2: INIZIO calcolo dell'identità sul bam read2..."
# python $scripts_folder/IdentityRevelations.py $bam $vcf IDENTITY_${output_name}_read2.tsv.gz read2
# echo "Step1.2: FINE calcolo dell'identità sul bam read2..."

# echo "Step2.2: INIZIO calcolo delle statistiche sull'identità del bam read2..."
# python $scripts_folder/IdentityRevelations_stats.py IDENTITY_${output_name}_read2.tsv.gz >> IDENTITY_${output_name}_read2.tsv.gz.stats
# echo "Step2.2: FINE calcolo delle statistiche sull'identità del bam read2..."

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Fine esecuzione: $(date)"
echo "Tempo totale di esecuzione: $execution_time secondi"