# Questo script estrae le reads da un file BAM tramite samtools,
# legge la stringa di qualità (campo QUAL di ogni allineamento),
# calcola la qualità media per ogni read (media dei punteggi Phred),
# e salva l'output in un file .tsv.gz con due colonne: read name e qualità media (valore intero).
# La qualità viene calcolata convertendo i caratteri ASCII della stringa QUAL in punteggi Phred (ord(char) - 33).

#nohup bash ../Compute_Q-measured.sh 20250513/Q_measured/alignment.rg.clipped.bam qualita_media_per_read_SALUS_WGS_no_capping.tsv.gz &

bam=$1

if [ -z "$bam" ] || [ -z "$2" ]; then
    echo "Usage: $0 <input.bam> <output.tsv.gz>"
    exit 1
fi
samtools view $1 | \
python3 -c "
import sys
for line in sys.stdin:
    fields = line.strip().split('\t')
    qual = fields[10]
    if qual:
        scores = [ord(c) - 33 for c in qual]
        mean_score = sum(scores) / len(scores)
        print(f'{fields[0]}\t{mean_score:.2f}')
" | gzip > $2
