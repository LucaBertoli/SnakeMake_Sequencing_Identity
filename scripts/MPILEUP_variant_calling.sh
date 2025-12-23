bam=$1

bcftools mpileup \
  -f /home/db/hg38_bwa2/hg38_chr1-22-XYM/Homo_sapiens_assembly38_noalt.fasta \
  -a AD,DP \
  -Q 0 \
  -q 0 \
  -B \
  --threads 10 \
  $bam \
| bcftools call \
  -m \
  --ploidy 1 \
  --threads 10 \
  -Oz \
  -o $(basename $bam .bam)_mpileup_permissive.vcf.gz

tabix $(basename $bam .bam)_mpileup_permissive.vcf.gz

/home/tools/bcftools-1.19/bcftools norm $(basename $bam .bam)_mpileup_permissive.vcf.gz -f /home/db/hg38_bwa2/hg38_chr1-22-XYM/Homo_sapiens_assembly38_noalt.fasta -m - both | bgzip > $(basename $bam .bam)_mpileup_permissive_split.vcf.gz
tabix $(basename $bam .bam)_mpileup_permissive_split.vcf.gz