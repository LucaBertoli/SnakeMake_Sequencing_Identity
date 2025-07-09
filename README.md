Questa repository contiene la collezione di script per il calcolo dell'identità su read sequenziate. 
Lo script principale à:
- Run_Identity_Workflow.py
questo script prende in input un BAM, un VCF ed un nome di file in output. Il bam deve essere generato allineando le read grezze senza eseguire trimming e senza eseguire deduplicazione/ricalibrazione. Il VCF deve essere normalizzato con bcftools norm, da utilizzare il GoldSeq GIAB se si usa NA12878 o il VCF chiamato sullo specifico campione. l'output name verrà utilizzato come basename dei file in output.