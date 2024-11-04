#!/bin/bash -l
#SBATCH -D /home/dgold/BLAST
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgold@ucdavis.edu 
#SBATCH -o /home/dgold/EXAMPLE-o%j.txt
#SBATCH -e /home/dgold/EXAMPLE-e%j.txt
#SBATCH -J EXAMPLE
#SBATCH -t 3:60:00 #60 minutes

cd /home/dgold/25_Oysters/1_Mapping
module load star

for i in ../0_Data/*.fq; do
  filename=$(basename "$i")
  STAR --runThreadN 16 --genomeDir ./genome --readFilesIn ./$i \
  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
  --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
  --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${filename%.clean.fq} ;\
done
