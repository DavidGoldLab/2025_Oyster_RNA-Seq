#!/bin/bash -l
#SBATCH -D /home/dgold/BLAST
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgold@ucdavis.edu 
#SBATCH -o /home/dgold/STAR-o%j.txt
#SBATCH -e /home/dgold/STAR-e%j.txt
#SBATCH -J STAR
#SBATCH -t 36:00:00 #60 minutes

cd /home/dgold/24_Oysters/1_Mapping

module load star

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ./genome \
--genomeFastaFiles ../0_Data/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna \
--sjdbGTFfile ../0_Data/genomic.gff \
--sjdbOverhang 100

for i in ../0_Data/*.fq; do
  filename=$(basename "$i")
  STAR --runThreadN 16 --genomeDir ./genome --readFilesIn ./$i \
  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
  --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
  --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${filename%.clean.fq} ;\
done
