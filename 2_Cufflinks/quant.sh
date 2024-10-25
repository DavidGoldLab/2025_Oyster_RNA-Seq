#!/bin/bash -l
#SBATCH -D /home/dgold/BLAST
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgold@ucdavis.edu 
#SBATCH -o /home/dgold/QUANT-o%j.txt
#SBATCH -e /home/dgold/QUANT-e%j.txt
#SBATCH -J QUANT
#SBATCH -t 48:00:00

cd /home/dgold/24_Oysters/2_Cufflinks
module load cufflinks # v.2.2.1

for i in ../1_Mapping/*bam; do \
 filename=$(basename "$i")
 cuffquant \
 Oyster_Final_Gene_Annotations.combined.gtf \
 --library-type fr-firststrand \
 -o CuffQuant_${filename%Aligned.sortedByCoord.out.bam} -p 16 $i ;
done

