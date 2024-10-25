#!/bin/bash -l
#SBATCH -D /home/dgold/BLAST
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgold@ucdavis.edu 
#SBATCH -o /home/dgold/CUFF-o%j.txt
#SBATCH -e /home/dgold/CUFF-e%j.txt
#SBATCH -J CUFF
#SBATCH -t 24:60:00 #60 minutes
cd /home/dgold/24_Oysters/2_Cufflinks
module load cufflinks # v.2.2.1
for i in ../1_Mapping/*Aligned.sortedByCoord.out.bam; do \
 filename=$(basename "$i")
 cufflinks --library-type fr-firststrand $i -p 16 -o ./${filename%Aligned.sortedByCoord.out.bam} ; \
done
