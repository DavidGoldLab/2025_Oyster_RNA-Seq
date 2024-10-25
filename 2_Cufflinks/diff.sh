#!/bin/bash -l
#SBATCH -D /home/dgold/BLAST
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgold@ucdavis.edu 
#SBATCH -o /home/dgold/DIFF-o%j.txt
#SBATCH -e /home/dgold/DIFF-e%j.txt
#SBATCH -J DIFF
#SBATCH -t 36:00:00 #60 minutes
cd /home/dgold/24_Oysters/2_Cufflinks
module load cufflinks # v.2.2.1
cuffdiff Oyster_Final_Gene_Annotations.combined.gtf \
C_38-1.cbx,C_38-2.cbx,C_38-3.cbx \
C_51-1.cbx,C_51-2.cbx,C_51-3.cbx \
E_38-1.cbx,E_38-2.cbx,E_38-3.cbx \
E_51-1.cbx,E_51-2.cbx,E_51-3.cbx \
-L E_38,E_51,C_38,C_51 \
-o cuffdiff -p 16 --library-type fr-firststrand
