#!/bin/bash -l
#SBATCH -D /home/dgold/
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgold@ucdavis.edu 
#SBATCH -o /home/dgold/QUAN-o%j.txt
#SBATCH -e /home/dgold/QUAN-e%j.txt
#SBATCH -J QUAN
#SBATCH -t 40:00:00 

cd /home/dgold/25_Oysters/2_Cufflinks
module load cufflinks # v.2.2.1

for i in ../1_Mapping/*bam; do \
 filename=$(basename "$i")
 cuffquant \
 Oyster_Final_Gene_Annotations.combined.gtf \
 --library-type fr-firststrand \
 -o CuffQuant_${filename%Aligned.sortedByCoord.out.bam} -p 16 $i ;
done

# Extract CBX files from subfolders and append sample names 
for i in */abundances.cxb; do\
 j=${i%/abundances.cxb}
 k=${j#CuffQuant_}
 cp $i $k.cbx ;
done

# Rename cxb files based on descriptions

mv SRR15846133.cbx C_51-1.cbx
mv SRR15846134.cbx C_51-2.cbx
mv SRR15846135.cbx C_51-3.cbx
mv SRR15846137.cbx C_38-1.cbx
mv SRR15846140.cbx E_51-1.cbx
mv SRR15846141.cbx E_51-2.cbx
mv SRR15846142.cbx E_51-3.cbx
mv SRR15846143.cbx E_38-1.cbx
mv SRR15846144.cbx E_38-2.cbx
mv SRR15846145.cbx E_38-3.cbx
mv SRR15846148.cbx C_38-2.cbx
mv SRR15846149.cbx C_38-3.cbx

#####################################
# Create count tables with Cuffnorm 
#####################################

cuffnorm Oyster_Final_Gene_Annotations.combined.gtf \
C_38-1.cbx,C_38-2.cbx,C_38-3.cbx \
C_51-1.cbx,C_51-2.cbx,C_51-3.cbx \
E_38-1.cbx,E_38-2.cbx,E_38-3.cbx \
E_51-1.cbx,E_51-2.cbx,E_51-3.cbx \
-L E_38,E_51,C_38,C_51 \
-o cuffnorm --library-type fr-firststrand

###################################################
# Calculate differential expression with Cuffdiff
###################################################

cuffdiff Oyster_Final_Gene_Annotations.combined.gtf \
C_38-1.cbx,C_38-2.cbx,C_38-3.cbx \
C_51-1.cbx,C_51-2.cbx,C_51-3.cbx \
E_38-1.cbx,E_38-2.cbx,E_38-3.cbx \
E_51-1.cbx,E_51-2.cbx,E_51-3.cbx \
-L E_38,E_51,C_38,C_51 \
-o cuffdiff -p 16 --library-type fr-firststrand
