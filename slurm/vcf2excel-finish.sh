#!/bin/bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -J vcf2excel2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:05:00
#SBATCH --mail-type=ALL
#SBATCH -p cclake

#SBATCH -o logs/job-%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#samples=`wc -l ${DISEASE_STATUS} | cut -f1 -d' '`
samples=`grep -vi 'omit' ${DISEASE_STATUS} | wc -l | tr -d ' '`
col_effect=`expr 2 \* $samples + 14`
col_seg=`expr 2 \* $samples + 4`

file=`ls *_filtered_and_sorted_for_excel.txt | head -1`
head -1 $file > header.txt

mkdir output
mv *_filtered_and_sorted_for_excel.txt *_vcf2excel_missing_effects.out *_vcf2excel_command_log.out output/


for i in `seq 1 5`; do
  (cat header.txt; awk -F"\t" -v c=$col_effect -v effect=$i '($c == effect)' output/*_filtered_and_sorted_for_excel.txt) > vcf2excel-effect${i}.txt
  LINECOUNT=`wc -l vcf2excel-effect${i}.txt | cut -f1 -d' '`
  if [[ $LINECOUNT == 1 ]]; then
    rm -f vcf2excel-effect${i}.txt
  fi
done
