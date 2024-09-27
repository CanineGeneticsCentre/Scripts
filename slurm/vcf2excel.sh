#!/bin/bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -J vcf2excel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --mail-type=ALL
#SBATCH -p cclake

#SBATCH -o logs/job-%A_%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

# Exported variables...
# DISEASE_STATUS, VCF_DIR, SCRIPTS, REF, CONFIG

CHR=`cut -f 1 -d':' ${FASTA}/${REF}-genomicsDB.intervals | cut -f 1 -d' ' | cut -d'_' -f 1 | uniq | head -${SLURM_ARRAY_TASK_ID} | tail -1`

if [[ ${#CHR} -lt 4 ]] ; then
  CHR="chr"${CHR}
fi

echo "perl ${SCRIPTS}/../perl/vcf2excel.pl --default --vcf ${VCF_DIR}/${REF}-${CHR}.ann.vcf.gz --status_file ${DISEASE_STATUS}"
perl ${SCRIPTS}/../perl/vcf2excel.pl --default --vcf ${VCF_DIR}/${REF}-${CHR}.ann.vcf.gz --status_file ${DISEASE_STATUS}
