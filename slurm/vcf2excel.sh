#!/bin/bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -J vcf2excel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --mail-type=ALL
#SBATCH -p skylake
#SBATCH --array 1-41

#SBATCH -o $HOME/rds/hpc-work/logs/job-%A_%a.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

# Exported variables...
# DISEASE_STATUS, FASTA, VCF_DIR, SCRIPTS, REF

CHR=`cut -f 1 -d':' ${FASTA}/${REF}-genomicsDB.intervals | cut -f 1 -d' ' | cut -d'_' -f 1 | sort -n | uniq | head -${SLURM_ARRAY_TASK_ID} | tail -1`

#perl /rds/project/rds-Qr3fy2NTCy0/Software/Git/Scripts/perl/vcf2excel.pl --default --vcf ${VCF_FILE} --status_file ${DISEASE_STATUS}
#perl /rds/project/rds-Qr3fy2NTCy0/Software/Git/Scripts/perl/vcf2excel.pl --default --vcf ${VCF}/vcf_chr/ens_WGS_219-chr\${SLURM_ARRAY_TASK_ID}.vcf.gz --status_file ${DISEASE_STATUS}

echo "perl ${SCRIPTS}/perl/vcf2excel.pl --default --vcf ${VCF_DIR}/${REF}-chr${CHR}.ann.vcf.gz --status_file ${DISEASE_STATUS}"
perl ${SCRIPTS}/perl/vcf2excel.pl --default --vcf ${VCF_DIR}/${REF}-chr${CHR}.ann.vcf.gz --status_file ${DISEASE_STATUS}