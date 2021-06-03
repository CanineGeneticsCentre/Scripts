#!/bin/bash

## Command to Run :  bash slurm_vcf2excel.sh <Disease Status File>

VCF_FILE='/rds/project/rds-Qr3fy2NTCy0/Data/ens_WGS_219_VEP.vcf.gz'
DISEASE_STATUS=$1
#SEG_SCORE=$2
#EFFECT_SCORE=1

[[ -z "$DISEASE_STATUS" ]] && { echo "ERROR: No disease status file provided for this run"; exit 1; }
if [ ! -e $VCF_FILE ]
then 
  echo "ERROR - Unable to find VCF file. Please check and try again - ${VCF_FILE}";
  exit 1;
fi

dos2unix ${DISEASE_STATUS}

sbatch <<EOT
#!/bin/bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -J vcf2excel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --mail-type=ALL
#SBATCH -p skylake

#SBATCH -o $HOME/rds/hpc-work/logs/job-%j.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

perl /rds/project/rds-Qr3fy2NTCy0/Software/Git/Scripts/perl/vcf2excel.pl --default --vcf ${VCF_FILE} --status_file ${DISEASE_STATUS}
EOT

