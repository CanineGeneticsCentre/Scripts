#!/bin/bash

## Command to Run :  bash slurm_vcf2excel.sh <Disease Status File>

VCF_FILE='/rds/project/rds-Qr3fy2NTCy0/Data/ens_WGS_219_VEP.vcf.gz'
DISEASE_STATUS=$1
#SEG_SCORE=$2
#EFFECT_SCORE=1

VCF=`dirname $VCF_FILE`

[[ -z "$DISEASE_STATUS" ]] && { echo "ERROR: No disease status file provided for this run"; exit 1; }
if [ ! -e $VCF_FILE ]; then 
  echo "ERROR - Unable to find VCF file. Please check and try again - ${VCF_FILE}";
  exit 1;
fi

count=`ls ${VCF}/vcf_chr/*.vcf.gz | wc -l`               # total number of VCF chr files available - should be 41!
if [ $count != 41 ]; then
  echo "ERROR - Unable to find chromosome specific VCF files. Please check and try again - ${VCF}";
  exit 1;
fi

dos2unix ${DISEASE_STATUS}

sbatch <<EOT
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

perl /rds/project/rds-Qr3fy2NTCy0/Software/Git/Scripts/perl/vcf2excel.pl --default --vcf ${VCF}/vcf_chr/ens_WGS_219-chr\${SLURM_ARRAY_TASK_ID}.vcf.gz --status_file ${DISEASE_STATUS}
EOT

