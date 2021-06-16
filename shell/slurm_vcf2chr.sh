#!/bin/bash

## Command to Run :  bash slurm_vcf2chr.sh <VCF_FILE>

#VCF_FILE='/rds/project/rds-Qr3fy2NTCy0/Data/ens_WGS_219_VEP.vcf.gz'
VCF_FILE=$1

[[ -z "$VCF_FILE" ]] && { echo "ERROR: No VCF file provided for this run"; exit 1; }
if [ ! -e $VCF_FILE ]
then 
  echo "ERROR - Unable to find VCF file. Please check and try again - ${VCF_FILE}";
  exit 1;
fi

VCF=`dirname $VCF_FILE`
mkdir vcf_chr; cd $VCF/vcf_chr

for CHR in {1..38}; do printf "chr$CHR\t1\t150000000\n$CHR\t1\t150000000" > $CHR.chr; printf "chr$CHR $CHR\n" >> chrs.list; done
printf "chrX\t1\t150000000\nX\t1\t150000000" > 39.chr; printf "chrX X\n" >> chrs.list;
printf "chrMT\t1\t150000000\nMT\t1\t150000000" > 40.chr; printf "chrMT MT\n" >> chrs.list;
tabix -l ${VCF_FILE} | tail -n +41 | awk -v OFS="\t" '{print $1, "1", "150000000"}' > 41.chr

sbatch <<EOT
#!/bin/bash

##SBATCH -A MELLERSH-SL3-CPU
#SBATCH -A KCGC-SL2-CPU
#SBATCH -J vcf2chr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --mail-type=ALL
#SBATCH -p skylake
#SBATCH --array 1-41

#SBATCH -o $HOME/rds/hpc-work/logs/job-%A_%a.out

module purge                                      # Removes all modules still loaded
module load rhel7/default-peta4                   # REQUIRED - loads the basic environment

module load bcftools-1.9-gcc-5.4.0-b2hdt5n        # bcftools
module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7    # bgzip/tabix

#VCF+="/vcf-chr/ens_WGS_219-chr"
#VCF+=\${SLURM_ARRAY_TASK_ID}
#VCF+=".vcf.gz"
#echo $VCF;

echo ens_WGS_219-chr\${SLURM_ARRAY_TASK_ID}.vcf.gz

bcftools view $VCF_FILE -R \${SLURM_ARRAY_TASK_ID}.chr | bcftools annotate --rename-chrs chrs.list | bcftools annotate --set-id +'%CHROM:%POS' | bgzip -c > ens_WGS_219-chr\${SLURM_ARRAY_TASK_ID}.vcf.gz
tabix -p vcf ens_WGS_219-chr\${SLURM_ARRAY_TASK_ID}.vcf.gz
rm -rf \${SLURM_ARRAY_TASK_ID}.chr


EOT

