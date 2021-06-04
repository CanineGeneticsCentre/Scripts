#!/bin/bash

## Command to Run :  bash slurm_vcf2chr.sh <VCF_FILE>

VCF_FILE='/rds/project/rds-Qr3fy2NTCy0/Data/ens_WGS_219_VEP.vcf.gz'

if [ ! -e $VCF_FILE ]
then 
  echo "ERROR - Unable to find VCF file. Please check and try again - ${VCF_FILE}";
  exit 1;
fi


for CHR in {1..38}; do printf "chr$CHR\n$CHR" > $CHR.chr; printf "chr$CHR $CHR\n" >> chrs.list; done
printf "chrX\nX" > 39.chr; printf "chrX X\n" >> chrs.list;
printf "chrMT\nMT" > 40.chr; printf "chrMT MT\n" >> chrs.list;
tabix -l ${VCF_FILE} | tail -n +41 > 41.chr

sbatch <<EOT
#!/bin/bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -J vcf2chr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:01:00
#SBATCH --mail-type=ALL
#SBATCH -p skylake
#SBATCH --array 1-41

#SBATCH -o $HOME/rds/hpc-work/logs/job-%A_%a.out

module load bcftools-1.9-gcc-5.4.0-b2hdt5n        # bcftools
module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7    # bgzip/tabix

VCF=`dirname $VCF_FILE`
VCF+="/vcf-chr/ens_WGS_219-chr\${SLURM_ARRAY_TASK_ID}.vcf.gz"

if [ ! -e $VCF ]
then
  echo "Need to create $VCF";
  # echo "chr$CHR $CHR" > chr$CHR.chrs
  # bcftools view $VCF_ALL --regions $CHR,chr$CHR | bcftools annotate --rename-chrs chrs.list | bcftools annotate --set-id +'%CHROM:%POS' | bgzip -c > $VCF
  # tabix -p vcf $VCF
  # rm -rf chrs.list
fi


EOT

