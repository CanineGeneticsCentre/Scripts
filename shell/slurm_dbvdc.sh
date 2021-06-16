#!/bin/bash

## Command to Run :  bash slurm_dbvdc.sh <SNPS FILE>

DBVDC='/rds/project/rds-Qr3fy2NTCy0/Data/Downloads/dbvdc.648.vars.ann.vcf.gz';
SNPS=$1

if [ ! -e $SNPS ]
then 
  echo "ERROR - Unable to find file of SNP positions to test - ${SNPS}";
  exit 1;
fi

sbatch <<EOT
#!/bin/bash

#SBATCH -A MELLERSH-SL3-CPU
##SBATCH -A KCGC-SL2-CPU
#SBATCH -J dbvdc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:00:30
#SBATCH --mail-type=END,FAIL
#SBATCH -p skylake

module purge                                      # Removes all modules still loaded
module load rhel7/default-peta4                   # REQUIRED - loads the basic environment

module load bcftools-1.9-gcc-5.4.0-b2hdt5n        # bcftools

bcftools filter -R ${SNPS} ${DBVDC} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > snps.out
bcftools query -l ${DBVDC} > samples.list

perl -ae 'print "CHR\tPOS\tREF\tALT"; while(<>){ chomp $_; print "\t".join("\t", $_."_A", $_."_B");} print "\n";' samples.list > dbvdc.648.out

perl /rds/project/rds-Qr3fy2NTCy0/Software/Git/Scripts/perl/check_dbvdc.pl ${SNPS}

rm -rf snps.out samples.list

EOT