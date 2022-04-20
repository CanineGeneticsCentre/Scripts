#!/bin/bash

## Command to Run :  bash slurm_vcf2excel.sh <Disease Status File> <CHR>


SCRIPTS=`dirname $0`

#VCF_FILE='/rds/project/rds-Qr3fy2NTCy0/Data/ens_WGS_219_VEP.vcf.gz'
DISEASE_STATUS=$1
#CHR=$2
#SEG_SCORE=$2
#EFFECT_SCORE=1

#VCF=`dirname $VCF_FILE`

#VCF='/rds/project/rds-Qr3fy2NTCy0/Data/VCF/'
GENOME='canfam3';
REF='cf3';
FASTA="/rds/project/rds-Qr3fy2NTCy0/Genomes/CanFam3.1/current"
VCF_DIR="/rds/project/rds-Qr3fy2NTCy0/Data/VCF/${GENOME}";

#printf "Which genome do you want to use?\n"
#printf "\t1. CanFam3 [default]\n"
#printf "\t2. CanFam4\n"
# Assign input value into a variable
#read answer

#if [[ -v $answer && $answer == "2" ]]; then
#    GENOME='canfam4';
#    REF='cf4';
#fi


[[ -z "$DISEASE_STATUS" ]] && { echo "ERROR: No disease status file provided for this run"; exit 1; }
#[[ -z "$CHR" ]] && { echo "ERROR: No chromosome speficied for this run"; exit 1; }

#VCF_FILE="${VCF_DIR}/${REF}-chr${CHR}.vcf.gz";

#if [ ! -e $VCF_FILE ]; then 
#  echo "ERROR - Unable to find VCF file. Please check and try again - ${VCF_FILE}";
#  exit 1;
#fi

#count=`ls ${VCF}/vcf_chr/*.vcf.gz | wc -l`               # total number of VCF chr files available - should be 41!
#count=`ls ${VCF}/vcf_chr/ens_WGS*.vcf.gz | wc -l`               # total number of VCF chr files available - should be 41!
count=`ls ${VCF_DIR}/${REF}-chr*.ann.vcf.gz | wc -l`               # total number of VCF chr files available - should be 41!
if [ $count != 41 ]; then
  echo "ERROR - Unable to find chromosome specific VCF files. Please check and try again - ${VCF}";
  exit 1;
fi

dos2unix ${DISEASE_STATUS}
echo sbatch --export=DISEASE_STATUS=${DISEASE_STATUS},FASTA=${FASTA},VCF_DIR=${VCF_DIR},SCRIPTS=${SCRIPTS},REF=${REF} ${SCRIPTS}/../slurm/vcf2excel.sh;
