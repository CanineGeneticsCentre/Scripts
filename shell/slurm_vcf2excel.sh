#!/bin/bash

## Command to Run :  bash slurm_vcf2excel.sh <Disease Status File> [<CHR>]


SCRIPTS=`dirname $0`

DISEASE_STATUS=$1
CHR=$2
[[ -z "$DISEASE_STATUS" ]] && { echo "ERROR: No disease status file provided for this run"; exit 1; }

GENOME='canfam4';
REF='cf4';

printf "Which genome do you want to use?\n"
printf "\t1. CanFam4 [default]\n"
printf "\t2. CanFam3\n"
# Assign input value into a variable
read answer

if [[ -n $answer && $answer == "2" ]]; then
    GENOME='canfam3';
    REF='cf3';
fi

FASTA="/rds/project/rds-9sJA7YGzZRc/Genomes/${GENOME}/current"
VCF_DIR="/rds/project/rds-9sJA7YGzZRc/VCF/${GENOME}";

DIR=`echo $RANDOM | md5sum | head -c 10`
mkdir -p $DIR/logs; cd $DIR
cp ../$DISEASE_STATUS .

count=`ls ${VCF_DIR}/${REF}-chr*.ann.vcf.gz | wc -l`               # total number of VCF chr files available - should be 41 for CF3, 42 for CF4!
if [ $count -lt 39 ]; then
  echo "ERROR - Unable to find chromosome specific VCF files. Please check and try again - ${VCF_DIR}";
  exit 1;
fi

if [ "$CHR" ]; then
  if [[ ${#CHR} -lt 4 ]] ; then
    CHR="chr"${CHR}
  fi
  LINE=$(grep -n "^${CHR}$" ${FASTA}/intervals/gdb-chromsomes.intervals | cut -f 1 -d":")
  ARRAY="${LINE}-${LINE}"
else
  ARRAY="1-${count}"
fi

dos2unix ${DISEASE_STATUS}
jid1=$(sbatch --array=${ARRAY} --export=DISEASE_STATUS=${DISEASE_STATUS},FASTA=${FASTA},VCF_DIR=${VCF_DIR},SCRIPTS=${SCRIPTS},REF=${REF} ${SCRIPTS}/../slurm/vcf2excel.sh);
echo $jid1;

sbatch --export=DISEASE_STATUS=${DISEASE_STATUS} --dependency=afterok:${jid1##* } ${SCRIPTS}/../slurm/vcf2excel-finish.sh
