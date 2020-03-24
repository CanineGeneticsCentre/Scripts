#!/bin/bash

## Command to Run :  bash ./scripts/slurm_gvcf2vcf.sh

export CF3="${HOME}/data/canfam3/"
export RESULTS="${HOME}/data"
#export SCRIPT_HOME=${HOME}/git/

uid=`date | md5sum | cut -c1-8`
RUN_NAME="gvcf2vcf_${uid}"

if [ -e ${RESULTS}/gVCF_batches/gvcf2vcf.list ]
then
	echo "submitting job..."
	jid1=$(sbatch -J gvcf2VCF /users/eschofield/scripts/gvcf2vcf/step1.sh)
else
	echo "No list of g.vcf files found named ${RESULTS}/gVCF_batches/gvcf2vcf.list"
fi
