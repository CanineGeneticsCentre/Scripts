#!/bin/bash

## Command to Run :  bash ./scripts/slurm_gvcfmerge.sh batch08

export CF3="${HOME}/data/canfam3/"
export RESULTS="${HOME}/data/samples"
#export SCRIPT_HOME=${HOME}/git/
BATCH=$1

uid=`date | md5sum | cut -c1-8`
RUN_NAME="gvcfmerge_${uid}"

if [ -e ${RESULTS}/${BATCH}.list ]
then
	jid1=$(sbatch -J g.vcf.merge /users/eschofield/scripts/gvcfmerge/step1.sh ${RUN_NAME} ${BATCH})
else
	echo "No list of files found named ${RESULTS}/${BATCH}.list"
fi
