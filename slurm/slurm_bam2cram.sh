#!/bin/bash

export CF3="${HOME}/data/canfam3/"
export RESULTS="${HOME}/data/samples"
#export SCRIPT_HOME=${HOME}/git/
SAMPLE=$1

uid=`date | md5sum | cut -c1-8`
RUN_NAME="bam2cram_${uid}"

jid1=$(sbatch -J ${SAMPLE}.bam2cram ${HOME}/scripts/bam2cram/bam2cram.sh ${RUN_NAME} ${SAMPLE})
echo $jid1;

