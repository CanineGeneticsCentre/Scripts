#!/bin/bash

export CF3="${HOME}/data/canfam3/"
export RESULTS="${HOME}/data/samples"
SAMPLE=$1

uid=`date | md5sum | cut -c1-8`
RUN_NAME="bam2gvcf_${uid}"

jid1=$(sbatch -J ${SAMPLE}.bam2gvcf ${HOME}/scripts/bam2gvcf/bam2gvcf.sh ${RUN_NAME} ${SAMPLE})
echo $jid1;

