#!/bin/bash

## Command to Run :  bash ./scripts/slurm_filterGVCF.sh <SAMPLE [BoT_32595]>

export CF3="${HOME}/data/canfam3/"
export RESULTS="${HOME}/data/samples"
#export SCRIPT_HOME=${HOME}/git/
SAMPLE=$1

jid1=$(sbatch -J ${SAMPLE}.gvcf2vcf ${HOME}/scripts/filterGVCF/gvcf2vcf.sh ${SAMPLE})
jid2=$(sbatch -J ${SAMPLE}.filterVCF --dependency=afterok:${jid1##* } ${HOME}/scripts/filterGVCF/filterVCF.sh ${SAMPLE})
jid3=$(sbatch -J ${SAMPLE}.finalStep --dependency=afterok:${jid2##* } ${HOME}/scripts/filterGVCF/finalStep.sh ${SAMPLE})
echo $jid3;