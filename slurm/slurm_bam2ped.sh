#!/bin/bash

export CF3="${HOME}/data/canfam3/"
export RESULTS="${HOME}/data/samples"
#export SCRIPT_HOME=${HOME}/git/
SAMPLE=$1

uid=`date | md5sum | cut -c1-8`
RUN_NAME="bam2ped_${uid}"

jid1=$(sbatch -J ${SAMPLE}.HaplotypeCaller /users/eschofield/scripts/bam2ped/step1.sh ${RUN_NAME} ${SAMPLE})
jid2=$(sbatch -J ${SAMPLE}.VariantFiltration --dependency=afterok:${jid1##* } /users/eschofield/scripts/bam2ped/step2.sh ${RUN_NAME} ${SAMPLE})
jid3=$(sbatch -J ${SAMPLE}.VariantAnnotator --dependency=afterok:${jid2##* } /users/eschofield/scripts/bam2ped/step3.sh ${RUN_NAME} ${SAMPLE})
jid4=$(sbatch -J ${SAMPLE}.VariantsToTable --dependency=afterok:${jid3##* } /users/eschofield/scripts/bam2ped/step4.sh ${RUN_NAME} ${SAMPLE})
jid5=$(sbatch -J ${SAMPLE}.TableToPed --dependency=afterok:${jid4##* } perl /users/eschofield/scripts/bam2ped/step5.pl ${RUN_NAME} ${SAMPLE})

# show dependencies in squeue output:
squeue -u $USER -o "%.8A  %.30j %.4C %.10m %.20E"