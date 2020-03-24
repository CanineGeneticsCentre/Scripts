#!/bin/bash

export DATA="${HOME}/data"
export GENOME="canfam4"
export SAMPLES="${DATA}/static/wgs_samples/"

SAMPLE=$1

uid=`date | md5sum | cut -c1-8`
RUN_NAME="${SAMPLE}_fastq2bam_${uid}"

mkdir -p "$RUN_NAME/logs"
cd $RUN_NAME

jid1=$(sbatch -J ${SAMPLE}.fastq2sam ${HOME}/scripts/fastq2bam/fastq2sam.sh ${SAMPLE})
jid2=$(sbatch -J ${SAMPLE}.sam2bam --dependency=afterok:${jid1##* } ${HOME}/scripts/fastq2bam/sam2bam.sh ${SAMPLE} ${RUN_NAME})
jid3=$(sbatch -J ${SAMPLE}.validateSam --dependency=afterok:${jid2##* } ${HOME}/scripts/fastq2bam/validateSam.sh ${SAMPLE})
jid4=$(sbatch -J ${SAMPLE}.markDuplicates --dependency=afterok:${jid3##* } ${HOME}/scripts/fastq2bam/markDuplicates.sh ${SAMPLE})
jid5=$(sbatch -J ${SAMPLE}.realignReads --dependency=afterok:${jid4##* } ${HOME}/scripts/fastq2bam/realignReads.sh ${SAMPLE})
jid6=$(sbatch -J ${SAMPLE}.realignIndels --dependency=afterok:${jid5##* } ${HOME}/scripts/fastq2bam/realignIndels.sh ${SAMPLE})
jid7=$(sbatch -J ${SAMPLE}.baseRecalibrator --dependency=afterok:${jid6##* } ${HOME}/scripts/fastq2bam/baseRecalibrator.sh ${SAMPLE})
jid8=$(sbatch -J ${SAMPLE}.printReads --dependency=afterok:${jid7##* } ${HOME}/scripts/fastq2bam/printReads.sh ${SAMPLE})
jid9=$(sbatch -J ${SAMPLE}.CollectMetrics --dependency=afterok:${jid8##* } ${HOME}/scripts/fastq2bam/collectMetrics.sh ${SAMPLE})
jid10=$(sbatch -J ${SAMPLE}.validateSam2 --dependency=afterok:${jid8##* } ${HOME}/scripts/fastq2bam/validateSam2.sh ${SAMPLE})
jid11=$(sbatch -J ${SAMPLE}.FinalStep --dependency=afterok:${jid9##* }:${jid10##* } ${HOME}/scripts/fastq2bam/finalStep.sh ${SAMPLE})

echo $jid11
echo logs/job-${jid11##* }.output
