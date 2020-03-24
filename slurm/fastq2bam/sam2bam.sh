#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBACTH --mem=16G
#SBATCH --nodes=1
#SBATCH -o logs/job-%j.out

module load apps/picard

SAMPLE=$1
RUN_NAME=$2

export PICARD_JAVA_TMPDIR="${HOME}/data/javatempdir"
export PICARD_JAVA_MEM_MX='16g'

time picard AddOrReplaceReadGroups INPUT=${SAMPLE}_aligned.sam OUTPUT=${SAMPLE}_aligned_sorted_rg.bam rgID=${RUN_NAME} LB=${GENOME} PL='ILLUMINA' PU=${SAMPLE} SM=${SAMPLE} SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

