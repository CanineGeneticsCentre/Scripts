#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBACTH --mem=4G
#SBATCH --nodes=1
#SBATCH -o logs/job-%j.out

module load apps/picard
module load apps/samtools

SAMPLE=$1

export PICARD_JAVA_TMPDIR="${HOME}/data/javatempdir"
export PICARD_JAVA_MEM_MX='4g'

time picard ValidateSamFile INPUT=${SAMPLE}_final.bam OUTPUT=validateSamFile_final.out MODE=SUMMARY MAX_OUTPUT=100 MAX_OPEN_TEMP_FILES=7900 VALIDATION_STRINGENCY=SILENT
time samtools flagstat ${SAMPLE}_final.bam > flagstat_final.out

