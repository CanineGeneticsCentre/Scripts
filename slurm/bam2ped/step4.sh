#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
##SBATCH -e outputs/job-%j.error
#SBATCH -o outputs/job-%j.output

module load apps/gatk

RUN_NAME=$1
SAMPLE=$2

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"
export GATK_JAVA_MEM_MX='32g'

OUTPUT_PATH="${RESULTS}/${SAMPLE}/${RUN_NAME}"
mkdir -p "$OUTPUT_PATH"

gatk -T VariantsToTable -R ${CF3}/ensembl/canfam3.fasta -F CHROM -F POS -F hd_chip.ID -F REF -F ALT -F FILTER -GF GT --showFiltered -V ${OUTPUT_PATH}/${SAMPLE}.g.vcf -o ${OUTPUT_PATH}/${SAMPLE}.table
