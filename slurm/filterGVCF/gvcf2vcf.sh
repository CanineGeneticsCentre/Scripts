#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
##SBATCH -e outputs/job-%j.error
#SBATCH -o outputs/job-%j.output

module load apps/gatk

SAMPLE=$1

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"

if test -f "${RESULTS}/${SAMPLE}/${SAMPLE}.gVCF.gz"; then
    mv ${RESULTS}/${SAMPLE}/${SAMPLE}.gVCF.gz ${RESULTS}/${SAMPLE}/${SAMPLE}.g.vcf.gz
    mv ${RESULTS}/${SAMPLE}/${SAMPLE}.gVCF.gz.tbi ${RESULTS}/${SAMPLE}/${SAMPLE}.g.vcf.gz.tbi
fi

gatk -T GenotypeGVCFs -R ${CF3}/ensembl/canfam3.fasta -V ${RESULTS}/${SAMPLE}/${SAMPLE}.g.vcf.gz -o ${RESULTS}/${SAMPLE}/${SAMPLE}.vcf.gz