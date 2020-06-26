#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBACTH --mem=32G
#SBATCH --ntasks=16
#SBATCH -o logs/job-%j.out

module load apps/gatk

SAMPLE=$1

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"
export GATK_JAVA_MEM_MX='32g'

time gatk -T BaseRecalibrator -nct 16 -R ${DATA}/${GENOME}/ensembl/${GENOME}.fasta -I ${SAMPLE}_clean.dedup.bam -knownSites ${DATA}/${GENOME}/${GENOME}_snps_all.vcf -o ${SAMPLE}.recal.grp -S LENIENT

