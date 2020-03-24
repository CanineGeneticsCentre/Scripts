#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
#SBACTH --mem=16G
#SBATCH -o outputs/job-%j.output
#SBATCH --nodes=1


module load apps/cramtools
module load apps/samtools

RUN_NAME=$1
SAMPLE=$2

export CRAMTOOLS_JAVA_TMPDIR="${HOME}/data/javatempdir"
export CRAMTOOLS_JAVA_MEM_MX='16g'

OUTPUT_PATH="${RESULTS}/${SAMPLE}"


#time cramtools cram -I ${OUTPUT_PATH}/${SAMPLE}\_final.bam -R ${CF3}/ensembl/canfam3.fasta -O  ${OUTPUT_PATH}/${SAMPLE}\_final.cram
rm -rf ${OUTPUT_PATH}/${SAMPLE}\_final.cram*
time samtools view -C ${OUTPUT_PATH}/${SAMPLE}\_final.bam -T ${CF3}/ensembl/canfam3.fasta -o ${OUTPUT_PATH}/${SAMPLE}.cram
time samtools index ${OUTPUT_PATH}/${SAMPLE}.cram

md5sum ${OUTPUT_PATH}/${SAMPLE}.cram
