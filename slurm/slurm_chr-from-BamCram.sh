#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
#SBACTH --mem=16G
#SBATCH --nodes=1
#SBATCH --job-name="chr.cram"
#SBATCH -o logs/job-%j.out
#SBATCH -e logs/job-%j.err

export DATA="${HOME}/data"
export GENOME="canfam4"

module load apps/samtools

SAMPLE=$1
CHR=$2

if [ -e ${SAMPLE}.cram ]; then
	samtools view -h -C -T ${DATA}/${GENOME}/ensembl/${GENOME}.fasta ${SAMPLE}.cram ${CHR} > ${SAMPLE}-chr${CHR}.cram
	samtools index ${SAMPLE}-chr${CHR}.cram
elif [ -e ${SAMPLE}.bam ]; then
	samtools view -h -b -T ${DATA}/${GENOME}/ensembl/${GENOME}.fasta ${SAMPLE}.bam ${CHR} > ${SAMPLE}-chr${CHR}.bam
	samtools index ${SAMPLE}-chr${CHR}.bam
else
	echo "ERROR - Neither BAM nor CRAM files found here for ${SAMPLE}";
fi