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

### if size of bam file >= 50M then delete corresponding sam file
bam_size=$(wc -c < ${SAMPLE}_dedup.bam)
if [ $bam_size -ge 50000000 ]; then
	rm -rf ${SAMPLE}_aligned_sorted_rg.bam
	rm -rf ${SAMPLE}_aligned_sorted_rg.bai
fi

time gatk -T RealignerTargetCreator -R ${DATA}/${GENOME}/ensembl/${GENOME}.fasta -I ${SAMPLE}_dedup.bam -nt 16 -o ${SAMPLE}_dedup.bam.intervals -mismatch 0.0 -S LENIENT --filter_mismatching_base_and_quals

