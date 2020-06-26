#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBACTH --mem=16G
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH -o logs/job-%j.out

module load apps/picard

SAMPLE=$1

export PICARD_JAVA_TMPDIR="${HOME}/data/javatempdir"
export PICARD_JAVA_MEM_MX='16g'

### if size of bam file >= 50M then delete corresponding sam file
bam_size=$(wc -c < ${SAMPLE}_aligned_sorted_rg.bam)
if [ $bam_size -ge 50000000 ]; then
	rm -rf ${SAMPLE}_aligned.sam
fi

time picard MarkDuplicates INPUT=${SAMPLE}_aligned_sorted_rg.bam OUTPUT=${SAMPLE}_dedup.bam VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true M=MarkDuplicates_metrics.out

