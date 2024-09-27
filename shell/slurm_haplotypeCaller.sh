#!/bin/bash

SAMPLE=$1
REF=$2
CHR=$3
POS=$4
FLANK=500

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; echo; echo "Usage: check_hc <SAMPLE> <REF> <CHR> <POS>"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; echo "Usage: check_hc <SAMPLE> <REF> <CHR> <POS>"; exit 1; }
[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run (e.g. chr23)"; echo "Usage: check_hc <SAMPLE> <REF> <CHR> <POS>"; exit 1; }
[[ -z "$POS" ]] && { echo "ERROR: No POSITION provided for this run (e.g. 81417279)"; echo "Usage: check_hc <SAMPLE> <REF> <CHR> <POS>"; exit 1; }

if [ "$CHR" ]; then
  if [[ ${#CHR} -lt 4 ]] ; then
    CHR="chr"${CHR}
  fi
fi

SCRIPTS=`dirname $0`
CFG="/rds/project/rds-Qr3fy2NTCy0/Genomes/${REF}.config"
source $CFG;

WGS_DIR="/rds/project/rds-9sJA7YGzZRc/Samples/${SAMPLE}";
BAM_FILE="${WGS_DIR}/${SAMPLE}-${REF}.bam";

if [ ! -e $BAM_FILE ]; then 
  echo "ERROR - Unable to find BAM file. Please check and try again - ${BAM_FILE}";
  exit 1;
fi

OUTPUT="${SAMPLE}_${CHR}-${POS}";

START="$((POS-FLANK))"
STOP="$((POS+FLANK))"
LOCATION="${CHR}:${START}-${STOP}"

sbatch <<EOT
#!/bin/bash

#SBATCH -A MELLERSH-SL3-CPU
#SBATCH -J HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 00:05:00
#SBATCH --mail-type=ALL
#SBATCH -p cclake

#SBATCH -o $HOME/rds/hpc-work/logs/job-%j.out

module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

module load $GATK
module load $TABIX

gatk HaplotypeCaller -R ${FASTA}/${GENOME}.fasta -I ${BAM_FILE} -bamout ${OUTPUT}.bam -O ${OUTPUT}.vcf -L ${LOCATION}

#bgzip ${OUTPUT}.vcf
#tabix -p vcf ${OUTPUT}.vcf
EOT
