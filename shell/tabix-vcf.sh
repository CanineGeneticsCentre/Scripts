#!/usr/bin/bash

module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7
module load bcftools-1.9-gcc-5.4.0-b2hdt5n

FILE=$1
POS=$2

[[ -z "$FILE" ]] && { echo "ERROR: No VCF FILE provided"; exit 1; }
[[ -z "$POS" ]] && { echo "ERROR: No POSITION provided"; exit 1; }

(bcftools view -h ${FILE} && tabix -p vcf ${FILE} ${POS}) | bgzip
