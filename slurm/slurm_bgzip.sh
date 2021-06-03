#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --mail-type=ALL
#SBATCH -p skylake

#SBATCH -o /home/es904/rds/hpc-work/logs/job-%j.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7    # bgzip/tabix

FILE=$1

bgzip $FILE
