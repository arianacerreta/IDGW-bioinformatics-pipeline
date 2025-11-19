#!/usr/bin/env bash
#SBATCH -p eight
#SBATCH -C "headnode"
#SBATCH -J zip_fastq
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=email@emailexample.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Load modules if needed
source /usr/modules/init/bash

# Variables
FASTQ="/path/individual_fastqs"
JOBS=$SLURM_CPUS_PER_TASK   # run 16 gzip jobs at once

cd "$FASTQ" || exit 1

echo "Compressing FASTQ files using $JOBS parallel gzip jobs..."

# Compress any remaining .fastq files in parallel
ls *.fastq | parallel -j "$JOBS" --will-cite '
    echo "Compressing {}..."
    gzip {}
'

# Final check
echo "Remaining uncompressed files:"
ls *.fastq 2>/dev/null || echo "All done!"
