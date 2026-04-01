#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J RG_sort
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Load samtools
source /usr/modules/init/bash
module load samtools

# Paths
OUTDIR="/path/to/outdir/RG_aligned_sorted"
BAMDIR="/path/to/bams" #location of aligned .bams

TOTAL_THREADS=$SLURM_CPUS_PER_TASK
JOBS=4
THREADS_PER_JOB=$(( TOTAL_THREADS / JOBS ))
export THREADS_PER_JOB OUTDIR

cd "$BAMDIR" || exit 1

ls *.bam | parallel --env THREADS_PER_JOB,OUTDIR -j "$JOBS" "
    SAMPLE=\$(basename {} .bam)

    # Extract fields
    i5=\$(echo \$SAMPLE | cut -d'_' -f3)           # 384i5
    idx=\$(echo \$SAMPLE | cut -d'_' -f4)          # 003

    plate=\$(echo \$SAMPLE | cut -d'_' -f5 | cut -d'.' -f1)  # P053

    # Extract lab ID
    lab_raw=\$(echo \$SAMPLE | sed -E 's/.*WolfComp_([^~]+).*/\\1/')
    sm=\$(echo \$lab_raw | sed 's/^WL-//')

    # RG consistent with old dataset
    rg=\"@RG\\tID:\${i5}_\${idx}_\${sm}\\tSM:\${sm}\\tLB:\${plate}\\tPL:ILLUMINA\\tPU:\${i5}_\${idx}\"

    # Add RG and sort
    samtools addreplacerg -r \"\$rg\" -O BAM {} \
    | samtools sort -@ \$THREADS_PER_JOB -o \"\$OUTDIR/\${SAMPLE}_RG_sort.bam\" -

    # Index if missing
    if [ ! -f \"\$OUTDIR/\${SAMPLE}_RG_sort.bai\" ]; then
        samtools index \"\$OUTDIR/\${SAMPLE}_RG_sort.bam\"
    fi

    echo \"Finished processing \$SAMPLE with RG: \$rg\"
"
