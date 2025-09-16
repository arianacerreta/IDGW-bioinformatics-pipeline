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
    field_count=\$(echo \"\$SAMPLE\" | awk -F'_' '{print NF}')

    case \$field_count in
        4)
            # i043_G12_GWAdapt9_PCRnegA9.bam
            index=\$(echo \$SAMPLE | cut -d'_' -f1)
            well=\$(echo \$SAMPLE | cut -d'_' -f2)
            lib=\$(echo \$SAMPLE | cut -d'_' -f3)
            sm=\$(echo \$SAMPLE | cut -d'_' -f4)
            rg=\"@RG\\tID:\${index}_\${well}_\${sm}\\tSM:\${sm}\\tLB:\${lib}\\tPL:ILLUMINA\\tPU:\${index}_\${well}\"
            ;;
        6)
            # i002_A01_GWAdapt11_PCR_Neg_AGW11.bam
            index=\$(echo \$SAMPLE | cut -d'_' -f1)
            well=\$(echo \$SAMPLE | cut -d'_' -f2)
            lib=\$(echo \$SAMPLE | cut -d'_' -f3)
            pcr=\$(echo \$SAMPLE | cut -d'_' -f4)
            neg=\$(echo \$SAMPLE | cut -d'_' -f5)
            sm=\$(echo \$SAMPLE | cut -d'_' -f6)
            rg=\"@RG\\tID:\${index}_\${well}_\${sm}\\tSM:\${pcr}\${neg}\${sm}\\tLB:\${lib}\\tPL:ILLUMINA\\tPU:\${index}_\${well}\"
            ;;
        *)
            # 5+ fields, normal samples
            index=\$(echo \$SAMPLE | cut -d'_' -f1)
            well=\$(echo \$SAMPLE | cut -d'_' -f2)
            lib=\$(echo \$SAMPLE | cut -d'_' -f3)
            sm=\$(echo \$SAMPLE | cut -d'_' -f4-)
            rg=\"@RG\\tID:\${index}_\${well}_\${sm}\\tSM:\${sm}\\tLB:\${lib}\\tPL:ILLUMINA\\tPU:\${index}_\${well}\"
            ;;
    esac

    # Add RG and sort
    samtools addreplacerg -r \"\$rg\" -O BAM {} \
    | samtools sort -@ \$THREADS_PER_JOB -o \"\$OUTDIR/\${SAMPLE}_RG_sort.bam\" -

     # Index if missing
    if [ ! -f \"\$OUTDIR/\${SAMPLE}_RG_sort.bai\" ]; then
        samtools index \"\$OUTDIR/\${SAMPLE}_RG_sort.bam\"
    fi

    echo \"Finished processing \$SAMPLE with RG: \$rg\"
"
