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
OUTDIR="/path/to/output/directory"
BAMDIR="/path/to/aligned/bams"

TOTAL_THREADS=$SLURM_CPUS_PER_TASK
JOBS=4
THREADS_PER_JOB=$(( TOTAL_THREADS / JOBS ))
export THREADS_PER_JOB OUTDIR

cd "$BAMDIR" || exit 1

ls *.bam | parallel --env THREADS_PER_JOB,OUTDIR -j "$JOBS" "
    SAMPLE=\$(basename {} .bam)
    NF=\$(echo \$SAMPLE | awk -F'_' '{print NF}')

    case \$NF in
        5)
            #i116_96_IDFGL0353_T9_NTCTray9.bam
            index=\$(echo \$SAMPLE | cut -d'_' -f1)
            well=\$(echo \$SAMPLE | cut -d'_' -f2)
            lib=\$(echo \$SAMPLE | cut -d'_' -f3)
            sm=\$(echo \$SAMPLE | cut -d'_' -f5)
            rg=\"@RG\\tID:\${index}_\${well}_\${sm}\\tSM:\${sm}\\tLB:\${lib}\\tPL:ILLUMINA\\tPU:\${index}_\${well}\"
            ;;
        6)
            #i328_2_IDFGL0419_T1_UI0742_C073772.bam
            index=\$(echo \$SAMPLE | cut -d'_' -f1)
            well=\$(echo \$SAMPLE | cut -d'_' -f2)
            lib=\$(echo \$SAMPLE | cut -d'_' -f3)
            sm=\$(echo \$SAMPLE | awk -F'_' '{print \$(NF-1)"_"\$NF}')
            rg=\"@RG\\tID:\${index}_\${well}_\${sm}\\tSM:\${sm}\\tLB:\${lib}\\tPL:ILLUMINA\\tPU:\${index}_\${well}\"
            ;;
        *)
            # all other fields, normal samples
            #i116_47_IDFGL0353_T9_CluVARI21C_0006_CHES.bam
            index=\$(echo \$SAMPLE | cut -d'_' -f1)
            well=\$(echo \$SAMPLE | cut -d'_' -f2)
            lib=\$(echo \$SAMPLE | cut -d'_' -f3)
            sm=\$(echo \$SAMPLE | awk -F'_' '{print \$(NF-2)"_"\$(NF-1)"_"\$NF}')
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
