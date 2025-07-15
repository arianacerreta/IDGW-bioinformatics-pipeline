#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J Clu_samtools_filter
#SBATCH --cpus-per-task=16
#SBATCH --mail-user email@emailexample.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#Script to filter .bam files based on MAPQ values

#Load the required modules
source /usr/modules/init/bash
module load samtools
module load parallel

#Define variables
NEWDIR="/path/to/samtools_filtered_BAMS"    #new directory you created to house outputs
MAPQ=20  #filters out everything <MAPQ
OUTDIR="${NEWDIR}/$(date +%y%m%d)_MAPQ${MAPQ}_pref_plateID_fastp" #new directory within NEWDIR with current date, pref=your genome prefix
BAMDIR="/path/to/bam_dir"
THREADS=16 #number of threads to use; should equal cpus-per-task in header

#Create output directory
mkdir -p "$OUTDIR"

#set working directory to where BAM files are located or exit if fails
cd "$BAMDIR" || exit 1

#Use parallel to filter BAM files
ls *.bam | parallel -j "$THREADS" '
    SAMPLE=$(basename {} .bam)
    samtools view -@ 1 -q '"$MAPQ"' -b -o '"$OUTDIR"'/${SAMPLE}_filtered.bam {} &&
    echo "Finished processing $SAMPLE"
'

