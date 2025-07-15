#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J Clu_bcftools_call
#SBATCH --cpus-per-task=32
#SBATCH --mail-user email@emailexample.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#Script to call all locations

#Load the required modules
source /usr/modules/init/bash
module load bcftools
module load samtools
module load htslib

#Define variables
REF="/path/to/genome/GCA_00000#####.fna"  #Path to reference genome (.fasta or .fna)
BED="/path/to/positions.bed"  #path to .bed file with chr,start,end in each row
NEWDIR="/path/to/bcftools_runs"    #new directory you created to house outputs
OUTDIR="${NEWDIR}/$(date +%y%m%d)_bcftools_pref_plateID_separate_minBQ20_maxDP5000_fastp"  #new directory within NEWDIR with current date; pref = your genome prefix
BAMLIST="/path/to/bam_list.txt" #This is where you can subset out any samples that you don't want to run
BAMDIR="/path/to/bamdir"     #path to directory with filtered BAM files to use
THREADS=32

maxDP=5000 #max depth to consider (saves time if lower)
minBQ=20 #Phred based quality score for single basepair

#Create output directories
mkdir -p "$OUTDIR"

cd "$BAMDIR" || exit 1      #set working directory to where BAM files are located

# Generate index files for each BAM file if not already present
for BAM in *.bam; do
    SAMPLE=$(basename "$BAM" .bam)  # Extract sample name from file

    # Check if index exists, if not, generate index
    if [ ! -f "$BAM.bai" ]; then
        echo "Indexing $BAM"
        samtools index "$BAM"
    else
        echo "Index for $BAM already exists, skipping indexing."
    fi
done

#Loop through all BAM files in the current directory

#for BAM in *.bam; do
 #       SAMPLE=$(basename "$BAM" .bam) #extract sample name from file
while read -r BAM; do
        SAMPLE=$(basename "$BAM" .bam)

        #Generate mpileup
        bcftools mpileup \
                --fasta-ref "$REF" \
                --regions-file "$BED" \
                --min-BQ "$minBQ" \
                --max-depth "$maxDP" \
                -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
                -Ou \
                --threads "$THREADS"\
                "$BAMDIR/$BAM" | \
        bcftools call \
                --multiallelic-caller \
                -Ov \
                --threads "$THREADS" \
                -o "$OUTDIR/${SAMPLE}_positions.vcf"

       echo "Finished processing $SAMPLE"
done < "$BAMLIST"
