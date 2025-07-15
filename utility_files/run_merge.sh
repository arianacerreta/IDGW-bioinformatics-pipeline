#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J Clu_bcftools_merge
#SBATCH --cpus-per-task=8
#SBATCH --mail-user email@emailexample.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#Script to normalize and then merge vcf files

#Load the required modules
source /usr/modules/init/bash
module load bcftools
module load htslib

#Define variables
REF="/path/to/genome/GCA_00000####.fna"  #Path to reference genome (.fasta or .fna)
DIR="/path/to/individual_vcfs" #path to all the bcfs you would like to merge
NORMDIR="${DIR}/normalized" #new directory
NEWDIR="${NORMDIR}/merged"

#make NORMDIR
mkdir -p "$NORMDIR"

#normalize VCFs
for VCF in "$DIR"/*.vcf; do
        SAMPLE=$(basename "$VCF" .vcf) #extract sample name from file

        bcftools norm -m +any -f "$REF" -o "$NORMDIR/${SAMPLE}_norm.vcf" "$VCF"

        #Compress VCF
        bgzip -f "$NORMDIR/${SAMPLE}_norm.vcf"

        #Index compressed VCF
        tabix -p vcf "$NORMDIR/${SAMPLE}_norm.vcf.gz"

        echo "Finished processing $SAMPLE"
done

#make NEWDIR and list in NEWDIR
mkdir -p "$NEWDIR"

# List VCF files and check if files exist
find "$NORMDIR" -name "*_norm.vcf.gz" > "$NEWDIR/vcf_list.txt"
if [ ! -s "$NEWDIR/vcf_list.txt" ]; then
    echo "No VCF files found. Exiting."
    exit 1
fi

#merge VCFs and save in unzipped VCF

bcftools merge -o "$NEWDIR/merged.vcf.gz" -Oz -l "$NEWDIR/vcf_list.txt"

# Index merged VCF
tabix -p vcf "$NEWDIR/merged.vcf.gz"

echo "VCF merge complete!"
