#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J Clu_bcftools_filter
#SBATCH --mail-user email@uemailexample.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#Script to filter VCF file

#Load the required modules
source /usr/modules/init/bash
module load bcftools
module load htslib

#Define variables
DP=20 #filters for a minimum depth of 20
DIR="/path/to/normalized/merged" #path to location of normalized and merged vcf
REMOVE="/path/to/remove_targets.txt" #text doc of positions to remove in chr, bp with no header
NEWDIR="${DIR}/$(date +%y%m%d)_filtered_DP${DP}"

mkdir -p ${NEWDIR}

echo "Filtering VCF read depth (DP< ${DP}), trimming ALT alleles, and removing specified sites..."

bcftools filter \
        -S . \
        -e "FORMAT/DP < ${DP}" \
        -Ou ${DIR}/merged.vcf.gz | \
bcftools view \
        --trim-alt-alleles \
        -T ^${REMOVE} \
        -Oz -o ${NEWDIR}/final_filtered.vcf.gz

#Index the filtered VCF
bcftools index ${NEWDIR}/final_filtered.vcf.gz

#unzip VCF
bgzip -d -c ${NEWDIR}/final_filtered.vcf.gz > ${NEWDIR}/final_filtered_unzipped.vcf

echo "Filtering complete. Output written to ${NEWDIR}"
