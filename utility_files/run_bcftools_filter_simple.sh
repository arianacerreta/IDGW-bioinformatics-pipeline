                                                                                            run_bcftools_filter.sh                                                                                                         #!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J Clu_bcftools_filter
#SBATCH --mail-user email@emailexample.com
#SBATCH --mail-type=BEGIN,END,FAIL

#Script to filter VCF file

#Load the required modules
source /usr/modules/init/bash
module load bcftools
module load htslib

#Define variables
DP=0
DIR="/path/to/normalized/merged" #this will be the folder from the previous step
NEWDIR="${DIR}/$(date +%y%m%d)_filtered_DP${DP}"

mkdir -p ${NEWDIR}

INPUT_VCF="${DIR}/merged.vcf.gz"
OUTPUT_VCF="${NEWDIR}/final_filtered.vcf.gz"

# Check if input file exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "ERROR: Input VCF not found: $INPUT_VCF"
    exit 1
fi

echo "Filtering VCF read depth (DP< ${DP}), trimming ALT alleles, and removing specified sites..."

bcftools filter -S . -e "FORMAT/DP < ${DP}" -Ou "$INPUT_VCF" | \
bcftools view --trim-alt-alleles -Oz -o "$OUTPUT_VCF"

#Index the filtered VCF
bcftools index "$OUTPUT_VCF"

#unzip VCF
bgzip -d -c "$OUTPUT_VCF" > "${NEWDIR}/final_filtered_unzipped.vcf"

echo "Filtering complete. Output written to ${NEWDIR}"
