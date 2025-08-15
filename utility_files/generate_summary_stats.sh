#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J summ_stats
#SBATCH --cpus-per-task=16
#SBATCH --mail-user email@emailexample.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#Script to calculate summary statistics for panel performance

#Load the required modules
source /usr/modules/init/bash
module load bcftools
module load samtools

#Define variables
MERGED="/path/to/merged.vcf.gz"
BAMDIR="/path/to/bam_dir"
TARGETS_BED="path/to/bedfile.bed" #BED file with aplicon target regions
THREADS=16 #number of threads to use; should equal cpus-per-task in header

#Create output directory
OUTDIR="path/to/output/directory/$(date +%y%m%d)_run_stats" #where you want your output directory created
mkdir -p "$OUTDIR"

#1 calculate and output depth per locus
echo "Calculating depth per locus..."
bcftools query -f '%CHROM\t%POS[\t%DP]\n' "$MERGED" > "$OUTDIR/depth_per_locus.tsv"

#2 calculate genotype completness per sample
echo "Calulating genotype completeness per sample..."
bcftools query -f '[%SAMPLE\t%GT\n]' "$MERGED" | \
awk '$2 != "./." { count[$1]++ } END { for (s in count) print s, count[s] }' \
> "$OUTDIR/completeness_per_sample.tsv"

# 3 mapping quality per sample
echo "Calculating mean mapping quality per sample..."

cd "$BAMDIR" || exit 1

echo -e "Sample\tMean_MAPQ_AllReads\tMean_MAPQ_Q30" > "$OUTDIR/mean_mapq.tsv"

ls *.bam | parallel -j "$THREADS" '
    SAMPLE=$(basename {} .bam)

   #mean MAPQ across all reads (no filtering)
   samtools view -F 0x904 {} | \
   awk "{sum+=\$5; n++} END {if (n>0) printf \"%0.2f\", sum/n; else printf \"NA\"}" > tmp_all.txt

    #mean MAPQ for reads with MAPQ >=30
    samtools view -q 30 -F 0x904 {} | \
    awk "{sum+=$5; n++} END {if (n>0) print \"%0.2f\", sum/n; else print printf \"NA\"}" > tmp_q30.txt

    echo -e "${SAMPLE}\t$(cat tmp_all.txt)\t$(cat tmp_q30.txt)"
 ' >> "$OUTDIR/mean_mapq.tsv"

rm -f tmp_all.txt tmp_q30.txt

# Percent reads mapping to target regions
echo "Calculating % reads on target..."
echo -e "Sample\tTotal_Reads\tReads_on_Target\tPercent_on_Target" > "$OUTDIR/reads_on_target.tsv"

ls *.bam | parallel -j "$THREADS" '
    SAMPLE=$(basename {} .bam)
    TOTAL=$(samtools view -c {})
    ON_TARGET=$(samtools view -c -L '"$TARGETS_BED"' {})
    if [ "$TOTAL" -gt 0 ]; then
      PCT=$(echo "scale=4; $ON_TARGET/$TOTAL*100" | bc)
    else
      PCT="NA"
    fi  
    echo -e "${SAMPLE}\t${TOTAL}\t${ON_TARGET}\t${PCT}"
' >> "$OUTDIR/reads_on_target.tsv"

echo "All summary statistics written to: $OUTDIR"
