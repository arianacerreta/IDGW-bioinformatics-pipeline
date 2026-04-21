#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J summ_stats_AdaptOptim
#SBATCH --cpus-per-task=16
#SBATCH --mail-user email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

#Script to calculate summary statistics for panel performance

#Load the required modules
source /usr/modules/init/bash
module load bcftools
module load samtools

#Define variables
MERGED="/path/unzipped.vcf"
BAMDIR1="/path/aligned/CanFam3_1_fastp_trimmed"
BAMDIR2="/path/aligned/CanFam3_1_fastp_trimmed"
TARGETS_BED="/path/info/adaptive_positions_canFam3_1.bed" #BED file with aplicon target regions
THREADS=16 #number of threads to use; should equal cpus-per-task in header

#Create output directory
OUTDIR="/path/summ_stats/$(date +%y%m%d)_CanFam3_1_details" #where you want your output directory created
mkdir -p "$OUTDIR"

#1 calculate and output depth per locus
echo "Calculating depth per locus..."
bcftools query -f '%CHROM\t%POS[\t%DP]\n' "$MERGED" > "$OUTDIR/depth_per_locus.tsv"

#2 calculate genotype completness per sample
echo "Calulating genotype completeness per sample..."
bcftools query -f '[%SAMPLE\t%GT\n]' "$MERGED" | \
awk '$2 != "./." { count[$1]++ } END { for (s in count) print s, count[s] }' \
> "$OUTDIR/completeness_per_sample.tsv"

# 3. Mapping quality per sample
echo "Calculating mean mapping quality per sample..."

# Output file for mean MAPQ
echo -e "Sample\tMean_MAPQ_AllReads\tMean_MAPQ_Q30" > "$OUTDIR/mean_mapq.tsv"

# Output file for histogram counts
hist_counts="$OUTDIR/mapq_histogram_counts.tsv"
> "$hist_counts"  # clear file

# Loop through BAMs in parallel
find "$BAMDIR1" "$BAMDIR2" -name "*.bam" | parallel -j "$THREADS" '
    SAMPLE=$(basename {} .bam)

    tmp_mean_all="tmp_mean_all_${SAMPLE}.txt"
    tmp_mean_q30="tmp_mean_q30_${SAMPLE}.txt"
    tmp_hist="tmp_hist_${SAMPLE}.tsv"

    # 1️⃣ean MAPQ across all aligned reads
    samtools view -F 0x904 {} | awk '\''{sum+=$5; n++} END {if(n>0) printf "%.2f\n", sum/n; else print "NA"}'\'' > $tmp_mean_all

    # 2️⃣ean MAPQ for reads with MAPQ >=30
    samtools view -q 30 -F 0x904 {} | awk '\''{sum+=$5; n++} END {if(n>0) printf "%.2f\n", sum/n; else print "NA"}'\'' > $tmp_mean_q30

    # 3️⃣istogram counts for all reads
    samtools view -F 0x904 {} | awk '\''{count[$5]++} END {for(val in count) print val, count[val]}'\'' > $tmp_hist

    # Append mean values to mean_mapq.tsv
    echo -e "${SAMPLE}\t$(cat $tmp_mean_all)\t$(cat $tmp_mean_q30)" >> '"$OUTDIR"'/mean_mapq.tsv

    # Append histogram counts to mapq_histogram_counts.tsv
    awk -v s="$SAMPLE" '\''{print s"\t"$1"\t"$2}'\'' $tmp_hist >> '"$hist_counts"'

    # Clean up temp files
    rm -f $tmp_mean_all $tmp_mean_q30 $tmp_hist
'
# Percent reads mapping to target regions
echo "Calculating % reads on target..."
echo -e "Sample\tTotal_Reads\tReads_on_Target\tPercent_on_Target" > "$OUTDIR/reads_on_target.tsv"

echo "DEBUG: listing BAMs found:"
find "$BAMDIR1" "$BAMDIR2" -name "*.bam" | head

find "$BAMDIR1" "$BAMDIR2" -name "*.bam" | parallel -j "$THREADS" '
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
