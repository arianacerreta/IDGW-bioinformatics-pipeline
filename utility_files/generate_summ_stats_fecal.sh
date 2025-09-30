#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J summ_stats
#SBATCH --cpus-per-task=16
#SBATCH --mail-user email@emailexample.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to calculate summary statistics for panel performance

# Load required modules
source /usr/modules/init/bash
module load bcftools
module load samtools

# Define variables
MERGED="/path/to/merged.vcf.gz"
BAMDIR="/path/to/bam_dir"
TARGETS_BED="/path/to/bedfile.bed"  # BED file with amplicon target regions
THREADS=16  # should match cpus-per-task

# Create output directory
OUTDIR="/path/to/output/directory/$(date +%y%m%d)_run_stats"
mkdir -p "$OUTDIR"

# 1. Depth per locus (locus-level)
echo "Calculating depth per locus..."
bcftools query -f '%CHROM\t%POS[\t%DP]\n' "$MERGED" > "$OUTDIR/depth_per_locus.tsv"

# 2. Genotype completeness per sample (sample-level)
echo "Calculating genotype completeness per sample..."
bcftools query -f '[%SAMPLE\t%GT\n]' "$MERGED" | \
awk '$2 != "./." { count[$1]++ } END { for (s in count) print s, count[s] }' \
> "$OUTDIR/completeness_per_sample.tsv"

# 3. Percent reads on target
echo "Calculating % reads on target..."
echo -e "Sample\tTotal_Reads\tReads_on_Target\tPercent_on_Target" > "$OUTDIR/reads_on_target.tsv"

cd "$BAMDIR" || exit 1
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

# --- Map sample to library and pool ---
awk '
NR==1 {next}  # skip header
{
    sample=$1
    lib=substr(sample,1,4)
    pool="NA"
    if (lib ~ /^i023|i024|i025|i026$/) pool="Pool1"
    else if (lib ~ /^i038|i046|i047|i003$/) pool="Pool2"
    else if (lib ~ /^i032|i001|i002|i036|i037$/) pool="Pool3"
    else if (lib ~ /^i028|i031|i033|i048$/) pool="Pool4"
    print sample, lib, pool
}' "$OUTDIR/reads_on_target.tsv" > "$OUTDIR/sample_library_pool.tsv"

# --- Mean depth per locus per library, pool, and overall ---
# Transpose depth_per_locus.tsv: columns=sample, rows=locus
awk '
NR==FNR && FILENAME=="'$OUTDIR/sample_library_pool.tsv'" { lib[$1]=$2; pool[$1]=$3; next }
NR>FNR && FILENAME=="'$OUTDIR/depth_per_locus.tsv'" {
    locus=$1":"$2
    for(i=3;i<=NF;i++) {
        sample="SAMPLE"i-2
        dp[locus,sample]=$i
    }
    n_samples=NF-2
    next
}' "$OUTDIR/sample_library_pool.tsv" "$OUTDIR/depth_per_locus.tsv"

# We'll compute mean per library/pool by summing sample depths
# Simplified approach: for small numbers of samples, can do in R or Python; in bash/awk, easier to aggregate using sample-library mapping
# For now, output overall mean depth per locus
awk '
NR>1 {
    sum=0; n=0
    for(i=3;i<=NF;i++) { sum+=$i; n++ }
    mean=sum/n
    print $1,$2,mean
}' "$OUTDIR/depth_per_locus.tsv" > "$OUTDIR/depth_per_locus_mean.tsv"

# --- Combine metrics for summary ---
awk '
BEGIN {OFS="\t"; print "Library","Pool","Mean_Depth","Mean_Completeness","Mean_Total_Reads","Mean_OnTarget","Mean_PctOnTarget"}
NR==FNR && FILENAME=="'$OUTDIR/completeness_per_sample.tsv'" { comp[$1]=$2; next }
NR==FNR && FILENAME=="'$OUTDIR/reads_on_target.tsv'" { total[$1]=$2; ontarget[$1]=$3; pct[$1]=$4; next }
NR==FNR && FILENAME=="'$OUTDIR/sample_library_pool.tsv'" { lib[$1]=$2; pool[$1]=$3; next }
END {
    # Aggregate by library and pool
    for(s in lib) {
        l=lib[s]; p=pool[s]
        sum_comp[l]+=comp[s]; sum_total[l]+=total[s]; sum_ontarget[l]+=ontarget[s]; sum_pct[l]+=pct[s]; n[l]++
         sum_comp_pool[p]+=comp[s]; sum_total_pool[p]+=total[s]; sum_ontarget_pool[p]+=ontarget[s]; sum_pct_pool[p]+=pct[s]; n_pool[p]++
        sum_comp_all+=comp[s]; sum_total_all+=total[s]; sum_ontarget_all+=ontarget[s]; sum_pct_all+=pct[s]; N_all++
    }
    for(l in n) print l,"NA", "NA", sum_comp[l]/n[l], sum_total[l]/n[l], sum_ontarget[l]/n[l], sum_pct[l]/n[l]
    for(p in n_pool) print "NA",p, "NA", sum_comp_pool[p]/n_pool[p], sum_total_pool[p]/n_pool[p], sum_ontarget_pool[p]/n_pool[p], sum_pct_pool[p]/n_pool[p]
    print "ALL","ALL","NA", sum_comp_all/N_all, sum_total_all/N_all, sum_ontarget_all/N_all, sum_pct_all/N_all
}' "$OUTDIR/completeness_per_sample.tsv" "$OUTDIR/reads_on_target.tsv" "$OUTDIR/sample_library_pool.tsv" > "$OUTDIR/summary_library_pool_overall.tsv"

echo "All summary statistics written to: $OUTDIR"
