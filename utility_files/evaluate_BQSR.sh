!/usr/bin/env bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J BQSR_eval
#SBATCH --cpus-per-task=4
#SBATCH --mail-user email@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

#Load the required modules
source /usr/modules/init/bash
module load gatk
module load R

#Variables and Paths
REF="/path/to/genome/genome.fna"   # The genome needs to have a.fai and .dict file associated with it
KNOWN="/path/to/known_sitesonly.sorted.CanFam3_1.vcf.gz"                 # has .tbi file associated with it
BAMDIR="/path/to/RG_aligned_sorted/bams"           # where *.bam live
RECALDIR="/path/to/recalibrated/directory"                    # where *.recal.table were were written
RECALBAMDIR="/path/to/recalibrated/directory/recal_bams" #path where the actual recalibrated .bams live
TABDIR="/path/to/new/dir/recalibration_tables_for_plots"     # new: where we’ll put before/after tables
PLOTDIR="/path/to/new/dir/recalibration_plots"               # new: where PDFs will go
THREADS=4                                            # number of parallel samples to process

export R_LIBS_USER=/mnt/ceph/user/R/x86_64-pc-linux-gnu-library/4.2
mkdir -p "$TABDIR" "$PLOTDIR"

# Build list of samples to process:
LIST="samples_to_plot.txt"
if [[ -f "$LIST" ]]; then
  mapfile -t SAMPLES < <(sed 's/.bam$//' "$LIST")
else
  mapfile -t SAMPLES < <(ls "$BAMDIR"/*.bam | head -n 8 | xargs -n1 -I{} basename {} .bam)
fi

# Sanity message
echo "Will generate plots for ${#SAMPLES[@]} samples:"
printf '  - %s\n' "${SAMPLES[@]}"

# Export vars for GNU parallel
export REF KNOWN BAMDIR RECALDIR TABDIR PLOTDIR

# Run per sample (in parallel)
printf '%s\n' "${SAMPLES[@]}" | parallel -j "$THREADS" '
  S={}

  BAM="$BAMDIR/${S}.bam"
  AFTER_TABLE="$RECALDIR/${S}.recal.table"

  # Check that input BAM exists
  if [[ ! -f "$BAM" ]]; then
    echo "[WARN] Missing BAM for $S: $BAM" >&2
    exit 0
  fi

  # Check that AFTER table exists
  if [[ ! -f "$AFTER_TABLE" ]]; then
    echo "[ERROR] AFTER table for $S does not exist: $AFTER_TABLE" >&2
    exit 1
 else
    echo "[OK] Using AFTER table for $S: $AFTER_TABLE"
  fi

  # BEFORE table from original BAM
  gatk BaseRecalibrator \
    -R "$REF" -I "$BAM" --known-sites "$KNOWN" \
    -O "$TABDIR/${S}.before.table" \
    > "$TABDIR/${S}.before.log" 2>&1

  # Generate plots
  gatk AnalyzeCovariates \
    -before "$TABDIR/${S}.before.table" \
    -after "$AFTER_TABLE" \
    -plots "$PLOTDIR/${S}.recalibration_plots.pdf" \
    > "$PLOTDIR/${S}.plots.log" 2>&1

  echo "Done: $S → $PLOTDIR/${S}.recalibration_plots.pdf"
'
