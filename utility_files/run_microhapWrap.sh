#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J microhapWrap_200_NTiss
#SBATCH --cpus-per-task=16
#SBATCH --mem=92G
#SBATCH --mail-user=email@email.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Safety settings
set -euo pipefail
trap 'echo "Error on line $LINENO"; exit 1' ERR

# Load samtools
source /usr/modules/init/bash
module purge
module load bowtie2
module load samtools
module load python/3.8.11        # or your python module

# variables and paths
mtype2="/path/Programs/microTyper/mtype2"
positions="/path/Genomes/wolfAmpRef_200/Clu200_pos_4.txt"
primers="/path/Genomes/primers_with_adapters.txt"
pres_abs="/path/Genomes/wolfAmpRef_5_OG_341loci/Clu341_pres_abs.txt"
fasta="/path/Genomes/wolfAmpRef_200/wolfAmpRef_200.fa"
FASTQ_DIR="/path/trimmed_fastp/min_length_50"
index_prefix="/path/Genomes/wolfAmpRef_200/wolfAmpRef_200"
bams_output_dir="/path/delomas-pipeline-scat-settings/bams"
geno_out="/path/delomas-pipeline-scat-settings/genos_OG_sry"
python_code="/path/scripts/microhapWrap.py"

mkdir -p "$bams_output_dir" "$geno_out"

cd "$FASTQ_DIR" || exit 1

# Debug info
echo "Running microhapWrap.py at $(date)"
echo "Using $(which python3)"
echo "CPU threads: $SLURM_CPUS_PER_TASK"

# Add mtype2 path so subprocess can find it
export PATH="/path/Programs/microTyper:$PATH"
which mtype2 || echo "⚠️  mtype2 not found in PATH"

# Align FASTQ files to BAM
echo "Aligning FASTQ files...(skipping existing BAMs)"
BAM_FILES=()

for fq in "$FASTQ_DIR"/*.fq.gz; do
    sample=$(basename "$fq" .fq.gz)
    bam="$bams_output_dir/${sample}.bam"

    #if BAM already exists and is non-empty, skip
    if [[ -s "$bam" ]]; then
        echo "✅ Skipping alignment for $sample — BAM already exists."
    else
    	echo "Aligning $fq -> $bam"
    	bowtie2 --local --no-unal -N 1 --rdg 5,5 --rfg 5,5 \
        	-x "$index_prefix" -U "$fq" 2> "$bams_output_dir/${sample}_bowtie2.log" | \
        	samtools view -F 276 -u | samtools sort -o "$bam"

    # Check BAM is non-empty
    	if [[ ! -s "$bam" ]]; then
        	echo "⚠️  BAM file $bam is empty! Check Bowtie2 logs: $Fbams_output_dir${sample}_bowtie2.log"
        	exit 1
    	fi
    fi
BAM_FILES+=("$bam")
done

echo "All FASTQ files aligned successfully or skipped (if already aligned)."

cd "$geno_out"

# Run the wrapper. 
echo "Running microhapWrap.py..."
python3 "$python_code" \
	-bam "${BAM_FILES[@]}" \
	-f "$FASTQ_DIR"/*.fq.gz \
	-p "$positions" \
	-r "$fasta" \
	-pa "$pres_abs" \
	--alignByLocus \
	-t "$SLURM_CPUS_PER_TASK"

echo "Run complete at $(date)"
