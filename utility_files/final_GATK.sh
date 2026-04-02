#!/usr/bin/env bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J GATK
#SBATCH --cpus-per-task=16
#SBATCH --mem 192G
#SBATCH --mail-user example@email.com
#SBATCH --mail-type=BEGIN,END,FAIL

#Load the required modules
source /usr/modules/init/bash
module load gatk/4.1.7.0

#Variables and Paths
REF="/path/to/genome/genome.fna"   #path ending in genome .fna
INDIVVCF="/path/to/individual/vcfs/indiv_vcfs_DPcap500"      #directory for individual vcfs; lol just realized this is never called, but should still be good for reference
BEDFILE="/path/to/bed/zzzz/200neutral.bed"               #path to where .bed file stored
FINALDIR="/path/to/directory/GATK/final_vcf" #new: will be made
MAPFILE="/path/to/map/in/indiv_vcfs_DPcap500/sample_map.txt" #path to the .txt made in the previous shell code

mkdir -p "$FINALDIR"

gatk --java-options "-Xmx96G" GenomicsDBImport \
  --genomicsdb-workspace-path "/path/to/GATK/genomics_db" \
  --batch-size 150 \
  -L "$BEDFILE" \
  --interval-padding 3 \
  --sample-name-map "$MAPFILE" \
  --merge-input-intervals true \
  --reader-threads 16

 echo "Combining GVCFs"

 gatk --java-options "-Xmx48G" GenotypeGVCFs \
  -R "$REF" \
  -V gendb:///path/to/GATK/genomics_db \
  -L "$BEDFILE" \
  --interval-padding 0 \
  --include-non-variant-sites true \
  -O "$FINALDIR/all.gatk.genotyped.vcf.gz"
