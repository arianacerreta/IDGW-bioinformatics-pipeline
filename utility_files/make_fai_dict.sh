#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH --cpus-per-task=1
#SBATCH -J fai_dict
#SBATCH --mail-user email@uidaho.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load gatk
module load samtools

cd path/to/reference_genomes/

gatk CreateSequenceDictionary \
  -R GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna \
  -O GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.dict

samtools faidx /path/to/GCA_000002285.4_Dog10K_Boxer_Tasha_genomic.fna
