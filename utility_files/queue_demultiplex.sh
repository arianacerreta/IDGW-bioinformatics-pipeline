#!/bin/bash
#SBATCH -p eight
#SBATCH -C "ceph"
#SBATCH -J demultiplexing
#SBATCH --cpus-per-task=16
#SBATCH --mail-user email@emailexample.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load required modules
source /usr/modules/init/bash
module load python

# Set variables and paths
work_dir=/workingdir/path
barcodes=$work_dir/path/to/demultiplex.csv
fastq=$work_dir/path/to/raw_data/zipped-data.fq.gz
out_dir=$work_dir/outdir/path

# Go to the script directory
cd $work_dir/scripts

# Run the demultiplexing script
./GTseq_BarcodeSplit_KML.py $barcodes $fastq $out_dir
