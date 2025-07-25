# IDGW-bioinformatics-pipeline
This is a repository for how to process raw genomic data from the Idaho Gray Wolf Project.

## 1. So the sequencing facility has sent you a link to your sequencing data...
This step assumes that you already have/know how to do the following:
- Access your institution's server
- Log into your server account through the console or an Ubuntu shell
- Basic commands in Bash such as `mkdir`, `cd`, and how to direct to a directory

Otherwise, I tried to include as many details for people new to using bash in command line. If you are already comfortable working in Unix and on servers, feel free to skip those details.

### Download the data
After logging into your server account, you will navigate into the appropriate directory (this is where the raw data will download) and copy this line of code that was sent to you from the sequencing facility and input it into your console. Example:

```wget -r -c -nH --cut-dirs=1 --no-parent --reject "index.html*" https://gc3fstorage.uoregon.edu/example/```

The files should begin downloading. Wait for this step to complete before moving on to the remaining steps.

## 2. Demultiplex the reads
This is how we go from all the raw sequences that are in one big fastq.gz file to fastq files for each well for a given plate. you will need to make a separate file that essentially outlines everything the computer will need to match up the well (e.g. A1) with the sample (e.g. Samp231) on a given plate (e.g. Plate3).

### Prepare your demultiplex file
1. Open excel or text editor of choice
2. Create headings in these exact order, spelling, and capitalization:

    - Sample, PlateID, i7_name, i7_sequence, i5_name, i5_sequence

3. Refer to the physical lab (or digital if they have been updated) PCR sheets and locate which i7 tagging primers were used for each library
4. If prepared in CGMEC at the University of Idaho, refer to [i7_primer_info.csv](utility_files/i7_primer_info.csv) under utility_files and choose the right barcode for the sequencing machine used:

    - **NovaSeq**: for i7 use the Reverse Complement
    - **MiSeq**: for i7 use the Reverse Complement
    - **NextSeq**: for i7 use the Reverse Complement
   
5. If prepared in CGMEC at the University of Idaho, refer to [i5_primer_info.csv](utility_files/i5_primer_info.csv) and fill in the i5 barcode and corresponding well:
   
    - **NovaSeq**: for i5 use the Reverse Complement
    - **MiSeq**: for i5 use the Barcode
    - **NextSeq**: for i5 use the Barcode
   
6. Save as .csv (referenced as demultiplex.csv in subsequent instructions). Example:
   | Sample | PlateID | i7_name | i7_sequence | i5_name | i5_sequence |
   |--------|---------|---------|-------------|---------|-------------|
   | Samp321 | Plate3 | i017 | CTCATC | A01 | CCGTTT |
   | Samp654 | Plate3 | i017 | CTCATC | B01 | AAGAGT |

### Run demultiplex code on the server
1. Download [GTseq_BarcodeSplit_KML.py](https://github.com/kiralong/gtseq_ref_align/tree/main/Main_pipeline) from Kira Long's [gtseq_ref_align](https://github.com/kiralong/gtseq_ref_align/tree/main) pipeline. This code was modified by Kira based on [Nate Campbell's demultiplexing script](https://github.com/GTseq/GTseek_utils/blob/Main/GTseq_BarcodeSplit_MP.py)
2. Read the details on Kira's page for further information. I have provided [queue_demultiplex.sh](utility_files/queue_demultiplex.sh) file for you to download and edit, but you could also copy the example code from Kira's page.
3. Place both the `GTseq_BarcodeSplit_KML.py` and `queue_demultiplex.sh` in your scripts directory on the server.
4. For the first time you are using the script, you will need to alter the permissions on the file:
   
    - In your scripts directory run: ```chmod 755 GTseq_BarcodeSplit.KML.py``` and ```chmod 755 queue_demultiplex.sh```

5. Beyond that, I would recommend that you don't edit `GTseq_BarcodeSplit_KML.py` unless you know what you are doing. Look for updates on Kira's github.
6. Copy over the demultiplex.csv to an appropriate folder on your server account.
7. Make an output directory called individual_fastqs (or something similar)
8. Edit `queue_demultiplex.sh`:

   ```nano queue_demultiplex.sh```
   
    - Edit line 2 for the appropriate partition
    - Edit line 6 with your email
    - Edit all paths in Lines 14-17
    - If needed edit Line 20 to direct to your scripts directory
    - Hit Ctrl X -> Y -> enter to save
      
9. Run the code on your server. If you have a SLURM scheduler (e.g. fortyfive at IIDS at University of Idaho), use:

    ```sbatch queue_demultiplex.sh```
    
*For reference, this step has taken ~45 minutes to demultiplex a NovaSeq run with 4 plates.*

### Zip all those individual fastqs
This will increase processing speed in subsequent steps
1. Go to the directory before your individuals_fastq directory and check to see that your individual_fastqs directory is in that folder using ```ls```
2. Use the code: ```gzip -r individual_fastqs/``` to recursively zip every file in the directory listed

*This will take a while (about ~1.5 hours). If you are feeling ambitious, you could write some code to zip the files in parallel, but I haven't done that yet. Once/if I do, I'll update this section.*

## 3. Trim Illumina adapter sequences and ployG tails
This is directly from [Kira Long's pipeline](https://github.com/kiralong/gtseq_ref_align), Step 2: Trim and remove adaptors. I have provided detailed instructions on how implement her shell code to run this step.

1. Open text editor of choice and make a name_map.tsv file with the name of each sample on a separate line. For example, Samp321.fastq.gz will become Samp321 in the list. A final format will look something like:

    Samp321

    Samp654

    Samp987

    Samp494

2. Make sure to save your .tsv with Unix line endings or this step will not work!
3. Transfer name_map.tsv to an appropriate location on your server.
4. Make directories such as ```work_dir/trimmed_fastp/min_length_50``` for your output files
5. Download [run_fastp.sh](https://github.com/kiralong/gtseq_ref_align/blob/main/Main_pipeline/run_fastp.sh) and place in your scripts directory
6. Edit ```run_fastp.sh```

    ```nano run_fastp.sh```
    - Edit header for IIDS server:
        - #SBATCH -p eight
        - #SBATCH -C "ceph"
        - #SBATCH -J Clu_fastp
        - #SBATCH --cpus-per-task=16
        - #SBATCH --mail-user `email@exampleemail.com`
        - #SBATCH --mail-type=BEGIN,END,FAIL
    - For IIDS server add ```source /usr/modules/init/bash``` in a separate line before ```module load fastp``` in line 9
    - Edit lines 12-14 with the appropriate directory paths
    - If your .tsv only has one column, delete ```cut -f 2 |``` from line 18

7. Save and give permissions to the shell code (```chmod 755 run_fastp.sh```)
8. Run ```sbatch run_fastp.sh```

*This step usually goes quickly (<10 min). Double-check the output directory to see if the last sample listed in the .tsv was trimmed. For some reason, this step sometimes leaves off the last sample. If it was left off, just make a .tsv with only that sample and rerun the `run_fastp.sh` but direct it to the .tsv with only one sample listed.*

## 4. Reference genome
This step is also directly based on [Kira Long's pipeline](https://github.com/kiralong/gtseq_ref_align), Step 3: Align to a reference genome. I have tried to write explicit instructions to elaborate on what Kira outlines in her pipeline to aid first-time users.

### Download your reference genome
*If you already have the reference genome you would like to use downloaded and a reference database created, then skip to **Align trimmed sequences**.*
This example is done with the current genome we are using for the Idaho Gray Wolf Project so that future researchers and students on the project can replicate the steps.

1. Go to [GenBank](https://www.ncbi.nlm.nih.gov/) and search "dog genome"
2. Select "Genomes" and find and click on "Dog10k_Boxer_Tasha"
3. On the row Submitted GenBank assembly find the three dots under "Actions" -> see more files on FTP -> right click on the fna.gz file -> copy link address
4. Back in the console in the directory you would like to store your genome, use the following code:

    ```wget link-from-step3```

    to download the genome.

5. Create a corresponding README (```nano README```) in the same folder with the following information:

    - Genome used: GCA_000002285.4, Dog10K_Boxer_Tasha, submitted GenBank assembly
    - Taxon: Canis lupus familiaris (dog)
    - Link: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002.285/GCA_000002285.4_Dog10K_Boxer_Tasha/
    - Notes: Downloaded the fna.gz

### Make a reference database for genome
1. Download [bwa_db.sh](https://github.com/kiralong/gtseq_ref_align/blob/main/Main_pipeline/bwa_db.sh)
2. Edit bwa_db.sh

    ```nano bwa_db.sh```
    - Edit header for IIDS server:
        - #SBATCH -p eight
        - #SBATCH -C "ceph"
        - #SBATCH --cpus-per-task=1
        - #SBATCH -J clu_bwa_db
        - #SBATCH --mail-user `email@emailexample.edu`
        - #SBATCH --mail-type=BEGIN,END,FAIL
    - Add ```source /usr/modules/init/bash``` after line 8 to work with IIDS servers
    - Edit lines 13-15 with your paths
    - *See Kira's github for explanation of ```$pref```*

3. Save and give your script permissions (```chmod 755 bwa_db.sh```)
4. Run ```sbatch bwa_db.sh``` in the console

*This step took ~3 hours for this genome but is highly dependent on genome size. If you would like to check on the run progress of this job, run ```squeue --me``` in the console.*

## 5. Align trimmed sequences 
This step is also directly based on [Kira Long's pipeline](https://github.com/kiralong/gtseq_ref_align), Step 3: Align to a reference genome. I have tried to write explicit instructions to elaborate on what Kira outlines in her pipeline to aid first-time users.

1. Make a directory for your aligned sequences. I recommend something like:

    ```~/aligned/Clu10kTash_fastp_trimmed```

2. Download [bwa_alignment.sh](https://github.com/kiralong/gtseq_ref_align/blob/main/Main_pipeline/bwa_alignment.sh) from Kira's github and put in your scripts directory.
3. Edit ```bwa_alignment.sh```

    ```nano bwa_alignment.sh```
    - Edit header for IIDS server:
        - #SBATCH -p eight
        - #SBATCH -C "ceph"
        - #SBATCH -J align_trimmed_clu_gtseq_data
        - #SBATCH --cpus-per-task=12    
        - #SBATCH --mail-user `email@emailexample.edu`
        - #SBATCH --mail-type=BEGIN,END,FAIL
    - Add ```source /usr/modules/init/bash``` after line 9 to work with IIDS servers
    - Edit lines 14-17 with your paths
    - If your .tsv only has one column, delete ```cut -f 2 |``` from line 20

4. Run ```sbatch bwa_alignment.sh``` in the console

*I have bumped cpus-per-task from 12 to 24 with no issues, but double-check server resources before doing so. This step has taken about 8 hours to run (with 12 cpus) for ~1000 samples.*

## 6. Genotyping
This is where this pipeline diverges from Kira's. If you are dealing with microhaplotypes, short indels, or positions that are highly likely to not be variable (i.e., diagnostic loci for species and you only have one species on your plates) continue with this pipeline which will use the software ```bcftools```, ```samtools```, and ```htslib```.

### Filter .bam files by mapping quality
1. If you have >1000 .bam files, run ```ulimit -n 5000``` to increase the amount of open files you can open at once. *5000 can be replaced with any reasonable number*
2. Make a new directory for your filtered .bam files. I recommend something like ```~/samtools_filtered_BAMS```
3. Download [run_samtools_filter.sh](utility_files/run_samtools_filter.sh) and put in scripts directory
4. Edit ```run_samtools_filter.sh```

    ```nano run_samtools_filter.sh```
    - Edit email in header
    - Edit variable and paths in lines 17-21

5. Save and give permissions, if needed (```chmod 755 run_samtools_filter.sh```)
6. Run ```sbatch run_samtools_filter.sh```

*A run of .bams from ~1000 samples took about 10 minutes. You can check the slurm record (```less slurm-#####.out```; to exit enter "q") to see the progress of the run. It will update with "Finished processing SAMPLE NAME" after each sample's .bam files have been filtered.*

### Call genotypes using .bed file and bcftools
1. Make sure your genome is unzipped. ```gzip -dk ./path/genome.fna.gz``` (You use "./" to orient to a subdirectory of your current directory)
2. Download [run_bcftools_call_separate.sh](utility_files/run_bcftools_call_separate.sh) and place in scripts directory
3. Create a bam_list.txt of the filtered .bam files you would like to genotype. This allows you to do a subset if needed. Don't forget to save with Unix line endings. One file name (e.g., i001_B04_GWAdapt7_UI1887_3E61_filtered.bam) per line.
4. Make a .bed file which lists all the chr, pos, pos that you would like to be called. Adaptive panel .bed files: [adaptive_positions_canFam3_1.bed](utility_files/adaptive_positions_canFam3_1.bed) and [adaptive_positions_Clu10kTash.bed](utility_files/adaptive_positions_Clu10kTash.bed). Neutral panel (modified subset from this [article](https://doi.org/10.1007/s12686-023-01301-x)): tempname_canFam3_1.bed and tempname_Clu10kTash.bed.
5. Double check your file limit with ```ulimit -n```. Increase your limit with ```ulimit -n 5000``` if needed.
6. Create a new directory for your bcftools runs. I recommend something like ```~/bcftools_runs```
7. Edit ```run_bcftools_call_separate.sh```

    ```nano run_bcftools_call_separate.sh```
    - Edit email in header
    - Edit variables and paths in lines 18-24
    - Edit your maxDP and minBQ in lines 26-27

8. Save and give permissions, if needed (```chmod 755 run_bcftools_call_separate.sh```).
9. Run ```sbatch run_bcftools_call_separate.sh```.

*This step most recently took 6 hours for ~1000 samples with CanFam3.1*

### Normalize and merge individual vcfs into once big vcf
1. Double check you have your file open limit high enough. ```ulimit -n 5000```
2. Download [run_merge.sh](utility_files/run_merge.sh) and place in your scripts directory
3. Edit ```run_merge.sh```

    ```nano run_merge.sh```
    - Edit email in header
    - Edit paths in lines 17-20

4. Save and give permissions, if needed (```chmod 755 run_merge.sh```)
5. Run ```sbatch run_merge.sh```

*This took about 15 minutes to run to merge ~1000 vcfs*

6. After the run has completed, direct to the merged directory and run ```module load htslib```
7. Unzip the vcf with ```bgzip -c -d merged.vcf.gz > unzipped.merged.vcf```
8. Inspect your unzipped.merged.vcf and see if all the loci and calls are expected. If you need to further filter your vcf by depth, to remove unnecessary loci (i.e., you had cast a wide net for an indel), or remove alternate alleles,nano move on to the other steps. Otherwise, this vcf can be used for subsequent analyses in R or other programs.

### Filter one more time (depth, remove unnecessary loci)
1. Download [run_bcftools_filter.sh](utility_files/run_bcftools_filter.sh) and place in scripts directory
2. If you have a list of loci to filter out, make a .txt with those positions as chr, bp and no header.
3. Edit ```run_bcftools_filter.sh```

    ```nano run_bcftools_filter.sh```
    - Edit email in header
    - Edit variables and paths in lines 16-19; DP is depth per locus for each individual
    - If you are not using remove_targets.txt with a list of targets to remove, then you should use this simpler shell code [run_bcftools_filter_simple.sh](utility_files/run_bcftools_filter_simple.sh)

4. Save and give permissions, if needed (```chmod 755 run_bcftools_filter.sh```)
5. Run ```sbatch run_bcftools_filter.sh```

*This should happen really fast (<2 min). This is currently the final step on the server. You can use the final_filtered_unzipped.vcf for R and other programs.*

# Final filtering and error rates in R

## 1. Download and install SNPfiltR[SNPfiltR](https://devonderaad.github.io/SNPfiltR/) (DeRaad 2022)
1. You can use this install_SNPfiltR.R[install_SNPfiltR.R](utility_files/R_code/install_SNPfiltR.R) code to help install this package. You can also refer to the package's github linked above.

