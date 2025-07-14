# IDGW-bioinformatics-pipeline
This is a repository for how to process raw genomic data from the Idaho Gray Wolf Project.

## So the sequencing facility has sent you a link to your sequencing data...
This step assumes that you already have/know how to do the following:
- Access your institution's server
- Log into your server account through the console or an Ubuntu shell
- Basic commands in Bash such as mkdir, cd, and how to direct to a directory

Otherwise, I tried to include as many details for people new to using bash in command line. If you are already comfortable working in Unix and on servers, feel free to skip those details.

### Download the data
After logging into your server account, you will navigate into the appropriate directory (this is where the raw data will download) and copy this line of code that was sent to you from the sequencing facility and input it into your console. Example:

```wget -r -c -nH --cut-dirs=1 --no-parent --reject "index.html*" https://gc3fstorage.uoregon.edu/example/```

The files should begin downloading. Wait for this step to complete before moving on to the remaining steps.

## Demultiplex the reads
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
3. Place both the GTseq_BarcodeSplit_KML.py and queue_demultiplex.sh in your scripts directory on the server.
4. For the first time you are using the script, you will need to alter the permissions on the file:
   
    - In your scripts directory run: ```chmod 755 GTseq_BarcodeSplit.KML.py``` and ```chmod 755 queue_demultiplex.sh```

5. Beyond that, I would recommend that you don't edit GTseq_BarcodeSplit_KML.py unless you know what you are doing. Look for updates on Kira's github.
6. Copy over the demultiplex.csv to an appropriate folder on your server account.
7. Make an output directory called individual_fastqs (or something similar)
8. Edit queue_demultiplex.sh:

   ```nano queue_demultiplex.sh```
   
    - Edit line 2 for the appropriate partition
    - Edit line 6 with your email
    - Edit all paths in Lines 14-17
    - If needed edit Line 20 to direct to your scripts directory
    - Hit Ctrl X -> Y -> enter to save
      
9. Run the code on your server. If you have a SLURM scheduler (e.g. fortyfive at IIDS at University of Idaho), use:

    ```sbatch queue_demultiplex.sh```
    
*For reference, this step has taken ~45 minutes to demultiplex a NovaSeq 6000 run with 4 plates.*

### Zip all those individual fastqs
This will increase processing speed in subsequent steps
1. Go to the directory before your individuals_fastq directory and check to see that your individual_fastqs directory is in that folder using ```ls```
2. Use the code: ```gzip -r individual_fastqs/``` to recursively zip every file in the directory listed

*This will take a while (about ~1.5 hours). If you are feeling ambitious, you could write some code to zip the files in parallel, but I haven't done that yet. Once/if I do, I'll update this section.*


