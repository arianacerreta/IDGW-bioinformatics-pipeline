# IDGW-bioinformatics-pipeline
This is a repository for how to process raw genomic data from the Idaho Gray Wolf Project.

## So the sequencing facility has sent you a link to your sequencing data...
This step assumes that you already have/know how to do the following:
- Access your institution's server
- Log into your server account through the console or an Ubuntu shell
- Basic commands in Bash such as mkdir, cd, and how to direct to a directory
### Download the data
After logging into your server account, you will navigate into the appropriate directory (this is where the raw data will download) and copy this line of code that was sent to you from the sequencing facility and input it into your console. Example:

```wget -r -c -nH --cut-dirs=1 --no-parent --reject "index.html*" https://gc3fstorage.uoregon.edu/example/```

The files should begin downloading. Wait for this step to complete before moving on to the remaining steps.

## Demultiplex the reads
This is how we go from all the raw sequences that are in one big fastq.gz file to fastq files for each well for a given plate. you will need to make a separate file that essentially outlines everything the computer will need to match up the well (e.g. A1) with the sample (e.g. Samp231) on a given plate (e.g. Plate3).

### Prepare your demultiplex file
1. Open excel or text editor of choice
2. Create headings in these exact order, spelling, and capitalization:
   a.	Sample, PlateID, i7_name, i7_sequence, i5_name, i5_sequence
3. Refer to the physical lab (or digital if they have been updated) PCR sheets and locate which i7 tagging primers were used for each library
4. If prepared in CGMEC at the University of Idaho, refer to [i7_primer_info.csv](utility_files/i7_primer_info.csv) under utility_files and choose the right barcode for the sequencing machine used:
   a. NovaSeq: for i7 use the Reverse Complement
   b. MiSeq: for i7 use the Reverse Complement
   c. NextSeq: for i7 use the Reverse Complement
6. If prepared in CGMEC at the University of Idaho, refer to [i5_primer_info.csv](utility_files/i5_primer_info.csv) and fill in the i5 barcode and corresponding well:
   a. NovaSeq: for i5 use the Reverse Complement
   b. MiSeq: for i5 use the Barcode
   c. NextSeq: for i5 use the Barcode
7. Save as .csv. Example:
   | Sample | PlateID | i7_name | i7_sequence | i5_name | i5_sequence |
   |--------|---------|---------|-------------|---------|-------------|
   | Samp321 | Plate3 | i017 | CTCATC | A01 | CCGTTT |
   | Samp654 | Plate3 | i017 | CTCATC | B01 | AAGAGT |

### Run demultiplex code on the server
1. Download [GTseq_BarcodeSplit_KML.py](https://github.com/kiralong/gtseq_ref_align/tree/main/Main_pipeline) from Kira Long's [gtseq_ref_align](https://github.com/kiralong/gtseq_ref_align/tree/main) pipeline. This code was modified by Kira based on [Nate Campbell's demultiplexing script](https://github.com/GTseq/GTseek_utils/blob/Main/GTseq_BarcodeSplit_MP.py)
2. Read the details on Kira's page for further information. I have provided queue_demultiplex.sh file for you to download and edit, but you could also copy the example code from Kira's page.


