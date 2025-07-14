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
4. If prepared in CGMEC, refer to "2024-08-22_i5_i7_Primer-Info_corrected.xlsx" 
