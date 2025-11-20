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

