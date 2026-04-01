## 6. Genotyping
### RG headers, sort, and index
Add RG headers to your .bam files and then sort and index them for subsequent steps using GATK. GATK will not work if you do not do this. This step will be highly dependent on your naming convention. I have included the shell code for the different sample pipelines we have used with the IDGW project thus far. GWAdapt ([RG_sort_index.sh](project-specific_files/RG_sort_index.sh)); Fecal Neutral 200 loci ([RG_sort_index_fecal.sh](project-specific_files/RG_sort_index_fecal.sh)); IDFGdata ([RG_sort_index_IDFGdata.sh](project-specific_files/RG_sort_index_IDFGdata.sh)); Neutral Tissue ([RG_sort_index_NeutralTissue.sh](project-specific_files/RG_sort_index_NeutralTissue.sh))
1. Download appropriate shell code and put in scripts directory.
2. Make a new directory for your .bams for this step. I recommed something like ```~/RG_aligned_sorted_Clu10kTash```
3. Edit correspond shell code for the appropriate samples.

   ```nano RG_sort_index.sh```
   - Edit email in header
   - Edit variables and paths in lines 14-15

4. Save and give permissions, if needed (```chmod 755 RG_sort_index.sh```)
5. Run ```sbatch RG_sort_index.sh```.

### Base recalibration
This recalibrates the base quality scores using known, confident SNPs. Since we run amplicon panels, these are the literature referenced SNPs for the GWAdapt Panel and the known variants fromt the Neutral Panel. The VCFs need to match the alignment you are using.

1. If you don't already have a .fai and .dict files associated with your reference genome, then do steps 1-3, otherwise skip to 4. Download [make_fai_dict.sh](utility_files/make_fai_dict.sh) and edit for the followling:

    ```nano make_fai_dict.sh```
     - Edit email in header
     - Edit variables and paths in lines 14 and 20
       
2. Save and give permissions, if needed (```chmod 755 make_fai_dict.sh```)
3. Run ```sbatch make_fai_dict.sh```)
4. Now that you have your .fai and .dict files in the same directory as your genome, download the corresponding VCF file and corresponding .tbi file. Make sure the .tbi file is in the same directory as the VCF. GWAdapt (CanFam3.1: [known77.sites.vcf.gz](project-specific_files/known77.sites.vcf.gz) and [known77.sites.vcf.gz.tbi](project-specific_files/known77.sites.vcf.gz.tbi)); Neutral 200 loci (Clu10kTash: [known_sites.sitesonly.Clu10kTash.vcf.gz](project-specific_files/known_sites.sitesonly.Clu10kTash.vcf.gz) and [known_sites.sitesonly.Clu10kTash.vcf.gz.tbi](project-specific_files/known_sites.sitesonly.Clu10kTash.vcf.gz.tbi); CanFam3.1: [known_sitesonly.sorted.CanFam3_1.vcf.gz](project-specific_files/known_sitesonly.sorted.CanFam3_1.vcf.gz) and [known_sitesonly.sorted.CanFam3_1.vcf.gz.tbi](project-specific_files/known_sitesonly.sorted.CanFam3_1.vcf.gz.tbi))
5. Download [run_BQSR.sh](utility_files/run_BQSR.sh) and make a new directory for your recalibrated .bams. I recommend making sure it is not nested in your previous folder.
6. Edit correspond shell code for the following:

   ```nano run_BQSR.sh```
   - Edit email in header
   - Edit variables and paths in lines 14-18

7. Save and give permissions, if needed (```chmod 755 run_BQSR.sh```)
8. Run ```sbatch run_BQSR.sh```

*For ~3000 samples, this took about 4.5 hours*

#### Optional:
#### Evaluate whether base recalibration worked

1. Download [evaluate_BQSR.sh](utility_files/evaluate_BQSR.sh) and edit the following:

   ```nano evaluate_BQSR.sh```
   - Edit email in header
   - Edit paths and variables in lines 15-22
   - Edit "user" or adjust path in line 24

2. Save and give permissions, if needed (```chmod 755 evaluate_BQSR.sh```)
3. Run ```sbatch evaluate_BQSR.sh```

This only looks at the first 8 samples so that you can inspect whether BQSR had any effect. With the small number of loci we are using, this step usually doesnt do much but won't hurt downstream processes and is recommended by GATK Best Practices when possible. For organisms without verified variants, you might not be able to use BQSR at all.

### Call variants with Haplotype Caller 
The following setion will call each of your recalibrate .bams separately using your reference genome. It will then combine them into a final VCF.

#### Call samples individually
1. Download [run_GATK.sh](utility_files/run_GATK.sh) and edit the following:

    ```nano run_GATK.sh```
   - Edit email in header
   - Edit paths and variables in lines 15-19
   - Optional: edit ```--cpus-per-task``` and ```--mem``` in header, ```THREADS``` (line 19), ```-Xmx##G``` (line 36), and ```--native-pair-hmm-threads``` (line 44)

The optional edits control how many threads and how much memory the job is requesting. I currently have it set a the max that I would do comfortably on the UI servers. These setting allowed ~1000 samples to run in ~1.5 hours. IMPORTANT: if several nodes are being used on the server (check with ```squeue``` and ```sinfo -s```) then you may want to lower your requests. The tradeoff is that your job will run slower, but won't be stuck in a queue waiting for resources.
 
2. Save and give permissions, if needed (```chmod 755 run_GATK.sh```)
3. Run ```sbatch run_GATK.sh```

*This took about 1 hr 15 min for ~1000 samples.*

#### Combine all individuals into "final" vcf
We will do more filtering in R later, so this is the final vcf produced on the server.

1. Download [final_GATK.sh](utility_files/final_GATK.sh) and edit the following:

   ```nano final_GATK.sh```
   - Edit email in header
   - Edit paths and variables
   - Edit the path in line 24. This should go to a non-existent folder (this function will make the directory) called "genomics_db." IMPORTANT: if this folder already exists, the code will fail. If you have to restart this code and the folder was created during the previous run, you should delete the old folder and start over
   - Edit line 36 to match line 24's path. Key things: 1) no quotation marks, 2) it has to start with ```gendb:///``` followed by the path, 3) make sure there are 3 "/" following the colon before the name of the first directory in your path.
   - Optional: edit ```--cpus-per-task``` and ```--mem``` in header, ```-Xmx##G``` (lines 23 & 34), ```--batch-size``` (line 25), and ```--reader-threads``` (line 30)
  
The optional edits control how many threads and how much memory the job is requesting. I currently have it set a the max that I would do comfortably on the UI servers. These setting allowed ~1000 samples to run in ~1 hour. IMPORTANT: if several nodes are being used on the server (check with ```squeue``` and ```sinfo -s```) then you may want to lower your requests. The tradeoff is that your job will run slower, but won't be stuck in a queue waiting for resources.

2. Save and give permissions, if needed (```chmod 755 final_GATK.sh```)
3. Run ```sbatch final_GATK.sh```
4. Optional (but necessary if importing to R): Unzip your final .vcf.gz with:

   - First, ```module load htslib```
   - Then, ```bgzip -c -d final.vcf.gz > unzipped.final.vcf```
     
*For ~1000 samples, this took about 1 hour to run final_GATK.sh.*

## 7. Generating panel summary stats
This code give syou some general information about how the runs did on completeness by sample, read depth by locus, and on target reads.

1. Download [generate_summary_stats.sh](utility_files/generate_summary_stats.sh) and edit the following:

    ```nano generate_summary_stats.sh```
   - Edit email in header
   - Edit paths and variables in lines 17-20 and line 23
  
2. Save and give permissions, if needed (```chmod 755 generate_summary_stats.sh```)
3. Run ```sbatch generate_summary_stats.sh```

*This should be fast. Definitely <5 min.*
