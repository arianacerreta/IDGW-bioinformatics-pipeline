## 6 Old Genotyping Pipeline using bcftools
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
