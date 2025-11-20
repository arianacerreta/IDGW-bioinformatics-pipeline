## 4. [Delomas et al. 2023](https://link.springer.com/article/10.1007/s12686-023-01301-x) method

1. Get mtype2 installed. Follow instructions on [Delomas' github](https://github.com/delomast/microTyper). For university servers, you will likely have to download it in your files and call it from there.
2. Download [microhapWrap.py](https://github.com/delomast/microhapWrap/blob/main/microhapWrap.py) from Delomas github
3. Verify you have a homebrew reference genome. This homebrew genome is a MUCH smaller .fasta style genome that we will align to. It consists of two lines per locus. The first line is denoted with >LOCUS_name. The second line is the locus sequence. Please see the .fa files in this repository for an example.

    - **Full set of neutral loci (341)**: wolfAmpRef_5.fa, with associated files -> wolfAmpRef_5.1.bt2, wolfAmpRef_5.2.bt2, wolfAmpRef_5.3.bt2, wolfAmpRef_5.4.bt2, wolfAmpRef_5.rev.1.bt2, wolfAmpRef_5.rev.2.bt2
    - **Fecal subset of neutral loci (200)**: wolfAmpRef_200.fa, with associated files-> wolfAmpRef_200.1.bt2, wolfAmpRef_200.2.bt2, wolfAmpRef_200.3.bt2, wolfAmpRef_200.4.bt2, wolfAmpRef_200.rev.1.bt2, wolfAmpRef_200.rev.2.bt2
    - If you need to generate the associated files for your unique .fa, use this code (at UI use a standalone server):
  
      ```module load bowtie2```
      
      ```bowtie2-build wolfAmpRef_5.fa wolfAmpRef_5```

   - The second argument is the prefix that will be used to make the .bt2 and rev.2.bt2 files.

4. You will need to know the location of these additional files:

   - Clu341_pos_4.txt or Clu200_pos_4.txt
   - Clu341_pres_abs.txt or Clu200_pres_abs.txt
   - primers_with_adapters.txt

5. Download [run_microhapWrap.sh](utility_files/run_microhapWrap.sh) and place in your scripts directory.
6. Edit ```run_microhapWrap.sh```

    ```nano run_microhapWrap.sh```

    - Edit header for IIDS server:
    
        - #SBATCH -p eight
        - #SBATCH -C "ceph"
        - #SBATCH --cpus-per-task=1
        - #SBATCH -J clu_bwa_db
        - #SBATCH --mail-user `email@emailexample.edu`
        - #SBATCH --mail-type=BEGIN,END,FAIL
        
    - Edit lines _ to _ with your paths
    - Edit bowtie alignment settings for fecal or tissue
        - Fecal:  
        
      ```
      bowtie2 --local --no-unal -N 1 --rdg 5,5 --rfg 5,5 \
          -x "$index_prefix" -U "$fq" 2> "$bams_output_dir/${sample}_bowtie2.log" | \
          samtools view -F 276 -u | samtools sort -o "$bam"
      ```

        - Tissue (same as in Delomas python code):

      ```
      bowtie2 --end-to-end --no-unal -N 1 --rdg 0,5 --rfg 0,5 --score-min L,0,-.76 \
          -x "$index_prefix" -U $fq 2> "$bams_output_dir/${sample}_bowtie2.log" | \
          samtools view -F 276 -u | samtools sort -o "$bam"
      ```

7. Save and give microhapWrap.py and run_microhapWrap.sh permissions (e.g., ```chmod 755 microhapWrap.py```)
8. If at UI, make sure you are back on the main cluster. Run ```sbatch microhapWrap.sh``` in the console.
