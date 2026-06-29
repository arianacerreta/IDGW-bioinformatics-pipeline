# IDGW-bioinformatics-pipeline
This is a repository for how to process raw genomic data from the Idaho Gray Wolf Project.

We have used several genotyping pipelines over the course of this project:
1. GATK
2. DELOMAS et al. 2023
3. BCFTOOLS

We no longer use the BCFTOOLS pipeline, but have left the instructions here for reference. We found that this pipeline works well for SNPs, but incorrectly genotyped short indels (<=3bp in length). If your panel includes indels (e.g., UI wolf adaptive panel; *Cerreta et al. in prep*), we recommend using the GATK pipeline.

If you are using the IDFG neutral wolf panel, we recommend the DELOMAS et al. 2023 pipeline. They have an entire GitHub repository for their method, but we outline how to use this pipeline with more explicit UI server specific instructions.

Follow the links to the appropriate steps below for each pipeline:

## GATK Pipeline

[Steps 1-3](STEPS1-3-ALL-PIPELINES.md) -> [Steps 4-5](STEP4-5-BCFTOOLS-and-GATK-PIPELINES.md) -> [Steps 6-7](STEP6-GATK-PIPELINE.md)

## DELOMAS et al. 2023 Pipeline

[Steps 1-3](STEPS1-3-ALL-PIPELINES.md) -> [Step 4](STEP4-DELOMAS-PIPELINE.md)

## BCFtools Pipeline

[Steps 1-3](STEPS1-3-ALL-PIPELINES.md) -> [Steps 4-5](STEP4-5-BCFTOOLS-and-GATK-PIPELINES.md) -> [Step 6](STEP6-BCFTOOLS-PIPELINE.md)

## [Common troubleshooting](HELP-JOBS-FAILING.md) if your jobs keep failing.

# R code
## Wolf Adaptive Panel
### Filtering and summaries

.vcf file from Step 6 of GATK Pipeline -> [Filtering by PCR negative (NTC)](utility_files/R_code/adaptive_methods/final_adapt_filter_with_PCRneg.R) -> [Replicate Error Calculations](utility_files/R_code/adaptive_methods/final_error_rates_GWAdapt.R) -> [Apply Filtering and Remove Replicates](utility_files/R_code/adaptive_methods/final_filter_remove-reps_GWAdapt.R) -> [Calculate Allele Frequencies](utility_files/R_code/adaptive_methods/final_freq_calc_GWAdapt.R)


To calculate summary stats from the adaptive panel:
[Summary stats](utility_files/R_code/adaptive_methods/summ_stats.R)

To evaluate k-locus concordance:
[K-locus code](utility_files/R_code/adaptive_methods/final_k-locus.R)

## SNP a Scat
### Filtering and matching

Use the microhap_genotypes.csv files output (Zenodo DOI: 10.5281/zenodo.20799436) from Step 4 of DELOMAS et al. 2023 Pipeline -> [Filtering and matching](utility_files/R_code/SNP-a-scat_R_code/final_matching_tissue_fecal.R)
This code should produce: 1) Figure S.4 (mismatch rates), 2) Figure 2.6, 3) data on heterozygosity, missingness, samples that overlap with microsatellites, matching, fecal capture

### Summary Results

Use zipped summary inputs (Zenodo DOI: 10.5281/zenodo.20799436) with [final_scatoptsummaries.R](utility_files/R_code/SNP-a-scat_R_code/final_scatoptsummaries.R)
This code should produce: 1) data in Tables 2-4; 2) data in Table S.1

Use qPCR_results.csv (Zenodo DOI: 10.5281/zenodo.20799436) and gen_data_fecal_2026-06-24.csv to summarize results based on qPCR concentration bins with [final_qpcr_bin_summary.R](utility_files/R_code/SNP-a-scat_R_code/final_qpcr_bin_summary.R)
This code will produce data from Table 5

### Simulation of PID and PIDsibs

Use populations.haps.vcf (Zenodo DOI: 10.5281/zenodo.20799436) and [simulation_200loci.R](utility_files/R_code/SNP-a-scat_R_code/simulation_200loci.R)
This code should produce Figure S.1

Use 70loci_freq_ScatTest2.csv (Zenodo DOI: 10.5281/zenodo.20799436) and [simulation_50loci.R](utility_files/R_code/SNP-a-scat_R_code/simulation_50loci.R)
This code should produce Figure S.2

### Determine Thresholds for Calling Sex Locus

Use pres_abs_genotypes.txt and ScatOpt4KnownSex.csv (Zenodo DOI: 10.5281/zenodo.20799436) and [final_sex_threshold.R](utility_files/R_code/SNP-a-scat_R_code/final_sex_threshold.R)
This should produce Figure S.3

## Dispersers

Use .vcf produced by [Filtering by PCR negative (NTC)](utility_files/R_code/adaptive_methods/final_adapt_filter_with_PCRneg.R)
