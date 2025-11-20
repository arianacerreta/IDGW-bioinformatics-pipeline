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



