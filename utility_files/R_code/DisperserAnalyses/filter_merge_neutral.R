#Written by Ariana Cerreta
#Updated: 29-June-2026

#code to filter and then merge both panels
#Panel 1: IDFGdata - contains 385 variants of putative neutral loci (dog +data needed to go with GWAdapt)
#Panel 2: NeutralTissue - contains 385 variants of putative neutral loci (wovles to go wiht GWAdapt)

library(vcfR)
library(tidyverse)
library(ggplot2)
library(SNPfiltR)

setwd("~/WOLVES/Analyses/GenomicValidation/")
#Neutral Tissue
vcf_neutraltissue<- read.vcfR("./outputs/CanFam3.1/filtered-by-neg_NeutralTissue_CanFam3.1_GATK.vcf") #vcf that has been filtered by PCR negatives; no longer includes PCR negatives
vcf_IDFG<- read.vcfR("./outputs/CanFam3.1/filtered-by-neg_IDFGdata_CanFam3.1_GATK.vcf") #vcf that has been filtered by PCR negatives; no longer includes PCR negatives

#filtering parameters (same filters as GWAdapt for consistency)
GQ<-30 #genotype quality of 30
DP<-as.numeric(20) #filter out depths less than 20
#keep all individuals and SNPs

hard_filter_rgq<-function (vcfR, depth = NULL, gq = NULL) 
{
  if (!inherits(vcfR, what = "vcfR")) {
    stop("specified vcfR object must be of class 'vcfR'")
  }
  if (!is.null(depth)) {
    if (is.numeric(depth) != "TRUE") {
      stop("specified depth cutoff must be numeric")
    }
    dp.matrix <- vcfR::extract.gt(vcfR, element = "DP", as.numeric = TRUE)
    if (sum(!is.na(dp.matrix)) < 0.5) {
      stop("genotype depth not specified in input vcf")
    }
    i <- round((sum(dp.matrix < depth, na.rm = TRUE)/sum(!is.na(dp.matrix))) * 
                 100, 2)
    message(i, "% of genotypes fall below a read depth of ", 
            depth, " and were converted to NA")
    dp.matrix[dp.matrix < depth] <- NA
    vcfR@gt[, -1][is.na(dp.matrix) == TRUE] <- NA
  }
  else {
    message("no depth cutoff provided, exploratory visualization will be generated.")
    dp.matrix <- vcfR::extract.gt(vcfR, element = "DP", as.numeric = TRUE)
    if (sum(!is.na(dp.matrix)) < 0.5) {
      message("genotype depth not specified in input vcf")
    }
    else {
      graphics::hist(dp.matrix, xlab = "genotype depth")
      graphics::abline(v = mean(dp.matrix, na.rm = TRUE), 
                       col = "red", lty = "dashed")
      graphics::hist(dp.matrix[dp.matrix < 25], xlab = "genotype depth")
    }
  }
  if (!is.null(gq)) {
    if (is.numeric(gq) != "TRUE") {
      stop("specified genotype quality cutoff must be numeric")
    }
    gq.matrix <- vcfR::extract.gt(vcfR, element = "GQ", as.numeric = TRUE)
    rgq.matrix<- vcfR::extract.gt(vcfR, element = "RGQ", as.numeric = TRUE)
    if (all(dim(gq.matrix) == dim(rgq.matrix))) {
      if (all(rownames(gq.matrix) == rownames(rgq.matrix))){
        if (all(colnames(gq.matrix) == colnames(rgq.matrix))){
          gq.matrix[is.na(gq.matrix)] <- rgq.matrix[is.na(gq.matrix)]
        }
      }
    }
    if (sum(!is.na(gq.matrix)) < 0.5) {
      stop("genotype quality not specified in input vcf")
    }
    j <- round((sum(gq.matrix < gq, na.rm = TRUE)/sum(!is.na(gq.matrix))) * 
                 100, 2)
    message(j, "% of genotypes fall below a genotype quality of ", 
            gq, " and were converted to NA")
    gq.matrix[gq.matrix < gq] <- NA
    vcfR@gt[, -1][is.na(gq.matrix) == TRUE] <- NA
  }
  else {
    message("no genotype quality cutoff provided, exploratory visualization will be generated.")
    gq.matrix <- vcfR::extract.gt(vcfR, element = "GQ", as.numeric = TRUE)
    if (sum(!is.na(gq.matrix)) < 0.5) {
      message("genotype quality not specified in input vcf")
    }
    else {
      graphics::hist(gq.matrix, xlab = "genotype quality")
      graphics::abline(v = mean(gq.matrix, na.rm = TRUE), 
                       col = "red", lty = "dashed")
    }
  }
  return(vcfR)
}

vcfR_neutraltissue <- hard_filter(vcfR = vcf_neutraltissue, depth = DP)%>%
  hard_filter_rgq(gq = GQ) 

vcfR_IDFG<-hard_filter(vcfR = vcf_IDFG, depth = DP)%>%
  hard_filter_rgq(gq = GQ)


gt_neutraltissue <- extract.gt(vcfR_neutraltissue, element = "GT", as.numeric = FALSE) #this no longer has PCR negatives
gt_df_neutral <- as.data.frame(gt_neutraltissue)

gt_IDFG <- extract.gt(vcfR_IDFG, element = "GT", as.numeric = FALSE) #this no longer has PCR negatives
gt_df_IDFG <- as.data.frame(gt_IDFG)

#fix column names
colnames(gt_df_neutral) <- str_extract(colnames(gt_df_neutral), "i[0-9A-Za-z_-]+(?=_RG_sort\\.recal)")
colnames(gt_df_IDFG) <- str_extract(colnames(gt_df_IDFG), "i[0-9A-Za-z_-]+(?=_RG_sort\\.recal)")

#append CHROM, POS
gt_df_neutral$CHROM <- getCHROM(vcfR_neutraltissue)
gt_df_neutral$POS <- as.numeric(getPOS(vcfR_neutraltissue))
gt_df_IDFG$CHROM <- getCHROM(vcfR_IDFG)
gt_df_IDFG$POS <- as.numeric(getPOS(vcfR_IDFG))

#append panel identifier
gt_df_neutral$PANEL <- rep("neutral", length(gt_df_neutral$CHROM))
gt_df_IDFG$PANEL <- rep("IDFG", length(gt_df_IDFG$CHROM))

#move CHROM and POS, and PANEL to front
gt_df_neutral <- gt_df_neutral %>% relocate(CHROM, POS, PANEL)
gt_df_IDFG <- gt_df_IDFG %>% relocate(CHROM, POS, PANEL)

#pivot to longform
gt_df_neutral <- pivot_longer(
  gt_df_neutral,
  cols = -c(CHROM, POS, PANEL),
  names_to = "SampleID",
  values_to = "GT"
)

gt_df_IDFG <- pivot_longer(
  gt_df_IDFG,
  cols = -c(CHROM, POS, PANEL),
  names_to = "SampleID",
  values_to = "GT"
)

gt_df_combined <- bind_rows(gt_df_neutral, gt_df_IDFG)

#extract depth
dp_neutral <- extract.gt(vcfR_neutraltissue, element = "DP", as.numeric = TRUE)
dp_neutral <- as.data.frame(dp_neutral)
dp_IDFG <- extract.gt(vcfR_IDFG, element = "DP", as.numeric = TRUE)
dp_IDFG <- as.data.frame(dp_IDFG)
colnames(dp_neutral) <- str_extract(colnames(dp_neutral), "i[0-9A-Za-z_-]+(?=_RG_sort\\.recal)")
colnames(dp_IDFG) <- str_extract(colnames(dp_IDFG), "i[0-9A-Za-z_-]+(?=_RG_sort\\.recal)")

# Calculate average read depth per sample
mean_depth_neutral <- as.data.frame(dp_neutral) %>%
  pivot_longer(cols = everything(), names_to = "SampleID", values_to = "Depth") %>%
  group_by(SampleID) %>%
  summarize(mean_depth = mean(Depth, na.rm = TRUE), .groups = "drop")

mean_depth_IDFG <- as.data.frame(dp_IDFG) %>%
  pivot_longer(cols = everything(), names_to = "SampleID", values_to = "Depth") %>%
  group_by(SampleID) %>%
  summarize(mean_depth = mean(Depth, na.rm = TRUE), .groups = "drop")


mean_depth_combined <- bind_rows(mean_depth_IDFG, mean_depth_neutral)

#define replicate pairs
sample_names_combined <- unique(gt_df_combined$SampleID)

#for IDFG

samples_df_combined <- tibble(SampleName = sample_names_combined) %>%
  mutate(
    # Extract UI_ID if present
    UI_ID = str_extract(SampleName, "UI\\d+"),
    
    # Everything after _T#_
    tail = str_extract(SampleName, "_T\\d+_.*") %>%
      str_replace("^_T\\d+_", ""),
    
    # If UI_ID is present, take everything after it; otherwise take full tail
    ReplicateID = if_else(!is.na(UI_ID),
                          str_replace(tail, paste0("^", UI_ID, "_"), ""),
                          tail)
  ) %>%
  select(-tail)%>%
  filter(!str_detect(SampleName, "CluVARI")) #remove dogs

rep_by_id <- samples_df_combined %>%
  group_by(ReplicateID) %>%
  filter(n() > 1) %>%  # only keep true replicates
  summarise(Pairs = list(as_tibble(
    t(combn(SampleName, 2)), .name_repair = ~ c("Sample1", "Sample2")
  )), .groups = "drop") %>%
  unnest(Pairs) %>%
  mutate(RepType = "Technical")

rep_by_ui_neutral <-samples_df_combined %>%
  group_by(UI_ID) %>%
  filter(n_distinct(ReplicateID) > 1, n() > 1) %>%  # more than one RepID = bio replicates
  summarise(Pairs = list(as_tibble(
    t(combn(SampleName, 2)), .name_repair = ~ c("Sample1", "Sample2")
  )), .groups = "drop") %>%
  unnest(Pairs) %>%
  mutate(RepType = "Biological")

replicate_map <- bind_rows(rep_by_id, rep_by_ui_neutral) %>%
  distinct(Sample1, Sample2, .keep_all = TRUE)

completeness_df_combined <- gt_df_combined %>%
  filter(!is.na(GT), GT != "./.") %>%
  group_by(SampleID) %>%
  summarize(non_missing_genotypes = n(), .groups = "drop") %>%
  left_join(mean_depth_combined, by = "SampleID")

rep_long_combined<- replicate_map %>%
  pivot_longer(cols = c(Sample1, Sample2), 
               names_to = "SampleCol", values_to = "SampleID") %>%
  left_join(completeness_df_combined, by = "SampleID")%>%
  mutate(
    MatchID=coalesce(ReplicateID,UI_ID)
  )

#use completeness and depth to choose the replicate to keep
# For each replicate group, pick highest completeness, then tie-break with highest mean depth
to_keep_combined <- rep_long_combined %>%
  group_by(MatchID) %>%
  arrange(desc(non_missing_genotypes), desc(mean_depth)) %>%
  slice(1) %>%
  pull(SampleID) 

#filter out the dropped replicates & dogs
gt_df_combined_noreps<-gt_df_combined %>%
  filter(!(SampleID %in% unlist(replicate_map[, c("Sample1", "Sample2")])) | 
           SampleID %in% to_keep_combined) %>%
  filter(!str_detect(SampleID, "CluVARI")) #remove dogs

#merging
# Add MatchID to neutral
gt_df_combined_noreps <- left_join(gt_df_combined_noreps,
                                  samples_df_combined %>% select(SampleName, ReplicateID, UI_ID),
                                  by = c("SampleID" = "SampleName"))

#which samples have both panel data
rep_ids_unqiue <- unique(gt_df_combined_noreps$UI_ID)
rep_ids_all <- unique(gt_df_combined_noreps$SampleID)

length(rep_ids_all) == length(rep_ids_unqiue)

length(unique(gt_df_combined_noreps$SampleID)) #1048
##remove UI1458 (same as UI0989 which has higher completeness)
##remove UI1459 (same as UI2455 which has higher completeness)
##remove UI1755 (same as UI1013 which has higher completeness)
##remove UI1401 controdicting information/genotype with usats
##remove UI2107 which is a dog

gt_df_final<-gt_df_combined_noreps %>%
  filter(!str_detect(SampleID, "UI1458")) %>%
  filter(!str_detect(SampleID, "UI1459")) %>%
  filter(!str_detect(SampleID, "UI1755")) %>%
  filter(!str_detect(SampleID, "UI1401")) %>%
  filter(!str_detect(SampleID, "UI2107")) 

length(unique(gt_df_final$SampleID)) #1043

#calculate MAF
#unphase gt_df_final
gt_df_final <- gt_df_final %>%
  mutate(
    GT_unphased = map_chr(GT, ~ {
      if (is.na(.x)) return(NA_character_)
      alleles <- str_split(.x, "[/|]")[[1]]
      paste(sort(alleles), collapse = "/")
    })
  )

#identify multiallelic loci and filter out
multiallelic_loci <- gt_df_final %>%
  filter(!is.na(GT_unphased)) %>%
  separate(GT_unphased, into = c("a1", "a2"), sep = "/", remove = FALSE) %>%
  mutate(a1 = as.integer(a1),
         a2 = as.integer(a2)) %>%
  group_by(CHROM, POS) %>%
  summarise(
    max_allele = max(c(a1, a2), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_allele > 1)

gt_df_biallelic <- gt_df_final %>%
  anti_join(multiallelic_loci, by = c("CHROM", "POS"))

#dummy amplicon name
gt_df_final2 <- gt_df_biallelic %>%
  arrange(CHROM, POS) %>%
  group_by(CHROM) %>%
  mutate(
    AmpliconIndex = cumsum(
      c(TRUE, diff(POS) > 100)
    )
  ) %>%
  ungroup()

gt_df_final2 <- gt_df_final2 %>%
  mutate(
    AmpliconID = paste(CHROM, AmpliconIndex, sep = "_")
  )

amplicon_lookup<- gt_df_final2 %>%
  select(CHROM, POS, AmpliconID) %>%
  distinct

#calculate maf to choose one SNP per amplicon
maf_df <- gt_df_final2 %>%
  filter(!is.na(GT_unphased)) %>%
  mutate(
    alt_alleles = case_when(
      GT_unphased == "0/0" ~ 0L,
      GT_unphased == "0/1" ~ 1L,
      GT_unphased == "1/1" ~ 2L,
      TRUE ~ NA_integer_
    )
  )

maf_df <- maf_df %>%
  group_by(AmpliconID, CHROM, POS) %>%
  summarise(
    n_samples = n_distinct(ReplicateID),
    n_called  = sum(!is.na(alt_alleles)),
    alt_count = sum(alt_alleles, na.rm = TRUE),
    total_alleles = 2 * n_called,
    alt_freq = alt_count / total_alleles,
    maf = pmin(alt_freq, 1 - alt_freq),
    .groups = "drop"
  )

maf_df <- maf_df %>%
  mutate(
    missing_prop = 1 - (n_called / n_samples)
  )

#choose best SNP
best_snp_per_amplicon <- maf_df %>%
  arrange(
    AmpliconID,
    desc(maf),
    missing_prop,
    desc(n_called)
  ) %>%
  group_by(AmpliconID) %>%
  slice(1) %>%
  ungroup()

best_snp_per_amplicon_no_zero <- best_snp_per_amplicon %>%
  filter(maf != 0) 

#prune original gt_df_final
gt_df_pruned <- gt_df_final2 %>%
  inner_join(
    best_snp_per_amplicon_no_zero %>%
      select(CHROM, POS),
    by = c("CHROM", "POS")
  )

# recalculate maf
maf_df_recalc <- gt_df_pruned %>%
  filter(!is.na(GT_unphased)) %>%
  mutate(
    alt_alleles = case_when(
      GT_unphased == "0/0" ~ 0L,
      GT_unphased == "0/1" ~ 1L,
      GT_unphased == "1/1" ~ 2L,
      TRUE ~ NA_integer_
    )
  )

maf_df_recalc <- maf_df_recalc %>%
  group_by(AmpliconID, CHROM, POS) %>%
  summarise(
    n_samples = n(),
    n_called  = sum(!is.na(alt_alleles)),
    alt_count = sum(alt_alleles, na.rm = TRUE),
    total_alleles = 2 * n_called,
    alt_freq = alt_count / total_alleles,
    maf = pmin(alt_freq, 1 - alt_freq),
    .groups = "drop"
  )

onlySNPs_maf_0.05<-maf_df_recalc %>%
  filter(maf >= 0.05) 

#filter final set of SNPs maf <0.05
gt_df_ready_for_merge <- gt_df_pruned %>%
  inner_join(
    onlySNPs_maf_0.05 %>%
      select(CHROM, POS),
    by = c("CHROM", "POS")
  )

saveRDS(gt_df_ready_for_merge, file = "./outputs/gt_df_ready_for_merge_neutral_CanFam3_1_GATK.rds")

#how to read back in when ready
gt_df_ready_neutral <- readRDS("./outputs/gt_df_ready_for_merge_neutral_CanFam3_1_GATK.rds")
#write.table(gt_df_combined_noreps, file = "./outputs/CanFam3.1/merged_neutral_IDFG_gt_CanFam3_1_GATK.txt", sep = "\t",
 #           row.names = FALSE)
