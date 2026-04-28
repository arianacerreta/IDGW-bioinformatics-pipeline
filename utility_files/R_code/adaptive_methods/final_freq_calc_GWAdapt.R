##### Final Frequency calculations for Adaptive Panel###

#libraries
library(tidyverse) #version 2.0.0: dplyr 1.2.1, forcats 1.0.1, lubridate 1.9.5, purrr 1.2.2, readr 2.2.0, stringr 1.6.0, tibble 3.3.1, tidyr 1.3.2
library(vcfR)

#inputs from .txt produced in final_filter_remove-reps_GWAdapt.R
genotypes<-read.delim2("./outputs/CanFam3.1/adaptive_canfam3_1_gt_GATK_4-28-26.txt")
#pull in minimally filtered vcf to pull REF and ALT
vcf<- read.vcfR("./outputs/filtered-by-neg_GWAdapt_CanFam3.1_GATK.vcf") #vcf without a DP filter by PCR negatives; no longer includes PCR negatives

#summary stats for paper
summary_genotype_completeness <- genotypes %>%
  group_by(SampleID) %>%
  summarise(
    total_loci = n(),
    genotyped_loci = sum(!is.na(GT) & GT != "./."),
    completeness = genotyped_loci / total_loci * 100
  ) 

genotype_completeness<- summary_genotype_completeness%>%
  summarise(avg_completeness = mean(completeness))

print(genotype_completeness)

#calculate genotype frequency
genotypes_freq <- genotypes %>%
  mutate(Locus = paste(CHROM, POS, sep = "_"))

###are there multiple alleles for any given locus
biallelic_multiallelic<- genotypes_freq %>%
  filter(!is.na(GT)) %>%
  # split genotype into two alleles
  separate(GT, into = c("A1", "A2"), sep = "[/|]", remove = FALSE) %>%
  # pivot to long format (one row per allele)
  pivot_longer(cols = c(A1, A2), values_to = "Allele") %>%
  # count unique alleles per locus
  group_by(Locus) %>%
  summarise(
    n_alleles = n_distinct(Allele, na.rm = TRUE),
    alleles = paste(sort(unique(Allele)), collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(
    locus_type = case_when(
      n_alleles == 1 ~ "invariant",
      n_alleles == 2 ~ "biallelic",
      n_alleles > 2 ~ "multiallelic"
    )
  )

alleles_per_locus <- genotypes_freq %>%
  filter(!is.na(GT)) %>%
  separate(GT, into = c("A1", "A2"), sep = "[/|]", remove = FALSE) %>%
  pivot_longer(cols = c(A1, A2), values_to = "Allele") %>%
  group_by(CHROM, POS, Locus) %>%
  summarise(
    allele_indices = list(sort(unique(as.integer(Allele)))),
    n_alleles = n_distinct(Allele),
    .groups = "drop"
  )

vcf_fix <- as.data.frame(vcf@fix) %>%
  transmute(
    CHROM,
    POS = as.numeric(POS),
    REF,
    ALT
  )

alleles_mapped <- alleles_per_locus %>%
  left_join(vcf_fix, by = c("CHROM", "POS")) %>%
  rowwise() %>%
  mutate(
    alt_vec = list(str_split(ALT, ",")[[1]]),
    all_alleles = list(c(REF, alt_vec)),
    kept_alleles = list(all_alleles[allele_indices + 1]),
    kept_alleles_str = paste(unlist(kept_alleles), collapse = ",")
  ) %>%
  ungroup()

allele_map <- alleles_mapped %>%
  mutate(
    new_index = map(allele_indices, ~ seq_along(.x) - 1)
  ) %>%
  unnest(c(allele_indices, new_index, kept_alleles)) %>%
  rename(
    old = allele_indices,
    new = new_index,
    allele = kept_alleles
  )

genotypes_reindexed <- genotypes_freq %>%
  filter(!is.na(GT)) %>%
  separate(GT, into = c("A1", "A2"), sep = "[/|]", remove = FALSE) %>%
  mutate(
    A1 = as.integer(A1),
    A2 = as.integer(A2)
  ) %>%
  # join mapping for each allele separately
  left_join(allele_map, by = c("CHROM", "POS", "Locus", "A1" = "old")) %>%
  rename(A1_new = new) %>%
  left_join(allele_map, by = c("CHROM", "POS", "Locus", "A2" = "old")) %>%
  rename(A2_new = new) %>%
  mutate(
    GT_reindexed = paste0(A1_new, "/", A2_new)
  )

allele_freqs <- genotypes_reindexed %>%
  pivot_longer(cols = c(A1_new, A2_new), values_to = "Allele") %>%
  group_by(Locus, Allele) %>%
  summarise(
    count = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    freq = count / sum(count)
  ) %>%
  ungroup()

allele_freqs <- allele_freqs %>%
  left_join(
    allele_map %>% select(Locus, new, allele),
    by = c("Locus", "Allele" = "new")
  )

geno_summary <- genotypes_reindexed %>%
  filter(!is.na(A1_new) & !is.na(A2_new)) %>%
  mutate(
    is_hom = A1_new == A2_new,
    is_het = A1_new != A2_new
  ) %>%
  group_by(Locus) %>%
  summarise(
    N = n(),
    n_hom = sum(is_hom),
    n_het = sum(is_het),
    prop_hom = n_hom / N,
    prop_het = n_het / N,
    .groups = "drop"
  )

allele_counts <- genotypes_reindexed %>%
  filter(!is.na(A1_new) & !is.na(A2_new)) %>%
  pivot_longer(cols = c(A1_new, A2_new), values_to = "allele") %>%
  group_by(Locus, allele) %>%
  summarise(
    MAC = n(),              # each row = one allele copy
    n_indiv = n_distinct(SampleID),
    .groups = "drop"
  )

allele_labels <- allele_map %>%
  group_by(Locus) %>%
  summarise(
    REF = allele[new == 0][1],
    ALT = paste(allele[new != 0], collapse = ","),
    n_alleles = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ALT = ifelse(ALT == "", ".", ALT)
  )

allele_flags <- genotypes_reindexed %>%
  filter(!is.na(A1_new) & !is.na(A2_new)) %>%
  mutate(
    is_hom = A1_new == A2_new
  ) %>%
  pivot_longer(cols = c(A1_new, A2_new), values_to = "allele") %>%
  group_by(Locus, allele) %>%
  summarise(
    MAC = n(),
    n_indiv = n_distinct(SampleID),
    n_hom_obs = sum(is_hom),
    n_het_obs = sum(!is_hom),
    .groups = "drop"
  ) %>%
  mutate(
    flag_MAC = MAC < 3,
    flag_single_individual = n_indiv == 1,
    flag_only_homozygous = n_het_obs == 0 & MAC > 0
  )

allele_flags <- allele_flags %>%
  mutate(
    flag_reason = case_when(
      flag_MAC & flag_single_individual & flag_only_homozygous ~
        "MAC_lt_3;single_individual;only_homozygous",
      
      flag_MAC & flag_single_individual ~
        "MAC_lt_3;single_individual",
      
      flag_MAC & flag_only_homozygous ~
        "MAC_lt_3;only_homozygous",
      
      flag_single_individual & flag_only_homozygous ~
        "single_individual;only_homozygous",
      
      flag_MAC ~ "MAC_lt_3",
      flag_single_individual ~ "single_individual",
      flag_only_homozygous ~ "only_homozygous",
      TRUE ~ NA_character_
    )
  )


final_table<- allele_freqs %>%
  left_join(allele_flags, by = c("Locus", "Allele" = "allele"))

write.csv(final_table, "./outputs/CanFam3.1/allele_freq_table_CanFam3_1_GATK_final_4-28-26.csv")

gt_summary <- genotypes_reindexed %>%
  group_by(Locus) %>%
  summarise(
    N = n(),
    
    prop_0_0 = mean(GT_reindexed == "0/0"),
    prop_0_1 = mean(GT_reindexed == "0/1"),
    prop_1_1 = mean(GT_reindexed == "1/1"),
    prop_0_2 = mean(GT_reindexed == "0/2"),
    prop_2_2 = mean(GT_reindexed == "2/2"),
    prop_3_3 = mean(GT_reindexed == "3/3"),
    
    .groups = "drop"
  )

ref_alt <- allele_map %>%
  group_by(Locus) %>%
  summarise(
    REF = allele[new == 0][1],
    ALT = paste(allele[new != 0], collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(
    ALT = ifelse(is.na(ALT) | ALT == "", ".", ALT)
  )

ref_freq <- allele_freqs %>%
  filter(Allele == 0) %>%
  transmute(
    Locus,
    REF_freq = freq
  )

final_locus_table <- gt_summary %>%
  left_join(ref_alt, by = "Locus") %>%
  left_join(ref_freq, by = "Locus") %>%
  mutate(
    ALT = ifelse(is.na(ALT), ".", ALT)
  )%>%
  select(Locus, REF, ALT, N, REF_freq, everything())

#of note: if a locus is invariant, the REF allele will reflect the only allele present even if it is not the REF allele on CanFam3.1; check the .vcf for the true REF allele

write.csv(final_locus_table, "./outputs/CanFam3.1/locus_table_CanFam3_1_GATK_final_4-28-26.csv")
