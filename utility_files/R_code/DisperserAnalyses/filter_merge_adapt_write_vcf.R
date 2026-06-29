#Written by Ariana Cerreta
#Updated: 29-June-2026
##Updates: new file which incorporates the 64 samples from AdaptOptim

#code to filter GWAdapt and then merge and save as vcf 
#Panel 1: GWAdapt - contains 282 amplicons of putative adaptive loci (330 unfiltered variants)
#Panel 2: Neutral - contains 200 amplicons of putative neutral loci (178 pruned variants; maf>0.05; 1 SNP per amplicon)

library(vcfR)
library(tidyverse)
library(ggplot2)
library(SNPfiltR)

vcf_adapt<- read.vcfR("../AdaptiveAnalyses/outputs/filtered-by-neg_GWAdapt_CanFam3.1_GATK.vcf") #vcf that has been filtered by PCR negatives; no longer includes PCR negatives

#read in metadata for AdaptOptim for UI_IDs
adaptoptim_meta<-read.csv("../AdaptiveAnalyses/inputs/AdaptOptim_WolfMeta_Age_Disp_Breed_Sex_Coat.csv")

#filtering parameters
GQ<-30 #genotype quality of 30
DP<-as.numeric(20) #filter out depths less than 20

#no SNP or sample filtering yet
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

vcfR_adapt<-hard_filter(vcfR = vcf_adapt, depth = DP)%>%
  hard_filter_rgq(gq = GQ) 

gt_adapt <- extract.gt(vcfR_adapt, element = "GT", as.numeric = FALSE) #this no longer has PCR negatives
gt_df_adapt <- as.data.frame(gt_adapt)


#fix column names
#colnames(gt_df_adapt) <- str_extract(colnames(gt_df_adapt), "i[0-9A-Za-z_-]+(?=_RG_sort\\.recal)")

#append CHROM, POS
gt_df_adapt$CHROM <- getCHROM(vcfR_adapt)
gt_df_adapt$POS <- as.numeric(getPOS(vcfR_adapt))

#append panel identifier
gt_df_adapt$PANEL <- rep("adapt", length(gt_df_adapt$CHROM))

#move CHROM and POS, and PANEL to front
gt_df_adapt <- gt_df_adapt %>% relocate(CHROM, POS, PANEL)

#pivot to longform
gt_df_adapt <- pivot_longer(
  gt_df_adapt,
  cols = -c(CHROM, POS, PANEL),
  names_to = "SampleID",
  values_to = "GT"
)

#extract depth
dp_adapt <- extract.gt(vcfR_adapt, element = "DP", as.numeric = TRUE)
dp_adapt <- as.data.frame(dp_adapt)

#colnames(dp_adapt) <- str_extract(colnames(dp_adapt), "i[0-9A-Za-z_-]+(?=_RG_sort\\.recal)")

mean_depth_adapt <- as.data.frame(dp_adapt) %>%
  pivot_longer(cols = everything(), names_to = "SampleID", values_to = "Depth") %>%
  group_by(SampleID) %>%
  summarize(mean_depth = mean(Depth, na.rm = TRUE), .groups = "drop")

#define replicate pairs
sample_names_adapt <- unique(gt_df_adapt$SampleID)

#extract data to make a pop map
popmap_long<-data.frame(id=colnames(vcfR_adapt@gt)[2:length(colnames(vcfR_adapt@gt))]) 

popmap<- popmap_long %>%
  mutate(i7 = str_split_i(id, "_", 1),
         well = str_split_i(id, "_", 2),
         platename = str_split_i(id, "_", 3),
         UI_ID = str_split_i(id, "_", 4),
         LabID = str_split_i(id, "_", 5)) 

rep_by_id_adapt <- popmap %>%
  group_by(LabID) %>%
  filter(n() > 1) %>%  # only keep true replicates
  summarise(Pairs = list(as_tibble(
    t(combn(id, 2)), .name_repair = ~ c("Sample1", "Sample2")
  )), .groups = "drop") %>%
  unnest(Pairs) %>%
  mutate(RepType = "Technical")

rep_by_ui_adapt <- popmap %>%
  filter(UI_ID != "NA") %>%   # <- drop NA UI_ID rows
  group_by(UI_ID) %>%
  filter(n_distinct(LabID) > 1, n() > 1) %>%  # more than one RepID = bio replicates
  summarise(Pairs = list(as_tibble(
    t(combn(id, 2)), .name_repair = ~ c("Sample1", "Sample2")
  )), .groups = "drop") %>%
  unnest(Pairs) %>%
  mutate(RepType = "Biological")

replicate_map_adapt <- bind_rows(rep_by_id_adapt, rep_by_ui_adapt) %>%
  distinct(Sample1, Sample2, .keep_all = TRUE)

completeness_df_adapt <- gt_df_adapt %>%
  filter(!is.na(GT), GT != "./.") %>%
  group_by(SampleID) %>%
  summarize(non_missing_genotypes = n(), .groups = "drop")%>%
  left_join(mean_depth_adapt, by = "SampleID")

rep_long_adapt<- replicate_map_adapt %>%
  pivot_longer(cols = c(Sample1, Sample2), 
               names_to = "SampleCol", values_to = "SampleID") %>%
  left_join(completeness_df_adapt, by = "SampleID")%>%
  mutate(
    MatchID=coalesce(LabID,UI_ID)
  )%>%
  filter(!str_detect(MatchID, "ZZAT")) %>% #remove ZZAT from being in there multiple times 
  filter(!(SampleID == "i048_D06_GWAdapt6_UI1755_YZL9" & SampleCol == "Sample2")) #YZL9 in twice


#use completeness and depth to choose the replicate to keep
# For each replicate group, pick highest completeness, then tie-break with highest mean depth
to_keep_adapt <- rep_long_adapt %>%
  group_by(MatchID) %>%
  arrange(desc(non_missing_genotypes), desc(mean_depth)) %>%
  slice(1) %>%
  pull(SampleID)

#filter out the dropped replicates
gt_df_adapt_noreps<-gt_df_adapt %>%
  filter(!(SampleID %in% unlist(replicate_map_adapt[, c("Sample1", "Sample2")])) | 
           SampleID %in% to_keep_adapt)
#merging
# Add ID
gt_df_adapt_noreps <- left_join(gt_df_adapt_noreps,
                                   popmap %>% select(id, LabID, UI_ID),
                                   by = c("SampleID" = "id")) 

gt_df_adapt_noreps<-gt_df_adapt_noreps %>%
  mutate(UI_ID = na_if(UI_ID, "NA"))

#which samples have both panel data
rep_ids_unqiue <- unique(na.omit(gt_df_adapt_noreps$UI_ID))
rep_ids_all <- unique(gt_df_adapt_noreps$SampleID)

length(rep_ids_all) == length(rep_ids_unqiue) #FALSE because several GTseek samples don't have UI_IDs

length(unique(gt_df_adapt_noreps$SampleID)) #1004
##remove UI1458 (same as UI0989 which was kept on neutral)
##remove UI1459 (same as UI2455 which was kept on neutral)
##remove UI1755 (same as UI1013 which was kept on neutral)
##remove UI1401 contradicting information/genotype with usats
##remove UI2107 which is a dog

gt_df_final<-gt_df_adapt_noreps %>%
  filter(!is.na(UI_ID))%>%
  filter(!str_detect(SampleID, "UI1458")) %>%
  filter(!str_detect(SampleID, "UI1459")) %>%
  filter(!str_detect(SampleID, "UI1755")) %>%
  filter(!str_detect(SampleID, "UI1401")) %>%
  filter(!str_detect(SampleID, "UI2107")) 

length(unique(gt_df_final$SampleID)) #951

#calculate MAF and identify pet loci

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

#read in pet locus metadata
pet_loci<-read.csv("../GenomicValidation/inputs/CanFam3.1/adaptive_loci_canfam3_1.csv", header = TRUE)

#add OG locus as amplicon identifier
gt_df_final2 <- gt_df_biallelic %>%
  left_join(
    pet_loci,
    by = c("CHROM","POS")
  )

#pet_loci should be able to serve as an amplicon lookup
#calculate maf
maf_df <- gt_df_final2 %>%
  #filter(!is.na(GT_unphased)) %>%
  mutate(
    alt_alleles = case_when(
      GT_unphased == "0/0" ~ 0L,
      GT_unphased == "0/1" ~ 1L,
      GT_unphased == "1/1" ~ 2L,
      TRUE ~ NA_integer_
    )
  )

maf_df <- maf_df %>%
  group_by(OG_name, CHROM, POS, Pet_SNP) %>%
  summarise(
    n_samples = n_distinct(LabID),
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
  group_by(OG_name) %>%
  group_modify(~ {
    
    df <- .x
    
    has_pet <- any(df$Pet_SNP, na.rm = TRUE)
    n_snps  <- nrow(df)
    
    if (has_pet) {
      
      pet_df <- df %>% filter(Pet_SNP)
      
      # Case 1: PET SNP exists with maf > 0
      if (any(pet_df$maf > 0, na.rm = TRUE)) {
        pet_df %>%
          arrange(desc(maf)) %>%
          slice(1)
        
        # Case 2: PET SNP exists but maf == 0
      } else {
        
        # Only one SNP in amplicon → keep it regardless
        if (n_snps == 1) {
          df
          
          # Multiple SNPs → choose highest MAF non-PET SNP
        } else {
          df %>%
            filter(!Pet_SNP) %>%
            arrange(desc(maf)) %>%
            slice(1)
        }
      }
      
      # Case 3: No PET SNPs at all
    } else {
      df %>%
        arrange(desc(maf)) %>%
        slice(1)
    }
  }) %>%
  ungroup()


#prune original gt_df_final
gt_df_pruned <- gt_df_final2 %>%
  inner_join(
    best_snp_per_amplicon%>%
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
  group_by(OG_name, CHROM, POS) %>%
  summarise(
    n_samples = n_distinct(LabID),
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

saveRDS(gt_df_ready_for_merge, file = paste0("./outputs/gt_df_ready_for_merge_GWAdapt_CanFam3_1_GATK", Sys.Date(),".rds"))

#clean environment
gt_df_ready_neutral <- readRDS("../GenomicValidation/outputs/gt_df_ready_for_merge_neutral_CanFam3_1_GATK.rds")

gt_df_ready_adapt <- readRDS("./outputs/gt_df_ready_for_merge_GWAdapt_CanFam3_1_GATK2026-06-29.rds")

#check shit
#unique UI_IDs in adaptive
length(unique(gt_df_ready_adapt$UI_ID)) #951
length(unique(gt_df_ready_adapt$LabID)) #951

#unique UI_IDs in neutral (has both IDFG and Uidaho200 loci)
length(unique(gt_df_ready_neutral$UI_ID)) #1043
length(unique(gt_df_ready_neutral$ReplicateID)) #1043


#rename neutral Replicate ID to LabID
gt_df_ready_neutral<- gt_df_ready_neutral %>% rename(LabID = ReplicateID)

#standardize format
# NOTE: UI1164 (adapt panel) and UI0698 (neutral panel) are the same biological individual (B517)
# Verified by matching microsatellite genotype; harmonized to UI0698 across panels
gt_df_adapt_stnd<-gt_df_ready_adapt %>%
  select(CHROM, POS, PANEL, SampleID, LabID, UI_ID, GT_unphased, OG_name)%>%
  rename(AmpliconID = OG_name)%>%
  mutate(
    UI_ID = if_else(
      PANEL == "adapt" & UI_ID == "UI1164",
      "UI0698",
      UI_ID
    )
  )

gt_df_neutral_stnd<-gt_df_ready_neutral %>%
  select(CHROM, POS, PANEL, SampleID, LabID, UI_ID, GT_unphased, AmpliconID)

#combine
gt_df_combined <- bind_rows(gt_df_neutral_stnd, gt_df_adapt_stnd)

#which samples have both panel data
rep_ids_neutral <- unique(gt_df_neutral_stnd$LabID) #1043
rep_ids_adapt <- unique(gt_df_adapt_stnd$LabID) #951

ui_ids_neutral <- unique(gt_df_neutral_stnd$UI_ID) #1043
ui_ids_adapt <- unique(gt_df_adapt_stnd$UI_ID) #951

shared_rep_ids <- intersect(rep_ids_neutral, rep_ids_adapt) #936
shared_rep_ids_op<-intersect(rep_ids_adapt,rep_ids_neutral) #936

shared_ui_ids <- intersect(ui_ids_neutral, ui_ids_adapt) #936
shared_ui_ids_op<-intersect(ui_ids_adapt,ui_ids_neutral) #936

#all gave the same of 936 in common
gt_df_shared <- gt_df_combined %>%
  filter(LabID %in% shared_rep_ids)

only_neutral <- setdiff(rep_ids_neutral, rep_ids_adapt) #107
only_adapt   <- setdiff(rep_ids_adapt, rep_ids_neutral) #15

only_neutral_ui <- setdiff(ui_ids_neutral, ui_ids_adapt) #107
only_adapt_ui   <- setdiff(ui_ids_adapt, ui_ids_neutral) #15

both_panels  <- shared_rep_ids #936

length(unique(gt_df_shared$LabID)) #936
length(unique(gt_df_shared$UI_ID)) #936

tibble(
  category = c("Only in Neutral Panel", "Only in Adapt Panel", "In Both Panels"),
  count = c(length(only_neutral), length(only_adapt), length(both_panels))
)

# recalculate maf because we dropped individuals exclusive to adpative and neutral
maf_df_recalc <- gt_df_shared %>%
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
    n_samples = n_distinct(LabID),
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
gt_df_ready_for_plink <- gt_df_shared %>%
  inner_join(
    onlySNPs_maf_0.05 %>%
      select(CHROM, POS),
    by = c("CHROM", "POS")
  )

gt_df_ready_for_plink %>%
  distinct(PANEL, LabID) %>%
  count(PANEL)

any(duplicated(
  gt_df_ready_for_plink %>% distinct(PANEL, CHROM, POS)
))

saveRDS(gt_df_ready_for_plink, file = paste0("./outputs/gt_df_ready_for_PLINK_", Sys.Date(), ".rds"))

vcf_adapt <- read.vcfR("./outputs/CanFam3.1/filtered-by-neg_GWAdapt_CanFam3.1_GATK.vcf")
vcf_neutral <- read.vcfR("../GenomicValidation/outputs//CanFam3.1/filtered-by-neg_NeutralTissue_CanFam3.1_GATK.vcf")

readRDS("../GenomicValidation//outputs/gt_df_ready_for_PLINK.rds")

extract_ref_alt <- function(vcf) {
  tibble(
    CHROM = getCHROM(vcf),
    POS   = as.integer(getPOS(vcf)),
    REF   = getREF(vcf),
    ALT   = getALT(vcf)
  ) %>%
    filter(!grepl(",", ALT)) %>%   # exclude multiallelic just in case
    distinct(CHROM, POS, .keep_all = TRUE)
}

refalt_raw_adapt <- data.frame(
  CHROM = getCHROM(vcf_adapt),
  POS   = getPOS(vcf_adapt),
  REF   = getREF(vcf_adapt),
  ALT   = getALT(vcf_adapt)
)

refalt_raw_adapt %>%
  count(CHROM, POS) %>%
  filter(n > 1)

refalt_raw_neutral <- data.frame(
  CHROM = getCHROM(vcf_neutral),
  POS   = getPOS(vcf_neutral),
  REF   = getREF(vcf_neutral),
  ALT   = getALT(vcf_neutral)
)

refalt_raw_neutral %>%
  count(CHROM, POS) %>%
  filter(n > 1)

refalt_fixed_adapt <- refalt_raw_adapt %>%
  mutate(
    ALT1 = str_split(ALT, ",", simplify = TRUE)[,1]
  ) %>%
  select(CHROM, POS, REF, ALT = ALT1)

refalt_fixed_neutral <- refalt_raw_neutral %>%
  mutate(
    ALT1 = str_split(ALT, ",", simplify = TRUE)[,1]
  ) %>%
  select(CHROM, POS, REF, ALT = ALT1)

refalt_all <- bind_rows(refalt_fixed_adapt, refalt_fixed_neutral) %>%
  distinct(CHROM, POS, .keep_all = TRUE)

gt_df_with_refalt <- gt_df_ready_for_plink %>%
  left_join(refalt_all, by = c("CHROM", "POS"))

missing_refalt <- gt_df_with_refalt %>%
  filter(is.na(REF) | is.na(ALT)) %>%
  distinct(CHROM, POS, PANEL, AmpliconID)

nrow(missing_refalt)

#only for trouble shooting
missing_refalt %>%
  mutate(
    in_adapt   = paste(CHROM, POS) %in%
      paste(refalt_adapt$CHROM, refalt_adapt$POS),
    in_neutral = paste(CHROM, POS) %in%
      paste(refalt_neutral$CHROM, refalt_neutral$POS)
  ) %>%
  count(in_adapt, in_neutral)

#save as vcf for PLINK
vcf_df <- gt_df_with_refalt %>%
  # create a unique variant ID: AmpliconID_CHROM_POS
  mutate(
    VAR_ID = paste0(AmpliconID, "_", CHROM, "_", POS),
    GT_vcf = ifelse(is.na(GT_unphased), "./.", GT_unphased)
  ) %>%
  select(CHROM, POS, VAR_ID, REF, ALT, GT_vcf, UI_ID) %>%
  pivot_wider(
    names_from = UI_ID,
    values_from = GT_vcf
  ) %>%
  arrange(CHROM, POS)

# Build VCF dataframe
vcf_out <- vcf_df %>%
  transmute(
    `#CHROM` = CHROM,
    POS = POS,
    ID = VAR_ID,
    REF = REF,
    ALT = ALT,
    QUAL = ".",
    FILTER = "PASS",
    INFO = ".",
    FORMAT = "GT",
    across(-c(CHROM, POS, VAR_ID, REF, ALT))
  )

vcf_file <- paste0("./outputs/gt_final_filtered", Sys.Date(), ".vcf")

# Write VCF header
vcf_meta <-  c(
  "##fileformat=VCFv4.2",
  "##source=gt_df_with_refalt",
  "##reference=CanFam3.1",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)
# Write meta-information
writeLines(vcf_meta, con = vcf_file)

write(
  paste(colnames(vcf_out), collapse = "\t"),
  file = vcf_file,
  append = TRUE
)

# Write data rows
write.table(
  vcf_out,
  file = vcf_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE
)

vcf_test <- read.vcfR("./outputs/gt_final_filtered2026-06-29.vcf")
