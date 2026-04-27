### Final Script: Filter with PCRneg for Adaptive Panel Analyses ###
### As recommended in Hervey et al. 2025 (https://onlinelibrary.wiley.com/doi/10.1002/ece3.71240), we required at least 2x depth as compared to the corresponding NTC for each plate
### e.g. if Locus1 for Sample1 for PlateA had a read depth of 3 and the NTC for Locus1 on Plate A had a read depth of 10, Locus1 was censored because 3 < 10*2
### e.g. if Locus2 for Sample1 for PlateA had a read depth of 30 and the NTC for Locus2 on Plate A had a read depth of 10, Locus2 was kept because 30 > 20

#libraries (base R version 4.2.3)
library(vcfR) #version 1.15.0
library(tidyverse) #version 2.0.0: dplyr 1.2.1, forcats 1.0.1, lubridate 1.9.5, purrr 1.2.2, readr 2.2.0, stringr 1.6.0, tibble 3.3.1, tidyr 1.3.2
library(ggplot2) #version 4.0.2
library(SNPfiltR) #version 1.0.1

#read in unzipped .vcf created in Step 6: Genotyping of the GATK Pipeline
vcf<- read.vcfR("../../AdaptiveOptim/GATK/CanFam3_1_final_vcf_GWAdapt_AdaptOptim/unzipped.all.gwadapt.adaptoptim.gatk.genotyped.vcf") #vcf without a DP filter from GATK

#read in metadata for AdaptOptim for UI_IDs
adaptoptim_meta<-read.csv("./inputs/AdaptOptim_WolfMeta_Age_Disp_Breed_Sex_Coat.csv")

#labels for saving
panel<- "GWAdapt"
genome<-"CanFam3.1"
pipeline<-"GATK"

#fix column names to remove path name
popmap <- data.frame(id = na.omit(colnames(vcf@gt) |> str_remove("_RG_sort\\.recal$") |> str_remove("_R[12]$") |> str_remove("~[0-9]+"))) %>%
  slice(-1)
popmap_long<-data.frame(id=colnames(vcf@gt)[2:length(colnames(vcf@gt))]) #as vamed in vcf

#reformat optimization names to match UI runs
popmap<- popmap %>%
  mutate(i7 = if_else(str_detect(id,"GTseek"), str_split_i(id, "_", 3), str_split_i(id, "_", 1)),
         well = if_else(str_detect(id,"GTseek"), str_split_i(id, "_", 4), str_split_i(id, "_", 2)),
         platename = if_else(str_detect(id,"GTseek"), str_split_i(id, "_", 1), str_split_i(id, "_", 3)),
         UI_ID = case_when(str_detect(id,"GTseek") ~ NA_character_ , 
                           str_detect(id,"PCR_Neg") ~ paste0(str_split_i(id, "_", 4), str_split_i(id, "_", 5)),
                           str_detect(id,"PCRneg") ~ "PCRNeg",
                           TRUE ~ str_split_i(id, "_", 4)),
         LabID = case_when(str_detect(id, "WL-") ~ str_split_i(id, "-", 2), 
                           str_detect(id,"GTseek") ~ str_split_i(id, "_", 6), 
                           str_detect(id, "PCR_Neg") ~  str_split_i(id, "_", 6),
                           str_detect(id, "PCRneg") ~ str_split_i(id, "_", 4),
                           TRUE ~ str_split_i(id, "_", 5))) %>%
  mutate(LabID = if_else(LabID == "C10807", "C010807", LabID))%>% #special case C10807 is C010807
  mutate(LabID = case_when(LabID == "AGW1" ~ "PCRnegA1", #special case fix PCR neg naming
                           LabID == "AGW2" ~ "PCRnegA2",
                           LabID == "AGW3" ~ "PCRnegA3",
                           LabID == "AGW4" ~ "PCRnegA4",
                           LabID == "AGW5" ~ "PCRnegA5",
                           LabID == "AGW6" ~ "PCRnegA6",
                           TRUE ~ LabID)) %>%
  left_join(adaptoptim_meta %>% select(LabID, UI_ID), by = "LabID", suffix = c("", "_meta")) %>%
  mutate(UI_ID = if_else(str_detect(id, "GTseek"),UI_ID_meta,UI_ID)) %>%
  select(-UI_ID_meta)%>%
  mutate(stand_id = paste(i7,well, platename, UI_ID,LabID, sep = "_"))

#make a pcr negative pop map
pcr_neg_popmap<-popmap_long %>% filter(grepl("PCR", id))
popmap_samples<-popmap_long %>% filter(!grepl("PCR",id))

#filter vcfs to have either only samples or only PCRnegatives
vcfR_pcr_neg <- vcf[samples = pcr_neg_popmap$id]
vcfR_0 <- vcf[samples = popmap_samples$id]
  
# get DP matrices
actual_dp <- extract.gt(vcfR_0, element = "DP", as.numeric = TRUE)
neg_dp<- extract.gt(vcfR_pcr_neg, element = "DP", as.numeric = TRUE)

# get locus keys
locus_key<-paste(vcfR_0@fix[, "CHROM"], vcfR_0@fix[, "POS"], sep = "_")
rownames(actual_dp)<- locus_key
rownames(neg_dp) <- locus_key

#rename columns (column names should match popmap$id)
colnames(actual_dp) <- colnames(actual_dp) |> str_remove("_RG_sort\\.recal$") |> str_remove("_R[12]$") |> str_remove("~[0-9]+")
colnames(neg_dp) <- colnames(neg_dp) |> str_remove("_RG_sort\\.recal$") |> str_remove("_R[12]$") |> str_remove("~[0-9]+")

#make a function to get the prefix
get_prefix<- function(x) str_extract(x, "i\\d{3}")

# === Build a lookup table for neg DPs ===
neg_depth_tbl <- as_tibble(neg_dp, rownames = "key") %>%
  pivot_longer(-key, names_to = "sample", values_to = "neg_dp") %>%
  mutate(prefix = get_prefix(sample)) %>%
  filter(!is.na(neg_dp))%>%
  group_by(key, prefix) %>%
  summarize(neg_dp = max(neg_dp), .groups = "drop")

#tibble of depths from samples 
actual_depth_tbl <- as_tibble(actual_dp, rownames = "key") %>%
  pivot_longer(-key, names_to = "sample", values_to = "actual_dp") %>%
  mutate(prefix = get_prefix(sample)) #GTseek Adapt Optim samples did not have a corresponding NTC to filter by; receives NA in prefix column

#join with neg lookup
depths_joined <- left_join(actual_depth_tbl, neg_depth_tbl, by = c("key", "prefix"))

#determine which to mask
depths_to_mask<- depths_joined %>%
  filter(!is.na(neg_dp) & actual_dp < 2*neg_dp)

# mask genotypes
actual_gt <- extract.gt(vcfR_0, element = "GT")
colnames(actual_gt) <- colnames(actual_dp) |> str_remove("_RG_sort\\.recal$") |> str_remove("_R[12]$") |> str_remove("~[0-9]+")

for (i in seq_len(nrow(depths_to_mask))) {
  row <- depths_to_mask[i, ]
  loc_idx <- which(rownames(actual_gt) == row$key)
  samp_idx <- which(colnames(actual_gt) == row$sample)
  vcfR_0@gt[loc_idx, samp_idx +1] <- sub("^[^:]+", "./.", vcfR_0@gt[loc_idx, samp_idx + 1])
}

name_key<- setNames(popmap$stand_id, popmap$id)
vcf_names_clean <-colnames(vcfR_0@gt)[-1] |>str_remove("_RG_sort\\.recal$") |>str_remove("_R[12]$") |>str_remove("~[0-9]+")
colnames(vcfR_0@gt)[-1] <- name_key[vcf_names_clean]
#write intermediate output
write.vcf(vcfR_0, file = paste0("./outputs/filtered-by-neg_", panel,"_", genome ,"_", pipeline, ".vcf"))
