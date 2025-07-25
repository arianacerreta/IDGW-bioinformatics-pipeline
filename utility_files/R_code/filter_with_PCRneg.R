#filtering by PCR negative depths
library(vcfR)
library(tidyverse)
library(ggplot2)
library(SNPfiltR)

#Neutral Tissue
vcf<- read.vcfR("./inputs/CanFam3.1/final_filtered_unzipped_CanFam3_1_DP0.vcf") #vcf without a DP filter from bcftools
panel<- "GWAdapt"
genome<-"CanFam3.1"
#fix column names to remove path name
popmap <- data.frame(id=na.omit(str_extract(colnames(vcf@gt), "i[0-9A-Za-z_-]+(?=_filtered\\.bam)")))# makes popmap object from ID's in vcf file
popmap_long<-data.frame(id=colnames(vcf@gt)[2:length(colnames(vcf@gt))]) 

pcr_neg_popmap<-popmap_long %>% filter(grepl("PCR", id))
popmap_samples<-popmap_long %>% filter(!grepl("PCR",id))

vcfR_pcr_neg <- vcf[samples = pcr_neg_popmap$id]
vcfR_0 <- vcf[samples = popmap_samples$id]

# get DP matrices
actual_dp <- extract.gt(vcfR_0, element = "DP", as.numeric = TRUE)
neg_dp<- extract.gt(vcfR_pcr_neg, element = "DP", as.numeric = TRUE)

# get locus keys
locus_key<-paste(vcfR_0@fix[, "CHROM"], vcfR_0@fix[, "POS"], sep = "_")
rownames(actual_dp)<- locus_key
rownames(neg_dp) <- locus_key

#rename columns
colnames(actual_dp) <- str_extract(colnames(actual_dp), "i[0-9A-Za-z_-]+(?=_filtered\\.bam)") 
colnames(neg_dp) <- str_extract(colnames(neg_dp), "i[0-9A-Za-z_-]+(?=_filtered\\.bam)") 

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
  mutate(prefix = get_prefix(sample))

#join with neg lookup
depths_joined <- left_join(actual_depth_tbl, neg_depth_tbl, by = c("key", "prefix"))

#determine which to mask
depths_to_mask<- depths_joined %>%
  filter(!is.na(neg_dp) & actual_dp < 2*neg_dp)

# mask genotypes
actual_gt <- extract.gt(vcfR_0, element = "GT")
colnames(actual_gt) <- str_extract(colnames(actual_gt), "i[0-9A-Za-z_-]+(?=_filtered\\.bam)") 

for (i in seq_len(nrow(depths_to_mask))) {
  row <- depths_to_mask[i, ]
  loc_idx <- which(rownames(actual_gt) == row$key)
  samp_idx <- which(colnames(actual_gt) == row$sample)
  vcfR_0@gt[loc_idx, samp_idx +1] <- sub("^[^:]+", "./.", vcfR_0@gt[loc_idx, samp_idx + 1])
}

#write intermediate output
write.vcf(vcfR_0, file = paste0("./outputs/filtered-by-neg_", panel, genome ,".vcf"))

#use this output for filtering and error calculations in error_rates.R
