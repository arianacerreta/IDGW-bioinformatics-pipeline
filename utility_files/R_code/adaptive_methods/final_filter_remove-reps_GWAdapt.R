#### filter GWAdapt and remove replicates  ####
###includes side quest to calculate mismatch between the only biological replicate###

#GWAdapt - contains 282 amplicons of putative adaptive loci (330 unfiltered variants)

#libraries (base R version 4.2.3)
library(vcfR) #version 1.15.0
library(tidyverse) #version 2.0.0: dplyr 1.2.1, forcats 1.0.1, lubridate 1.9.5, purrr 1.2.2, readr 2.2.0, stringr 1.6.0, tibble 3.3.1, tidyr 1.3.2
library(ggplot2) #version 4.0.2
library(SNPfiltR) #version 1.0.1

vcf_adapt<-read.vcfR("./outputs/filtered-by-neg_GWAdapt_CanFam3.1_GATK.vcf")

#filtering parameters
GQ<-30 #genotype quality of 30
DP<-as.numeric(20) #filter out depths less than 20
SNP<-0.5 #50% of samples have locus called
samp<-0.5 #missingness rate of sample, i.e. 50% complete gentoype

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
  hard_filter_rgq(gq = GQ) %>%
  missing_by_snp(cutoff = SNP) %>%
  missing_by_sample(cutoff = samp)

gt_adapt <- extract.gt(vcfR_adapt, element = "GT", as.numeric = FALSE) #this no longer has PCR negatives
gt_df_adapt <- as.data.frame(gt_adapt)

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

mean_depth_adapt <- as.data.frame(dp_adapt) %>%
  pivot_longer(cols = everything(), names_to = "SampleID", values_to = "Depth") %>%
  group_by(SampleID) %>%
  summarize(mean_depth = mean(Depth, na.rm = TRUE), .groups = "drop")

#define replicate pairs
sample_names_adapt <- unique(gt_df_adapt$SampleID)

#for GWAdapts
samples_df_adapt <- tibble(SampleName = sample_names_adapt) %>%
  mutate(
    ReplicateID = str_split_i(SampleName, "_",5)  # grab the last underscore-separated bit
  ) %>%
  mutate(
    UI_ID = str_split_i(SampleName, "_", 4)
  )

rep_by_id_adapt <- samples_df_adapt %>%
  group_by(ReplicateID) %>%
  filter(n() > 1) %>%  # only keep true replicates
  summarise(Pairs = list(as_tibble(
    t(combn(SampleName, 2)), .name_repair = ~ c("Sample1", "Sample2")
  )), .groups = "drop") %>%
  unnest(Pairs) %>%
  mutate(RepType = "Technical")

rep_by_ui_adapt <-samples_df_adapt %>%
  filter(UI_ID != "NA") %>%   # <- drop NA UI_ID rows
  group_by(UI_ID) %>%
  filter(n_distinct(ReplicateID) > 1, n() > 1) %>%  # more than one RepID = bio replicates
  summarise(Pairs = list(as_tibble(
    t(combn(SampleName, 2)), .name_repair = ~ c("Sample1", "Sample2")
  )), .groups = "drop") %>%
  unnest(Pairs) %>%
  mutate(RepType = "Biological")

#### Side Quest to Calculate mismatch rate between biological replicates ###
{
  score_gt_match<-function(gt1,gt2){
  #Handle missing genotypes
  if (is.na(gt1) & is.na(gt2)) return("BothNA")
  else if (is.na(gt1) | is.na(gt2)) return("OneNA")
  
  #split alleles
  a1 <- unlist(strsplit(gt1, "[/|]"))
  a2 <- unlist(strsplit(gt2, "[/|]"))
  
  #sort alleles for unphased comparison
  a1<- sort(a1)
  a2<- sort(a2)
  
  #score
  if (identical(a1,a2)) return("Full Match")
  else if (any(a1 %in% a2)) return("Partial Match")
  else return("Mismatch")
}

paired_gt_df <- rep_by_ui_adapt %>%
  left_join(gt_df_adapt,
            by = c("Sample1" = "SampleID"),
            relationship = "many-to-many") %>%
  rename(GT1 = GT,
         CHROM = CHROM,
         POS = POS) %>%
  left_join(gt_df_adapt, by = c("Sample2" = "SampleID", "CHROM", "POS")) %>%
  rename(GT2 = GT)

gt_pairs_mismatch <- paired_gt_df %>%
  rowwise() %>%
  mutate(MatchType = score_gt_match(GT1, GT2)) %>%
  ungroup()

#summarize error rates
summary_by_pair <- gt_pairs_mismatch %>%
  group_by(Sample1, Sample2) %>%
  summarize(
    N = n(),
    Full_Match = sum(MatchType == "Full Match", na.rm = TRUE),
    Partial_Match = sum(MatchType == "Partial Match", na.rm = TRUE),
    Mismatch = sum(MatchType == "Mismatch", na.rm = TRUE),
    Missing_Full = sum(MatchType == "BothNA", na.rm = TRUE),
    Missing_Partial = sum(MatchType == "OneNA", na.rm = TRUE),
    Full_Match_Rate = Full_Match / (Full_Match + Partial_Match + Mismatch),
    Mismatch_Rate = (Mismatch + Partial_Match) / (Full_Match + Partial_Match + Mismatch),
    Full_Missingness = Missing_Full / N,
    Partial_Missingness = Missing_Partial / N,
    Inclusive_Missingness = (Missing_Full + Missing_Partial) / N
  )
#pair with i027_C06_GWAdapt11_UI1755_ZZAT and i048_D06_GWAdapt6_UI1755_YZL9 reported in manuscript
}
##### end side quest #####
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
    MatchID=coalesce(ReplicateID,UI_ID)
  )%>%
  filter(!str_detect(SampleID, "ZZAT"))

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
                                samples_df_adapt %>% select(SampleName, ReplicateID, UI_ID),
                                by = c("SampleID" = "SampleName"))

rep_ids_unique <- unique(gt_df_adapt_noreps$ReplicateID)
rep_ids_all <- unique(gt_df_adapt_noreps$SampleID)

length(rep_ids_all) == length(rep_ids_unique)

length(unique(gt_df_adapt_noreps$SampleID)) #871

##remove UI2107 which is a dog

gt_df_final<-gt_df_adapt_noreps %>%
  filter(!str_detect(SampleID, "UI2107")) 

length(unique(gt_df_final$SampleID)) #870
write.table(gt_df_adapt_noreps, file = "./outputs/CanFam3.1/adaptive_canfam3_1_gt_GATK_4-28-26.txt", sep = "\t",
            row.names = FALSE)
