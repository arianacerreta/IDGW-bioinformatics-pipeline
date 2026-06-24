#Tissue/Scat Match Code
#Developed from code written by Elise Stacy modified by Ariana Cerreta

library(adegenet)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(reshape2)

# Plotting
library(ggplot2)
library(gridExtra)

# For organizing matching data
library(igraph)

#Read in Files
{
  IDFG_haps <- read.csv("./inputs/IDFG_data/microhap_genotypes.csv", na.strings = "") #source from microhaplotypes_genotypes_raw.zip from https://github.com/kboudinot/IDGW_structure
  scatopt4 <- read.csv("./inputs/scatopt4/microhap_genotypes.csv", na.strings = "") #use microhap_genotypes_scatopt4.csv from Zenodo
  finalscat <- read.csv("./inputs/finalscat/microhap_genotypes.csv", na.strings = "")#use microhap_genotypes_final_scat.csv from Zenodo
  
  #Combine fecal with IDFG to call alleles at the same time
  both<- scatopt4 %>% rbind(finalscat) %>% mutate(Type = "fecal")
  df_fecal_raw <- both  # rename to your actual object
  
  IDFG_haps<- IDFG_haps %>% mutate(Type = "IDFG")
  df_IDFG_raw<-IDFG_haps
  
  df_IDFG_fecal_raw<- df_fecal_raw %>% rbind(df_IDFG_raw)
  
  #clear up space and don't get confused
  rm(both, df_fecal_raw, df_IDFG_raw, finalscat, IDFG_haps, scatopt4)
}
# Clean raw data
snp_df <- df_IDFG_fecal_raw %>%
  mutate(Indiv = Indiv |>              
           basename() |>                      # remove full path
           sub("\\.1\\.bam$", "", x = _)      # remove trailing ".1.bam"
  )

# Identify all unique alleles per locus
{
  alleles_by_locus <- snp_df %>% select(Locus, Allele1, Allele2) %>%
    pivot_longer(cols = c(Allele1, Allele2), values_to = "Allele") %>%
    filter(!is.na(Allele)) %>% distinct(Locus, Allele) %>%      
    group_by(Locus) %>% mutate(Allele_Index = dense_rank(Allele)) %>%
    ungroup()
  
  # Attach numeric indices and make ordered genotype codes
  df_encoded <- snp_df %>% left_join(alleles_by_locus, by = c("Locus", "Allele1" = "Allele")) %>%
    rename(A1_Code = Allele_Index) %>% left_join(alleles_by_locus, by = c("Locus", "Allele2" = "Allele")) %>%
    rename(A2_Code = Allele_Index) %>%
    mutate(
      Genotype_Code = ifelse(
        is.na(A1_Code) | is.na(A2_Code),
        NA_character_,
        sprintf("%d%d", A1_Code, A2_Code)  # keep A1/A2 order exactly
      )
    )
  }

# Filter to only loci on fecal panel (because IDFG data has all 341)
loci_200panel<-read.table("../../Genomes/wolfAmpRef_200/Clu200_pos_4.txt", header = TRUE)
locinames<-unique(loci_200panel$locus) #there are 199 because one amplicon is the sex locus SRY

#only keep loci that are in locinames for df_encoded
df_encoded<-df_encoded %>% filter(Locus %in% locinames)

#now we need to separate back into IDFG and fecal to clean-up metadata (i.e. LabID, etc.) since ID conventions were different
df_encoded_fecal<- df_encoded %>% filter(Type == "fecal")

df_encoded_IDFG<- df_encoded %>% filter(Type == "IDFG")

# For df_encoded_fecal: 1) extract metadata, 2) filter by PCR negative depth, 3) remove PCR negatives
#extract metadata
{df_encoded_extract_fecal<- df_encoded_fecal %>%
    mutate(
      parts = str_split(Indiv, "_"),
      # NEW: remove "redo", "redo1", "redoXYZ", etc.
      parts = map(parts, ~ .x[!str_detect(.x, "^redo")])
    ) %>%
    mutate(
      Platei7 = map_chr(parts, ~ .x[1]),
      Well    = map_chr(parts, ~ .x[2]),
      LabID   = map_chr(parts, ~ .x[length(.x)]),
      ui_candidate = map_chr(parts, ~ if (length(.x) >= 2) .x[length(.x) - 1] else NA_character_),
      ui_present   = map_lgl(ui_candidate, ~ !is.na(.x) && str_starts(.x, "UI")),
      start_idx = map_int(parts, ~ if (length(.x) >= 3 && .x[3] == "IDGW") 4L else 3L),
      end_idx = map_int(parts, ~ {
        n <- length(.x)
        if (n < 3) return(NA_integer_)
        if (n >= 2 && str_starts(.x[n-1], "UI")) n - 2L else n - 1L
      }),
      ExperimentID = map2_chr(parts, start_idx, ~ {
        p <- .x
        s <- .y
        n <- length(p)
        # compute end dynamically again to ensure consistency:
        e <- if (n >= 2 && str_starts(p[n-1], "UI")) n - 2L else n - 1L
        if (is.na(s) || e < s) NA_character_ else paste(p[s:e], collapse = "_")
      }),
      UI_ID = map2_chr(ui_candidate, ui_present, ~ if (.y) .x else NA_character_)
    ) %>%
    select(-parts, -ui_candidate, -ui_present, -start_idx, -end_idx) %>%
    select(Indiv, LabID, Platei7, ExperimentID, Well, UI_ID, everything())
}
#filter by PCR negative depths
{
##add in depths so that you can filter by depth
df_encoded_depth_fecal<-df_encoded_extract_fecal %>%
  mutate(
    TotalDepth = case_when(
      is.na(Allele1) | is.na(Allele2) |
        is.na(Allele1_count) | is.na(Allele2_count) ~ NA_real_,  # missing data
      Allele1 == Allele2 ~ Allele1_count,                      # homozygous
      TRUE ~ Allele1_count + Allele2_count                     # heterozygous
    ))%>%
  select(LabID, Platei7, ExperimentID, Well, Locus, Allele1, Allele2, Allele1_count, Allele2_count,
         TotalDepth, everything())

#filter by PCRneg depths
neg_depths<- df_encoded_depth_fecal %>%
  filter(str_detect(LabID, "PCR")) %>% mutate(TotalDepth = replace_na(TotalDepth, 0)) %>%   # convert missing to 0
  group_by(Platei7, Locus) %>%
  summarize(
    neg_max = max(TotalDepth, na.rm = TRUE), #since multiple NTC per plate, chose largest
    .groups = "drop"
  )

df_filtered_fecal <- df_encoded_depth_fecal %>% left_join(neg_depths, by = c("Platei7", "Locus"))

df_filtered_fecal <- df_filtered_fecal %>%
  mutate(
    TotalDepth_filtered = case_when(
      # keep PCR negatives untouched
      str_detect(LabID, "PCR") ~ TotalDepth,
      # if no negative exists for this plate+locus, keep original
      is.na(neg_max) ~ TotalDepth,
      # censor if depth < 2 × negative depth
      TotalDepth < 2 * neg_max ~ NA_real_,
      TRUE ~ TotalDepth),
    # censor genotype whenever the filtered depth is NA just in case
    Genotype_Code = if_else(
      is.na(TotalDepth_filtered),
      NA_character_,
      Genotype_Code))

genotype_matrix_fecal <- df_filtered_fecal %>%
  select(Indiv, LabID, Platei7, ExperimentID, Well, Type, Locus, Genotype_Code) %>%
  pivot_wider(
    names_from = Locus,
    values_from = Genotype_Code
  )

#identify technical replicates
rep_map <- genotype_matrix_fecal %>%
  distinct(LabID, Well, .keep_all = FALSE) %>%  # one row per LabID x Well
  group_by(LabID) %>%
  mutate(Replicate = row_number()) %>%
  ungroup()

genotype_matrix_fecal <- genotype_matrix_fecal %>%
  left_join(rep_map, by = c("LabID", "Well")) %>%
  select(LabID, Platei7, ExperimentID, Well, Replicate, everything())
}
#remove PCR negatives
{
#since we already used negative for quality control of genotypes, filter out  PCR negatives
genotype_matrix_fecal<- genotype_matrix_fecal%>% filter(!str_detect(LabID, "PCR"))
}

#number of unique LabIDs
length(unique(genotype_matrix_fecal$LabID)) #790

# Scoring matching microhaplotype genotypes to calculate mismatch rate
#make a function
score_gt_match<-function(gt1,gt2){
  #Handle missing genotypes
  if (is.na(gt1) & is.na(gt2)) return("BothNA")
  else if (is.na(gt1) | is.na(gt2)) return("OneNA")
  
  # Convert to character
  gt1 <- as.character(gt1)
  gt2 <- as.character(gt2)
  
  #split alleles
  a1 <- strsplit(gt1, "")[[1]]
  a2 <- strsplit(gt2, "")[[1]]
  
  #sort alleles for unphased comparison
  a1<- sort(a1)
  a2<- sort(a2)
  
  #score
  if (identical(a1,a2)) return("Full Match")
  else if (any(a1 %in% a2)) return("Partial Match")
  else return("Mismatch")
}

format_paired_gt_df <- function(genotype_matrix) {
  # Identify replicate groups
  replicate_map <- genotype_matrix %>%
    select(Indiv, LabID) %>%
    group_by(LabID) %>%
    filter(n() >= 2) %>%
    summarise(pairs = list(as.data.frame(t(combn(Indiv, 2)))), .groups = "drop") %>%
    unnest(pairs) %>%
    rename(Sample1 = V1, Sample2 = V2)
  
  # Pivot genotype matrix to long format
  gt_long <- genotype_matrix %>%
    select(Indiv, starts_with("Clu_")) %>%   # adjust if your locus names differ
    pivot_longer(cols = -Indiv,
                 names_to = "Locus",
                 values_to = "GT")
  
  # Join genotypes for sample pairs
  paired_gt_df <- replicate_map %>%
    left_join(gt_long, by = c("Sample1" = "Indiv"), relationship = "many-to-many") %>%
    rename(GT1 = GT) %>%
    left_join(gt_long, by = c("Sample2" = "Indiv", "Locus"), relationship = "many-to-many") %>%
    rename(GT2 = GT)
  
  # Apply scoring
  paired_gt_df <- paired_gt_df %>%
    mutate(MatchScore = mapply(score_gt_match, GT1, GT2))
  
  return(paired_gt_df)
}

# Process data and create summaries in one function
calculate_gt_pairs_errors <- function(paired_gt_df) {
  # Create gt_pairs_mismatch dataframe
  gt_pairs_mismatch <- paired_gt_df %>%
    rowwise() %>%
    mutate(MatchType = score_gt_match(GT1, GT2)) %>%
    ungroup()
  
  # Calculate summary_by_pair
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
  
  # Calculate summary_by_locus
  summary_by_locus <- gt_pairs_mismatch %>%
    group_by(Locus) %>%
    summarize(
      N=n(),
      Full_Match = sum(MatchType == "Full Match", na.rm = TRUE),
      Mismatch = sum(MatchType == "Mismatch", na.rm = TRUE),
      Partial_Match = sum(MatchType == "Partial Match", na.rm = TRUE),
      Missing_Full = sum(MatchType == "BothNA", na.rm = TRUE),
      Missing_Partial = sum(MatchType == "OneNA", na.rm = TRUE),
      Full_Match_Rate = Full_Match / (Full_Match + Partial_Match + Mismatch),
      Both_Mismatch_Rate = (Partial_Match + Mismatch) / (Full_Match + Partial_Match + Mismatch),
      Full_Missingness = Missing_Full / N,
      Partial_Missingness = Missing_Partial / N,
      Inclusive_Missingness = (Missing_Full + Missing_Partial) / N
    )
  
  # Return both summaries
  list(summary_by_pair = summary_by_pair,
       summary_by_locus = summary_by_locus)
}

# Initialize storage data frames
pair_errors <- list()
locus_errors <- list()

depth_thresholds <- seq(10, 20, by = 5)
missing_thresholds <- seq(0.3, 0.50, by = 0.05)

# Pivot df_heat to wide format to match genotype_matrix_noneg
## based on  df_filtered_fecal which still has PCR negatives
depth_wide_fecal <- df_filtered_fecal %>%
  select(Indiv, Locus, TotalDepth) %>%
  pivot_wider(names_from = Locus, values_from = TotalDepth, names_prefix = "Depth_")

#number of rows in depth_wide_fecal and genotype_matrix_fecal 
#(should be larger than genotype_matrix_fecal because of PCR negatives)
nrow(depth_wide_fecal) #1656
nrow(genotype_matrix_fecal) #1620

# genotype_matrix_fecal already has IndivID as key
genotype_with_depth <- genotype_matrix_fecal %>%left_join(depth_wide_fecal, by = "Indiv")  

#check to see PCR negatives removed (i.e., matches genotype_matrix_fecal)
nrow(genotype_with_depth) #1620

#Loop through depth values by increments of 5
for (depth in depth_thresholds) {
  # Loop through missingness by increments of 0.05 (0.4 threshold would allow up to 40% missingness in a sample)
  for (missing in missing_thresholds) {
    
    genotype_filtered <- genotype_with_depth %>%
      mutate(across(starts_with("Clu_"), 
                    ~ {
                      depth_col <- paste0("Depth_", cur_column())
                      # If depth is below threshold, set genotype to NA
                      ifelse(get(depth_col) < depth, NA, .x)
                    }))
    
    genotype_filtered <- genotype_filtered %>%
      rowwise() %>%
      mutate(
        total_loci = length(c_across(starts_with("Clu_"))),        # total loci
        loci_missing = sum(is.na(c_across(starts_with("Clu_")))),  # count missing
        prop_missing = loci_missing / total_loci
      ) %>%
      ungroup() %>%
      filter(prop_missing <= missing)  # missing threshold
    
    
    # Format paired genotype data frame
    popmap_samples <- genotype_filtered %>% select(LabID, Platei7, Well)
    
    paired_gt_df <- format_paired_gt_df(genotype_filtered)
    
    # Calculate error rate
    error_rate <- calculate_gt_pairs_errors(paired_gt_df)
    
    # Add current Depth and GQ values to error summaries
    pair_error <- cbind(error_rate$summary_by_pair, Depth = depth, Missing = missing)
    locus_error <- cbind(error_rate$summary_by_locus, Depth = depth, Missing = missing)
    
    # Accumulate error data
    pair_errors <- rbind(pair_errors, pair_error)
    locus_errors <- rbind(locus_errors, locus_error)
  }
}

#summarise mismatch results
pair_errors$ID_1 <- sub(".*_", "", pair_errors$Sample1)
pair_errors$ID_2 <- sub(".*_", "", pair_errors$Sample2)

pair_err_summ <- pair_errors %>% # Summarize mismatch data
  group_by(Depth, Missing) %>%
  summarise(
    N_Samples = n(), # Replicates are paired, so N_samples is # of replicate pairs retained under the filtering scheme
    Avg_Mismatch = round(mean(Mismatch_Rate, na.rm = TRUE), 4),
    Max_Mismatch = round(max(Mismatch_Rate, na.rm = TRUE), 4),
    Min_Locus_Overlap = min(((1 - Inclusive_Missingness) * N), na.rm = TRUE),
    Avg_Missing = round(mean(Missing_Full, na.rm = TRUE), 1)
  )

# Plot error rates
mismatch_plot<-ggplot(pair_err_summ, aes(x=Missing, y = Depth, fill = Avg_Mismatch)) +
  geom_tile() +
  scale_fill_gradient2(low = "lightblue", mid ="purple", high = "red", midpoint = 0.003) + # Adjust to the midpoints of your error distribution
  geom_text(aes(label = paste("N =", N_Samples, "pairs")), color = "black", size = 3) +
  labs(title = "Sample replicate pairs retained and average mismatch", subtitle = "across depth and missingness filtering", fill = "Avg. Mismatch") +
  xlab("Sample Missingness")
mismatch_plot

ggsave("./outputs/fecal_mismatch_filtering.png", plot = mismatch_plot, width = 8, height = 6, dpi = 300)


pair_errors_15_0.5<- pair_errors %>% 
  filter(Depth == 15, Missing == 0.50)

mean(pair_errors_15_0.5$Mismatch_Rate) #0.004264144
max(pair_errors_15_0.5$Mismatch_Rate) #0.1882353

#Apply filtering scheme and calculate heterozygosity
## We selected depth = 15, missingess = 0.5

# Apply filtering scheme so that you can look at individual heterozygosity
##by depth = 15
genotype_filtered_15_0.5_fecal <- genotype_with_depth %>%
  mutate(across(starts_with("Clu_"), 
                ~ {
                  depth_col <- paste0("Depth_", cur_column())
                  # If depth is below threshold, set genotype to NA
                  ifelse(get(depth_col) < 15, NA, .x)
                }))

#by samplemissingness = 0.5
genotype_filtered_15_0.5_fecal <- genotype_filtered_15_0.5_fecal %>%
  rowwise() %>%
  mutate(
    total_loci = length(c_across(starts_with("Clu_"))),        # total loci
    loci_missing = sum(is.na(c_across(starts_with("Clu_")))),  # count missing
    prop_missing = loci_missing / total_loci
  ) %>%
  ungroup() %>%
  filter(prop_missing <= 0.5)

mean(genotype_filtered_15_0.5_fecal$prop_missing) #0.07835678

# Convert matrix object to adegenet's genlind format for microhaplotype
gen_data_fecal<-genotype_filtered_15_0.5_fecal %>%
  select(-LabID, -Platei7, -Well, -ExperimentID, -Replicate,-total_loci, -loci_missing, -prop_missing,-Type,
         -starts_with("Depth_"))

ind.names<-gen_data_fecal$Indiv
gen_data_fecal$Indiv<- NULL

#updated with filtering scheme; only 1000 samples
genind_fecal <- df2genind(
  gen_data_fecal,
  ind.names = ind.names,
  ploidy = 2,
  sep = "",
  NA.char = "NA"
)

# Convert to matrix
geno <- tab(genind_fecal, NA.method = "asis")

# split allele columns by locus
locus_per_col <- sub("\\.[0-9]+$", "", colnames(geno))

locus_cols <- split(seq_len(ncol(geno)), locus_per_col)

het_matrix <- sapply(locus_cols, function(cols){
  rowSums(geno[, cols, drop = FALSE] == 1, na.rm = TRUE) == 2
})

non_missing <- sapply(locus_cols, function(cols){
  rowSums(!is.na(geno[, cols, drop = FALSE])) > 0
})

num_het <- rowSums(het_matrix, na.rm = TRUE)

num_genotyped <- rowSums(non_missing)

total_loci <- ncol(het_matrix)

missing_rate <- 1 - (num_genotyped / total_loci)

#per individual #het loci/#total called loci
Ho <- rowSums(het_matrix) / rowSums(non_missing)

qc_table <- data.frame(
  Individual = rownames(geno),
  Ho = Ho,
  Heterozygous_Loci = num_het,
  Genotyped_Loci = num_genotyped,
  Missing_Rate = missing_rate
)

qc_table <- qc_table[order(qc_table$Ho), ] 

summary(Ho)
#Min.     1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.008696 0.456853 0.484849 0.483376 0.517767 0.775701 
sd(Ho) #0.06355988

# High het and high replicate mismatch is a good indication of contamination
# also sort out Ho <0.2, we had some that had abnormally low Ho based on panel design

#already filtered for missingness (0.5) and depth (15): genotype_filtered_15_0.5_fecal

High_het <- qc_table %>% filter(Ho > 0.65 | Ho<0.2)

High_rep_error <- pair_errors %>% 
  filter(Depth == 15, Missing == 0.50) %>%
  filter(Mismatch_Rate > 0.05)

Remove_samples_fecal <- unique(c(High_het$Individual, High_rep_error$Sample1, High_rep_error$Sample2))

#filter out samples that failed by het and mismatch

genotype_with_depth_fecal<- genotype_filtered_15_0.5_fecal %>%
  filter(!(Indiv %in% Remove_samples_fecal))

best_reps<- genotype_with_depth_fecal %>%
  #compute missingness and mean depth per row
  rowwise() %>%
  mutate(
    total_loci = length(c_across(starts_with("Clu_"))),        # total loci
    loci_missing = sum(is.na(c_across(starts_with("Clu_")))),  # count missing
    prop_missing = loci_missing / total_loci,
    
    mean_depth = mean(c_across(starts_with("Depth_Clu_")), na.rm = TRUE)) %>%
  ungroup() %>%
  select(prop_missing, mean_depth, everything())

kept_best<- best_reps %>%
  # For each LabID, keep the row with:
  # (1) lowest missingness
  # (2) highest mean_depth (tie-breaker)
  group_by(LabID) %>%
  arrange(prop_missing, desc(mean_depth)) %>%
  slice(1) %>%     # keep the best row
  ungroup()

#double check everything still meets minimum filtering criteria
genotype_filtered_final_fecal <- kept_best %>%
  mutate(across(starts_with("Clu_"), 
                ~ {
                  depth_col <- paste0("Depth_", cur_column())
                  # If depth is below threshold, set genotype to NA
                  ifelse(get(depth_col) < 15, NA, .x)
                }))


genotype_filtered_final_fecal <- genotype_filtered_final_fecal %>%
  rowwise() %>%
  mutate(
    total_loci = length(c_across(starts_with("Clu_"))),        # total loci
    loci_missing = sum(is.na(c_across(starts_with("Clu_")))),  # count missing
    prop_missing = loci_missing / total_loci
  ) %>%
  ungroup() %>%
  filter(prop_missing <= 0.5)

#individual identification success rate
length(unique(genotype_filtered_final_fecal$LabID))/length(unique(genotype_matrix_fecal$LabID))

###export and format with diagnostic loci for database
#all loci except for sry

clean_fecal_delomas_pipeline<- genotype_filtered_final_fecal%>%
  select(-prop_missing,-mean_depth, -Platei7, -ExperimentID, 
         -Well, -Replicate, -starts_with("Depth_"), -total_loci, -loci_missing)%>%
  pivot_longer(
    cols= starts_with("Clu"),
    names_to = "Locus",
    values_to = "GT"
  )

write.csv(clean_fecal_delomas_pipeline, file = paste0("./outputs/fecal_all_genos_long_", Sys.Date(), ".csv"))  
write.csv(alleles_by_locus, file = paste0("./outputs/microhap_allele_key_",Sys.Date(), ".csv"))
