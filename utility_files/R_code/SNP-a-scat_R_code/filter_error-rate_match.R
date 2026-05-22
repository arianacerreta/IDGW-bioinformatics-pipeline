---
title: "Tissue_Scat_match"
author: "Ariana Cerreta"
date: "2025-12-08"
output: html_document
#Some of this code was developed from code written by Elise Stacy
---
```{r}
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
```
Read in Files
```{r}
IDFG_haps <- read.csv("./inputs/IDFG_data/microhap_genotypes_IDFG_haps.csv", na.strings = "")

scatopt4 <- read.csv("./inputs/scatopt4/microhap_genotypes.csv", na.strings = "")

finalscat <- read.csv("./inputs/finalscat/microhap_genotypes_final_scat.csv", na.strings = "")

```
Combine fecal with IDFG to call alleles at the same time
```{r}
both<- scatopt4 %>%
  rbind(finalscat)%>%
  mutate(Type = "fecal")

df_fecal_raw <- both  # rename to your actual object

IDFG_haps<- IDFG_haps %>%
  mutate(Type = "IDFG")

df_IDFG_raw<-IDFG_haps

df_IDFG_fecal_raw<- df_fecal_raw %>%
  rbind(df_IDFG_raw)

#clear upspace and don't get confused
rm(both, df_fecal_raw, df_IDFG_raw, finalscat, IDFG_haps, scatopt4)

# Clean raw data
snp_df <- df_IDFG_fecal_raw %>%
  mutate(Indiv = Indiv |>              
      basename() |>                      # remove full path
      sub("\\.1\\.bam$", "", x = _)      # remove trailing ".1.bam"
  )

# Identify all unique alleles per locus
alleles_by_locus <- snp_df %>%
  select(Locus, Allele1, Allele2) %>%
  pivot_longer(cols = c(Allele1, Allele2), values_to = "Allele") %>%
  filter(!is.na(Allele)) %>%      
  distinct(Locus, Allele) %>%      
  group_by(Locus) %>%
  mutate(Allele_Index = dense_rank(Allele)) %>%
  ungroup()

# Attach numeric indices and make ordered genotype codes
df_encoded <- snp_df %>%
  left_join(alleles_by_locus,
            by = c("Locus", "Allele1" = "Allele")) %>%
  rename(A1_Code = Allele_Index) %>%
  left_join(alleles_by_locus,
            by = c("Locus", "Allele2" = "Allele")) %>%
  rename(A2_Code = Allele_Index) %>%
  mutate(
    Genotype_Code = ifelse(
      is.na(A1_Code) | is.na(A2_Code),
      NA_character_,
      sprintf("%d%d", A1_Code, A2_Code)  # keep A1/A2 order exactly
    )
  )
```
Filter to only loci on fecal panel
```{r}
#first which loci are we removing
loci_200panel<-read.table("../../Genomes/wolfAmpRef_200/Clu200_pos_4.txt", header = TRUE)
locinames<-unique(loci_200panel$locus)

#only keep loci that are in locinames for df_encoded
df_encoded<-df_encoded %>%
  filter(Locus %in% locinames)
```
Clean up data for fecal
```{r}
###now we need to separate to clean-up metadata (i.e. LabID, etc.)
df_encoded_fecal<- df_encoded %>%
  filter(Type == "fecal")

df_encoded_IDFG<- df_encoded %>%
  filter(Type == "IDFG")

# For df_encoded_fecal: 1) extract metadata, 2) filter by PCR negative depth, 3) remove PCR negatives
df_encoded_extract_fecal<- df_encoded_fecal %>%
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


#add in depths so that you can filter by depth
df_encoded_depth_fecal<-df_encoded_extract_fecal %>%
  mutate(
    TotalDepth = case_when(
      is.na(Allele1) | is.na(Allele2) |
        is.na(Allele1_count) | is.na(Allele2_count) ~ NA_real_,  # missing data
      
      Allele1 == Allele2 ~ Allele1_count,                      # homozygous
      
      TRUE ~ Allele1_count + Allele2_count                     # heterozygous
    )
  )%>%
  select(LabID, Platei7, ExperimentID, Well, Locus, Allele1, Allele2, Allele1_count, Allele2_count,
         TotalDepth, everything())

#filter by PCRneg depths
neg_depths<- df_encoded_depth_fecal %>%
  filter(str_detect(LabID, "PCR")) %>%
  mutate(TotalDepth = replace_na(TotalDepth, 0)) %>%   # convert missing to 0
  group_by(Platei7, Locus) %>%
  summarize(
    neg_max = max(TotalDepth, na.rm = TRUE),
    .groups = "drop"
  )

df_filtered_fecal <- df_encoded_depth_fecal %>%
  left_join(neg_depths, by = c("Platei7", "Locus"))

df_filtered_fecal <- df_filtered_fecal %>%
  mutate(
    TotalDepth_filtered = case_when(
      # keep PCR negatives untouched
      str_detect(LabID, "PCR") ~ TotalDepth,
      
      # if no negative exists for this plate+locus, keep original
      is.na(neg_max) ~ TotalDepth,
      
      # censor if depth < 2 × negative depth
      TotalDepth < 2 * neg_max ~ NA_real_,
      
      TRUE ~ TotalDepth
    ),
    
    # NEW: censor genotype whenever the filtered depth is NA
    Genotype_Code = if_else(
      is.na(TotalDepth_filtered),
      NA_character_,
      Genotype_Code
    )
  )

genotype_matrix_fecal <- df_filtered_fecal %>%
  select(Indiv, LabID, Platei7, ExperimentID, Well, Type, Locus, Genotype_Code) %>%
  pivot_wider(
    names_from = Locus,
    values_from = Genotype_Code
  )

rep_map <- genotype_matrix_fecal %>%
  distinct(LabID, Well, .keep_all = FALSE) %>%  # one row per LabID x Well
  group_by(LabID) %>%
  mutate(Replicate = row_number()) %>%
  ungroup()

genotype_matrix_fecal <- genotype_matrix_fecal %>%
  left_join(rep_map, by = c("LabID", "Well")) %>%
  select(LabID, Platei7, ExperimentID, Well, Replicate, everything())

#since we already used negative for quality control of genotypes, filter out  PCR negatives
genotype_matrix_fecal<- genotype_matrix_fecal%>% #used to be called genotype_matrix_noneg
  filter(!str_detect(LabID, "PCR")) 

#number of unique LabIDs
length(unique(genotype_matrix_fecal$LabID))
```
Plot depth for fecal
```{r}
df_heat_fecal <- df_filtered_fecal %>%
  select(Indiv, ExperimentID, Locus, TotalDepth)%>%
  mutate(depth_bin = cut(
    TotalDepth,
    breaks = c(-Inf, 10, 25, 50, 100,Inf),
    labels = c("0–10", "11–25", "26–50", "51–100", ">100"),
    right = TRUE
  ))

depth_colors <- c(
  "0–10"     = "red",
  "11–25"    = "orange",
  "26–50"    = "yellow",
  "51–100"   = "green",
  ">100"     = "darkgreen"
)

ggplot(df_heat_fecal, aes(x = Indiv, y = Locus, fill = depth_bin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Depth (binned)",
    values = depth_colors,
    drop = FALSE  # keeps all bins in legend even if not present
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),      # hide x labels
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),      # hide y labels
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Individuals",
    y = "Loci",
    fill = "Total Depth"
  )
```
## Visualize sample depth and missingness by plate
## Loop through plotting heatmap and histogram for each plate
```{r}
#for microhaplotypes

plots <- map(unique(df_heat_fecal$ExperimentID), function(exp) {

  df_sub <- df_heat_fecal %>% filter(ExperimentID == exp)

  p_heat <- ggplot(df_sub, aes(x = Indiv, y = Locus, fill = depth_bin)) +
    geom_tile() +
    scale_fill_manual(name = "Avg Depth (binned)", values = depth_colors, drop = FALSE) +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(title = paste("Heatmap — Experiment", exp))

  p_hist <- ggplot(df_sub, aes(x = TotalDepth)) +
    geom_histogram(bins = 40) +
    theme_bw() +
    labs(title = paste("Depth Histogram — Experiment", exp))

  list(
    experiment = exp,
    heatmap = p_heat,
    histogram = p_hist
  )
})

for (i in 1:length(unique(df_heat_fecal$ExperimentID))) {
  cat("### Experiment:", plots[[i]]$experiment, "\n")
  print(plots[[i]]$heatmap)
  print(plots[[i]]$histogram)
}

```
## Plot summary stats by plate
```{r}
#these still include PCRneg because developed from df_filtered
sample_df <- df_heat_fecal %>%
  group_by(Indiv, ExperimentID) %>%
  summarise(avg_depth= mean(TotalDepth, na.rm = TRUE),
            loci_GT= sum(is.finite(TotalDepth)),
            total_loci= n_distinct(Locus),
            prop_GT= loci_GT/total_loci,
            missing= sum(is.na(TotalDepth))/total_loci
            )%>%
  ungroup()%>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 10, 25, 50, 100,Inf),
    labels = c("0–10", "11–25", "26–50", "51–100", ">100"),
    right = TRUE
  ))

ggplot(sample_df, aes(x = ExperimentID, y = missing)) +
  geom_jitter(aes(color = depth_bin), width = 0.2, size = 2, alpha = 0.7) +
  scale_fill_manual(
    name = "Avg Depth (binned)",
    values = depth_colors,
    drop = FALSE)+  # keeps all bins in legend even if not present
  labs(
    x = "Plate",
    y = "Missingness (proportion)",
    color = "Mean depth",
    title = "Missingness per plate, colored by sample mean depth"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
```
## Scoring matching genotypes microhaplotypes to calculate mismatch rate
```{r}
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
depth_wide_fecal <- df_heat_fecal %>%
  select(Indiv, Locus, TotalDepth) %>%
  pivot_wider(names_from = Locus, values_from = TotalDepth, names_prefix = "Depth_")

# genotype_matrix_noneg already has IndivID as key
genotype_with_depth <- genotype_matrix_fecal %>%
  left_join(depth_wide_fecal, by = "Indiv")  

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
```
Summarise mismatch results
```{r}
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
```
We selected depth = 15, missingess = 0.5
```{r}

# Apply filtering scheme so that you can look at individual heterozygosity
genotype_filtered_15_0.5_fecal <- genotype_with_depth %>%
      mutate(across(starts_with("Clu_"), 
                ~ {
                  depth_col <- paste0("Depth_", cur_column())
                  # If depth is below threshold, set genotype to NA
                  ifelse(get(depth_col) < 15, NA, .x)
                }))

    
  genotype_filtered_15_0.5_fecal <- genotype_filtered_15_0.5_fecal %>%
      rowwise() %>%
      mutate(
        total_loci = length(c_across(starts_with("Clu_"))),        # total loci
        loci_missing = sum(is.na(c_across(starts_with("Clu_")))),  # count missing
        prop_missing = loci_missing / total_loci
      ) %>%
      ungroup() %>%
      filter(prop_missing <= 0.5)

    # Convert matrix object to adegenet's genlind format for microhaplotype

  gen_data_fecal<-genotype_filtered_15_0.5_fecal %>%
    select(-LabID, -Platei7, -Well, -ExperimentID, -Replicate,-total_loci, -loci_missing, -prop_missing,-Type,
           -starts_with("Depth_"))
 
ind.names<-gen_data_fecal$Indiv
gen_data_fecal$Indiv<- NULL

#updated with filtering scheme only 1000 samples
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

Ho <- rowSums(het_matrix) / rowSums(non_missing)

qc_table <- data.frame(
  Individual = rownames(geno),
  Ho = Ho,
  Heterozygous_Loci = num_het,
  Genotyped_Loci = num_genotyped,
  Missing_Rate = missing_rate
)

qc_table <- qc_table[order(qc_table$Ho), ] 

####adegenet way####
# split genind into one genind per locus
loci_list <- seploc(genind_fecal)

# determine heterozygosity per locus per individual
het_matrix <- sapply(loci_list, function(loc){

  g <- tab(loc, NA.method = "asis")

  # heterozygous if two allele columns == 1
  rowSums(g == 1, na.rm = TRUE) == 2
})

# determine loci with data
non_missing <- sapply(loci_list, function(loc){

  g <- tab(loc, NA.method = "asis")

  rowSums(!is.na(g)) > 0
})

# observed heterozygosity per individual
Ho_ade <- rowSums(het_matrix) / rowSums(non_missing)

# High het and high replicate mismatch is a good indication of contamination
#also sort out Ho <0.2, we had some that had abnormally low
# Would need to run mixed sample experiments to validate a true cutoff
# A cluster of high het in wells on a plate are a good indication of sample mixing in the lab
# Looking at a sample plate map and replicate error and heterozygosity by sample can help sort this out
```

```{r}
#already filtered for missingness (0.5) and depth (15)
genotype_filtered_15_0.5_fecal

High_het <- qc_table %>% filter(Ho > 0.65 | Ho<0.2)
High_rep_error <- pair_errors %>% 
  filter(Depth == 15, Missing == 0.50) %>%
  filter(Mismatch_Rate > 0.05)

Remove_samples_fecal <- unique(c(High_het$Individual, High_rep_error$Sample1, High_rep_error$Sample2))


#filter out samples that failed by het and mismatch
genotype_with_depth_fecal<- genotype_with_depth %>%
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

write.csv(clean_fecal_delomas_pipeline, file = "./outputs/fecal_all_genos_long_5-22-26.csv")  
write.csv(alleles_by_locus, file = "./outputs/microhap_allele_key_5-22-26.csv")

####use genotype_filtered_final_fecal
  #remove diagnostic loci for matching among biological replicates and harvest samples
  diag_loc<-read.csv("./inputs/neutral_diagnostic_loci_new.csv", header = TRUE)

diag_loc<-diag_loc %>%
  mutate(plain_locus = str_extract(OG_locus_ID, ".*(?=_[^_]*$)"))

remove_loci<- diag_loc$plain_locus

locinames<-locinames[!locinames %in% remove_loci]  
  
gen_data_fecal<-genotype_filtered_final_fecal %>%
    select(-LabID, -Platei7, -Well, -ExperimentID, -Replicate,-total_loci, -loci_missing,
           -prop_missing,-Type,-mean_depth, -starts_with("Depth_"))%>%
    select(Indiv, any_of(locinames))
    
write.csv(gen_data_fecal, file = "./outputs/gen_data_fecal_19Mar2026.csv")

ind.names_fecal<-gen_data_fecal$Indiv
gen_data_fecal$Indiv<- NULL

#updated final
genind_fecal <- df2genind(
  gen_data_fecal,
  ind.names = ind.names_fecal,
  ploidy = 2,
  sep = "",
  NA.char = "NA"
)

#Filtering IDFG data
  # Long to wide: one row per individual
genotype_matrix_IDFG <- df_encoded_IDFG %>%
  select(Indiv, Locus,Type, Genotype_Code) %>%
  pivot_wider(
    names_from = Locus,
    values_from = Genotype_Code
  ) %>%
  mutate(
    Indiv_long  = Indiv,
    Indiv_short = str_extract(Indiv, "Clu.*|NTC.*"),
    Indiv_short = if_else(is.na(Indiv_short), Indiv, Indiv_short)
  ) %>%
  select(Indiv, Indiv_long, Indiv_short, everything())  

genotype_missing <- genotype_matrix_IDFG %>% 
  filter(if_all(5:203, is.na))   #rows where all genotypes are NA

genotype_missing <- genotype_missing$Indiv #199 individuals when we only have 200 loci

#remove individuals with missing genotype and add totaldepth for depth filtering    
snp_df_filt_IDFG <- df_encoded_IDFG %>%
      filter(!Indiv %in% genotype_missing) %>% 
  mutate( #add total depth information
    TotalDepth = case_when(
      is.na(Allele1) | is.na(Allele2) |
        is.na(Allele1_count) | is.na(Allele2_count) ~ NA_real_,  # missing data
      
      Allele1 == Allele2 ~ Allele1_count,                      # homozygous
      
      TRUE ~ Allele1_count + Allele2_count                     # heterozygous
    )
  )%>%
  select(Indiv, Locus, Type, Allele1, Allele2, Allele1_count, Allele2_count,
         TotalDepth, everything())

#replace genotype matrix
genotype_matrix_IDFG <- snp_df_filt_IDFG %>%
  select(Indiv, Locus, Type, Genotype_Code) %>%
  pivot_wider(
    names_from = Locus,
    values_from = Genotype_Code
  ) %>%
  mutate(
    Indiv_long  = Indiv,
    Indiv_short = str_extract(Indiv, "Clu.*|NTC.*"),
    Indiv_short = if_else(is.na(Indiv_short), Indiv, Indiv_short)
  ) %>%
  select(Indiv, Indiv_long, Indiv_short, everything()) 

# Pivot snp_df_filt to wide format to match genotype_matrix
depth_wide_IDFG <- snp_df_filt_IDFG %>%
  select(Indiv, Locus, TotalDepth) %>%
  pivot_wider(names_from = Locus, values_from = TotalDepth, names_prefix = "Depth_")

# join with genotype_matrix_IDFG to filter by depth
genotype_with_depth_IDFG <- genotype_matrix_IDFG %>%
  left_join(depth_wide_IDFG, by = "Indiv")  

# Apply filtering scheme from scat
genotype_filtered_IDFG_final <- genotype_with_depth_IDFG %>%
      mutate(across(starts_with("Clu_"), 
                ~ {
                  depth_col <- paste0("Depth_", cur_column())
                  # If depth is below threshold, set genotype to NA
                  ifelse(get(depth_col) < 15, NA, .x)
                }))

    
  genotype_filtered_IDFG_final <- genotype_filtered_IDFG_final %>%
      rowwise() %>%
      mutate(
        total_loci = length(c_across(starts_with("Clu_"))),        # total loci
        loci_missing = sum(is.na(c_across(starts_with("Clu_")))),  # count missing
        prop_missing = loci_missing / total_loci
      ) %>%
      ungroup() %>%
      filter(prop_missing <= 0.5)

gen_data_IDFG<-genotype_filtered_IDFG_final %>%
    select(-Indiv_long, -Indiv_short,-total_loci, -loci_missing, -prop_missing,-Type,
           -starts_with("Depth_"))%>%
    select(Indiv, any_of(locinames))

ind.names_IDFG<-gen_data_IDFG$Indiv
gen_data_IDFG$Indiv<- NULL

genind_IDFG <- df2genind(
  gen_data_IDFG,
  ind.names = ind.names_IDFG,
  ploidy = 2,
  sep = "",
  NA.char = "NA"
)

#merge with fecal genind
combined <- repool(genind_fecal, genind_IDFG)

```
Matching among all samples (scat and harvest)
# Assume your data is in a genind object called `gi`
```{r}
gi<-combined
n_ind <- nInd(gi)

#extract allele count matrix
X <- gi@tab
loci <- rep(locNames(gi), gi@loc.n.all)   # locus name for each allele column
alleles <- colnames(gi@tab)               # e.g. "Clu_10.AC"

# Initialize result matrices
allele_matches <- matrix(0, n_ind, n_ind)
shared_loci <- matrix(0, n_ind, n_ind)
rownames(allele_matches) <- indNames(gi)
colnames(allele_matches) <- indNames(gi)
rownames(shared_loci) <- indNames(gi)
colnames(shared_loci) <- indNames(gi)

# Helper: extract alleles for one individual at one locus

get_alleles <- function(i, loc) { 
  cols <- which(loci == loc) 
  counts <- X[i, cols] 
  # If ANY allele count is NA → locus is missing return NA flag
  if (any(is.na(counts))) return(NA) 
  # otherwise return alleles for that individual at that locus 
  rep(alleles[cols], counts) }

# Loop through pairs of individuals: old without progress prints
for (i in 1:(n_ind - 1)) {
  cat("Starting individual", i, "of", n_ind, "\n")
  
  for (j in (i + 1):n_ind) {
    #cat("  Pair:", i, "-", j, "\n")
    
    for (loc in locNames(gi)) {

      a1 <- get_alleles(i, loc)
      a2 <- get_alleles(j, loc)

      # skip locus if either is NA (missing)
      if (is.na(a1)[1] || is.na(a2)[1]) next 

      shared_loci[i, j] <- shared_loci[i, j] + 1

      # genotype matching rule:
      # each allele is listed once per copy (diploid → 2)
      # sort alleles to ignore order
      a1_sorted <- sort(a1)
      a2_sorted <- sort(a2)
      
      # get 2 it perfect match, get 1 if half match, get 0 if no match
      if (identical(a1_sorted, a2_sorted)) {
          match<- 2
          
          } else if(length(intersect(a1_sorted, a2_sorted)) > 0){
            match<-1
          } else {
            match <-0
          }
      
      allele_matches[i, j] <- allele_matches[i, j] + match
    }

    allele_matches[j, i] <- allele_matches[i, j]
    shared_loci[j, i] <- shared_loci[i, j]
    
    #cat("  Finished pair", i, "-", j, 
     #   " | shared loci =", shared_loci[i,j], 
      #  " | allele matches =", allele_matches[i,j], "\n")
  }
  cat("Finished individual", i, "of", n_ind, "\n")
}


##load in data
#shared_loci<-readRDS("./outputs/pairwise/shared_loci.rds")
#allele_matches<-readRDS("./outputs/pairwise/allele_matches.rds")

# Convert to long format
shared_long  <- melt(shared_loci, varnames=c("ind1", "ind2"), value.name="loci_shared") %>%
    mutate(
    ind1 = as.character(ind1),
    ind2 = as.character(ind2)
  ) %>%
  filter(ind1 > ind2)

allele_long  <- melt(allele_matches, varnames=c("ind1", "ind2"), value.name="allele_matches") %>%
  mutate(
    ind1 = as.character(ind1),
    ind2 = as.character(ind2)
  ) %>%
  filter(ind1 > ind2)

matching_df <- left_join(shared_long, allele_long)
matching_df$prop_mismatching <- 1 - (matching_df$allele_matches / (2* matching_df$loci_shared)) 

p1 <- ggplot(matching_df, aes(x=prop_mismatching)) + 
  geom_histogram(binwidth=0.01, fill="blue", alpha=0.7) +
  theme_minimal() +
  geom_vline(aes(xintercept = 0.015), color = "black") +
  labs(title="Distribution of Proportion Mismatched Alleles between Samples", x="Proportion Mismatched Alleles", y="Frequency")

matching_tail_df <- matching_df %>% filter(prop_mismatching < 0.2)

p2 <- ggplot(matching_tail_df, aes(x=prop_mismatching)) + 
  geom_histogram(binwidth=0.005, fill="blue", alpha=0.7) +
  geom_vline(aes(xintercept = 0.015), color = "black") +
  theme_minimal() +
  labs(x="Proportion Mismatched Alleles", y="Frequency")

grid.arrange(p1, p2, nrow = 2)
```
# Matches
## Filter and match
Only allow matches between samples with at least X number genotyped loci shared between them and a mismatch rate of X
```{r}
#import known matches from microsatellites
known_matches<-read.csv("./inputs/tissue-scat-samples_for-validation-and-error-calculations.csv", header = TRUE)

known_matches_long<- known_matches %>%
  pivot_longer(cols = everything(), names_to = "Scat_num", values_to = "Sample") %>%
  filter(Sample != "") %>%                      # remove blank strings
  separate(Sample, into = c("UI_ID", "Scat_ID"), sep = "_") %>%
  group_by(UI_ID) %>%
  mutate(Replicate = row_number()) %>%            # arbitrary count per individual
  ungroup() %>%
  select(UI_ID, Scat_ID, Replicate)

matching_df_labID<-matching_df %>%
  mutate(
    IndivID1 = str_extract(ind1, "[^/]+$"),
    IndivID1 = sub("\\..*$", "", IndivID1),
    parts1 = str_split(IndivID1, "_"),
    parts1 = map(parts1, ~ .x[!str_detect(.x, "^redo")]))%>%
   mutate(
    IndivID2 = str_extract(ind2, "[^/]+$"),
    IndivID2 = sub("\\..*$", "", IndivID2),
    parts2 = str_split(IndivID2, "_"),
    parts2 = map(parts2, ~ .x[!str_detect(.x, "^redo")]))%>%
  mutate(
    LabID1   = map_chr(parts1, ~ .x[length(.x)]),
    LabID2   = map_chr(parts2, ~ .x[length(.x)]))%>%
  select(-parts1, -parts2, -IndivID2, -IndivID1)

known_pairs<- known_matches_long%>%
  select(UI_ID, LabID = Scat_ID) %>%              # rename for clarity
  group_by(UI_ID) %>%
  summarise(LabID_pairs = list(combn(LabID, 2, simplify = FALSE))) %>%
  unnest(LabID_pairs) %>%
  mutate(
    LabID1 = map_chr(LabID_pairs, 1),
    LabID2 = map_chr(LabID_pairs, 2)
  ) %>%
  select(-LabID_pairs)

# Now join with matching_df to find empirical pairs
empirical_matches <- matching_df_labID %>%
  inner_join(known_pairs, by = c("LabID1", "LabID2"))

empirical_matches_rev <- matching_df_labID %>%
  inner_join(known_pairs, by = c("LabID1" = "LabID2", "LabID2" = "LabID1"))

empirical_matches_all <- bind_rows(empirical_matches, empirical_matches_rev)

# Now you can look at the distributions
summary(empirical_matches_all$prop_mismatching)
summary(empirical_matches_all$loci_shared)

# Optional: plot distributions to help pick cutoffs
ggplot(empirical_matches_all, aes(x = prop_mismatching)) + 
  geom_histogram(binwidth = 0.005) +
  theme_minimal()

ggplot(empirical_matches_all, aes(x = loci_shared)) + 
  geom_histogram(binwidth = 1) +
  theme_minimal()

matching_w_enough_loci <- matching_df %>% filter(loci_shared >= 50) # Only allow comparisons with at least X number loci shared (based on PID)

matching_df %>% filter(loci_shared >= 50)

min(matching_df$loci_shared)

#max prop mismatch was 0.013 so I'm going to set it at a conservative 0.05
matches_propmismatch <- matching_w_enough_loci %>% filter(prop_mismatching < 0.05) # A match must have mismatch of less than X% 
```
## Cluster matching samples to individual IDs with igraph
```{r}
df <- matches_propmismatch[,1:2]

# Step 1: Create a graph from the pairwise matches
g <- graph_from_data_frame(df, directed = FALSE)

# Step 2: Find connected components
components <- components(g)

# Step 3: Create a data frame with groups
clusters_df <- data.frame(
  individual = names(components$membership),
  group = components$membership
)

# Optional: group them together for viewing
grouped <- clusters_df %>%
  group_by(group) %>%
  summarise(members = paste(sort(individual), collapse = ", "))

# Samples matching other samples
head(grouped)

write.csv(grouped, "./outputs/pairwise/grouped_individuals_20Mar2026.csv")
write.csv(matches_propmismatch, "./outputs/pairwise/matches_promismatch_20Mar2026.csv")
```
## Format matching samples with samples that didn't have a match
```{r}
# Add all non-matching samples
all_samples <- data.frame(individual = as.character(unique(c(matching_w_enough_loci$ind1, matching_w_enough_loci$ind2))))


grouped_long<- grouped %>%
  separate_rows(members, sep = ",") %>%
  mutate(
    members = str_trim(members),
    group = as.character(group)   # 👈 force to character
    )

#label sample types: tissue or scat; IDGW/ScatOpt4 = scat; CluIDFG = tissue; CluVAR = dog so exclude
grouped_long <- grouped_long %>%
  mutate(
    sample_type = case_when(
      str_detect(members, "CluVAR") ~ "exclude",
      str_detect(members, "CluIDFG") ~ "tissue",
      str_detect(members, "IDGW|ScatOpt4") ~ "scat",
      TRUE ~ "other"
    )
  )

all_samples2 <- all_samples %>%
  rename(members = individual) %>%
  mutate(
    sample_type = case_when(
      str_detect(members, "CluVAR") ~ "exclude",
      str_detect(members, "CluIDFG") ~ "tissue",
      str_detect(members, "IDGW|ScatOpt4") ~ "scat",
      TRUE ~ "other"
    )
  )

ungrouped <- all_samples2 %>%
  filter(!members %in% grouped_long$members) %>%
  mutate(group = paste0("singleton_", row_number()))%>%
  select(group, members, sample_type)

combined_unique <- bind_rows(
  grouped_long,
  ungrouped
)

combined_unique <- combined_unique %>%
  filter(sample_type != "exclude")%>%
  filter(sample_type != "other")

individuals<-combined_unique %>%
  group_by(group) %>%
  summarise(
    has_scat = any(sample_type == "scat"),
    has_tissue = any(sample_type == "tissue"),
    n_scat_samples = sum(sample_type == "scat"),
    n_tissue_samples = sum(sample_type == "tissue"),
    .groups = "drop"
  )
#How many unique individuals with scat
n_scat_individuals <- individuals %>%
  filter(has_scat) %>%
  nrow()

n_scat_individuals_df <- individuals %>%
  filter(has_scat)
#check:
sum(n_scat_individuals_df$has_scat)
#[1] 264
sum(n_scat_individuals_df$has_tissue)
#[1] 71

#how many scat individuals match and tissue sample
n_scat_with_tissue<- individuals %>%
  filter(has_scat & has_tissue) %>%
  nrow()

###scat recaptures
scat_summary <- combined_unique %>%
  filter(sample_type == "scat") %>%
  group_by(group) %>%
  summarise(
    n_scat = n(),
    recaptures = n_scat - 1,
    .groups = "drop"
  )
mean_recaptures <- scat_summary %>%
  summarise(mean_recaptures = mean(recaptures),
            sd=sd(recaptures))

scat_count<- scat_summary %>%
  count(n_scat)%>%
  arrange(n_scat)

recap_counts<-scat_summary %>%
  count(recaptures) %>%
  arrange(recaptures)

prop_singletons<-scat_summary %>%
  summarise(
    prop_singletons = mean(n_scat == 1)
  )


capture_plot<-ggplot(scat_count, aes(x = n_scat, y = n)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, max(scat_count$n_scat), by = 1)) +
  labs(
    title = "Capture Frequency Distribution",
    x = "Number of Scats",
    y = "Number of Individuals"
  ) +
  theme_minimal(base_size = 16)+
  theme(
    # ❌ remove vertical gridlines (x direction)
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # ✅ keep horizontal gridlines (y direction)
    panel.grid.major.y = element_line(),
    panel.grid.minor.y = element_line(),
    
    # 🔤 ensure axis text is at least size 12
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

######old code########
df <- left_join(all_samples, clusters_df) %>% arrange(group)

# Step 1: Find the maximum group number
max_group <- max(df$group, na.rm = TRUE)

# Step 2: Assign new groups to NAs
df_final <- df %>%
  mutate(group = if_else(
    is.na(group),
    max_group + cumsum(is.na(group)),
    as.numeric(group)  # ensure group stays numeric
  ))

snp_inds_df <- df_final %>%
  group_by(group) %>%
    mutate(Individual_snp = sprintf("IDGW_%03d", group)) %>% 
    arrange(group)
```
## Plot distribution of # of matches
```{r}
matches_dist_df <- data.frame(table(snp_inds_df$Individual_snp))

colnames(matches_dist_df) <- c("Individual", "Samples Matching")

ggplot(matches_dist_df, aes(x = `Samples Matching`)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Matches",
       x = "Number of Samples Matching One Individual",
       y = "Count") +
  theme_minimal()
```                     
