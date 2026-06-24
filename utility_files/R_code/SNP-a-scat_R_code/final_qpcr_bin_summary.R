#comparison with qPCR bins and final scat
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

#read in qPCR results
qPCR<-read.csv("./inputs/qPCR_results.csv")

#read in starting Delomas pipeline microhaps
scatopt4 <- read.csv("./inputs/scatopt4/microhap_genotypes.csv", na.strings = "") #use microhap_genotypes_scatopt4.csv from Zenodo
finalscat <- read.csv("./inputs/finalscat/microhap_genotypes.csv", na.strings = "")#use microhap_genotypes_final_scat.csv from Zenodo

#combine and extract metadata from sample names
both<- scatopt4 %>%
  rbind(finalscat)
{
both<- both %>%
  mutate(
    IndivID = str_extract(Indiv, "[^/]+$"),
    IndivID = sub("\\..*$", "", IndivID),
    parts = str_split(IndivID, "_"),
    # NEW: remove "redo", "redo1", "redoXYZ", etc.
    parts = map(parts, ~ .x[!str_detect(.x, "^redo")])
  ) %>%
  mutate(
    Platei7 = map_chr(parts, ~ .x[1]),
    Well    = map_chr(parts, ~ .x[2]),
    LabID   = map_chr(parts, ~ .x[length(.x)]),
    
    # detect UI candidate (second from last) and whether it is a UI id
    ui_candidate = map_chr(parts, ~ if (length(.x) >= 2) .x[length(.x) - 1] else NA_character_),
    ui_present   = map_lgl(ui_candidate, ~ !is.na(.x) && str_starts(.x, "UI")),
    
    # start index: 3 normally, but if 3rd element == "IDGW" start at 4
    start_idx = map_int(parts, ~ if (length(.x) >= 3 && .x[3] == "IDGW") 4L else 3L),
    
    # end index: normally second-to-last (length-1), but if UI present we stop before UI (length-2)
    end_idx = map_int(parts, ~ {
      n <- length(.x)
      if (n < 3) return(NA_integer_)
      if (n >= 2 && str_starts(.x[n-1], "UI")) n - 2L else n - 1L
    }),
    
    # build ExperimentID safely (if end < start, produce NA)
    ExperimentID = map2_chr(parts, start_idx, ~ {
      p <- .x
      s <- .y
      n <- length(p)
      # compute end dynamically again to ensure consistency:
      e <- if (n >= 2 && str_starts(p[n-1], "UI")) n - 2L else n - 1L
      if (is.na(s) || e < s) NA_character_ else paste(p[s:e], collapse = "_")
    }),
    
    # UI_ID only if ui_present
    UI_ID = map2_chr(ui_candidate, ui_present, ~ if (.y) .x else NA_character_)
  ) %>%
  select(-parts, -ui_candidate, -ui_present, -start_idx, -end_idx, -Indiv) %>%
  select(IndivID, LabID, Platei7, ExperimentID, Well, UI_ID, everything())
  }

#filter out PCRnegatives
both<- both %>%
  filter(!str_detect(LabID, "PCR")) 

#clean the qPCR nomenclature
qPCR_clean<-qPCR %>%
  mutate(
    SampleClean = str_trim(str_remove(Sample.Name, "\\s*\\([^)]*\\)"))
  )

qPCR_unique<-unique(qPCR_clean$SampleClean)
length(qPCR_unique) #2283 samples were quantified on qPCR

sequenced_unique<-unique(both$LabID)
length(sequenced_unique) #790 samples were sequenced

same_id<-intersect(qPCR_unique,sequenced_unique)

#samples >0.08
qPCR_passed<-qPCR_clean%>%
  filter(Quantity.Mean > 0.08)

unique_passed_qpcr<-unique(qPCR_passed$SampleClean)
length(unique_passed_qpcr) #787

length(unique_passed_qpcr)/length(qPCR_unique) #0.3447219

passed_only <- setdiff(unique_passed_qpcr, sequenced_unique)

sequenced_only <- setdiff(sequenced_unique, unique_passed_qpcr)

diff_seq_qPCR_passed<-setdiff(unique_passed_qpcr,sequenced_unique)

both_passed_and_sequenced <- intersect(unique_passed_qpcr, sequenced_unique) 
length(both_passed_and_sequenced)#779

sequenced_only_in_qPCR <- sequenced_only[sequenced_only %in% qPCR_unique]

sequenced_only_not_in_qPCR <- setdiff(sequenced_only, qPCR_unique)

passed_not_sequenced <- setdiff(unique_passed_qpcr, sequenced_unique)
length(passed_not_sequenced)   # 8

sequenced_failed_qpcr <- setdiff(sequenced_unique, unique_passed_qpcr)
length(sequenced_failed_qpcr)  # 11

sequenced_failed_qpcr_df <- qPCR_clean %>%
  filter(SampleClean %in% sequenced_failed_qpcr) %>%
  select(SampleClean, Quantity.Mean) %>%
  mutate(category = "Sequenced_but_failed_qPCR")

passed_not_sequenced_df <- qPCR_clean %>%
  filter(SampleClean %in% passed_not_sequenced) %>%
  select(SampleClean, Quantity.Mean) %>%
  mutate(category = "Passed_but_not_sequenced")

bind_df <- bind_rows(sequenced_failed_qpcr_df, passed_not_sequenced_df)

# View results
bind_df

#final dataset
final_samples<-read.csv("./outputs/gen_data_fecal_2026-06-24.csv")

#extract info
{
final_samples<-final_samples%>%
  mutate(
    IndivID = str_extract(Indiv, "[^/]+$"),
    IndivID = sub("\\..*$", "", IndivID),
    parts = str_split(IndivID, "_"),
    # NEW: remove "redo", "redo1", "redoXYZ", etc.
    parts = map(parts, ~ .x[!str_detect(.x, "^redo")])
  ) %>%
  mutate(
    Platei7 = map_chr(parts, ~ .x[1]),
    Well    = map_chr(parts, ~ .x[2]),
    LabID   = map_chr(parts, ~ .x[length(.x)]),
    
    # detect UI candidate (second from last) and whether it is a UI id
    ui_candidate = map_chr(parts, ~ if (length(.x) >= 2) .x[length(.x) - 1] else NA_character_),
    ui_present   = map_lgl(ui_candidate, ~ !is.na(.x) && str_starts(.x, "UI")),
    
    # start index: 3 normally, but if 3rd element == "IDGW" start at 4
    start_idx = map_int(parts, ~ if (length(.x) >= 3 && .x[3] == "IDGW") 4L else 3L),
    
    # end index: normally second-to-last (length-1), but if UI present we stop before UI (length-2)
    end_idx = map_int(parts, ~ {
      n <- length(.x)
      if (n < 3) return(NA_integer_)
      if (n >= 2 && str_starts(.x[n-1], "UI")) n - 2L else n - 1L
    }),
    
    # build ExperimentID safely (if end < start, produce NA)
    ExperimentID = map2_chr(parts, start_idx, ~ {
      p <- .x
      s <- .y
      n <- length(p)
      # compute end dynamically again to ensure consistency:
      e <- if (n >= 2 && str_starts(p[n-1], "UI")) n - 2L else n - 1L
      if (is.na(s) || e < s) NA_character_ else paste(p[s:e], collapse = "_")
    }),
    
    # UI_ID only if ui_present
    UI_ID = map2_chr(ui_candidate, ui_present, ~ if (.y) .x else NA_character_)
  ) %>%
  select(-parts, -ui_candidate, -ui_present, -start_idx, -end_idx, -Indiv) %>%
  select(IndivID, LabID, Platei7, ExperimentID, Well, UI_ID, everything())
}

final_sample_names<-unique(final_samples$LabID)

#are any of the failed but sequenced samples in the final dataset
# Check which sequenced_failed_qpcr samples are in the final dataset
in_final <- intersect(sequenced_failed_qpcr, final_sample_names)
length(in_final) #2

#calculate genotype success for plotting and stuff
# Identify loci columns
loci_cols <- grep("^Clu_", colnames(final_samples), value = TRUE)

# Calculate % genotype success per individual
gt_success <- final_samples %>%
  rowwise() %>%
  mutate(
    n_loci = length(loci_cols),
    n_called = sum(!is.na(c_across(all_of(loci_cols)))),
    pct_success = n_called / n_loci * 100
  ) %>%
  ungroup()

# Keep the highest Quantity.Mean per sample
qPCR_max <- qPCR_clean %>%
  group_by(SampleClean) %>%
  summarise(Quantity.Mean = max(Quantity.Mean, na.rm = TRUE)) %>%
  ungroup()

bin_levels <- c("FAILED","LOW","MID1","MID2","HIGH")

qPCR_max <- qPCR_max %>%
  mutate(qpcr_bin = case_when(
    Quantity.Mean < 0.08 ~ "FAILED",
    Quantity.Mean >= 0.08 & Quantity.Mean < 0.2 ~ "LOW",
    Quantity.Mean >= 0.2 & Quantity.Mean < 0.4 ~ "MID1",
    Quantity.Mean >= 0.4 & Quantity.Mean < 0.9 ~ "MID2",
    Quantity.Mean >= 0.9 ~ "HIGH",
    TRUE ~ NA_character_
  )) %>%
  mutate(qpcr_bin = factor(qpcr_bin, levels = bin_levels))

qPCR_max_sequenced <- qPCR_max %>%
  filter(SampleClean %in% sequenced_unique)

# final_samples has LabID column
final_with_qpcr <- qPCR_max_sequenced %>%
  mutate(in_final = ifelse(SampleClean %in% final_samples$LabID, "Yes", "No"))

summary_stacked <- final_with_qpcr %>%
  group_by(qpcr_bin, in_final) %>%
  summarise(n = n()) %>%
  ungroup()

summary_stacked <- summary_stacked %>%
  mutate(qpcr_bin = factor(qpcr_bin, levels = bin_levels))

gt_with_qpcr <- gt_success %>%
  left_join(qPCR_max_sequenced %>% select(SampleClean, qpcr_bin),
            by = c("LabID" = "SampleClean")) %>%
  mutate(in_final = ifelse(LabID %in% final_samples$LabID, "Yes", "No"))

# Only include samples that made it to the final dataset
avg_gt_per_bin <- gt_with_qpcr %>%
  filter(LabID %in% final_samples$LabID) %>%
  group_by(qpcr_bin) %>%
  summarise(avg_gt_success = mean(pct_success, na.rm = TRUE))%>%
  rename(qpcr_bin= qpcr_bin)

summary_stacked <- summary_stacked %>%
  left_join(avg_gt_per_bin, by = "qpcr_bin")

totals <- summary_stacked %>%
  group_by(qpcr_bin) %>%
  summarise(total_n = sum(n)) %>%
  left_join(avg_gt_per_bin, by = "qpcr_bin")

sum(totals$total_n) #790

summary_bin_2 <- gt_with_qpcr %>%
  group_by(qpcr_bin) %>%
  summarise(
    n_samples = n(),
    avg_gt_success = mean(pct_success, na.rm = TRUE),
    median_gt_success = median(pct_success, na.rm = TRUE)
  ) %>%
  arrange(qpcr_bin)

sum(summary_bin_2$n_samples) #472

final_summary<- totals %>%
  left_join(summary_bin_2, by= "qpcr_bin")%>%
  rename(retained = n_samples) %>%
  mutate(perc.retained = (retained/total_n)*100)

#total perc. retained
sum(summary_bin_2$n_samples)/sum(totals$total_n) #0.5974684
#overall gt_success
mean(gt_success$pct_success) #93.59
