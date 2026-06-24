library(tidyverse)

#ScatOpt2 (Optimization Srategy 2)
summ_stats<-as_tibble(read.table("../../ScatTrials/Scat_Trial_2/delomas-pipeline/genos_final_sry/summary_stats_fixed.txt", header = TRUE))
{
#reformat
summ_stats<- summ_stats %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(LabID = str_extract(IndivID, "(?<=_)[^_]+(?=\\.1$)")) %>%
  mutate(Platei7 = str_extract(IndivID, "^[^_]+")) %>%
  mutate(ExperimentID = str_split(IndivID, "_", simplify = TRUE)[,3])


#In Table 2
#avg. GTsuccess  per experiment and % on target (of trimmed fqs)
exp_summ<- summ_stats %>%
  group_by(ExperimentID) %>%
  summarise(
    avg_GT = mean(GenotypeSuccess, na.rm = TRUE), # avg. genotypesuccess
    on_target = sum(CountAlign, na.rm = TRUE) / sum(TotalReads, na.rm = TRUE), # proportion aligned
    mean_on_targ = mean(PropAlign, na.rm = TRUE), # mean on targ
    sd_on_targ = sd(PropAlign, na.rm = TRUE) #sd on targ
    ) %>%
  ungroup()

#TableS.1, ScatOpt2
sum(summ_stats$CountAlign) #number aligned after trimming
sum(summ_stats$TotalReads) #number reads after trimming, min lgenth 50bp
sum(summ_stats$CountAlign)/560260956 #prop. read aligned of total reads 560260956 from original .fastq from sequencing facility
sum(summ_stats$TotalReads)/560260956 # prop. reads remain after trim; total reads 560260956 from original .fastq from sequencing facility
sum(summ_stats$CountAlign)/sum(summ_stats$TotalReads) #proportion aligned of remaining timmed reads

#depth per locus
indiv_depth_locus<-as_tibble(read.table("../../ScatTrials/Scat_Trial_2/delomas-pipeline/genos_final_sry/align_by_ind_locus.txt", header = TRUE))
loci<-as_tibble(read.table("../../Genomes/wolfAmpRef_334/Clu334_pos_4.txt", header = TRUE))

all_indivs<- unique(indiv_depth_locus$Indiv)
all_loci<-unique(loci$locus)
all_loci<-c(all_loci,"Clu_sry")

full_grid<-tidyr::expand_grid(
  Indiv = all_indivs,
  Locus = all_loci
)

indiv_depth_full<-full_grid %>%
  left_join(indiv_depth_locus, by = c("Indiv", "Locus")) %>%
  mutate(CountAlign = replace_na(CountAlign,0))

indiv_depth_full<- indiv_depth_full %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(LabID = str_extract(IndivID, "(?<=_)[^_]+(?=\\.1$)")) %>%
  mutate(Platei7 = str_extract(IndivID, "^[^_]+")) %>%
  mutate(ExperimentID = str_split(IndivID, "_", simplify = TRUE)[,3])

locus_depth<- indiv_depth_full %>%
  group_by(ExperimentID, Locus) %>%
  summarise(avg_depth = mean(CountAlign, na.rm = TRUE))%>%
  ungroup()

locus_depth <- locus_depth %>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 10, 25, 50, 100,Inf),
    labels = c("0–10", "11–25", "26–50", "51–100", ">100"),
    right = TRUE
  ))

#for each indiv
indiv_depth <- indiv_depth_full %>%
  group_by(ExperimentID, LabID) %>%
  summarise(avg_depth= mean(CountAlign, na.rm = TRUE))%>%
  ungroup()

indiv_depth <- indiv_depth %>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 10, 25, 50, 100,Inf),
    labels = c("0–10", "11–25", "26–50", "51–100", ">100"),
    right = TRUE
  ))

lab_order <- indiv_depth %>%
  group_by(LabID) %>%
  summarise(mean_depth = mean(avg_depth, na.rm = TRUE)) %>%
  arrange(mean_depth) %>%   # lowest mean at top of factor levels
  pull(LabID)

indiv_depth <- indiv_depth %>%
  mutate(LabID = factor(LabID, levels = lab_order))

#Table 2 Depth
#for each experiment
experiment_summary <- indiv_depth_full %>%
  group_by(ExperimentID) %>%
  summarise(
    mean_depth_all = mean(CountAlign, na.rm = TRUE),
    median_depth   = median(CountAlign, na.rm = TRUE),
    sd_depth       = sd(CountAlign, na.rm = TRUE),
    n_indiv        = n_distinct(Indiv),
    n_loci         = n_distinct(Locus)
  )

}

#ScatOpt3 (Optimization Strategy 3 and 4)

summ_stats_3_3<-as_tibble(read.table("../../ScatTrials/ScatOpt3June2024/delomas-pipeline/genos_final_sry/summary_stats_fixed.txt", header = TRUE))

summ_stats_3_6exp<-as_tibble(read.table("../../ScatTrials/ScatOpt3June2024/delomas-pipeline-70/genos_final_sry/summary_stats_fixed.txt", header = TRUE))

summ_stats_opt3<- rbind(summ_stats_3_3,summ_stats_3_6exp)
{
#reformat
summ_stats_opt3<- summ_stats_opt3 %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(LabID = str_extract(IndivID, "(?<=_)[^_]+(?=\\.1$)")) %>%
  mutate(Platei7 = str_extract(IndivID, "^[^_]+")) %>%
  mutate(ExperimentID = str_split(IndivID, "_", simplify = TRUE)[,3])

#Table 3
#avg. GTsuccess  per experiment and % on target (of trimmed fqs)
exp_summ<- summ_stats_opt3 %>%
  group_by(ExperimentID) %>%
  summarise(
    avg_GT = mean(GenotypeSuccess, na.rm = TRUE), # avg. genotypesuccess
    on_target = sum(CountAlign, na.rm = TRUE) / sum(TotalReads, na.rm = TRUE) # proportion aligned
  ) %>%
  ungroup()

#Table S1
sum(summ_stats_opt3$CountAlign) #number aligned after trimming
sum(summ_stats_opt3$TotalReads) #number reads after trimming, min legenth 50bp
sum(summ_stats_opt3$CountAlign)/552954523 #prop. read aligned of total reads 552954523 from original .fastq from sequencing facility
sum(summ_stats_opt3$TotalReads)/552954523 # prop. reads remain after trim; total reads 552954523 from original .fastq from sequencing facility
sum(summ_stats_opt3$CountAlign)/sum(summ_stats_opt3$TotalReads) #proportion aligned of remaining trimmed reads

#only Optimization Strategy 3 (3Exp)
sum(summ_stats_3_3$CountAlign)
sum(summ_stats_3_3$TotalReads)
sum(summ_stats_3_3$CountAlign)/sum(summ_stats_3_3$TotalReads)

#only Optimization Strategy 4 (1-6Exp)
sum(summ_stats_3_6exp$CountAlign)
sum(summ_stats_3_6exp$TotalReads)
sum(summ_stats_3_6exp$CountAlign)/sum(summ_stats_3_6exp$TotalReads)

#depth per locus
##s8NEw, s8new2, s9low with 334 loci
indiv_depth_locus<-as_tibble(read.table("../../ScatTrials/ScatOpt3June2024/delomas-pipeline/genos_final_sry/align_by_ind_locus.txt", header = TRUE))
loci334<-as_tibble(read.table("../../Genomes/wolfAmpRef_334/Clu334_pos_4.txt", header = TRUE))

all_indivs<- unique(indiv_depth_locus$Indiv)
all_loci<-unique(loci334$locus)
all_loci<-c(all_loci,"Clu_sry")

full_grid<-tidyr::expand_grid(
  Indiv = all_indivs,
  Locus = all_loci
)

indiv_depth_full<-full_grid %>%
  left_join(indiv_depth_locus, by = c("Indiv", "Locus")) %>%
  mutate(CountAlign = replace_na(CountAlign,0))

indiv_depth_full<- indiv_depth_full %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(LabID = str_extract(IndivID, "(?<=_)[^_]+(?=\\.1$)")) %>%
  mutate(Platei7 = str_extract(IndivID, "^[^_]+")) %>%
  mutate(ExperimentID = str_split(IndivID, "_", simplify = TRUE)[,3])

locus_depth<- indiv_depth_full %>%
  group_by(ExperimentID, Locus) %>%
  summarise(avg_depth = mean(CountAlign, na.rm = TRUE))%>%
  ungroup()

locus_depth <- locus_depth %>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 10, 25, 50, 100,Inf),
    labels = c("0–10", "11–25", "26–50", "51–100", ">100"),
    right = TRUE
  ))

#for each indiv
indiv_depth <- indiv_depth_full %>%
  group_by(ExperimentID, LabID) %>%
  summarise(avg_depth= mean(CountAlign, na.rm = TRUE))%>%
  ungroup()

indiv_depth <- indiv_depth %>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 10, 25, 50, 100,Inf),
    labels = c("0–10", "11–25", "26–50", "51–100", ">100"),
    right = TRUE
  ))

lab_order <- indiv_depth %>%
  group_by(LabID) %>%
  summarise(mean_depth = mean(avg_depth, na.rm = TRUE)) %>%
  arrange(mean_depth) %>%   # lowest mean at top of factor levels
  pull(LabID)

indiv_depth <- indiv_depth %>%
  mutate(LabID = factor(LabID, levels = lab_order))

#Table 3
#for each experiment
experiment_summary <- indiv_depth_full %>%
  group_by(ExperimentID) %>%
  summarise(
    mean_depth_all = mean(CountAlign, na.rm = TRUE),
    median_depth   = median(CountAlign, na.rm = TRUE),
    sd_depth       = sd(CountAlign, na.rm = TRUE),
    n_indiv        = n_distinct(Indiv),
    n_loci         = n_distinct(Locus)
  )

#depth per locus
##six experiments 70 loci
indiv_depth_locus<-as_tibble(read.table("../../ScatTrials/ScatOpt3June2024/delomas-pipeline-70/genos_final_sry/align_by_ind_locus.txt", header = TRUE))
loci70<-as_tibble(read.table("../../Genomes/wolfAmpRef_70/Clu70_pos_4.txt", header = TRUE))

all_indivs<- unique(indiv_depth_locus$Indiv)
all_loci<-unique(loci70$locus)
all_loci<-c(all_loci,"Clu_sry")

full_grid<-tidyr::expand_grid(
  Indiv = all_indivs,
  Locus = all_loci
)

indiv_depth_full<-full_grid %>%
  left_join(indiv_depth_locus, by = c("Indiv", "Locus")) %>%
  mutate(CountAlign = replace_na(CountAlign,0))

indiv_depth_full<- indiv_depth_full %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(LabID = str_extract(IndivID, "(?<=_)[^_]+(?=\\.1$)")) %>%
  mutate(Platei7 = str_extract(IndivID, "^[^_]+")) %>%
  mutate(ExperimentID = str_split(IndivID, "_", simplify = TRUE)[,3])

locus_depth<- indiv_depth_full %>%
  group_by(ExperimentID, Locus) %>%
  summarise(avg_depth = mean(CountAlign, na.rm = TRUE))%>%
  ungroup()

locus_depth <- locus_depth %>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 25, 50, 100, 500 ,Inf),
    labels = c("0–25", "26–50", "51–100", "101-500", ">500"),
    right = TRUE
  ))

#for each indiv
indiv_depth <- indiv_depth_full %>%
  group_by(ExperimentID, LabID) %>%
  summarise(avg_depth= mean(CountAlign, na.rm = TRUE))%>%
  ungroup()

indiv_depth <- indiv_depth %>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 25, 50, 100, 500 ,Inf),
    labels = c("0–25", "26–50", "51–100", "101-500", ">500"),
    right = TRUE
  ))

lab_order <- indiv_depth %>%
  group_by(LabID) %>%
  summarise(mean_depth = mean(avg_depth, na.rm = TRUE)) %>%
  arrange(mean_depth) %>%   # lowest mean at top of factor levels
  pull(LabID)

indiv_depth <- indiv_depth %>%
  mutate(LabID = factor(LabID, levels = lab_order))

#Table 3
#for each experiment
experiment_summary <- indiv_depth_full %>%
  group_by(ExperimentID) %>%
  summarise(
    mean_depth_all = mean(CountAlign, na.rm = TRUE),
    median_depth   = median(CountAlign, na.rm = TRUE),
    sd_depth       = sd(CountAlign, na.rm = TRUE),
    n_indiv        = n_distinct(Indiv),
    n_loci         = n_distinct(Locus)
  )
}

# ScatOpt4 (Optimization Strategy 5) & FinalScat
summ_stats_4<-as_tibble(read.table("../../ScatTrials/ScatOpt4Oct2024/delomas-pipeline/genos_final_sry/summary_stats_fixed.txt", header = TRUE))
summ_stats_final<-as_tibble(read.table("../../FinalScat/delomas-pipeline/genos_final_sry/summary_stats_fixed.txt", header = TRUE))

summ_stats<-rbind(summ_stats_4, summ_stats_final)

#reformat
summ_stats<- summ_stats %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(LabID = str_extract(IndivID, "(?<=_)[^_]+(?=\\.1$)")) %>%
  mutate(Platei7 = str_extract(IndivID, "^[^_]+")) %>%
  mutate(ExperimentID = str_split(IndivID, "_", simplify = TRUE)[,3])

summ_stats<- summ_stats %>%
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


#avg. GTsuccess  per experiment and % on target (of trimmed fqs)
#remove PCRneg
summ_stats_samps<- summ_stats %>%
  filter(!str_detect(LabID, "PCR")) 

#Table 4
exp_summ<- summ_stats_samps %>%
  group_by(ExperimentID) %>%
  summarise(
    avg_GT = round(mean(GenotypeSuccess, na.rm = TRUE), 2), # avg. genotypesuccess
    on_target = round(sum(CountAlign, na.rm = TRUE) / sum(TotalReads, na.rm = TRUE), 2) # proportion aligned
  ) %>%
  ungroup()

# Scat Opt4
summ_stats_scatopt4<- summ_stats %>%
  filter(ExperimentID == "ScatOpt4")

#Table S.1
sum(summ_stats_scatopt4$CountAlign) #number aligned after trimming
sum(summ_stats_scatopt4$TotalReads) #number reads after trimming, min lgenth 50bp
sum(summ_stats_scatopt4$CountAlign)/27994396 #prop. read aligned of total reads 27994396 from original .fastq from sequencing facility
sum(summ_stats_scatopt4$TotalReads)/27994396 # prop. reads remain after trim; total reads 27994396 from original .fastq from sequencing facility
sum(summ_stats_scatopt4$CountAlign)/sum(summ_stats_scatopt4$TotalReads) #proportion aligned of remaining timmed reads

#Final Scat
#8596 s-s4: 503664617
summ_stats_8596 <- summ_stats %>%
  filter(ExperimentID %in% c("S1_MID1", "S2_MID1", "S3_MID2", "S4_MID2"))

sum(summ_stats_8596$CountAlign) #number aligned after trimming
sum(summ_stats_8596$TotalReads) #number reads after trimming, min lgenth 50bp
sum(summ_stats_8596$CountAlign)/503664617 #prop. read aligned of total reads  from original .fastq from sequencing facility
sum(summ_stats_8596$TotalReads)/503664617 # prop. reads remain after trim; total reads  from original .fastq from sequencing facility
sum(summ_stats_8596$CountAlign)/sum(summ_stats_8596$TotalReads) #proportion aligned of remaining timmed reads

#9327: Adams_Waits_IDGW_S5_6_8_12_070925: 465664952
summ_stats_9327 <- summ_stats %>%
  filter(ExperimentID %in% c("S5_MID2", "S6_HIGH" , "S8_HIGH","S12_LOW"))

sum(summ_stats_9327$CountAlign) #number aligned after trimming
sum(summ_stats_9327$TotalReads) #number reads after trimming, min lgenth 50bp
sum(summ_stats_9327$CountAlign)/465664952 #prop. read aligned of total reads  from original .fastq from sequencing facility
sum(summ_stats_9327$TotalReads)/465664952 # prop. reads remain after trim; total reads  from original .fastq from sequencing facility
sum(summ_stats_9327$CountAlign)/sum(summ_stats_9327$TotalReads) #proportion aligned of remaining timmed reads


#9328: Adams_Waits_IDGW_S7_13_14_16_17_070925: 502755114

summ_stats_9328 <- summ_stats %>%
  filter(ExperimentID %in% c("S7_HIGH","S13_LOW"  ,"S14_LOW" ,"S16_MID1_MID2","S17_MID2_HIGH"))

sum(summ_stats_9328$CountAlign) #number aligned after trimming
sum(summ_stats_9328$TotalReads) #number reads after trimming, min lgenth 50bp
sum(summ_stats_9328$CountAlign)/502755114 #prop. read aligned of total reads  from original .fastq from sequencing facility
sum(summ_stats_9328$TotalReads)/502755114 # prop. reads remain after trim; total reads  from original .fastq from sequencing facility
sum(summ_stats_9328$CountAlign)/sum(summ_stats_9328$TotalReads) #proportion aligned of remaining timmed reads

#9329:Adams_Waits_IDGW_S9-11_S15_070925: 467129496

summ_stats_9329 <- summ_stats %>%
  filter(ExperimentID %in% c("S9_LOW", "S10_LOW" ,"S11_LOW" ,"S15_LOW"))

sum(summ_stats_9329$CountAlign) #number aligned after trimming
sum(summ_stats_9329$TotalReads) #number reads after trimming, min lgenth 50bp
sum(summ_stats_9329$CountAlign)/467129496 #prop. read aligned of total reads  from original .fastq from sequencing facility
sum(summ_stats_9329$TotalReads)/467129496 # prop. reads remain after trim; total reads  from original .fastq from sequencing facility
sum(summ_stats_9329$CountAlign)/sum(summ_stats_9329$TotalReads) #proportion aligned of remaining timmed reads

#depth per locus
indiv_depth_locus<-as_tibble(read.table("../../ScatTrials/ScatOpt4Oct2024/delomas-pipeline/genos_final_sry/align_by_ind_locus.txt", header = TRUE))
indiv_depth_locus_final<-as_tibble(read.table("../../FinalScat/delomas-pipeline/genos_final_sry/align_by_ind_locus.txt", header = TRUE))

indiv_depth_locus<- rbind(indiv_depth_locus, indiv_depth_locus_final)
loci<-as_tibble(read.table("../../Genomes/wolfAmpRef_200/Clu200_pos_4.txt", header = TRUE))

all_indivs<- unique(indiv_depth_locus$Indiv)
all_loci<-unique(loci$locus)
all_loci<-c(all_loci,"Clu_sry")

full_grid<-tidyr::expand_grid(
  Indiv = all_indivs,
  Locus = all_loci
)

indiv_depth_full<-full_grid %>%
  left_join(indiv_depth_locus, by = c("Indiv", "Locus")) %>%
  mutate(CountAlign = replace_na(CountAlign,0))

indiv_depth_full<- indiv_depth_full %>%
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

locus_depth<- indiv_depth_full %>%
  group_by(ExperimentID, Locus) %>%
  summarise(avg_depth = mean(CountAlign, na.rm = TRUE))%>%
  ungroup()

locus_depth <- locus_depth %>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 10, 25, 50, 100,Inf),
    labels = c("0–10", "11–25", "26–50", "51–100", ">100"),
    right = TRUE
  ))

#for each indiv
indiv_depth <- indiv_depth_full %>%
  group_by(LabID, ExperimentID) %>%
  summarise(avg_depth= mean(CountAlign, na.rm = TRUE))%>%
  ungroup()

indiv_depth <- indiv_depth %>%
  mutate(depth_bin = cut(
    avg_depth,
    breaks = c(-Inf, 10, 25, 50, 100,Inf),
    labels = c("0–10", "11–25", "26–50", "51–100", ">100"),
    right = TRUE
  ))

lab_order <- indiv_depth %>%
  group_by(LabID) %>%
  summarise(mean_depth = mean(avg_depth, na.rm = TRUE)) %>%
  arrange(mean_depth) %>%   # lowest mean at top of factor levels
  pull(LabID)

indiv_depth <- indiv_depth %>%
  mutate(LabID = factor(LabID, levels = lab_order)) 

#Table 4
#for each experiment
experiment_summary <- indiv_depth_full %>%
  filter(!str_detect(LabID, "PCR")) %>%
  group_by(ExperimentID) %>%
  summarise(
    mean_depth_all = mean(CountAlign, na.rm = TRUE),
    median_depth   = median(CountAlign, na.rm = TRUE),
    sd_depth       = sd(CountAlign, na.rm = TRUE),
    n_indiv        = n_distinct(IndivID),
    n_loci         = n_distinct(Locus)
  )
