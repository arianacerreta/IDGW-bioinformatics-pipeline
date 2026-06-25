#final set presence absence thresholds

library(tidyverse)

### ScatOpt4
#NOTE Delomas output .txt does not include NA for uncalled sexes and the read.table will balk.
#you will need to edit the .txt to include NA for those columns for this code to work
pres_abs<-as_tibble(read.table("./inputs/pres_abs_genotypes.txt", header = TRUE))

known_sex<-as_tibble(read.csv("./inputs/ScatOpt4KnownSex.csv", header = TRUE))

known_sex_noNA<-known_sex%>%
  filter(scat_sex != "NA")%>%
  filter(scat_sex != "U")

#reformat
pres_abs_formatted<- pres_abs %>%
  mutate(IndivID = str_extract(Indiv, "[^/]+$")) %>%
  mutate(LabID = str_extract(IndivID, "(?<=_)[^_]+(?=\\.1$)")) %>%
  mutate(Platei7 = str_extract(IndivID, "^[^_]+")) %>%
  mutate(ExperimentID = str_split(IndivID, "_", simplify = TRUE)[,3]) %>%
  select(-Indiv) %>%
  select(LabID, everything()) %>%
  left_join(known_sex, by = c("LabID" = "Sample")) %>%
  select(LabID, Genotype, scat_sex, tissue_field_otherfecal_sex, everything())

pres_abs_known<-pres_abs_formatted%>%
  filter(scat_sex != "NA")%>%
  filter(scat_sex != "U")

#set INF to OtherAligned Reads; this essentially gives values as if there was 1 read for the SRY locus
pres_abs_known2 <- pres_abs_known %>%
  mutate(
    Ratio = ifelse(is.infinite(Ratio), OtherAlignedReads, Ratio)
  )

pres_abs_summ<- pres_abs_known2 %>%
  group_by(scat_sex) %>%
  summarise(
    minRatio = min(Ratio, na.rm = TRUE), # minimum ratio
    maxRatio = max(Ratio, na.rm = TRUE), # max ratio
    meanRatio = mean(Ratio, na.rm = TRUE), #mean
    minReads = min(OtherAlignedReads, na.rm = TRUE),
    maxReads = max(OtherAlignedReads, na.rm = TRUE
    )
  ) %>%
  ungroup()

ratio_summary <- pres_abs_known2 %>%
  group_by(scat_sex) %>%
  summarize(
    min = min(Ratio),
    q05 = quantile(Ratio, 0.05),
    q10 = quantile(Ratio, 0.10),
    median = median(Ratio),
    q90 = quantile(Ratio, 0.90),
    q95 = quantile(Ratio, 0.95),
    max = max(Ratio)
  )

#evaluate
thresholds <- expand_grid(
  MaxRatioPres = seq(0, 22600, by = 100),
  MinRatioAbs  = seq(0, 3600, by = 100)
) %>% 
  filter(MinRatioAbs > MaxRatioPres)

results <- thresholds %>%
  mutate(
    misclass = map2_int(MaxRatioPres, MinRatioAbs,
                        ~ sum({
                          call <- case_when(
                            pres_abs_known2$Ratio < .x ~ "M",
                            pres_abs_known2$Ratio > .y ~ "F",
                            TRUE ~ "missing"
                          )
                          call != pres_abs_known2$scat_sex & call != "missing"
                        })),
    missing = map2_int(MaxRatioPres, MinRatioAbs,
                       ~ sum({
                         call <- case_when(
                           pres_abs_known2$Ratio < .x ~ "M",
                           pres_abs_known2$Ratio > .y ~ "F",
                           TRUE ~ "missing"
                         )
                         call == "missing"
                       }))
  )

results2 <- results %>%
  mutate(
    n = nrow(pres_abs_known2),
    accuracy = 1 - misclass / n,
    missing_rate = missing / n,
    # Only evaluate accuracy on non-missing samples
    effective_accuracy = 1 - misclass / (n - missing),
    # Combined score — you can tune this (just an example)
    score = effective_accuracy * (1 - missing_rate)
  )

#first misclassification rank, missingness second
min_mis <- min(results2$misclass)
results_mis <- results2 %>%
  filter(misclass == min_mis)

min_missing <- min(results_mis$missing)

results_mis_miss <- results_mis %>%
  filter(missing == min_missing)

best_combo <- results_mis_miss %>%
  arrange(MaxRatioPres, desc(MinRatioAbs)) %>%
  slice(1)

MRA<-best_combo$MinRatioAbs
MRP<-best_combo$MaxRatioPres

pres_abs_test <- pres_abs_known %>%
  mutate(
    call = case_when(
      Ratio < MRP ~ "M",
      Ratio > MRA ~ "F",
      TRUE ~ "missing"
    )
  )%>%
  select(LabID, Genotype,scat_sex, call, ExperimentID, Ratio, LocusReads, OtherAlignedReads)
table(pres_abs_test$call, pres_abs_test$scat_sex)

#visualize

ggplot(pres_abs_known2, aes(x=Ratio , fill = scat_sex))+
  geom_histogram(alpha = 0.4, bins = 45)+
  labs(fill = "Microsatellite Sex", y = "Count")+
  geom_vline(xintercept = MRA, color = "#F8766D", size = 1)+
  geom_vline(xintercept = MRP, color = "#00BFC4", size = 1)+
  geom_segment(aes(x = 190, y = 6.2, xend = 130, yend = 6.2),
               arrow = arrow(length = unit(0.3,"cm"),type = "closed"), size= 1, color = "#00BFC4"
  ) +
  geom_segment(aes(x = 3700, y = 6.2, xend = 5400, yend = 6.2),
               arrow = arrow(length = unit(0.3,"cm"),type = "closed"), size= 1, color = "#F8766D"
  ) +
  annotate("text", x= 100, y=6.2, label = "Male", color = "#00BFC4", size = 5)+
  annotate("text", x= 8000, y=6.2, label = "Female", color = "#F8766D", size = 5)+
  scale_x_log10()+
  theme_bw()+
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))
