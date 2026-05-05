### k-locus analysis ###

#libraries
library(tidyverse) #version 2.0.0: dplyr 1.2.1, forcats 1.0.1, lubridate 1.9.5, purrr 1.2.2, readr 2.2.0, stringr 1.6.0, tibble 3.3.1, tidyr 1.3.2
library(vcfR)

#inputs from .txt produced in final_filter_remove-reps_GWAdapt.R
genotypes<-read.delim2("./path/CanFam3.1/adaptive_canfam3_1_gt_GATK_4-28-26.txt")

#metadata
AdaptOptim_meta<-read.csv("./inputs/AdaptOptim_WolfMeta_Age_Disp_Breed_Sex_Coat.csv")

#same metadata from dispersers analysis
Adapt_meta<-read.csv("../Dispersers/inputs/breeder_disperser_genpop_reference.csv")

#Known coat color key
#fix weird naming conventions
check <- AdaptOptim_meta %>%
  filter(!is.na(Coat))%>%
  select(UI_ID, Coat_optim = Coat) %>%
  full_join(
    Adapt_meta %>% select(UI_ID, Coat_meta = Coat),
    by = "UI_ID"
  )

check <- check %>%
  mutate(match = Coat_optim == Coat_meta)

coat_key<- check %>%
  filter(Coat_optim == Coat_meta | is.na(Coat_optim) | is.na(Coat_meta)) %>%
  mutate(Coat = coalesce(Coat_optim, Coat_meta)) %>%
  select(UI_ID, Coat)%>%
  mutate(Coat = if_else(Coat == "gray", "Gray", Coat))
 
  

#filter to just the k-locus; ALT (1) is the deletion based on .vcf
black_expected <- genotypes %>%
  filter(CHROM == "CM000016.3", POS == 58965448) %>% ####canFam3.1
  mutate(Genotype_Black = case_when(
    GT %in% c("0/1", "1/0", "1/1", "0|1","1|0","1|1") ~ TRUE,
    GT %in% c("0/0", "0|0") ~ FALSE,
    TRUE ~ NA  # handles missing or malformed GT
  )) %>%
  left_join(coat_key, by = "UI_ID") %>%
  mutate(Coat = if_else(is.na(Coat), "Unknown", Coat))


#has a genotype call at this locus
n_with_call_at_locus<- black_expected %>%
  filter(!is.na(Genotype_Black)) %>%
  nrow()

#849 individuals have a coat color call

#has a field coat color recorded and genotype expected to be black
n_black_geno_no_unknown_field<-black_expected %>%
  filter(!is.na(Genotype_Black)) %>%
  filter(Genotype_Black == TRUE) %>%
  filter(!str_to_lower(Coat) == "unknown")%>%
  nrow()

# 100 wolves with black genotype and a field coat color call

#wolves with a genotype indicating they should be black including unknowns
n_black_by_geno <-black_expected %>% 
  filter(Genotype_Black == TRUE) %>%
  nrow()

#196

#coat = black & genotype = black
n_correctly_black <- black_expected %>%
  filter(Genotype_Black == TRUE, str_to_lower(Coat) == "black") %>%
  nrow()

#81

#all wolves with genotypes and how many are expected to be black based on genotype
(n_black_by_geno/n_with_call_at_locus)*100 #23.09%

#no unknowns, of wolves that were genotypically black and have coat information, what proportion were called black in field and by geno
(n_correctly_black/n_black_geno_no_unknown_field)*100 #81%

# False positives: recorded as Black in field, but genotype says NOT black
n_black_from_field <- black_expected %>%
  filter(str_to_lower(Coat) == "black") %>%
  filter(!is.na(Genotype_Black)) #83

n_black_by_geno_from_field <-n_black_from_field %>%
  filter(Genotype_Black == TRUE) %>%
  nrow() #81

n_notblack_by_geno_from_field <-n_black_from_field %>%
  filter(Genotype_Black == FALSE) %>%
  nrow() #2

# Count how many came from each actual coat color category (excluding "Unknown")
black_mismatch_counts <- black_expected %>%
  filter(Genotype_Black == FALSE,
         str_to_lower(Coat) != "unknown") %>%
  filter(str_to_lower(Coat) == "black") %>%
  group_by(Coat) %>%
  summarise(n = n())


# Define non-black coat colors of interest
non_black_colors <- c("gray", "white", "brown", "black/tan", "gray/tan", "brown/gray", "tan")

# Count how many were called black incorrectly from those categories
false_positive_detail <- black_expected %>%
  filter(Genotype_Black == TRUE,
         str_to_lower(Coat) %in% non_black_colors) %>%
  mutate(Coat = str_to_title(Coat)) %>%  # Capitalize for display
  group_by(Coat) %>%
  summarise(False_Black_Calls = n()) %>%
  arrange(desc(False_Black_Calls))


# Collapse Coat into "black", "not_black", or NA (for unknowns)
black_summary <- black_expected %>%
  mutate(
    Coat_obs = case_when(
      str_to_lower(Coat) == "black" ~ "black",
      str_to_lower(Coat) == "unknown" ~ NA_character_,
      is.na(Coat) ~ NA_character_,
      TRUE ~ "not_black"
    ),
    Genotype_exp = case_when(
      Genotype_Black == TRUE ~ "black",
      Genotype_Black == FALSE ~ "not_black",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Coat_obs), !is.na(Genotype_exp)) %>%  # drop missing
  count(Genotype_exp, Coat_obs) %>%
  pivot_wider(names_from = Coat_obs,
              values_from = n,
              values_fill = 0)

# Add row totals
black_summary <- black_summary %>%
  mutate(Total = black + not_black)

# Add column totals
black_summary <- bind_rows(
  black_summary,
  summarise(black_summary, 
            Genotype_exp = "Total",
            black = sum(black),
            not_black = sum(not_black),
            Total = sum(Total))
)

black_summary

# Extract values from the confusion matrix
a <- black_summary$black[black_summary$Genotype_exp == "black"]       # True Black (TP)
b <- black_summary$not_black[black_summary$Genotype_exp == "black"]   # False Black (FP)
c <- black_summary$black[black_summary$Genotype_exp == "not_black"]   # False Not-Black (FN)
d <- black_summary$not_black[black_summary$Genotype_exp == "not_black"] # True Not-Black (TN)

N <- a + b + c + d

# Accuracy and diagnostic metrics
metrics <- tibble(
  Accuracy = (a + d) / N,
  Sensitivity = a / (a + c),   # P(Observed black | Genotype black)
  Specificity = d / (b + d),   # P(Observed not_black | Genotype not_black)
  Precision   = a / (a + b),   # P(Genotype black | Observed black)
  False_Positive_Rate = b / (b + d),
  False_Negative_Rate = c / (a + c)
)

metrics

chi<-chisq.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))

chi$observed
chi$expected
