### Final GWAdapt summ_stats ###
library(tidyverse) #version 2.0.0: dplyr 1.2.1, forcats 1.0.1, lubridate 1.9.5, purrr 1.2.2, readr 2.2.0, stringr 1.6.0, tibble 3.3.1, tidyr 1.3.2

#% ontarget stats
perc_ontarg <- read_tsv("./path/summ_stats_260421_CanFam3_1_allGWAdapt_allAdaptOptim_GATK/reads_on_target.tsv")

sum(perc_ontarg$Total_Reads) #total number of raw reads
# 634760495

sum(perc_ontarg$Reads_on_Target) #total number of reads aligned
# 423073425

#percentage
(sum(perc_ontarg$Reads_on_Target)/sum(perc_ontarg$Total_Reads))*100
#66.65088

hist_ontarget<-ggplot(perc_ontarg)+
  geom_histogram(aes(Percent_on_Target), na.rm = TRUE, fill = "gray", color = "black")+
  xlim(0,100)+
  geom_vline(xintercept = (sum(perc_ontarg$Reads_on_Target) / 
                             sum(perc_ontarg$Total_Reads))* 100, 
             color = "blue", linetype = "dashed", size = 1)+
  theme_minimal()+
  xlab("Percent on Target")+
  ylab("Count")
