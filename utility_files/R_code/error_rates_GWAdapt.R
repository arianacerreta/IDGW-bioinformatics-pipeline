### Code for error rates for GWAdapt
library(vcfR)
library(tidyverse)
library(ggplot2)
library(SNPfiltR)

#Load library
vcf<- read.vcfR("./outputs/CanFam3.1/filtered-by-neg_GWAdaptCanFam3.1.vcf") #vcf without a DP filter by PCR negatives; no longer includes PCR negatives

#panel ID (important if you are looking at several panels)
panel<-"GWAdapt"
genome<-"CanFam3.1"

#extract data to make a pop map
popmap_long<-data.frame(id=colnames(vcf@gt)[2:length(colnames(vcf@gt))]) 
vcfR <- vcf
## Scoring matching genotypes
#make a function
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

#run through filters of missingness
{filters<- seq(from = 0.1, to = 0.9, by = 0.1)
  depths<-c(30,20,10,5)
  snps<-seq(from = 0.70, to = 0.95, by = 0.05)
  miss_mismatch<-as.data.frame(matrix(ncol = 7, nrow =0))
  colnames(miss_mismatch)<-c("Depth","SNPmissing","MissingnessFilter","Sample1","Sample2","Missingness","Mismatch")
  locus_mismatch<-as.data.frame(matrix(ncol = 7, nrow =0))
  colnames(locus_mismatch)<-c("Depth","SNPmissing","MissingnessFilter","CHROM","POS","Missingness","Mismatch")
}

#tell the for loop which panel (this helps with naming conventions)
GWAdapt<- TRUE

for (j in depths) {
  vcfR_dp <- hard_filter(vcfR = vcfR, depth = j)
  
  for (h in snps) {
    vcfSNPs<- missing_by_snp(vcfR = vcfR_dp, cutoff = h)
    
    for (i in filters) {
      vcftemp <- missing_by_sample(vcfR = vcfSNPs, cutoff = i)
      gt <- extract.gt(vcftemp, element = "GT", as.numeric = FALSE) #this no longer has PCR negatives
      gt_df <- as.data.frame(gt)
      #fix column names
      colnames(gt_df) <- str_extract(colnames(gt_df), "i[0-9A-Za-z_-]+(?=_filtered\\.bam)")
      #append CHROM, POS
      gt_df$CHROM <- getCHROM(vcftemp)
      gt_df$POS <- as.numeric(getPOS(vcftemp))
      #move CHROM and POS to front
      gt_df <- gt_df %>% relocate(CHROM, POS)
      #pivot to longform
      gt_df <- pivot_longer(
        gt_df,
        cols = -c(CHROM, POS),
        names_to = "SampleID",
        values_to = "GT"
      )
      
      #define replicate pairs
      sample_names <- unique(gt_df$SampleID)
      
      #for GWAdapts
      if (GWAdapt == TRUE) {
        samples_df <- tibble(SampleName = sample_names) %>%
          mutate(
            temp = str_extract(SampleName, "GWAdapt\\d+_.*") %>%
              str_replace("GWAdapt\\d+_", ""),
            ReplicateID = str_extract(temp, "[^_]+$")  # grab the last underscore-separated bit
          )
      }
      #for neutral panel
      else {
        samples_df <- tibble(SampleName = sample_names) %>%
          mutate(
            temp = str_extract(SampleName, "N200_T\\d+_.*") %>%
              str_replace("N200_T\\d+_", ""),
            ReplicateID = str_extract(temp, "[^_]+$")  # grab the last underscore-separated bit
          )
      }
      replicate_map <- samples_df %>%
        group_by(ReplicateID) %>%
        filter(n() > 1) %>%  # only keep true replicates
        summarise(Pairs = list(as_tibble(
          t(combn(SampleName, 2)), .name_repair = ~ c("Sample1", "Sample2")
        )), .groups = "drop") %>%
        unnest(Pairs)
      
      #add short sampleID column to gt_df for matc replicates
      gt_df <- gt_df %>%
        mutate(
          ShortID = str_extract(SampleID, "[^_]+_[^_]+$"),
          ReplicateID = str_extract(ShortID, "[^_]+$")  # grab the last underscore-separated bit
        )
      # gt_df has: CHROM, POS, SampleID (full name), GT, ShortID, ReplicateID
      # replicate_map has: ReplicateID (short), Sample1 (full name), Sample2 (full name)
      
      # First, join original sample genotype
      paired_gt_df <- replicate_map %>%
        left_join(gt_df,
                  by = c("Sample1" = "SampleID"),
                  relationship = "many-to-many") %>%
        rename(GT1 = GT,
               CHROM = CHROM,
               POS = POS) %>%
        left_join(gt_df, by = c("Sample2" = "SampleID", "CHROM", "POS")) %>%
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
      
      #summarize_by_locus
      summary_by_locus <- gt_pairs_mismatch %>%
        group_by(CHROM, POS) %>%
        summarize(
          N = n(),
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
      
      summary_by_pair_clean <- summary_by_pair %>%
        filter(!is.na(Mismatch_Rate))
      
      summary_by_locus_clean <- summary_by_locus %>%
        filter(!is.na(Both_Mismatch_Rate))
      #store all the data
      temp_matrix <- matrix(
        data = NA,
        nrow = length(summary_by_pair_clean$Sample1),
        ncol = 7
      )
      temp_matrix[, 1] <- rep(j, n = length(summary_by_pair_clean$Sample1))
      temp_matrix[, 2] <- rep(h, n = length(summary_by_pair_clean$Sample1))
      temp_matrix[, 3] <- rep(i, n = length(summary_by_pair_clean$Sample1))
      temp_matrix[, 4] <- summary_by_pair_clean$Sample1
      temp_matrix[, 5] <- summary_by_pair_clean$Sample2
      temp_matrix[, 6] <- summary_by_pair_clean$Inclusive_Missingness
      temp_matrix[, 7] <- summary_by_pair_clean$Mismatch_Rate
      colnames(temp_matrix) <- c("Depth",
                                 "SNPmissing",
                                 "MissingnessFilter",
                                 "Sample1",
                                 "Sample2",
                                 "Missingness",
                                 "Mismatch")
      
      #for loci too
      temp_matrix2 <- matrix(
        data = NA,
        nrow = length(summary_by_locus_clean$CHROM),
        ncol = 7
      )
      temp_matrix2[, 1] <- rep(j, n = length(summary_by_locus_clean$CHROM))
      temp_matrix2[, 2] <- rep(h, n = length(summary_by_pair_clean$Sample1))
      temp_matrix2[, 3] <- rep(i, n = length(summary_by_locus_clean$CHROM))
      temp_matrix2[, 4] <- summary_by_locus_clean$CHROM
      temp_matrix2[, 5] <- summary_by_locus_clean$POS
      temp_matrix2[, 6] <- summary_by_locus_clean$Inclusive_Missingness
      temp_matrix2[, 7] <- summary_by_locus_clean$Both_Mismatch_Rate
      colnames(temp_matrix2) <- c("Depth",
                                  "SNPmissing",
                                  "MissingnessFilter",
                                  "CHROM",
                                  "POS",
                                  "Missingness",
                                  "Mismatch")
      
      #append to big matrix
      miss_mismatch <- rbind(miss_mismatch, temp_matrix)
      locus_mismatch <- rbind(locus_mismatch, temp_matrix2)
    }
  }
}

####to see how many samples and loci would be retained for each filter scheme###
retained<-as.data.frame(matrix(ncol = 5, nrow =0))
colnames(retained)<-c("Depth","LociFilter" ,"MissingnessFilter", "Samples", "Loci")

for (j in depths) {
  vcfR_dp <- hard_filter(vcfR = vcfR, depth = j)
  for (h in snps) {
    vcfSNPs <- missing_by_snp(vcfR = vcfR_dp, cutoff = h)
    for (i in filters) {
      vcftemp <- missing_by_sample(vcfR = vcfSNPs, cutoff = i)
      temp_matrix <- matrix(data = NA,
                            nrow = 1,
                            ncol = 5)
      colnames(temp_matrix) <- c("Depth",
                                 "LociFilter" ,
                                 "MissingnessFilter",
                                 "Samples",
                                 "Loci")
      temp_matrix[, 1] <- j
      temp_matrix[, 2] <- h
      temp_matrix[, 3] <- i
      temp_matrix[, 4] <- ncol(vcftemp@gt) - 1
      temp_matrix[, 5] <- nrow(vcftemp@fix)
      retained <- rbind(retained, temp_matrix)
    }
  }
}

###make the data table in the right format
retained$Samples<-as.numeric(retained$Samples)
retained$Loci<-as.numeric(retained$Loci)
retained$Depth<-as.numeric(retained$Depth)
retained$LociFilter<-as.numeric(retained$LociFilter)
retained$MissingnessFilter<-as.numeric(retained$MissingnessFilter)

####PLOTTING WITH THREE FILTERING STEPS####
#summarize data
miss_mismatch$Mismatch <- as.numeric(miss_mismatch$Mismatch)
miss_mismatch$Missingness <- as.numeric(miss_mismatch$Missingness)

summary_tbl <- miss_mismatch %>%
  group_by(Depth, SNPmissing, MissingnessFilter) %>%
  summarize(
    mean_mismatch = mean(Mismatch, na.rm = TRUE),
    mean_missingness = mean(Missingness, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
summary_tbl$Depth <- as.numeric(as.character(summary_tbl$Depth))
summary_tbl$SNPmissing <- as.numeric(as.character(summary_tbl$SNPmissing))
summary_tbl$MissingnessFilter <- as.numeric(as.character(summary_tbl$MissingnessFilter))

summary_tbl_binned <- summary_tbl %>%
  mutate(
    mismatch_bin = case_when(
      mean_mismatch < 0.005 ~ "< 0.005",
      mean_mismatch < 0.01 ~ "0.005 - 0.01",
      mean_mismatch < 0.015 ~ "0.01 - 0.015",
      mean_mismatch < 0.02 ~ "0.015 - 0.02",
      mean_mismatch >= 0.02 ~ "> 0.02",
      TRUE ~ NA_character_
    ),
    mismatch_bin = factor(
      mismatch_bin,
      levels = c("< 0.005", "0.005 - 0.01", "0.01 - 0.015", "0.015 - 0.02", "> 0.02")
    )
  ) %>%
  mutate(Depth = factor(
    Depth,
    levels = c(5, 10, 20, 30),
    labels = c("Read Depth 5", "Read Depth 10", "Read Depth 20", "Read Depth 30")
  ))

mismatch_colors <- c(
  "< 0.005" = "green2",
  "0.005 - 0.01" = "greenyellow",
  "0.01 - 0.015" = "yellow",
  "0.015 - 0.02" = "orange",
  "> 0.02" = "orangered"
)
#plot mismatch
facet_plot<-ggplot(summary_tbl_binned,
                   aes(
                     x = as.factor(SNPmissing),
                     y = as.factor(MissingnessFilter),
                     fill = mismatch_bin
                   )) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), color = "black", size = 4) +
  facet_wrap(~ Depth) +
  scale_fill_manual(values = mismatch_colors, name = "Mismatch Rate") +
  labs(
    x = "Locus Completeness Threshold",
    y = "Missingness Filter by Sample",
    title = "Mismatch Rate Across Filtering Parameters",
    subtitle = "Tile labels show replicate pair count (n)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(paste0("./outputs/CanFam3.1/missmatch_missingness_DP_facet_", panel, "_", genome, ".pdf"), 
       plot = facet_plot, width = 9, height = 5, units = "in")
#save as png or pdf

#retained loci
loci_colors <- c(
  "300+" = "green2",
  "250-300" = "greenyellow",
  "200-250" = "yellow",
  "150-200" = "orange",
  "<150" = "orangered"
)
retained_binned <- retained %>%
  mutate(loci_bin = case_when(
    Loci < 150 ~ "<150",
    Loci < 200 ~ "150-200",
    Loci < 250 ~ "200-250",
    Loci < 300 ~ "250-300",
    Loci >= 300 ~ "300+",
    TRUE ~ NA_character_
  )) %>%
  mutate(Depth = factor(
    Depth,
    levels = c(5, 10, 20, 30),
    labels = c("Read Depth 5", "Read Depth 10", "Read Depth 20", "Read Depth 30")
  ))

retained_plot<-ggplot(retained_binned,
                      aes(
                        x = as.factor(LociFilter),
                        y = as.factor(MissingnessFilter),
                        fill = loci_bin
                      )) +
  geom_tile(color = "white") +
  geom_text(aes(label = Samples), color = "black", size = 4) +
  facet_wrap( ~ Depth) +
  scale_fill_manual(values = loci_colors, name = "Loci") +
  labs(
    x = "Locus Completeness Threshold",
    y = "Missingness Filter by Sample",
    title = "Samples and Loci Retained Across Filtering Parameters",
    subtitle = "Tile labels show sample count (n)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(paste0("./outputs/CanFam3.1/retained_missingness_DP_facet_", panel,"_",genome, ".png"), 
       plot = retained_plot, width = 9, height = 5, units = "in")
