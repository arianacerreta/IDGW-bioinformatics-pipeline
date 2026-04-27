##### Final code for calculating error rates for the GWAdapt panel #####

#libraries (base R version 4.2.3)
library(vcfR) #version 1.15.0
library(tidyverse) #version 2.0.0: dplyr 1.2.1, forcats 1.0.1, lubridate 1.9.5, purrr 1.2.2, readr 2.2.0, stringr 1.6.0, tibble 3.3.1, tidyr 1.3.2
library(ggplot2)#version 4.0.2
library(SNPfiltR) #version 1.0.1

#Load vcf file generated in final_filter_with_PCRneg.R
vcf<- read.vcfR("./outputs/filtered-by-neg_GWAdapt_CanFam3.1_GATK.vcf") #vcf without a DP filter by PCR negatives; no longer includes PCR negatives

#panel ID (important if you are looking at several panels)
panel<-"GWAdapt"
genome<-"CanFam3.1"
pipeline<-"GATK"

#extract data to make a pop map
popmap_long<-data.frame(id=colnames(vcf@gt)[2:length(colnames(vcf@gt))]) 
vcfR <- vcf

popmap<- popmap_long %>%
  mutate(
    sample_id=str_split_i(id, "_", 5)
  )

length(unique(popmap$sample_id)) #number of uniquely ID lab samples

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
{filters<- seq(from = 0.2, to = 0.9, by = 0.1) #sample completeness
  depths<-c(20,10,5)
  snps<-seq(from = 0.80, to = 0.50, by = -0.05) #locus completeness
  quality<-c(30, 20, 10)
  miss_mismatch<-as.data.frame(matrix(ncol = 8, nrow =0))
  colnames(miss_mismatch)<-c("Depth","SNPmissing","MissingnessFilter","Quality","Sample1","Sample2","Missingness","Mismatch")
  locus_mismatch<-as.data.frame(matrix(ncol = 8, nrow =0))
  colnames(locus_mismatch)<-c("Depth","SNPmissing","MissingnessFilter","Quality","CHROM","POS","Missingness","Mismatch")
}

#altered function
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

#reminder: naming convention is i7/plateID_well_platename_UIID_sample
for (j in depths) {
  vcfR_dp <- hard_filter(vcfR = vcfR, depth = j)
  
  for (k in quality) {
    vcfR_gq<- hard_filter_rgq(vcfR= vcfR_dp, gq = k)
    
    for (h in snps) {
      vcfSNPs<- missing_by_snp(vcfR = vcfR_gq, cutoff = h)
      
      for (i in filters) {
        vcftemp <- missing_by_sample(vcfR = vcfSNPs, cutoff = i)
        gt <- extract.gt(vcftemp, element = "GT", as.numeric = FALSE) #this no longer has PCR negatives
        gt_df <- as.data.frame(gt)
        
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
        samples_df<- tibble(SampleName = sample_names) %>%
          mutate(
            ReplicateID = str_split_i(SampleName, "_", 5)  # grab the last underscore-separated bit
          )
    
        replicate_map <- samples_df %>%
          group_by(ReplicateID) %>%
          filter(n() > 1) %>%  # only keep true replicates
          summarise(Pairs = list(as_tibble(
            t(combn(SampleName, 2)), .name_repair = ~ c("Sample1", "Sample2")
          )), .groups = "drop") %>%
          unnest(Pairs)
        
        #add short sampleID column to gt_df for match replicates
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
          ncol = 8
        )
        temp_matrix[, 1] <- rep(j, n = length(summary_by_pair_clean$Sample1))
        temp_matrix[, 2] <- rep(h, n = length(summary_by_pair_clean$Sample1))
        temp_matrix[, 3] <- rep(i, n = length(summary_by_pair_clean$Sample1))
        temp_matrix[, 4] <- rep(k, n = length(summary_by_pair_clean$Sample1))
        temp_matrix[, 5] <- summary_by_pair_clean$Sample1
        temp_matrix[, 6] <- summary_by_pair_clean$Sample2
        temp_matrix[, 7] <- summary_by_pair_clean$Inclusive_Missingness
        temp_matrix[, 8] <- summary_by_pair_clean$Mismatch_Rate
        colnames(temp_matrix) <- c("Depth",
                                   "SNPmissing",
                                   "MissingnessFilter",
                                   "Quality",
                                   "Sample1",
                                   "Sample2",
                                   "Missingness",
                                   "Mismatch")
        
        #for loci too
        temp_matrix2 <- matrix(
          data = NA,
          nrow = length(summary_by_locus_clean$CHROM),
          ncol = 8
        )
        temp_matrix2[, 1] <- rep(j, n = length(summary_by_locus_clean$CHROM))
        temp_matrix2[, 2] <- rep(h, n = length(summary_by_locus_clean$CHROM))
        temp_matrix2[, 3] <- rep(i, n = length(summary_by_locus_clean$CHROM))
        temp_matrix2[, 4] <- rep(k, n = length(summary_by_locus_clean$CHROM))
        temp_matrix2[, 5] <- summary_by_locus_clean$CHROM
        temp_matrix2[, 6] <- summary_by_locus_clean$POS
        temp_matrix2[, 7] <- summary_by_locus_clean$Inclusive_Missingness
        temp_matrix2[, 8] <- summary_by_locus_clean$Both_Mismatch_Rate
        colnames(temp_matrix2) <- c("Depth",
                                    "SNPmissing",
                                    "MissingnessFilter",
                                    "Quality",
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
}

####to see how many samples and loci would be retained for each filter scheme###
retained<-as.data.frame(matrix(ncol = 6, nrow =0))
colnames(retained)<-c("Depth","LociFilter" ,"MissingnessFilter", "Quality" ,"Samples", "Loci")
for (j in depths) {
  vcfR_dp <- hard_filter(vcfR = vcfR, depth = j)
  for (k in quality) {
    vcfR_gq<- hard_filter_rgq(vcfR= vcfR_dp, gq = k)
    
    for (h in snps) {
      vcfSNPs <- missing_by_snp(vcfR = vcfR_gq, cutoff = h)
      for (i in filters) {
        vcftemp <- missing_by_sample(vcfR = vcfSNPs, cutoff = i)
        temp_matrix <- matrix(data = NA,
                              nrow = 1,
                              ncol = 6)
        colnames(temp_matrix) <- c("Depth",
                                   "LociFilter" ,
                                   "MissingnessFilter",
                                   "Quality",
                                   "Samples",
                                   "Loci")
        temp_matrix[, 1] <- j
        temp_matrix[, 2] <- h
        temp_matrix[, 3] <- i
        temp_matrix[, 4] <- k
        temp_matrix[, 5] <- ncol(vcftemp@gt) - 1
        temp_matrix[, 6] <- nrow(vcftemp@fix)
        retained <- rbind(retained, temp_matrix)
      }
    }
  }
}

###make the data table in the right format
retained$Samples<-as.numeric(retained$Samples)
retained$Loci<-as.numeric(retained$Loci)
retained$Depth<-as.numeric(retained$Depth)
retained$Quality<-as.numeric(retained$Quality)
retained$LociFilter<-as.numeric(retained$LociFilter)
retained$MissingnessFilter<-as.numeric(retained$MissingnessFilter)

####PLOTTING WITH THREE FILTERING STEPS####
#summarize data
miss_mismatch$Mismatch <- as.numeric(miss_mismatch$Mismatch)
miss_mismatch$Missingness <- as.numeric(miss_mismatch$Missingness)

summary_tbl <- miss_mismatch %>%
  group_by(Depth, SNPmissing, MissingnessFilter, Quality) %>%
  summarize(
    mean_mismatch = mean(Mismatch, na.rm = TRUE),
    mean_missingness = mean(Missingness, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

summary_tbl$Quality <- as.numeric(as.character(summary_tbl$Quality))
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
  ) )%>%
  mutate(Quality = factor(
    Quality,
    levels = c(10,20,30),
    labels = c("GQ 10", "GQ 20", "GQ 30")
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
  facet_grid(Depth ~ Quality) +
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

ggsave(paste0("./outputs/CanFam3.1/missmatch_missingness_DP_facet_", panel, "_", genome,"_",pipeline, ".png"), 
       plot = facet_plot, width = 9, height = 5, units = "in")
