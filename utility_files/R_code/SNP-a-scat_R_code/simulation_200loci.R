#Updated by Ariana Cerreta 4-7-2026

###load in libraries###
library(tidyverse)
library(CKMRsim)
library(vcfR)
library(ggplot2)
library(gridExtra)

#####edited functions#####
##edit to CKMRsim::reindex_markers function for Issue#9 on github (https://github.com/eriqande/CKMRsim/issues/9)
reindex_markers<- function(M){
  M %>% dplyr::ungroup() %>% dplyr::arrange(Chrom, Pos, desc(Freq)) %>% 
    #dplyr::group_by(Chrom) %>% 
    dplyr::mutate(locidx = as.integer(factor(Locus, levels = unique(Locus)))) %>% 
    dplyr::group_by(Chrom, Locus) %>% 
    dplyr::mutate(alleidx = as.integer(factor(Allele, levels = unique(Allele))), newfreq = Freq/sum(Freq)) %>% 
    dplyr::select(-AlleIdx, -LocIdx, -Freq) %>% 
    rename(Freq = newfreq, AlleIdx = alleidx, LocIdx = locidx) %>% 
    dplyr::ungroup()
}

#####read in data#####
hap.vcf<-read.vcfR("inputs/Test2_i004_whitelist-loci-only_DP20/populations.haps.vcf")

#####extract and format data####
gt<-extract.gt(hap.vcf, return.alleles = FALSE, convertNA = TRUE)

#convert to dataframe and mutate rownames to a column value called Locus
gt_df<-as.data.frame(gt)
gt_df<- gt_df %>%
  rownames_to_column(var = "Locus")

#extract sample names
sample_names<-colnames(gt) #have to refer to origianal gt file to avoid having 'Locus' included as a sample name

#extract CHROM and POS from VCF
chrom_pos<- as.data.frame(hap.vcf@fix[,c("CHROM", "POS", "ID")])
chrom_pos$POS<- as.numeric(chrom_pos$POS) #convert POS to numeric

# extra numerical part from CHROM and convert to interger
chrom_pos$Chrom <- as.integer(str_extract(chrom_pos$CHROM,"\\d+"))

#add locus as a unique identifier for merging
chrom_pos$Locus <- chrom_pos$ID

#calculate GT%
loci<- as.numeric(length(gt_df$Locus))
num_samps<-ncol(hap.vcf@gt)-1
freq_GT_by_sample<-matrix(data=NA, nrow = num_samps, ncol = 2)
colnames(freq_GT_by_sample)<- c("sample","GTfreq")

for (i in 2:ncol(hap.vcf@gt)){
  num_na<-as.numeric(sum(is.na(gt_df[,i]))) #how many NAs are in an individuals col
  GTfreq<- (loci-num_na)/loci #calculate the frequency of present loci (aka %GT)
  name<-colnames(gt_df)[i] #which sample are we calculating for
  freq_GT_by_sample[i-1,1]<-name
  freq_GT_by_sample[i-1,2]<-GTfreq
}

all_colnames<-colnames(gt_df)
censor_samples<-freq_GT_by_sample[freq_GT_by_sample[,"GTfreq"]< 0.5, "sample"]
samples_keep<- setdiff(all_colnames, censor_samples)
gt_df<-gt_df[,samples_keep]

#calculate loci coverage
loci_cov<-matrix(data = NA, nrow =nrow(gt_df), ncol= 2)
colnames(loci_cov)<-c("Locus","Cov")

for (i in 1:nrow(gt_df)){
  num_na<-as.numeric(sum(is.na(gt_df[i,])))
  samps<-ncol(gt_df)-1
  cov<-(samps-num_na)/samps
  loci_name<-gt_df[i,1]
  loci_cov[i,1]<-loci_name
  loci_cov[i,2]<-cov
}

loci_names<-gt_df[,1]
censor_loci<-loci_cov[loci_cov[,"Cov"]<0.5, "Locus"]
loci_keep<- setdiff(loci_names,censor_loci)
gt_df<- gt_df %>%
  filter(Locus %in% loci_keep)

#long form genotype tibble
long_geno<- gt_df %>%
  left_join(chrom_pos, by = "Locus") %>%
  pivot_longer(cols = -c(Locus,POS,ID, Chrom, CHROM), names_to = "Indiv", values_to = "Genotype") %>%
  mutate(gene_copy_1 = str_split(Genotype, "/", simplify = TRUE)[,1],
         gene_copy_2 = str_split(Genotype, "/", simplify = TRUE)[,2]) %>%
  select(-c(Genotype,ID, CHROM)) %>%
  pivot_longer(cols = starts_with("gene_copy"), names_to = "gene_copy", values_to = "Allele") %>%
  mutate(gene_copy = ifelse(gene_copy == "gene_copy_1", "1", "2")) %>%
  mutate(Allele = ifelse(Allele == "" | Allele == "NA", NA, Allele))

#save as tibble
long_geno<- as_tibble(long_geno)

#calculate allele frequencies
allele_freqs <- long_geno %>%
  count(Locus, Allele) %>%
  group_by(Locus) %>%
  mutate(Freq = n/sum(n),
         Chrom = as.character(long_geno$Chrom[match(Locus,long_geno$Locus)]),
         Pos = long_geno$POS[match(Locus,long_geno$Locus)]) %>%
  ungroup() %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele))

#implement updated function
afreqs_ready <- reindex_markers(allele_freqs)

#2024.23.04 - adding on to DJM (Dana Morin)'s script to import .csv and format allele frequencies and run simulation-ALC

al.freq.ls<-NULL

wide_afreqs<-matrix(ncol = 1+max(afreqs_ready$AlleIdx), nrow= length(unique(afreqs_ready$LocIdx)))
colnames(wide_afreqs)<-c("Marker",1:max(afreqs_ready$AlleIdx))
wide_afreqs[,1]<-unique(afreqs_ready$Locus)

afreq_wide<- afreqs_ready %>% 
  select(Locus, AlleIdx, Freq)%>%
  pivot_wider(
    names_from = AlleIdx,
    values_from =Freq,
    names_prefix = "Marker"
  ) %>%
  as.data.frame()

#formatting your afreqs_wide to work with Dana's code

for (i in 1:length(afreq_wide$Locus)){
  addon<-unname(unlist(c(afreq_wide[i,-1])))#accommodates any number of frequencies per locus
  addon<-addon[!is.na(addon)]
  addon<-list(addon)
  al.freq.ls<-(c(al.freq.ls,addon))
  
}
#####SIMULATION code to randomly choose X loci, calculate PID, PIDsibs, repeat#####
#adapted from Dana Morin's original code

N<-1000 #number of simulations
x<-200 #max number of loci to simulate through

PID.pop.tbl<-matrix(data = NA, nrow = x, ncol = 3) #matrix for output of sim for PID
PID.pop.sib.tbl<-matrix(data = NA, nrow = x, ncol = 3) #matrix for output of sim for PIDsib

for (j in 1:x){ #for each count of loci
  
  PIpop.loci<-NULL
  PIpop.sib.loci<-NULL
  
  for(h in 1:N){ #do N simulations
    samp<-sample(al.freq.ls, j, replace = FALSE) #sample j alleles from # available
    PIpop <- NULL #clear temp obj
    PIpop.sib <- NULL #clear temp obj
    for(i in 1:length(samp)){#waits et al 2001
      
      r1 <- sum(samp[[i]]^4)
      r2 <- outer(X=samp[[i]],Y=samp[[i]],FUN=function(X,Y)(2*X*Y)^2)
      r <- sum(r2[lower.tri(r2)])+r1 #has correction to match eqn 1 in waits et al 2001
      s <- 0.25 + (0.5 * sum(samp[[i]]^2)) + (0.5*sum(samp[[i]]^2)^2) -
        (0.25*sum(samp[[i]]^4))
      
      PIpop <- c(PIpop,r)
      PIpop.sib <- c(PIpop.sib,s)
    }
    PIpop.loci<-c(PIpop.loci,prod(PIpop))
    PIpop.sib.loci<-c(PIpop.sib.loci, prod(PIpop.sib))
  }
  PID.pop.tbl[j,1]<-j
  PID.pop.tbl[j,2]<-mean(PIpop.loci)
  PID.pop.tbl[j,3]<-sd(PIpop.loci)
  PID.pop.sib.tbl[j,1]<-j
  PID.pop.sib.tbl[j,2]<-mean(PIpop.sib.loci)
  PID.pop.sib.tbl[j,3]<-sd(PIpop.sib.loci)
}

PID.pop.sib.tbl<-as.data.frame(PID.pop.sib.tbl)
colnames(PID.pop.sib.tbl)<-c("Loci", "Avg.PIDsib", "SD.PIDsib")
PID.pop.tbl<-as.data.frame((PID.pop.tbl))
colnames(PID.pop.tbl)<-c("Loci", "Avg.PID", "SD.PIDsib")

#####plotting####
#ggplot graphing 1-2-2025 ALC
tempPID<-PID.pop.tbl
tempPID$Type<-rep("PID", length(tempPID$Loci))
colnames(tempPID)<-c("Loci", "Avg.Probability", "SD","Type")

tempPIDsib<-PID.pop.sib.tbl
tempPIDsib$Type<-rep("PIDsibs", length(tempPIDsib$Loci))
colnames(tempPIDsib)<-c("Loci", "Avg.Probability", "SD" ,"Type")

PIDcombo<-rbind(tempPID,tempPIDsib)

plot1<-ggplot(PIDcombo, aes(x=Loci, y=Avg.Probability ,color =Type, linetype = Type))+
  geom_line(size=1.5)+
  coord_cartesian(xlim = c(0,200), ylim = c(0,0.01))+
  scale_color_manual(values=c("gray68","black"))+
  scale_linetype_manual(values=c("solid","dashed"))+
  theme(text=element_text(size = 16, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(size = 1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.width = unit(2,'cm'),
        axis.text = element_text(size = 14),
        legend.position = "top",
        title = element_text(size=14))

zoom_plot<-ggplot(PIDcombo, aes(x=Loci, y=Avg.Probability ,color =Type, linetype = Type))+
  geom_line(size=1.5)+
  coord_cartesian(xlim = c(0,50), ylim = c(0,0.005))+
  scale_color_manual(values=c("gray68","black"))+
  scale_linetype_manual(values=c("solid","dashed"))+
  theme(text=element_text(size = 16, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(size = 14),
        legend.position = "none",
        title = element_text(size=14))

grid.arrange(plot1,zoom_plot, ncol=1)
