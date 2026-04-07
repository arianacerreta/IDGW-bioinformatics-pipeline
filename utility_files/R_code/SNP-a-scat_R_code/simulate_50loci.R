#DJM
#20240312 - modifying old script from Dec 2015

#ALC
#2024.23.04 - adding on to DJM's script to import .csv and format allele frequencies

#need al.freq.ls, which is a list of vectors for allele frequencies for each locus
#replace al.freq.ls with actual allele frequencies for each locus, but for example:

#al.freq.ls.example<-list(c(0.3,0.6, 0.1),c(0.3,0.7),c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), 
#                 c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), 
#                 c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), 
#                 c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), c(0.3,0.7), c(0.3,0.7),
#                 c(0.3,0.7), c(0.3,0.7), c(0.3,0.7))

#read in your own CSV file
#format: each row has a locus, first column is name of locus, columns 1-(max alleles/per locus) contain frequencies
#for microhaps, I treated each unique combination as its own allele and calculated frequencies using
#the populations.hapstats.tsv file that populations in stacks will output for you

library(ggplot2)
ST2_freq<-read.csv("inputs/70loci_freq_ScatTest2.csv", header = TRUE)
al.freq.ls<-NULL

#formatting your CSV to work with Dana's code
for (i in 1:length(ST2_freq$Marker)){
  addon<-unname(unlist(c(ST2_freq[i,-1])))#accommodates any number of frequencies per locus
    #addon<-c(ST2_freq[i,2],ST2_freq[i,3],ST2_freq[i,4],ST2_freq[i,5], ST2_freq[i,6]) #old code which had max 5 freq per allele
  addon<-addon[!is.na(addon)]
  addon<-list(addon)
  al.freq.ls<-(c(al.freq.ls,addon)) 
}

#ALC 2024.23.04
#### code to randomly choose X loci, calculate PID, PIDsib, repeat

N<-1000 #number of simulations
x<-50 #max number of loci to simulate through

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
colnames(PID.pop.sib.tbl)<-c("Loci", "Avg.PIDsib", "SD")
PID.pop.tbl<-as.data.frame((PID.pop.tbl))
colnames(PID.pop.tbl)<-c("Loci", "Avg.PID", "SD")

#ggplot graphing 1-2-2025 ALC
tempPID<-PID.pop.tbl
tempPID$Type<-rep("PID", length(tempPID$Loci))
colnames(tempPID)<-c("Loci", "Probability", "SD","Type")

tempPIDsib<-PID.pop.sib.tbl
tempPIDsib$Type<-rep("PIDsibs", length(tempPIDsib$Loci))
colnames(tempPIDsib)<-c("Loci", "Probability", "SD","Type")

PIDcombo<-rbind(tempPID,tempPIDsib)

ggplot(PIDcombo, aes(x=Loci, y=Probability ,color =Type, linetype = Type))+
  geom_line(size=1.5)+
  coord_cartesian(xlim = c(0,50), ylim = c(0,0.01))+
  scale_color_manual(values=c("gray68","black"))+
  scale_linetype_manual(values=c("solid","dashed"))+
  geom_hline(yintercept = 0.001, color = "red", size = 1)+
  theme(text=element_text(size = 16, family = "sans"),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(size = 1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.width = unit(2,'cm'),
        axis.text = element_text(size = 14),
        legend.position = c(0.85,0.9),
        title = element_text(size=14))                    
