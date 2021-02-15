#########################################
#           Unphased Together           #
#########################################
## Load Libraries ###############################################################################################################
library(ggplot2)
library(xlsx)         # Export in Excel
library(tidyverse)    # data manipulation
library(cluster)      # clustering algorithms
library(stringr)

#### Whole Genome data #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/unphased_together/whole_genome")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_no_hp.txt")

## Data processing ##############################################################################################################


for (file in input_files){
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # Define gene from input file
  gene = strsplit(file, split = "_")[[1]][2]
  
  # Process dataframe for plottin
  data_all = data.frame(count=c(data$RepeatHMM, data$STRique, data$Tandem_Genotypes),
                        program=rep(c("RepeatHMM","STRique","Tandem-Genotypes"),each=nrow(data)),
                        color=rep(c("red","blue","orange"),each=nrow(data)))

  
  # Create plot for sample and gene
  ggplot(data_all, aes(x=program, y=count)) +
    geom_jitter(size=2, color = data_all$color) +
    ylab("Repeat Count") +
    xlab("")+
    scale_x_discrete(breaks = waiver(), labels=c(paste0("RepeatHMM\n(n = ", length(na.omit(data$RepeatHMM)),")"),
                                                 paste0("STRique\n(n = ", length(na.omit(data$STRique)),")"),
                                                 paste0("Tandem Genotypes\n(n = ", length(na.omit(data$Tandem_Genotypes)),")"))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=14)) +
    geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "grey") +
    coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(sample," ",gene), subtitle = paste0("Expected values: 18/18")) +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey") #+
    #geom_hline(yintercept=191, linetype="dashed", size = 0.25, color = "darkgrey")
  
  ggsave(path = "plots", filename = paste0(sample,"_",gene,"_unphased_hline.pdf"))
  
  
}







#### Targeted data #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/unphased_together/targeted")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_no_hp.txt")

## Data processing ##############################################################################################################


for (file in input_files){
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # Define gene from input file
  gene = strsplit(file, split = "_")[[1]][2]
  
  # Process dataframe for plottin
  data_all = data.frame(count=c(data$RepeatHMM, data$STRique, data$Tandem_Genotypes),
                        program=rep(c("RepeatHMM","STRique","Tandem-Genotypes"),each=nrow(data)),
                        color=rep(c("red","blue","orange"),each=nrow(data)))
  
  
  # Create plot for sample and gene
  ggplot(data_all, aes(x=program, y=count)) +
    geom_jitter(size=2, color = data_all$color) +
    ylab("Repeat Count") +
    xlab("")+
    scale_x_discrete(breaks = waiver(), labels=c(paste0("RepeatHMM\n(n = ", length(na.omit(data$RepeatHMM)),")"),
                                                 paste0("STRique\n(n = ", length(na.omit(data$STRique)),")"),
                                                 paste0("Tandem Genotypes\n(n = ", length(na.omit(data$Tandem_Genotypes)),")"))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=14)) +
    geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "pink") +
    coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(sample," ",gene), subtitle = paste0("Expected values: 20/183-199")) +
    geom_hline(yintercept=20, linetype="dashed", size = 0.25, color = "darkgrey") +
    geom_hline(yintercept=191, linetype="dashed", size = 0.25, color = "darkgrey")
  
  
  ggsave(path = "plots", filename = paste0(sample,"_",gene,"_unphased_hline.pdf"))
  
}



















#### Combined data #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/unphased_together/combined")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_no_hp.txt")

## Data processing ##############################################################################################################


for (file in input_files){
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # Define gene from input file
  gene = strsplit(file, split = "_")[[1]][2]
  
  # Process dataframe for plottin
  data_all = data.frame(count=c(data$RepeatHMM, data$STRique, data$Tandem_Genotypes),
                        program=rep(c("RepeatHMM","STRique","Tandem-Genotypes"),each=nrow(data)),
                        color=rep(c("red","blue","orange"),each=nrow(data)))
  
  # Create plot for sample and gene
  ggplot(data_all, aes(x=program, y=count)) +
    geom_jitter(size=2, color = data_all$color) +
    ylab("Repeat Count") +
    xlab("")+
    scale_x_discrete(breaks = waiver(), labels=c(paste0("RepeatHMM\n(n = ", length(na.omit(data$RepeatHMM)),")"),
                                                 paste0("STRique\n(n = ", length(na.omit(data$STRique)),")"),
                                                 paste0("Tandem Genotypes\n(n = ", length(na.omit(data$Tandem_Genotypes)),")"))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=14)) +
    geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "grey") +
    coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(sample," ",gene), subtitle = paste0("Expected values: 18/18")) +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey") +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey")
  
  ggsave(path = "plots", filename = paste0(sample,"_",gene,"_unphased_hline.pdf"))
  
}


