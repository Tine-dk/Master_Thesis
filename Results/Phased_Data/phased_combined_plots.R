### Phased Combined Plots #############
## Load Libraries ###############################################################################################################
library(ggplot2)
library(xlsx)         # Export in Excel
library(tidyverse)    # data manipulation
library(cluster)      # clustering algorithms
library(stringr)




#### Tandem_Genotypes #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/phased_combined/tandem_genotypes/")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_combined.txt")

## Data processing ##############################################################################################################


for (file in input_files){
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # Define gene from input file
  gene = strsplit(file, split = "_")[[1]][2]
  
  # Process dataframe for plottin
  data_all = data.frame(count=c(data$Illumina, data$LongShot),haplotype=c(data$Haplotype), method=rep(c("illumina","longshot"),each=nrow(data)))
  data_all = na.omit(data_all)
  
  # Create plot for sample and gene
  ggplot(data_all, aes(x=method, y=count,color=haplotype,fill=haplotype)) +
    geom_jitter(size=2) +
    ylab("Repeat Count") +
    xlab("") +
    scale_x_discrete(breaks = waiver(), labels=c(paste0("Illumina\n(hap1 n = ", nrow(filter(data_all, method == "illumina" & haplotype == "hap1")),")",
                                                        " \n(hap2 n = ", nrow(filter(data_all, method == "illumina" & haplotype == "hap2")),")"),
                                                 paste0("LongShot\n(hap1 n = ", nrow(filter(data_all, method == "longshot" & haplotype == "hap1")),")",
                                                        " \n(hap2 n = ", nrow(filter(data_all, method == "longshot" & haplotype == "hap2")),")"))) +
    # xticks(paste0("RepeatHMM\nn = ", length(na.omit(data$RepeatHMM)))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=14)) +
    geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "grey") +
    coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(sample," ",gene), subtitle = paste0("Expected values: 18/18")) +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey")
  
  ggsave(path = "plots", filename = paste0(sample,"_",gene,"_tandemgenotypes_phased_hline.pdf"))
}



#### RepeatHMM #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/phased_combined/RepeatHMM/")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_combined.txt")

## Data processing ##############################################################################################################


for (file in input_files){
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # Define gene from input file
  gene = strsplit(file, split = "_")[[1]][2]
  
  # Process dataframe for plottin
  data_all = data.frame(count=c(data$Illumina, data$LongShot),haplotype=c(data$Haplotype), method=rep(c("illumina","longshot"),each=nrow(data)))
  data_all = na.omit(data_all)
  
  # Create plot for sample and gene
  ggplot(data_all, aes(x=method, y=count,color=haplotype,fill=haplotype)) +
    geom_jitter(size=2) +
    ylab("Repeat Count") +
    xlab("") +
    scale_x_discrete(breaks = waiver(), labels=c(paste0("Illumina\n(hap1 n = ", nrow(filter(data_all, method == "illumina" & haplotype == "hap1")),")",
                                                        " \n(hap2 n =", nrow(filter(data_all, method == "illumina" & haplotype == "hap2")),")"),
                                                 paste0("LongShot\n(hap1 n = ", nrow(filter(data_all, method == "longshot" & haplotype == "hap1")),")",
                                                        " \n(hap2 n =", nrow(filter(data_all, method == "longshot" & haplotype == "hap2")),")"))) +
    # xticks(paste0("RepeatHMM\nn = ", length(na.omit(data$RepeatHMM)))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=14)) +
    geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "grey") +
    coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(sample," ",gene), subtitle = paste0("Expected values: 18/18")) +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey") +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey")
  
  ggsave(path = "plots", filename = paste0(sample,"_",gene,"_repeathmm_phased_hline.pdf"))
}






#### STRique #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/phased_combined/STRique/")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_results.txt")

## Data processing ##############################################################################################################

for (file in input_files){
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  for (gene in c("HTT","fmr1")){
    # subset input into gene specific table and perform Threshold Count Filtering
    data_subset = data %>% filter(target == gene, score_prefix > 4 & score_suffix > 4)
    
    # Get illumina-guided haplotype results
    if(gene == "fmr1"){
      gene = "FMR1"
    }
    
    illumina_guided_files <- list.files("../../SNV_analyse/Output_haplotype_phasing/", pattern = paste0(sample,"_",gene,"_reads"),full.names = T)
    if (!identical(illumina_guided_files, character(0))){
      illumina_data = read.table(illumina_guided_files, header = T)
      names(illumina_data) = c("ID","illumina_haplotype", "illumina_origin")
      
      data_subset = merge(data_subset, illumina_data, by = "ID", all.x = T)
    }

    # Get longshot haplotype results
    chr = ifelse(gene == "HTT", "chr4", "chrX")
    longshot_files <- list.files("../../longshot/readnames/combined/", pattern = paste0(sample,"_combined_longshot_hp_",chr,"_readnames"),full.names = T)
    if (!identical(longshot_files, character(0))){
      longshot_data = read.table(longshot_files, header = F, fill = NA)
      if(ncol(longshot_data) > 1){
        colnames(longshot_data) <- c("ID", "longshot_haplotype")
        longshot_data = longshot_data %>% filter(longshot_haplotype != "")
        longshot_data$longshot_haplotype[longshot_data$longshot_haplotype=="HP:i:1"] <- "hap1"
        longshot_data$longshot_haplotype[longshot_data$longshot_haplotype=="HP:i:2"] <- "hap2"
        
        data_subset = merge(data_subset, longshot_data, by = "ID", all.x = T)
      }
    }
    
    # Process dataframe for plottin
    data_all = data.frame(count=rep(c(data_subset$count), 2), haplotype=c(data_subset$illumina_haplotype, data_subset$longshot_haplotype), method=rep(c("illumina","longshot"),each=nrow(data_subset)))
    data_all = na.omit(data_all)
    
    # Plot the data
    ggplot(data_all, aes(x=method, y=count,color=haplotype,fill=haplotype)) +
      geom_jitter(size=2) +
      ylab("Repeat Count") +
      xlab("") +
      scale_x_discrete(breaks = waiver(), labels=c(paste0("Illumina\n(hap1 n = ", nrow(filter(data_all, method == "illumina" & haplotype == "hap1")),")",
                                                          " \n(hap2 n =", nrow(filter(data_all, method == "illumina" & haplotype == "hap2")),")"),
                                                   paste0("LongShot\n(hap1 n = ", nrow(filter(data_all, method == "longshot" & haplotype == "hap1")),")",
                                                          " \n(hap2 n = ", nrow(filter(data_all, method == "longshot" & haplotype == "hap2")),")"))) +
      theme(axis.text=element_text(size=14), 
            axis.title=element_text(size=14)) +
      geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "grey") +
      coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(sample," ",gene), subtitle = paste0("Expected values: 20/183-199")) +
      geom_hline(yintercept=20, linetype="dashed", size = 0.25, color = "darkgrey") +
      geom_hline(yintercept=191, linetype="dashed", size = 0.25, color = "darkgrey")
    
      
    ggsave(path = "plots", filename = paste0(sample,"_",gene,"_STRique_phased_hline.pdf"))
    
    }
}



#### TRiCoLOR #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/phased_combined/tricolor/")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_combined.txt")

## Data processing ##############################################################################################################


for (file in input_files){
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # Define gene from input file
  gene = strsplit(file, split = "_")[[1]][2]
  
  # Process dataframe for plottin
  data_all = data.frame(count=c(data$Illumina, data$LongShot),haplotype=c(data$Haplotype), method=rep(c("illumina","longshot"),each=nrow(data)))
  data_all = na.omit(data_all)
  
  # Create plot for sample and gene
  ggplot(data_all, aes(x=method, y=count,color=haplotype,fill=haplotype)) +
    geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.75, size = 4) +
    ylab("Repeat Count") +
    xlab("") +
    scale_x_discrete(breaks = waiver(), labels=c(paste0("Illumina\n(hap1 n = ", nrow(filter(data_all, method == "illumina" & haplotype == "hap1")),")",
                                                        " \n(hap2 n =", nrow(filter(data_all, method == "illumina" & haplotype == "hap2")),")"),
                                                 paste0("LongShot\n(hap1 n = ", nrow(filter(data_all, method == "longshot" & haplotype == "hap1")),")",
                                                        " \n(hap2 n =", nrow(filter(data_all, method == "longshot" & haplotype == "hap2")),")"))) +
    # xticks(paste0("RepeatHMM\nn = ", length(na.omit(data$RepeatHMM)))) +
    theme(axis.text=element_text(size=14), 
          axis.title=element_text(size=14)) +
    geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "grey") +
    coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(sample," ",gene), subtitle = paste0("Expected values: 18/18")) +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey") +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey")
  
  ggsave(path = "plots", filename = paste0(sample,"_",gene,"_tricolor_phased_hline.pdf"))
}





#### Program Means Illumina #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/phased_combined/comparison/")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_combined_illumina.txt")

## Data processing ##############################################################################################################


for (file in input_files){
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # Define gene from input file
  gene = strsplit(file, split = "_")[[1]][2]
  
  # Process dataframe for plottin
  data_all = data.frame(count=data$Count_mean, haplotype=c(data$Haplotype), program=rep(c("RepeatHMM","STRique", "Tandem-Genotypes", "TRiCoLOR"),each=nrow(data)))
  data_all = na.omit(data_all)
  
  # Create plot for sample and gene
  ggplot(data, aes(x=Program, y=Count_mean,color=Haplotype,fill=Haplotype)) +
    geom_point(position=position_jitter(h=0.1, w=0.1),
               shape = 21, alpha = 0.75, size = 4) +
    ylab("Repeat Count") +
    xlab("") +
    scale_x_discrete(breaks = waiver(), labels=c("RepeatHMM","STRique", "Tandem-Genotypes", "TRiCoLOR")) +
    theme(axis.text=element_text(size=8), 
          axis.title=element_text(size=12)) +
    geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "grey") +
    coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(sample," ",gene," Illumina HP"), subtitle = paste0("Expected values: 18/18")) +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey") +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey")
  
  ggsave(path = "plots", filename = paste0(sample,"_",gene,"_comparison_Illumina_phased_hline.pdf"))
  
}



#### Program Means LongSHot #########################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/phased_combined/comparison/")
## Read tables ##############################################################################################################
input_files = list.files(pattern = "_combined_longshot.txt")

## Data processing ##############################################################################################################


for (file in input_files){
  
  # read input
  data = read.table(file, header = T, na.strings = "NA")
  
  # Define sample from input file
  sample = strsplit(file, split = "_")[[1]][1]
  
  # Define gene from input file
  gene = strsplit(file, split = "_")[[1]][2]
  
  # Process dataframe for plottin
  data_all = data.frame(count=data$Count_mean, haplotype=c(data$Haplotype), program=rep(c("RepeatHMM","STRique", "Tandem-Genotypes", "TRiCoLOR"),each=nrow(data)))
  data_all = na.omit(data_all)
  
  # Create plot for sample and gene
  ggplot(data, aes(x=Program, y=Count_mean,color=Haplotype,fill=Haplotype)) +
    geom_point(position=position_jitter(h=0.1, w=0.1),
               shape = 21, alpha = 0.75, size = 4) +
    ylab("Repeat Count") +
    xlab("") +
    scale_x_discrete(breaks = waiver(), labels=c("RepeatHMM","STRique", "Tandem-Genotypes", "TRiCoLOR")) +
    theme(axis.text=element_text(size=8), 
          axis.title=element_text(size=12)) +
    geom_violin(width=0.1, color="grey", alpha=0.2, lwd = 0.4, fill = "grey") +
    coord_cartesian(ylim = c(0,max(na.omit(data_all$count)))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(sample," ",gene," LongShot HP"), subtitle = paste0("Expected values: 18/18")) +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey") +
    geom_hline(yintercept=18, linetype="dashed", size = 0.25, color = "darkgrey")
  
  ggsave(path = "plots", filename = paste0(sample,"_",gene,"_comparison_longshot_phased_hline.pdf"))
  
}
