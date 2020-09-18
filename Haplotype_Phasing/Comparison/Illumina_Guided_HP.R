##########################
### Illumina Guided HP ###
##########################

#################################################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/SNV_analyse")

## Load Libraries ###############################################################################################################
require(ggplot2)       # Logos Plot
require(ggseqlogo)     # Logos Plot
#library(kmer)
#require(xlsx)
#library(readxl)
#library(stringdist)    # Clustering of Strings i.e. the SNV k-mers


### Loading the table ###########################################################################################################
GM05538_HTT_window40000 <- read.table("G87-GM05538_HTT_bases_window40000.txt", header = TRUE, sep = ",")
GM06891_HTT_window40000 <- read.table("G87-GM06891_HTT_bases_window40000.txt", header = TRUE, sep = ",")
GM07541_HTT_window40000 <- read.table("G87-GM07541_HTT_bases_window40000.txt", header = TRUE, sep = ",")
GM07861_HTT_window40000 <- read.table("G87-GM07861_HTT_bases_window40000.txt", header = TRUE, sep = ",")
GM20239_HTT_window40000 <- read.table("G87-GM20239_HTT_bases_window40000.txt", header = TRUE, sep = ",")
GM07541_FMR1_window40000 <- read.table("G87-GM07541_FMR1_bases_window40000.txt", header = TRUE, sep = ",")
GM20239_FMR1_window40000 <- read.table("G87-GM20239_FMR1_bases_window40000.txt", header = TRUE, sep = ",")
GM04284_FMR1_window100000 <- read.table("G87-GM04284_FMR1_bases_window100000_AF_0_20_to_0_80.txt", header = TRUE, sep = ",")
GM05538_FMR1_window100000 <- read.table("G87-GM05538_FMR1_bases_window100000_AF_0_20_to_0_80.txt", header = TRUE, sep = ",")
GM07861_FMR1_window100000 <- read.table("G87-GM07861_FMR1_bases_window100000_AF_0_20_to_0_80.txt", header = FALSE, sep = ",") #If header is set to true all T-nucleotides are set to TRUE. It's an error
GM04284_HTT_window120000 <- read.table("G87-GM04284_HTT_bases_window120000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",")
GM06891_FMR1_window200000 <- read.table("G87-GM06891_FMR1_bases_window200000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",")

### Raw sequence logo containing all data (both alleles) ########################################################################
## Using ggseqlogo and character vectors

## Create tables with zeros instead of NA and deletions (-)
# GM04284 HTT Window 120000
GM04284_HTT_window120000[is.na(GM04284_HTT_window120000)] <- 0
GM04284_HTT_window120000[GM04284_HTT_window120000 == "-"] <- 0
# GM05538 HTT Window 40000
GM05538_HTT_window40000[is.na(GM05538_HTT_window40000)] <- 0
GM05538_HTT_window40000[GM05538_HTT_window40000 == "-"] <- 0
# GM06891 HTT window 40000
GM06891_HTT_window40000[is.na(GM06891_HTT_window40000)] <- 0
GM06891_HTT_window40000[GM06891_HTT_window40000 == "-"] <- 0
# GM07541 HTT window 40000
GM07541_HTT_window40000[is.na(GM07541_HTT_window40000)] <- 0
GM07541_HTT_window40000[GM07541_HTT_window40000 == "-"] <- 0
# GM07861 HTT window 40000
GM07861_HTT_window40000[is.na(GM07861_HTT_window40000)] <- 0
GM07861_HTT_window40000[GM07861_HTT_window40000 == "-"] <- 0
# GM20239 HTT window 40000
GM20239_HTT_window40000[is.na(GM20239_HTT_window40000)] <- 0
GM20239_HTT_window40000[GM20239_HTT_window40000 == "-"] <- 0
# GM04284 FMR1 window 100000
GM04284_FMR1_window100000[is.na(GM04284_FMR1_window100000)] <- 0
GM04284_FMR1_window100000[GM04284_FMR1_window100000 == "-"] <- 0
# GM05538 FMR1 window 100000
GM05538_FMR1_window100000[is.na(GM05538_FMR1_window100000)] <- 0
GM05538_FMR1_window100000[GM05538_FMR1_window100000 == "-"] <- 0
# GM06891 FMR1 window 200000
GM06891_FMR1_window200000[is.na(GM06891_FMR1_window200000)] <- 0
GM06891_FMR1_window200000[GM06891_FMR1_window200000 == "-"] <- 0
# GM07541 FMR1 window 40000
GM07541_FMR1_window40000[is.na(GM07541_FMR1_window40000)] <- 0
GM07541_FMR1_window40000[GM07541_FMR1_window40000 == "-"] <- 0
# GM07861 FMR1 window 100000
GM07861_FMR1_window100000[is.na(GM07861_FMR1_window100000)] <- 0
GM07861_FMR1_window100000[GM07861_FMR1_window100000 == "-"] <- 0
# GM20239 FMR1 window 40000
GM20239_FMR1_window40000[is.na(GM20239_FMR1_window40000)] <- 0
GM20239_FMR1_window40000[GM20239_FMR1_window40000 == "-"] <- 0

## Make a character vector for the sequences to be added to #######################################################################
# GM0284 HTT window 120000
Character_vector_GM04284_HTT_window120000 = character(nrow(GM04284_HTT_window120000))
# GM05538 HTT window 40000
Character_vector_GM05538_HTT_window40000 = character(nrow(GM05538_HTT_window40000))
# GM06891 HTT window 40000
Character_vector_GM06891_HTT_window40000 = character(nrow(GM06891_HTT_window40000))
# GM07541 HTT window 40000
Character_vector_GM07541_HTT_window40000 = character(nrow(GM07541_HTT_window40000))
# GM07861 HTT window 40000
Character_vector_GM07861_HTT_window40000 = character(nrow(GM07861_HTT_window40000))
# GM20239 HTT window 40000
Character_vector_GM20239_HTT_window40000 = character(nrow(GM20239_HTT_window40000))
# GM04284 FMR1 window 100000
Character_vector_GM04284_FMR1_window100000 = character(nrow(GM04284_FMR1_window100000))
# GM05538 FMR1 window 100000
Character_vector_GM05538_FMR1_window100000 = character(nrow(GM05538_FMR1_window100000))
# GM06891 FMR1 window 200000
Character_vector_GM06891_FMR1_window200000 = character(nrow(GM06891_FMR1_window200000))
# GM07541 FMR1 window 40000
Character_vector_GM07541_FMR1_window40000 = character(nrow(GM07541_FMR1_window40000))
# GM07861 FMR1 window 100000
Character_vector_GM07861_FMR1_window100000 = character(nrow(GM07861_FMR1_window100000))
# GM20239 FMR1 window 40000
Character_vector_GM20239_FMR1_window40000 = character(nrow(GM20239_FMR1_window40000))

## Add the sequences to the character vector ######################################################################################
# GM04284 HTT
i <- 1
while(i < (nrow(GM04284_HTT_window120000)+1)) {
  Character_vector_GM04284_HTT_window120000[i] <- paste(c(as.character(as.character(GM04284_HTT_window120000[i, 2:ncol(GM04284_HTT_window120000)]))),collapse = "")
  i <- i + 1
}
# GM05538 HTT
i <- 1
while(i < (nrow(GM05538_HTT_window40000)+1)) {
  Character_vector_GM05538_HTT_window40000[i] <- paste(c(as.character(as.character(GM05538_HTT_window40000[i, 2:ncol(GM05538_HTT_window40000)]))),collapse = "")
  i <- i + 1
}
# GM06891 HTT
i <- 1
while(i < (nrow(GM06891_HTT_window40000)+1)) {
  Character_vector_GM06891_HTT_window40000[i] <- paste(c(as.character(as.character(GM06891_HTT_window40000[i, 2:ncol(GM06891_HTT_window40000)]))),collapse = "")
  i <- i + 1
}
# GM07541 HTT
i <- 1
while(i < (nrow(GM07541_HTT_window40000)+1)) {
  Character_vector_GM07541_HTT_window40000[i] <- paste(c(as.character(as.character(GM07541_HTT_window40000[i, 2:ncol(GM07541_HTT_window40000)]))),collapse = "")
  i <- i + 1
}
# GM07861 HTT
i <- 1
while(i < (nrow(GM07861_HTT_window40000)+1)) {
  Character_vector_GM07861_HTT_window40000[i] <- paste(c(as.character(as.character(GM07861_HTT_window40000[i, 2:ncol(GM07861_HTT_window40000)]))),collapse = "")
  i <- i + 1
}
# GM20239 HTT
i <- 1
while(i < (nrow(GM20239_HTT_window40000)+1)) {
  Character_vector_GM20239_HTT_window40000[i] <- paste(c(as.character(as.character(GM20239_HTT_window40000[i, 2:ncol(GM20239_HTT_window40000)]))),collapse = "")
  i <- i + 1
}
# GM04284 FMR1
i <- 1
while(i < (nrow(GM04284_FMR1_window100000)+1)) {
  Character_vector_GM04284_FMR1_window100000[i] <- paste(c(as.character(as.character(GM04284_FMR1_window100000[i, 2:ncol(GM04284_FMR1_window100000)]))),collapse = "")
  i <- i + 1
}
# GM05528 FMR1
i <- 1
while(i < (nrow(GM05538_FMR1_window100000)+1)) {
  Character_vector_GM05538_FMR1_window100000[i] <- paste(c(as.character(as.character(GM05538_FMR1_window100000[i, 2:ncol(GM05538_FMR1_window100000)]))),collapse = "")
  i <- i + 1
}
# GM06891 FMR1
i <- 1
while(i < (nrow(GM06891_FMR1_window200000)+1)) {
  Character_vector_GM06891_FMR1_window200000[i] <- paste(c(as.character(as.character(GM06891_FMR1_window200000[i, 2:ncol(GM06891_FMR1_window200000)]))),collapse = "")
  i <- i + 1
}
# GM07541 FMR1
i <- 1
while(i < (nrow(GM07541_FMR1_window40000)+1)) {
  Character_vector_GM07541_FMR1_window40000[i] <- paste(c(as.character(as.character(GM07541_FMR1_window40000[i, 2:ncol(GM07541_FMR1_window40000)]))),collapse = "")
  i <- i + 1
}
# GM07861 FMR1
i <- 1
while(i < (nrow(GM07861_FMR1_window100000)+1)) {
  Character_vector_GM07861_FMR1_window100000[i] <- paste(c(as.character(as.character(GM07861_FMR1_window100000[i, 2:ncol(GM07861_FMR1_window100000)]))),collapse = "")
  i <- i + 1
}
# GM20239 FMR1
i <- 1
while(i < (nrow(GM20239_FMR1_window40000)+1)) {
  Character_vector_GM20239_FMR1_window40000[i] <- paste(c(as.character(as.character(GM20239_FMR1_window40000[i, 2:ncol(GM20239_FMR1_window40000)]))),collapse = "")
  i <- i + 1
}



## Create the sequence logo #############################################################################################################################################
# GM04284 HTT
ggseqlogo(Character_vector_GM04284_HTT_window120000, method = 'bits', seq_type = 'dna')
# GM05538 HTT
ggseqlogo(Character_vector_GM05538_HTT_window40000, method = 'bits', seq_type = 'dna')
# GM06891 HTT
ggseqlogo(Character_vector_GM06891_HTT_window40000, method = 'bits', seq_type = 'dna')
# GM07541 HTT
ggseqlogo(Character_vector_GM07541_HTT_window40000, method = 'bits', seq_type = 'dna')
# GM07861 HTT
ggseqlogo(Character_vector_GM07861_HTT_window40000, method = 'bits', seq_type = 'dna')
# GM20239 HTT
ggseqlogo(Character_vector_GM20239_HTT_window40000, method = 'bits', seq_type = 'dna')
# GM04284 FMR1
ggseqlogo(Character_vector_GM04284_FMR1_window100000, method = 'bits', seq_type = 'dna')
# GM05538 FMR1
ggseqlogo(Character_vector_GM05538_FMR1_window100000, method = 'bits', seq_type = 'dna')
# GM06891 FMR1
ggseqlogo(Character_vector_GM06891_FMR1_window200000, method = 'bits', seq_type = 'dna')
# GM07541 FMR1
ggseqlogo(Character_vector_GM07541_FMR1_window40000, method = 'bits', seq_type = 'dna')
# GM07861 FMR1
ggseqlogo(Character_vector_GM07861_FMR1_window100000[2:10], method = 'bits', seq_type = 'dna')
# GM20239 FMR1
ggseqlogo(Character_vector_GM20239_FMR1_window40000, method = 'bits', seq_type = 'dna')


## Modify the sequence logo #############################################################################################################################################
# GM05538 HTT
ggseqlogo(Character_vector_GM05538_HTT_window40000, method = 'bits', seq_type = 'dna') + 
  annotate('text', x=12.5, y=1.15, label='Total Data') + 
  annotate('text', x=22.5, y=1.15, label='GM05538') + 
  annotate('text', x=32.5, y=1.15, label='Window 40000') + 
  annotate('rect', xmin = 2.5, xmax = 4.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 5.5, xmax = 8.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 11.5, xmax = 12.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 17.5, xmax = 18.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 20.5, xmax = 21.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 23.5, xmax = 26.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 27.5, xmax = 35.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 36.5, xmax = 38.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 39.5, xmax = 41.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 42.5, xmax = 45.5, ymin = -0.025, ymax = 1.075, alpha = .1, col='black', fill='yellow')

# GM06891 HTT
#   1:6, 8:17, 20:23, 26:31, 33:35, 37:41, 43:48, 51:57,  59:61, 64:66, 69:70
ggseqlogo(Character_vector_GM06891_HTT_window40000, method = 'bits', seq_type = 'dna') + 
  annotate('text', x=15.5, y=1.15, label='Total Data') + 
  annotate('text', x=35.5, y=1.15, label='GM06891') + 
  annotate('text', x=55.5, y=1.15, label='Window 40000') + 
  annotate('rect', xmin = 0.5, xmax = 6.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 7.5, xmax = 17.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 19.5, xmax = 23.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 25.5, xmax = 31.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 32.5, xmax = 35.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 36.5, xmax = 41.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 42.5, xmax = 48.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 50.5, xmax = 57.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 58.5, xmax = 61.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 63.5, xmax = 66.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 68.5, xmax = 70.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow')


# GM07541 HTT
#   4, 6:7, 9, 11:12, 14:16, 18, 20:21, 23:24, 26, 32, 34:35, 40, 43:53, 55:57, 59:60, 62:66, 68:69, 71
ggseqlogo(Character_vector_GM07541_HTT_window40000, method = 'bits', seq_type = 'dna') + 
  annotate('text', x=15.5, y=1.15, label='Total Data') + 
  annotate('text', x=35.5, y=1.15, label='GM07541') + 
  annotate('text', x=55.5, y=1.15, label='Window 40000') + 
  annotate('rect', xmin = 3.5, xmax = 4.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 5.5, xmax = 7.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 8.5, xmax = 9.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 10.5, xmax = 12.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 13.5, xmax = 16.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 17.5, xmax = 18.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 19.5, xmax = 21.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 22.5, xmax = 24.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 25.5, xmax = 26.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 31.5, xmax = 32.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 33.5, xmax = 35.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 39.5, xmax = 40.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 42.5, xmax = 53.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 54.5, xmax = 57.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 58.5, xmax = 60.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 61.5, xmax = 66.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 67.5, xmax = 69.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 70.5, xmax = 71.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow')


# GM06891 HTT
#   1:6, 8:17, 20:23, 26:31, 33:35, 37:41, 43:48, 51:57,  59:61, 64:66, 69:70
ggseqlogo(Character_vector_GM06891_HTT_window40000, method = 'bits', seq_type = 'dna') + 
  annotate('text', x=15.5, y=1.15, label='Total Data') + 
  annotate('text', x=35.5, y=1.15, label='GM06891') + 
  annotate('text', x=55.5, y=1.15, label='Window 40000') + 
  annotate('rect', xmin = 0.5, xmax = 6.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 7.5, xmax = 17.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 19.5, xmax = 23.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 25.5, xmax = 31.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 32.5, xmax = 35.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 36.5, xmax = 41.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 42.5, xmax = 48.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 50.5, xmax = 57.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 58.5, xmax = 61.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 63.5, xmax = 66.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 68.5, xmax = 70.5, ymin = -0.025, ymax = 1.5, alpha = .1, col='black', fill='yellow')


# GM20239 HTT
ggseqlogo(Character_vector_GM20239_HTT_window40000, method = 'bits', seq_type = 'dna') + 
  annotate('text', x=15.5, y=1.5, label='Total Data') + 
  annotate('text', x=35.5, y=1.5, label='GM20239') + 
  annotate('text', x=55.5, y=1.5, label='Window 40000') + 
  annotate('rect', xmin = 4.5, xmax = 5.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 6.5, xmax = 7.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 8.5, xmax = 13.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 14.5, xmax = 15.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 16.5, xmax = 17.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 18.5, xmax = 19.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 20.5, xmax = 23.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 24.5, xmax = 26.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 29.5, xmax = 34.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 35.5, xmax = 36.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 37.5, xmax = 38.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 39.5, xmax = 49.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 50.5, xmax = 54.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 56.5, xmax = 64.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') 

# GM20239 FMR1
ggseqlogo(Character_vector_GM20239_FMR1_window40000, method = 'bits', seq_type = 'dna') + 
  annotate('text', x=45.5, y=1.75 , label='Total Data') + 
  annotate('text', x=65.5, y=1.75, label='GM20239') + 
  annotate('text', x=85.5, y=1.75, label='Window 40000') + 
  annotate('rect', xmin = 2.5, xmax = 3.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 4.5, xmax = 7.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 8.5, xmax = 9.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 10.5, xmax = 11.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 12.5, xmax = 22.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 25.5, xmax = 28.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 29.5, xmax = 39.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 40.5, xmax = 43.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 44.5, xmax = 48.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 49.5, xmax = 59.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 60.5, xmax = 61.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 62.5, xmax = 64.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 67.5, xmax = 73.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 74.5, xmax = 79.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 80.5, xmax = 83.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 84.5, xmax = 85.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 87.5, xmax = 91.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 94.5, xmax = 99.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 102.5, xmax = 104.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 105.5, xmax = 108.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 109.5, xmax = 111.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 113.5, xmax = 116.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 117.5, xmax = 118.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 119.5, xmax = 121.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') +
  annotate('rect', xmin = 122.5, xmax = 132.5, ymin = -0.025, ymax = 2, alpha = .1, col='black', fill='yellow') 


## Remove SNV's from the table that are not good enough for HP ##########################################################################################################
# GM05538 HTT
#   Positions kept for further processing: 3:4, 6:8, 12, 18, 21, 24:26, 28:35, 37:38, 40:41, 43:45
kept_SNVs <- c(3:4, 6:8, 12, 18, 21, 24:26, 28:35, 37:38, 40:41, 43:45)
GM05538_HTT_window40000_subset_0.75_bits <- subset(GM05538_HTT_window40000[, kept_SNVs+1]) 
Character_vector_GM05538_HTT_window40000_subset_0.75_bits = character(nrow(GM05538_HTT_window40000))

z <- ncol(GM05538_HTT_window40000_subset_0.75_bits)
i <- 1
while(i < (nrow(GM05538_HTT_window40000)+1)) {
  Character_vector_GM05538_HTT_window40000_subset_0.75_bits[i] <- paste(c(as.character(as.character(GM05538_HTT_window40000_subset_0.75_bits[i, 1:z]))),collapse = "")
  i <- i + 1
}

ggseqlogo(Character_vector_GM05538_HTT_window40000_subset_0.75_bits, method = 'bits', seq_type = 'dna')


# GM06891 HTT
#   Positions kept for further processing: 1:6, 8:17, 20:23, 26:31, 33:35, 37:41, 43:48, 51:57,  59:61, 64:66, 69:70.
kept_SNVs <- c(1:6, 8:17, 20:23, 26:31, 33:35, 37:41, 43:48, 51:57,  59:61, 64:66, 69:70)
GM06891_HTT_window40000_subset_0.75_bits <- subset(GM06891_HTT_window40000[, kept_SNVs+1]) 
Character_vector_GM06891_HTT_window40000_subset_0.75_bits = character(nrow(GM06891_HTT_window40000))
z <- ncol(GM06891_HTT_window40000_subset_0.75_bits)
i <- 1
while(i < (nrow(GM06891_HTT_window40000)+1)) {
  Character_vector_GM06891_HTT_window40000_subset_0.75_bits[i] <- paste(c(as.character(as.character(GM06891_HTT_window40000_subset_0.75_bits[i, 1:z]))),collapse = "")
  i <- i + 1
}

ggseqlogo(Character_vector_GM06891_HTT_window40000_subset_0.75_bits, method = 'bits', seq_type = 'dna')


# GM07541 HTT
#   Positions kept for further processing: 1:6, 8:17, 20:23, 26:31, 33:35, 37:41, 43:48, 51:57,  59:61, 64:66, 69:70.
kept_SNVs <- c(4, 6:7, 9, 11:12, 14:16, 18, 20:21, 23:24, 26, 32, 34:35, 40, 43:53, 55:57, 59:60,
               62:66, 68:69, 71)
GM07541_HTT_window40000_subset_0.75_bits <- subset(GM07541_HTT_window40000[, kept_SNVs+1]) 
Character_vector_GM07541_HTT_window40000_subset_0.75_bits = character(nrow(GM07541_HTT_window40000))
z <- ncol(GM07541_HTT_window40000_subset_0.75_bits)
i <- 1
while(i < (nrow(GM07541_HTT_window40000)+1)) {
  Character_vector_GM07541_HTT_window40000_subset_0.75_bits[i] <- paste(c(as.character(as.character(GM07541_HTT_window40000_subset_0.75_bits[i, 1:z]))),collapse = "")
  i <- i + 1
}

ggseqlogo(Character_vector_GM07541_HTT_window40000_subset_0.75_bits, method = 'bits', seq_type = 'dna')

# GM07861 HTT
#   Positions kept for further processing: 1:5
kept_SNVs <- c(1:4)
GM07861_HTT_window40000_subset <- subset(GM07861_HTT_window40000[, kept_SNVs+1]) 
Character_vector_GM07861_HTT_window40000_subset = character(nrow(GM07861_HTT_window40000))
z <- ncol(GM07861_HTT_window40000_subset)
i <- 1
while(i < (nrow(GM07861_HTT_window40000)+1)) {
  Character_vector_GM07861_HTT_window40000_subset[i] <- paste(c(as.character(as.character(GM07861_HTT_window40000_subset[i, 1:z]))),collapse = "")
  i <- i + 1
}

ggseqlogo(Character_vector_GM07861_HTT_window40000_subset, method = 'bits', seq_type = 'dna')



# GM20239 HTT
#   Positions kept for further processing: 5, 7, 9:13, 15, 17, 19, 21:23, 25.26, 30:34, 36, 38, 40:49, 51:54, 57:64
kept_SNVs <- c(5, 7, 9:13, 15, 17, 19, 21:23, 25.26, 30:34, 36, 38, 40:49, 51:54, 57:64)
GM20239_HTT_window40000_subset_0.75_bits <- subset(GM20239_HTT_window40000[, kept_SNVs+1]) 
Character_vector_GM20239_HTT_window40000_subset_0.75_bits = character(nrow(GM20239_HTT_window40000))
z <- ncol(GM20239_HTT_window40000_subset_0.75_bits)
i <- 1
while(i < (nrow(GM20239_HTT_window40000)+1)) {
  Character_vector_GM20239_HTT_window40000_subset_0.75_bits[i] <- paste(c(as.character(as.character(GM20239_HTT_window40000_subset_0.75_bits[i, 1:z]))),collapse = "")
  i <- i + 1
}

ggseqlogo(Character_vector_GM20239_HTT_window40000_subset_0.75_bits, method = 'bits', seq_type = 'dna')


# GM05538 FMR1
#   Positions kept for further processing: 1:2
kept_SNVs <- c(1:2)
GM05538_FMR1_window100000_subset <- subset(GM05538_FMR1_window100000[, kept_SNVs+1]) 
Character_vector_GM05538_FMR1_window100000_subset = character(nrow(GM05538_FMR1_window100000))
z <- ncol(GM05538_FMR1_window100000_subset)
i <- 1
while(i < (nrow(GM05538_FMR1_window100000)+1)) {
  Character_vector_GM05538_FMR1_window100000_subset[i] <- paste(c(as.character(as.character(GM05538_FMR1_window100000_subset[i, 1:z]))),collapse = "")
  i <- i + 1
}

ggseqlogo(Character_vector_GM05538_FMR1_window100000_subset, method = 'bits', seq_type = 'dna')


# GM07541 FMR1
#   Positions kept for further processing: 1:6
kept_SNVs <- c(1:6)
GM07541_FMR1_window40000_subset <- subset(GM07541_FMR1_window40000[, kept_SNVs+1]) 
Character_vector_GM07541_FMR1_window40000_subset = character(nrow(GM07541_FMR1_window40000))
z <- ncol(GM07541_FMR1_window40000_subset)
i <- 1
while(i < (nrow(GM07541_FMR1_window40000)+1)) {
  Character_vector_GM07541_FMR1_window40000_subset[i] <- paste(c(as.character(as.character(GM07541_FMR1_window40000_subset[i, 1:z]))),collapse = "")
  i <- i + 1
}

ggseqlogo(Character_vector_GM07541_FMR1_window40000_subset, method = 'bits', seq_type = 'dna')


# GM20239 FMR1
#   Positions kept for further processing: 3, 5:7, 9, 11, 13:22, 26:28, 30:39, 41:43, 45:48, 50:59, 61, 63:64, 68:73, 
#   75:79, 81:83, 85, 88:91, 95:99, 103:104, 106:108, 110:111, 114:116, 118, 120:121, 123:132
kept_SNVs <- c(3, 5:7, 9, 11, 13:22, 26:28, 30:39, 41:43, 45:48, 50:59, 61, 63:64, 68:73, 
               75:79, 81:83, 85, 88:91, 95:99, 103:104, 106:108, 110:111, 114:116, 118, 120:121, 123:132)
GM20239_FMR1_window40000_subset_0.75_bits <- subset(GM20239_FMR1_window40000[, kept_SNVs+1])
Character_vector_GM20239_FMR1_window40000_subset_0.75_bits = character(nrow(GM20239_FMR1_window40000))

z <- ncol(GM20239_FMR1_window40000_subset_0.75_bits)
i <- 1
while(i < (nrow(GM20239_FMR1_window40000)+1)){
  Character_vector_GM20239_FMR1_window40000_subset_0.75_bits[i] <- paste(c(as.character(GM20239_FMR1_window40000_subset_0.75_bits[i, 1:z])),collapse = "")
  i <- i + 1
}

ggseqlogo(Character_vector_GM20239_FMR1_window40000_subset_0.75_bits, method = 'bits', seq_type = 'dna')



## Making the consensus sequence for the two alleles
##   Open the table in excel and go through each column and assign the reads one by one to an allele which
##   are given two different colors for visual clarity.
##   When allt he columns have been analysed, the row numbers assigned the two alleles should be given below.
##   Then new logo plots are made.
# GM04284 HTT
write.csv(GM04284_HTT_window120000, file = "GM04284_HTT_window120000_subset.txt")
allele_1_rows_GM04284_HTT <- c(1, 3, 7)
Character_vector_GM04284_HTT_window120000_allele_1 <- Character_vector_GM04284_HTT_window120000[allele_1_rows_GM04284_HTT]
Character_vector_GM04284_HTT_window120000_allele_1
allele_2_rows_GM04284_HTT <- c(2, 4:5, 8)
Character_vector_GM04284_HTT_window120000_allele_2 <- Character_vector_GM04284_HTT_window120000[allele_2_rows_GM04284_HTT]
Character_vector_GM04284_HTT_window120000_allele_2
ggseqlogo(Character_vector_GM04284_HTT_window120000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM04284_HTT_window120000_allele_2, method = 'bits', seq_type = 'dna')

# GM05538 HTT
write.csv(GM05538_HTT_window40000_subset_0.75_bits, file = "GM05538_HTT_window40000_subset_0.75_bits.txt")
allele_1_rows_GM05538_HTT <- c(1:6, 8, 12, 15, 18:19, 21, 27:29, 31, 33, 39, 47:49, 51:52, 57:59,
                               64, 66, 69, 73:75, 78, 80, 82, 84:86, 88:89, 92:94, 97, 100, 102:106,
                               108, 111:112, 116, 119, 121, 123:124, 127, 129:131, 135, 137)
Character_vector_GM05538_HTT_window40000_allele_1 <- Character_vector_GM05538_HTT_window40000[allele_1_rows_GM05538_HTT]
Character_vector_GM05538_HTT_window40000_allele_1
allele_2_rows_GM05538_HTT <- c(7, 11, 16:17, 20, 22:23, 30, 36, 38, 40:42, 44:45, 50, 53, 55:56,
                               60:62, 65, 67:68, 71:77, 79, 81, 83, 91, 95:96, 98:99, 101, 107,
                               109:110, 113:115, 117:118, 120, 122, 125:126, 128, 132)
Character_vector_GM05538_HTT_window40000_allele_2 <- Character_vector_GM05538_HTT_window40000[allele_2_rows_GM05538_HTT]
Character_vector_GM05538_HTT_window40000_allele_2
ggseqlogo(Character_vector_GM05538_HTT_window40000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM05538_HTT_window40000_allele_2, method = 'bits', seq_type = 'dna')

# GM06891 HTT
write.csv(GM06891_HTT_window40000_subset_0.75_bits, file = "GM06891_HTT_window40000_subset_0_75_bits.txt")
allele_1_rows_GM06891_HTT <- c(1, 4, 6:7, 10:12, 14:15, 17, 20:21, 23:24, 27, 29, 31, 33:34, 37, 39:43,
                               47, 49:50, 52:54, 56:57, 60:62, 64:68, 71, 74:76, 80:81, 84:86, 88:89,
                               94, 98:99, 103, 106, 109:110, 115:116, 118, 121:125, 131:132, 136:137,
                               140, 147, 154, 158:159, 163, 166, 168, 170:172, 175:176, 178, 184:186,
                               189:193, 196, 200:203, 205, 209:211, 213:214, 216:217, 222, 225, 228:229,
                               242, 246:248, 251:252, 258, 265:268: 270:273, 276:278, 280, 282:284,
                               286:288, 292, 298:299, 301, 305:306, 311:313)
Character_vector_GM06891_HTT_window40000_allele_1 <- Character_vector_GM06891_HTT_window40000[allele_1_rows_GM06891_HTT]
Character_vector_GM06891_HTT_window40000_allele_1
allele_2_rows_GM06891_HTT <- c(2:3, 5, 8:9, 13, 16, 18:19, 22, 25:26, 32, 35, 38, 44:46, 48, 51, 58:59,
                               63, 69:70, 72, 77:79, 82:83, 87, 90:92, 95:97, 100:102, 104:105, 107:108,
                               111:114, 117, 119:120, 126:130, 133:135, 138:139, 141:146, 148:151, 153,
                               155:157, 160:162, 164:165, 167, 173:174, 177, 179:181, 187:188, 194:195,
                               197:199, 204, 206:208, 212, 215, 218:221, 226, 231:233, 237:238, 243:245,
                               249:250, 253:257, 259, 261:264, 269, 274, 279, 281, 285, 289:291, 295:297,
                               300, 304, 308:310)
Character_vector_GM06891_HTT_window40000_allele_2 <- Character_vector_GM06891_HTT_window40000[allele_2_rows_GM06891_HTT]
Character_vector_GM06891_HTT_window40000_allele_2
ggseqlogo(Character_vector_GM06891_HTT_window40000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM06891_HTT_window40000_allele_2, method = 'bits', seq_type = 'dna')


# GM07541 HTT
write.csv(GM07541_HTT_window40000_subset_0.75_bits, file = "GM07541_HTT_window40000_subset_0_75_bits.txt")
allele_1_rows_GM07541_HTT <- c(1:8, 12:14, 17, 20:21, 25, 28:29, 32, 35, 37:39, 42, 45:49, 53, 57:58,
                               60, 64:66, 72, 76, 78:80, 85, 91, 94, 98, 102, 109:110, 112, 115, 117,
                               119, 121:123, 125:126, 129:131, 133:136, 139:140, 142:144, 147:148, 150:151,
                               153, 156, 160)
Character_vector_GM07541_HTT_window40000_allele_1 <- Character_vector_GM07541_HTT_window40000[allele_1_rows_GM07541_HTT]
Character_vector_GM07541_HTT_window40000_allele_1
allele_2_rows_GM07541_HTT <- c(10:11, 15:16, 18:19, 23:24, 26:27, 30:31, 33:34, 36, 40:41, 43:44,
                               50, 54:56, 59, 63, 67:71, 74, 81, 83:84, 86:90, 92:93, 95:97,
                               99:101, 103, 106:108, 111, 113:114, 116, 118, 120, 124, 132,
                               137:138, 141, 145:146, 149, 152, 154:155, 157:159)
Character_vector_GM07541_HTT_window40000_allele_2 <- Character_vector_GM07541_HTT_window40000[allele_2_rows_GM07541_HTT]
Character_vector_GM07541_HTT_window40000_allele_2
ggseqlogo(Character_vector_GM07541_HTT_window40000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM07541_HTT_window40000_allele_2, method = 'bits', seq_type = 'dna')

# GM07861 HTT
write.csv(GM07861_HTT_window40000_subset, file = "GM07861_HTT_window40000_subset.txt")
allele_1_rows_GM07861_HTT <- c(1, 5, 9:10, 14:15, 17:18, 20, 27:32, 36, 38, 40:41)
Character_vector_GM07861_HTT_window40000_allele_1 <- Character_vector_GM07861_HTT_window40000[allele_1_rows_GM07861_HTT]
Character_vector_GM07861_HTT_window40000_allele_1
allele_2_rows_GM07861_HTT <- c(2:4, 6:8, 11:13, 16, 19, 21:26, 33:35, 37, 39, 42:43, 45)
Character_vector_GM07861_HTT_window40000_allele_2 <- Character_vector_GM07861_HTT_window40000[allele_2_rows_GM07861_HTT]
Character_vector_GM07861_HTT_window40000_allele_2
ggseqlogo(Character_vector_GM07861_HTT_window40000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM07861_HTT_window40000_allele_2, method = 'bits', seq_type = 'dna')

# GM20239 HTT
write.csv(GM20239_HTT_window40000_subset_0.75_bits, file = "GM20239_HTT_window40000_subset_0_75_bits.txt")
allele_1_rows_GM20239_HTT <- c(1, 3, 6, 8:9, 14:15, 19, 24, 29:30, 32:33, 35, 38:39, 41, 43, 47,
                               49:50, 52, 54:58, 61:66, 68:70, 73:74, 77, 79, 84:85, 87, 91:92,
                               94:98, 100, 105, 113:114, 116, 118, 122:124, 129:135, 138:139,
                               146:147, 149:150, 157, 159:160, 162, 165, 168, 173:174, 176:177,
                               180, 186, 189, 191, 193:197, 199:202, 204:209, 211, 216, 220, 222,
                               224:225, 228:229, 231:232, 236, 238:240, 242, 245, 251, 253, 256,
                               258, 260:261, 264:266, 268, 272, 274, 294, 298:299, 301:303, 307:311,
                               313:317, 319, 323:325, 327:328, 332, 334)
Character_vector_GM20239_HTT_window40000_allele_1 <- Character_vector_GM20239_HTT_window40000[allele_1_rows_GM20239_HTT]
Character_vector_GM20239_HTT_window40000_allele_1
allele_2_rows_GM20239_HTT <- c(4:5, 7, 11:13, 16, 18, 22, 25:27, 31, 34, 36:37, 40, 42, 45:46,
                               48, 51, 53, 59:60, 67, 71:72, 78, 80:81, 83, 86, 88:90, 99, 101,
                               106, 111, 117, 119, 121, 125:127, 136:137, 140:145, 148, 151:156,
                               158, 161, 163, 171, 175, 179, 182, 184:185, 187:188, 190, 198,
                               203, 210, 212, 214:215, 218:219, 221, 223, 227, 230, 233:235,
                               237, 243:244, 247:250, 252, 254:255, 257, 259, 262:263, 267, 269:271,
                               278, 282, 293, 300, 304:306, 312, 318, 320:322, 326, 329:331,
                               336, 338:339)
Character_vector_GM20239_HTT_window40000_allele_2 <- Character_vector_GM20239_HTT_window40000[allele_2_rows_GM20239_HTT]
Character_vector_GM20239_HTT_window40000_allele_2
ggseqlogo(Character_vector_GM20239_HTT_window40000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM20239_HTT_window40000_allele_2, method = 'bits', seq_type = 'dna')

# GM05538 FMR1
write.csv(GM05538_FMR1_window100000_subset, file = "GM05538_FMR1_window100000_subset.txt")
allele_1_rows_GM05538_FMR1 <- c(16:17, 27)
Character_vector_GM05538_FMR1_window100000_allele_1 <- Character_vector_GM05538_FMR1_window100000[allele_1_rows_GM05538_FMR1]
Character_vector_GM07541_FMR1_window40000_allele_1
allele_2_rows_GM05538_FMR1 <- c(19:24, 26, 28)
Character_vector_GM05538_FMR1_window100000_allele_2 <- Character_vector_GM05538_FMR1_window100000[allele_2_rows_GM05538_FMR1]
Character_vector_GM05538_FMR1_window100000_allele_2
ggseqlogo(Character_vector_GM05538_FMR1_window100000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM05538_FMR1_window100000_allele_2, method = 'bits', seq_type = 'dna')

# GM07541 FMR1
write.csv(GM07541_FMR1_window40000_subset, file = "GM07541_FMR1_window40000_subset.txt")
allele_1_rows_GM07541_FMR1 <- c(1:3, 5:11, 17, 20, 22, 25, 28:29, 31, 34, 36, 40, 42, 48, 50, 59, 61:62, 64, 66, 69)
Character_vector_GM07541_FMR1_window40000_allele_1 <- Character_vector_GM07541_FMR1_window40000[allele_1_rows_GM07541_FMR1]
Character_vector_GM07541_FMR1_window40000_allele_1
allele_2_rows_GM07541_FMR1 <- c(4, 12:16, 18:19, 21, 23:24, 26:27, 30, 32, 35, 37:39, 43:47, 51:58, 60, 63, 65, 67:68)
Character_vector_GM07541_FMR1_window40000_allele_2 <- Character_vector_GM07541_FMR1_window40000[allele_2_rows_GM07541_FMR1]
Character_vector_GM07541_FMR1_window40000_allele_2
ggseqlogo(Character_vector_GM07541_FMR1_window40000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM07541_FMR1_window40000_allele_2, method = 'bits', seq_type = 'dna')

# GM20239 FMR1
write.csv(GM20239_FMR1_window40000_subset_0.75_bits, file = "GM20239_FMR1_window40000_subset_0_75_bits.txt")
allele_1_rows_GM20239_FMR1 <- c(2, 7, 13, 16, 26, 35, 37:38, 40:41, 43:44, 46:51, 55:56, 62, 71:73, 81,
                                85:89, 92, 95:97, 100:101, 103:104, 106:108, 111, 115:117, 121:122, 124,
                                126, 130, 133:135, 137:138, 140:143, 145, 147:149, 152:156, 158, 160,
                                162, 167:171, 174:175, 177, 180, 182:183, 185:186, 188:192, 195, 197:199,
                                203:205, 209, 215, 227, 229:233, 236:238, 242, 245, 247:249, 251, 254:255,
                                257, 260:262, 268:270, 272, 276:277, 283, 286, 289:290, 292, 294:297,
                                302:303, 305, 307:308, 310, 313, 315:318)
Character_vector_GM20239_FMR1_window40000_allele_1 <- Character_vector_GM20239_FMR1_window40000[allele_1_rows_GM20239_FMR1]
Character_vector_GM20239_FMR1_window40000_allele_1
allele_2_rows_GM20239_FMR1 <- c(3:5, 9, 12, 19, 21, 29:30, 36, 45, 53, 57:58, 61, 69:70, 74:80, 82:84,
                                90:91, 93:94, 98:99, 102, 105, 112, 114, 118:120, 123, 125, 127:129, 131:132,
                                136, 139, 144, 146, 150:151, 157, 159, 161, 163:166, 172:173, 176, 178:179,
                                181, 184, 187, 193:194, 196, 200:202, 206:208, 210:214, 216, 221, 223,
                                225:226, 228, 239:241, 243, 246, 250, 253, 258:259, 263:264, 266:267,
                                171, 273, 275, 278:281, 284:285, 288, 291, 293, 298:301, 304, 306, 309, 311:312,
                                314, 319)
Character_vector_GM20239_FMR1_window40000_allele_2 <- Character_vector_GM20239_FMR1_window40000[allele_2_rows_GM20239_FMR1]
Character_vector_GM20239_FMR1_window40000_allele_2
ggseqlogo(Character_vector_GM20239_FMR1_window40000_allele_1, method = 'bits', seq_type = 'dna')
ggseqlogo(Character_vector_GM20239_FMR1_window40000_allele_2, method = 'bits', seq_type = 'dna')



## Create a new matrix with the read ID and the allele assignment.
# GM04284 HTT
colum_names <- c("Allele_Assignment")
GM04284_HTT_window120000_HP <- matrix(nrow = nrow(GM04284_HTT_window120000), ncol = 1)
colnames(GM04284_HTT_window120000_HP) <- colum_names
rownames(GM04284_HTT_window120000_HP) <- GM04284_HTT_window120000[,1]
GM04284_HTT_window120000_HP[allele_1_rows_GM04284_HTT, 1] <- 1
GM04284_HTT_window120000_HP[allele_2_rows_GM04284_HTT, 1] <- 2
write.csv(GM04284_HTT_window120000_HP, file = "GM04284_HTT_window120000_HP.txt")
table(GM04284_HTT_window120000_HP)



# GM05538 HTT
colum_names <- c("Allele_Assignment")
GM05538_HTT_window40000_HP <- matrix(nrow = nrow(GM05538_HTT_window40000), ncol = 1)
colnames(GM05538_HTT_window40000_HP) <- colum_names
rownames(GM05538_HTT_window40000_HP) <- GM05538_HTT_window40000[,1]
GM05538_HTT_window40000_HP[allele_1_rows_GM05538_HTT, 1] <- 1
GM05538_HTT_window40000_HP[allele_2_rows_GM05538_HTT, 1] <- 2
write.csv(GM05538_HTT_window40000_HP, file = "GM05538_HTT_window40000_HP_2.txt")
table(GM05538_HTT_window40000_HP)



# GM06891 HTT
colum_names <- c("Allele_Assignment")
GM06891_HTT_window40000_HP <- matrix(nrow = nrow(GM06891_HTT_window40000), ncol = 1)
colnames(GM06891_HTT_window40000_HP) <- colum_names
rownames(GM06891_HTT_window40000_HP) <- GM06891_HTT_window40000[,1]
GM06891_HTT_window40000_HP[allele_1_rows_GM06891_HTT, 1] <- 1
GM06891_HTT_window40000_HP[allele_2_rows_GM06891_HTT, 1] <- 2
write.csv(GM06891_HTT_window40000_HP, file = "GM06891_HTT_window40000_HP.txt")
table(GM06891_HTT_window40000_HP)


# GM07541 HTT
colum_names <- c("Allele_Assignment")
GM07541_HTT_window40000_HP <- matrix(nrow = nrow(GM07541_HTT_window40000), ncol = 1)
colnames(GM07541_HTT_window40000_HP) <- colum_names
rownames(GM07541_HTT_window40000_HP) <- GM07541_HTT_window40000[,1]
GM07541_HTT_window40000_HP[allele_1_rows_GM07541_HTT, 1] <- 1
GM07541_HTT_window40000_HP[allele_2_rows_GM07541_HTT, 1] <- 2
write.csv(GM07541_HTT_window40000_HP, file = "GM07541_HTT_window40000_HP.txt")
table(GM07541_HTT_window40000_HP)

# GM07861 HTT
colum_names <- c("Allele_Assignment")
GM07861_HTT_window40000_HP <- matrix(nrow = nrow(GM07861_HTT_window40000), ncol = 1)
colnames(GM07861_HTT_window40000_HP) <- colum_names
rownames(GM07861_HTT_window40000_HP) <- GM07861_HTT_window40000[,1]
GM07861_HTT_window40000_HP[allele_1_rows_GM07861_HTT, 1] <- 1
GM07861_HTT_window40000_HP[allele_2_rows_GM07861_HTT, 1] <- 2
write.csv(GM07861_HTT_window40000_HP, file = "GM07861_HTT_window40000_HP.txt")
table(GM07861_HTT_window40000_HP)

# GM20239 HTT
colum_names <- c("Allele_Assignment")
GM20239_HTT_window40000_HP <- matrix(nrow = nrow(GM20239_HTT_window40000), ncol = 1)
colnames(GM20239_HTT_window40000_HP) <- colum_names
rownames(GM20239_HTT_window40000_HP) <- GM20239_HTT_window40000[,1]
GM20239_HTT_window40000_HP[allele_1_rows_GM20239_HTT, 1] <- 1
GM20239_HTT_window40000_HP[allele_2_rows_GM20239_HTT, 1] <- 2
write.csv(GM20239_HTT_window40000_HP, file = "GM20239_HTT_window40000_HP.txt")
table(GM20239_HTT_window40000_HP)

# GM05538 FMR1
colum_names <- c("Allele_Assignment")
GM05538_FMR1_window100000_HP <- matrix(nrow = nrow(GM05538_FMR1_window100000), ncol = 1)
colnames(GM05538_FMR1_window100000_HP) <- colum_names
rownames(GM05538_FMR1_window100000_HP) <- GM05538_FMR1_window100000[,1]
GM05538_FMR1_window100000_HP[allele_1_rows_GM05538_FMR1, 1] <- 1
GM05538_FMR1_window100000_HP[allele_2_rows_GM05538_FMR1, 1] <- 2
write.csv(GM05538_FMR1_window100000_HP, file = "GM05538_FMR1_window100000_HP.txt")
table(GM05538_FMR1_window100000_HP)

# GM07541 FMR1
colum_names <- c("Allele_Assignment")
GM07541_FMR1_window40000_HP <- matrix(nrow = nrow(GM07541_FMR1_window40000), ncol = 1)
colnames(GM07541_FMR1_window40000_HP) <- colum_names
rownames(GM07541_FMR1_window40000_HP) <- GM07541_FMR1_window40000[,1]
GM07541_FMR1_window40000_HP[allele_1_rows_GM07541_FMR1, 1] <- 1
GM07541_FMR1_window40000_HP[allele_2_rows_GM07541_FMR1, 1] <- 2
write.csv(GM07541_FMR1_window40000_HP, file = "GM07541_FMR1_window40000_HP.txt")
table(GM07541_FMR1_window40000_HP)

# GM20239 FMR1
colum_names <- c("Allele_Assignment")
GM20239_FMR1_window40000_HP <- matrix(nrow = nrow(GM20239_FMR1_window40000), ncol = 1)
colnames(GM20239_FMR1_window40000_HP) <- colum_names
rownames(GM20239_FMR1_window40000_HP) <- GM20239_FMR1_window40000[,1]
GM20239_FMR1_window40000_HP[allele_1_rows_GM20239_FMR1, 1] <- 1
GM20239_FMR1_window40000_HP[allele_2_rows_GM20239_FMR1, 1] <- 2
write.csv(GM20239_FMR1_window40000_HP, file = "GM20239_FMR1_window40000_HP.txt")
table(GM20239_FMR1_window40000_HP)



##################################################################################################################
###################################  NOT COMPLETED AND NOT WORKING  ##############################################
##################################################################################################################
### K-mer Classification
## Making two matrices: one for the k-mers and one for the classification 
#GM05538_HTT_window40000_strings <- matrix(, ncol = 43, nrow = 198 )
#GM05538_HTT_window40000_class <- matrix(, ncol = 15, nrow = 198 )
#read_ids <- rownames(GM05538_HTT_window40000)
#rownames(GM05538_HTT_window40000_strings) <- read_ids
#rownames(GM05538_HTT_window40000_class) <- read_ids


## Filling the k-mer table 
#i <- 1
#j <- 1
#s <- 2
#e <- 4
#while (i <= nrow(GM05538_HTT_window40000_strings)){
#  if (j <= ncol(GM05538_HTT_window40000_strings)) {
#    X <-  subset(GM05538_HTT_window40000[i, s:e])
#    X <- apply(X, 1, paste, collapse = ",")
#    GM05538_HTT_window40000_strings[i,j] <- X
#    j <- j + 1
#    s <- s + 1
#    e <- e + 1
#  } else {
#    i <- i + 1
#    j <- 1
#    s <- 2
#    e <- 4
#  }
#}


## Filling of the classification table
# NAs and deletions ("-") are given as an x
# I want to cluster every column isolated and add the cluster assignments to the "class" table

# 1. Create a matrix that calculates the distance between 
#   the standard event types and the event levels from the data.

#dist_matrix <- stringdistmatrix(tolower(evl),ET, useNames=TRUE ,method = "dl")


########### This part is not completed and does not solve the HP ##########################
# Call the stringdistmatrix function and request 4 groups
#uniquemodels <- unique(as.character(GM05538_HTT_window40000_strings[,2]))
#distancemodels <- stringdistmatrix(uniquemodels, uniquemodels, method = "soundex")
#rownames(distancemodels) <- uniquemodels
#hc <- hclust(as.dist(distancemodels))

# Visualize the dendogram
#plot(hc)
#rect.hclust(hc, k = 2)

# Visualize the grouping
#dfClust <- data.frame(uniquemodels, cutree(hc, k = 2))
#names(dfClust) <- c('modelname', 'cluster')
#plot(table(dfClust$cluster))
#print(paste('Average number of models per cluster:', mean(table(dfClust$cluster))))

# Lets look at the top groups and see what the algorithm did:
#t <- table(dfClust$cluster)
#t <- cbind(t,t/length(dfClust$cluster))
#t <- t[order(t[,2],decreasing = TRUE),]
#p <- data.frame(factorName = rownames(t), binCount = t[,1], percentFound = t[,2])
#dfClust <- merge(x = dfClust, y = p, by.x = 'cluster', by.y = 'factorName', all.x = T)
#names(dfClust) <- c('cluster', 'modelname')
#head(dfClust[c('cluster', 'modelname')], 200)


#########################################################################################################
