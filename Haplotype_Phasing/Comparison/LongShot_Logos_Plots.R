#### This script creates logos plots for LongShot ########

#################################################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/logos_plots/longshot_snv_df")

## Load Libraries ###############################################################################################################
require(ggplot2)
require(ggseqlogo)

### Loading the table ###########################################################################################################
GM04284_HTT <- read.table("longshot_GM04284_HTT_het_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM04284_HTT)[1] <- c("read_id")

GM05538_HTT <- read.table("longshot_GM05538_HTT_het_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM05538_HTT)[1] <- c("read_id")

GM06891_HTT <- read.table("longshot_GM06891_HTT_het_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM06891_HTT)[1] <- c("read_id")

GM07541_HTT <- read.table("longshot_GM07541_HTT_het_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM07541_HTT)[1] <- c("read_id")

GM07861_HTT <- read.table("longshot_GM07861_HTT_het_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM07861_HTT)[1] <- c("read_id")

GM20239_HTT <- read.table("longshot_GM20239_HTT_het_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM20239_HTT)[1] <- c("read_id")

GM04284_FMR1 <- read.table("longshot_GM04284_FMR1_hom_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM04284_FMR1)[1] <- c("read_id")

GM05538_FMR1 <- read.table("longshot_GM05538_FMR1_hom_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM05538_FMR1)[1] <- c("read_id")

GM06891_FMR1 <- read.table("longshot_GM06891_FMR1_hom_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM06891_FMR1)[1] <- c("read_id")

GM07541_FMR1 <- read.table("longshot_GM07541_FMR1_het_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM07541_FMR1)[1] <- c("read_id")

GM07861_FMR1 <- read.table("longshot_GM07861_FMR1_hom_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM07861_FMR1)[1] <- c("read_id")

GM20239_FMR1 <- read.table("longshot_GM20239_FMR1_het_bases_window40000_DP10.txt", header = TRUE, sep = ",")
colnames(GM20239_FMR1)[1] <- c("read_id")


### Raw sequence logo containing all data (both alleles) ########################################################################
## Using ggseqlogo and character vectors

## Create tables with zeros instead of NA and deletions (-)
GM04284_HTT[is.na(GM04284_HTT)] <- 0
GM04284_HTT[GM04284_HTT == "-"] <- 0

GM05538_HTT[is.na(GM05538_HTT)] <- 0
GM05538_HTT[GM05538_HTT == "-"] <- 0

GM06891_HTT[is.na(GM06891_HTT)] <- 0
GM06891_HTT[GM06891_HTT == "-"] <- 0

GM07541_HTT[is.na(GM07541_HTT)] <- 0
GM07541_HTT[GM07541_HTT == "-"] <- 0

GM07861_HTT[is.na(GM07861_HTT)] <- 0
GM07861_HTT[GM07861_HTT == "-"] <- 0

GM20239_HTT[is.na(GM20239_HTT)] <- 0
GM20239_HTT[GM20239_HTT == "-"] <- 0

GM04284_FMR1[is.na(GM04284_FMR1)] <- 0
GM04284_FMR1[GM04284_FMR1 == "-"] <- 0

GM05538_FMR1[is.na(GM05538_FMR1)] <- 0
GM05538_FMR1[GM05538_FMR1 == "-"] <- 0

GM06891_FMR1[is.na(GM06891_FMR1)] <- 0
GM06891_FMR1[GM06891_FMR1 == "-"] <- 0

GM07541_FMR1[is.na(GM07541_FMR1)] <- 0
GM07541_FMR1[GM07541_FMR1 == "-"] <- 0

GM07861_FMR1[is.na(GM07861_FMR1)] <- 0
GM07861_FMR1[GM07861_FMR1 == "-"] <- 0

GM20239_FMR1[is.na(GM20239_FMR1)] <- 0
GM20239_FMR1[GM20239_FMR1 == "-"] <- 0


## Make a character vector for the sequences to be added to 
Character_vector_GM04284_HTT = character(nrow(GM04284_HTT))

Character_vector_GM05538_HTT = character(nrow(GM05538_HTT))

Character_vector_GM06891_HTT = character(nrow(GM06891_HTT))

Character_vector_GM07541_HTT = character(nrow(GM07541_HTT))

Character_vector_GM07861_HTT = character(nrow(GM07861_HTT))

Character_vector_GM20239_HTT = character(nrow(GM20239_HTT))

Character_vector_GM04284_FMR1 = character(nrow(GM04284_FMR1))

Character_vector_GM05538_FMR1 = character(nrow(GM05538_FMR1))

Character_vector_GM06891_FMR1 = character(nrow(GM06891_FMR1))

Character_vector_GM07541_FMR1 = character(nrow(GM07541_FMR1))

Character_vector_GM07861_FMR1 = character(nrow(GM07861_FMR1))

Character_vector_GM20239_FMR1 = character(nrow(GM20239_FMR1))


## Add the sequences to the character vector
# GM04284 HTT
i <- 1
while(i < (nrow(GM04284_HTT)+1)) {
  Character_vector_GM04284_HTT[i] <- paste(c(as.character(as.character(GM04284_HTT[i, 2:ncol(GM04284_HTT)]))),collapse = "")
  i <- i + 1
}
# GM05538 HTT
i <- 1
while(i < (nrow(GM05538_HTT)+1)) {
  Character_vector_GM05538_HTT[i] <- paste(c(as.character(as.character(GM05538_HTT[i, 2:ncol(GM05538_HTT)]))),collapse = "")
  i <- i + 1
}
# GM06891 HTT
i <- 1
while(i < (nrow(GM06891_HTT)+1)) {
  Character_vector_GM06891_HTT[i] <- paste(c(as.character(as.character(GM06891_HTT[i, 2:ncol(GM06891_HTT)]))),collapse = "")
  i <- i + 1
}
# GM07541 HTT
i <- 1
while(i < (nrow(GM07541_HTT)+1)) {
  Character_vector_GM07541_HTT[i] <- paste(c(as.character(as.character(GM07541_HTT[i, 2:ncol(GM07541_HTT)]))),collapse = "")
  i <- i + 1
}
# GM07861 HTT
i <- 1
while(i < (nrow(GM07861_HTT)+1)) {
  Character_vector_GM07861_HTT[i] <- paste(c(as.character(as.character(GM07861_HTT[i, 2:ncol(GM07861_HTT)]))),collapse = "")
  i <- i + 1
}
# GM20239 HTT
i <- 1
while(i < (nrow(GM20239_HTT)+1)) {
  Character_vector_GM20239_HTT[i] <- paste(c(as.character(as.character(GM20239_HTT[i, 2:ncol(GM20239_HTT)]))),collapse = "")
  i <- i + 1
}
# GM04284 FMR1
i <- 1
while(i < (nrow(GM04284_FMR1)+1)) {
  Character_vector_GM04284_FMR1[i] <- paste(c(as.character(as.character(GM04284_FMR1[i, 2:ncol(GM04284_FMR1)]))),collapse = "")
  i <- i + 1
}
# GM05528 FMR1
i <- 1
while(i < (nrow(GM05538_FMR1)+1)) {
  Character_vector_GM05538_FMR1[i] <- paste(c(as.character(as.character(GM05538_FMR1[i, 2:ncol(GM05538_FMR1)]))),collapse = "")
  i <- i + 1
}
# GM06891 FMR1
i <- 1
while(i < (nrow(GM06891_FMR1)+1)) {
  Character_vector_GM06891_FMR1[i] <- paste(c(as.character(as.character(GM06891_FMR1[i, 2:ncol(GM06891_FMR1)]))),collapse = "")
  i <- i + 1
}
# GM07541 FMR1
i <- 1
while(i < (nrow(GM07541_FMR1)+1)) {
  Character_vector_GM07541_FMR1[i] <- paste(c(as.character(as.character(GM07541_FMR1[i, 2:ncol(GM07541_FMR1)]))),collapse = "")
  i <- i + 1
}
# GM07861 FMR1
i <- 1
while(i < (nrow(GM07861_FMR1)+1)) {
  Character_vector_GM07861_FMR1[i] <- paste(c(as.character(as.character(GM07861_FMR1[i, 2:ncol(GM07861_FMR1)]))),collapse = "")
  i <- i + 1
}
# GM20239 FMR1
i <- 1
while(i < (nrow(GM20239_FMR1)+1)) {
  Character_vector_GM20239_FMR1[i] <- paste(c(as.character(as.character(GM20239_FMR1[i, 2:ncol(GM20239_FMR1)]))),collapse = "")
  i <- i + 1
}


## Create the sequence logo
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/logos_plots/longshot_logos_plots")
# GM04284 HTT
jpeg("GM04284_HTT_longshot.jpg")
ggseqlogo(Character_vector_GM04284_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM05538 HTT
jpeg("GM05538_HTT_longshot.jpg")
ggseqlogo(Character_vector_GM05538_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM06891 HTT
jpeg("GM06891_HTT_longshot.jpg")
ggseqlogo(Character_vector_GM06891_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM07541 HTT
jpeg("GM07541_HTT_longshot.jpg")
ggseqlogo(Character_vector_GM07541_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM07861 HTT
jpeg("GM07861_HTT_longshot.jpg")
ggseqlogo(Character_vector_GM07861_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM20239 HTT
jpeg("GM20239_HTT_longshot.jpg")
ggseqlogo(Character_vector_GM20239_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM04284 FMR1
jpeg("GM04284_FMR1_longshot.jpg")
ggseqlogo(Character_vector_GM04284_FMR1, method = 'bits', seq_type = 'dna')
dev.off()

# GM05538 FMR1
jpeg("GM05538_FMR1_longshot.jpg")
ggseqlogo(Character_vector_GM05538_FMR1, method = 'bits', seq_type = 'dna')
dev.off()

# GM06891 FMR1
jpeg("GM06891_FMR1_longshot.jpg")
ggseqlogo(Character_vector_GM06891_FMR1, method = 'bits', seq_type = 'dna')
dev.off()

# GM07541 FMR1
jpeg("GM07541_FMR1_longshot.jpg")
ggseqlogo(Character_vector_GM07541_FMR1, method = 'bits', seq_type = 'dna')
dev.off()

# GM07861 FMR1
jpeg("GM07861_FMR1_longshot.jpg")
ggseqlogo(Character_vector_GM07861_FMR1[2:10], method = 'bits', seq_type = 'dna')
dev.off()

# GM20239 FMR1
jpeg("GM20239_FMR1_longshot.jpg")
ggseqlogo(Character_vector_GM20239_FMR1, method = 'bits', seq_type = 'dna')
dev.off()


### Making the consensus sequence for the two alleles ###############################################################################
## (Subset the data frames according to haplotypes then make new logos plots)

## Load Read-Names for the alleles
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/logos_plots/longshot_readnames")
# check the text files for missing values

GM04284_HTT_readnames <- read.table("GM04284_longshot_hp_chr4_readnames.txt", header = F)
colnames(GM04284_HTT_readnames) <- c("read_id", "hp_tag")

GM05538_HTT_readnames <- read.table("GM05538_longshot_hp_chr4_readnames.txt", header = F)
colnames(GM05538_HTT_readnames) <- c("read_id", "hp_tag")

GM06891_HTT_readnames <- read.table("GM06891_longshot_hp_chr4_readnames.txt", header = F)
colnames(GM06891_HTT_readnames) <- c("read_id","hp_tag")

GM07541_HTT_readnames <- read.table("GM07541_longshot_hp_chr4_readnames.txt", header = F)
colnames(GM07541_HTT_readnames) <- c("read_id","hp_tag")

GM07861_HTT_readnames <- read.table("GM07861_longshot_hp_chr4_readnames.txt", header = F)
colnames(GM07861_HTT_readnames) <- c("read_id","hp_tag")

GM20239_HTT_readnames <- read.table("GM20239_longshot_hp_chr4_readnames.txt", header = F)
colnames(GM20239_HTT_readnames) <- c("read_id","hp_tag")

GM04284_FMR1_readnames <- read.table("GM04284_longshot_hp_chrX_readnames.txt", header = F)
colnames(GM04284_FMR1_readnames) <- c("read_id","hp_tag")

GM05538_FMR1_readnames <- read.table("GM05538_longshot_hp_chrX_readnames.txt", header = F)
colnames(GM05538_FMR1_readnames) <- c("read_id","hp_tag")

GM06891_FMR1_readnames <- read.table("GM06891_longshot_hp_chrX_readnames.txt", header = F)
colnames(GM06891_FMR1_readnames) <- c("read_id","hp_tag")

GM07541_FMR1_readnames <- read.table("GM07541_longshot_hp_chrX_readnames.txt", header = F)
colnames(GM07541_FMR1_readnames) <- c("read_id","hp_tag")

GM07861_FMR1_readnames <- read.table("GM07861_longshot_hp_chrX_readnames.txt", header = F)
colnames(GM07861_FMR1_readnames) <- c("read_id","hp_tag")

GM20239_FMR1_readnames <- read.table("GM20239_longshot_hp_chrX_readnames.txt", header = F)
colnames(GM20239_FMR1_readnames) <- c("read_id","hp_tag")


## Subset Longshot data frames into haplotypes
GM04284_HTT_longshot_hp1 <- subset(GM04284_HTT_readnames[GM04284_HTT_readnames$hp_tag == "HP:i:1",])
GM04284_HTT_longshot_hp2 <- subset(GM04284_HTT_readnames[GM04284_HTT_readnames$hp_tag == "HP:i:2",])
GM04284_FMR1_longshot_hp1 <- subset(GM04284_FMR1_readnames[GM04284_FMR1_readnames$hp_tag == "HP:i:1",])
GM04284_FMR1_longshot_hp2 <- subset(GM04284_FMR1_readnames[GM04284_FMR1_readnames$hp_tag == "HP:i:2",])

GM05538_HTT_longshot_hp1 <- subset(GM05538_HTT_readnames[GM05538_HTT_readnames$hp_tag == "HP:i:1",])
GM05538_HTT_longshot_hp2 <- subset(GM05538_HTT_readnames[GM05538_HTT_readnames$hp_tag == "HP:i:2",])
GM05538_FMR1_longshot_hp1 <- subset(GM05538_FMR1_readnames[GM05538_FMR1_readnames$hp_tag == "HP:i:1",])
GM05538_FMR1_longshot_hp2 <- subset(GM05538_FMR1_readnames[GM05538_FMR1_readnames$hp_tag == "HP:i:2",]) 

GM06891_HTT_longshot_hp1 <- subset(GM06891_HTT_readnames[GM06891_HTT_readnames$hp_tag == "HP:i:1",])
GM06891_HTT_longshot_hp2 <- subset(GM06891_HTT_readnames[GM06891_HTT_readnames$hp_tag == "HP:i:2",])
GM06891_FMR1_longshot_hp1 <- subset(GM06891_FMR1_readnames[GM06891_FMR1_readnames$hp_tag == "HP:i:1",])
GM06891_FMR1_longshot_hp2 <- subset(GM06891_FMR1_readnames[GM06891_FMR1_readnames$hp_tag == "HP:i:2",]) 

GM07541_HTT_longshot_hp1 <- subset(GM07541_HTT_readnames[GM07541_HTT_readnames$hp_tag == "HP:i:1",])
GM07541_HTT_longshot_hp2 <- subset(GM07541_HTT_readnames[GM07541_HTT_readnames$hp_tag == "HP:i:2",])
GM07541_FMR1_longshot_hp1 <- subset(GM07541_FMR1_readnames[GM07541_FMR1_readnames$hp_tag == "HP:i:1",])
GM07541_FMR1_longshot_hp2 <- subset(GM07541_FMR1_readnames[GM07541_FMR1_readnames$hp_tag == "HP:i:2",]) 

GM07861_HTT_longshot_hp1 <- subset(GM07861_HTT_readnames[GM07861_HTT_readnames$hp_tag == "HP:i:1",])
GM07861_HTT_longshot_hp2 <- subset(GM07861_HTT_readnames[GM07861_HTT_readnames$hp_tag == "HP:i:2",])
GM07861_FMR1_longshot_hp1 <- subset(GM07861_FMR1_readnames[GM07861_FMR1_readnames$hp_tag == "HP:i:1",])
GM07861_FMR1_longshot_hp2 <- subset(GM07861_FMR1_readnames[GM07861_FMR1_readnames$hp_tag == "HP:i:2",]) 

GM20239_HTT_longshot_hp1 <- subset(GM20239_HTT_readnames[GM20239_HTT_readnames$hp_tag == "HP:i:1",])
GM20239_HTT_longshot_hp2 <- subset(GM20239_HTT_readnames[GM20239_HTT_readnames$hp_tag == "HP:i:2",])
GM20239_FMR1_longshot_hp1 <- subset(GM20239_FMR1_readnames[GM20239_FMR1_readnames$hp_tag == "HP:i:1",])
GM20239_FMR1_longshot_hp2 <- subset(GM20239_FMR1_readnames[GM20239_FMR1_readnames$hp_tag == "HP:i:2",]) 



### Create subset of the data frames for allele logos plots
GM04284_HTT_allele_1 <- merge(GM04284_HTT, GM04284_HTT_longshot_hp1)
GM04284_HTT_allele_2 <- merge(GM04284_HTT, GM04284_HTT_longshot_hp2)
GM04284_FMR1_allele_1 <- merge(GM04284_FMR1, GM04284_FMR1_longshot_hp1)
GM04284_FMR1_allele_2 <- merge(GM04284_FMR1, GM04284_FMR1_longshot_hp2)

GM05538_HTT_allele_1 <- merge(GM05538_HTT, GM05538_HTT_longshot_hp1)
GM05538_HTT_allele_2 <- merge(GM05538_HTT, GM05538_HTT_longshot_hp2)
GM05538_FMR1_allele_1 <- merge(GM05538_FMR1, GM05538_FMR1_longshot_hp1)
GM05538_FMR1_allele_2 <- merge(GM05538_FMR1, GM05538_FMR1_longshot_hp2)

GM06891_HTT_allele_1 <- merge(GM06891_HTT, GM06891_HTT_longshot_hp1)
GM06891_HTT_allele_2 <- merge(GM06891_HTT, GM06891_HTT_longshot_hp2)
GM06891_FMR1_allele_1 <- merge(GM06891_FMR1, GM06891_FMR1_longshot_hp1)
GM06891_FMR1_allele_2 <- merge(GM06891_FMR1, GM06891_FMR1_longshot_hp2)

GM07541_HTT_allele_1 <- merge(GM07541_HTT, GM07541_HTT_longshot_hp1)
GM07541_HTT_allele_2 <- merge(GM07541_HTT, GM07541_HTT_longshot_hp2)
GM07541_FMR1_allele_1 <- merge(GM07541_FMR1, GM07541_FMR1_longshot_hp1)
GM07541_FMR1_allele_2 <- merge(GM07541_FMR1, GM07541_FMR1_longshot_hp2)

GM07861_HTT_allele_1 <- merge(GM07861_HTT, GM07861_HTT_longshot_hp1)
GM07861_HTT_allele_2 <- merge(GM07861_HTT, GM07861_HTT_longshot_hp2)
GM07861_FMR1_allele_1 <- merge(GM07861_FMR1, GM07861_FMR1_longshot_hp1)
GM07861_FMR1_allele_2 <- merge(GM07861_FMR1, GM07861_FMR1_longshot_hp2)

GM20239_HTT_allele_1 <- merge(GM20239_HTT, GM20239_HTT_longshot_hp1)
GM20239_HTT_allele_2 <- merge(GM20239_HTT, GM20239_HTT_longshot_hp2)
GM20239_FMR1_allele_1 <- merge(GM20239_FMR1, GM20239_FMR1_longshot_hp1)
GM20239_FMR1_allele_2 <- merge(GM20239_FMR1, GM20239_FMR1_longshot_hp2)



## Make a character vector for the sequences to be added to for the allele logos plots
Character_vector_GM04284_HTT_allele_1 = character(nrow(GM04284_HTT_allele_1))
Character_vector_GM04284_HTT_allele_2 = character(nrow(GM04284_HTT_allele_2))
Character_vector_GM04284_FMR1_allele_1 = character(nrow(GM04284_FMR1_allele_1))
Character_vector_GM04284_FMR1_allele_2 = character(nrow(GM04284_FMR1_allele_2))

Character_vector_GM05538_HTT_allele_1 = character(nrow(GM05538_HTT_allele_1))
Character_vector_GM05538_HTT_allele_2 = character(nrow(GM05538_HTT_allele_2))
Character_vector_GM05538_FMR1_allele_1 = character(nrow(GM05538_FMR1_allele_1))
Character_vector_GM05538_FMR1_allele_2 = character(nrow(GM05538_FMR1_allele_2))

Character_vector_GM06891_HTT_allele_1 = character(nrow(GM06891_HTT_allele_1))
Character_vector_GM06891_HTT_allele_2 = character(nrow(GM06891_HTT_allele_2))
Character_vector_GM06891_FMR1_allele_1 = character(nrow(GM06891_FMR1_allele_1))
Character_vector_GM06891_FMR1_allele_2 = character(nrow(GM06891_FMR1_allele_2))

Character_vector_GM07541_HTT_allele_1 = character(nrow(GM07541_HTT_allele_1))
Character_vector_GM07541_HTT_allele_2 = character(nrow(GM07541_HTT_allele_2))
Character_vector_GM07541_FMR1_allele_1 = character(nrow(GM07541_FMR1_allele_1))
Character_vector_GM07541_FMR1_allele_2 = character(nrow(GM07541_FMR1_allele_2))

Character_vector_GM07861_HTT_allele_1 = character(nrow(GM07861_HTT_allele_1))
Character_vector_GM07861_HTT_allele_2 = character(nrow(GM07861_HTT_allele_2))
Character_vector_GM07861_FMR1_allele_1 = character(nrow(GM07861_FMR1_allele_1))
Character_vector_GM07861_FMR1_allele_2 = character(nrow(GM07861_FMR1_allele_2))

Character_vector_GM20239_HTT_allele_1 = character(nrow(GM20239_HTT_allele_1))
Character_vector_GM20239_HTT_allele_2 = character(nrow(GM20239_HTT_allele_2))
Character_vector_GM20239_FMR1_allele_1 = character(nrow(GM20239_FMR1_allele_1))
Character_vector_GM20239_FMR1_allele_2 = character(nrow(GM20239_FMR1_allele_2))



## Add the sequences to the character vector
# GM04284 HTT
i <- 1
while(i < (nrow(GM04284_HTT_allele_1)+1)) {
  Character_vector_GM04284_HTT_allele_1[i] <- paste(c(as.character(as.character(GM04284_HTT_allele_1[i, 2:ncol(GM04284_HTT_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM04284_HTT_allele_2)+1)) {
  Character_vector_GM04284_HTT_allele_2[i] <- paste(c(as.character(as.character(GM04284_HTT_allele_2[i, 2:ncol(GM04284_HTT_allele_2)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM04284_FMR1_allele_1)+1)) {
  Character_vector_GM04284_FMR1_allele_1[i] <- paste(c(as.character(as.character(GM04284_FMR1_allele_1[i, 2:ncol(GM04284_FMR1_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM04284_FMR1_allele_2)+1)) {
  Character_vector_GM04284_FMR1_allele_2[i] <- paste(c(as.character(as.character(GM04284_FMR1_allele_2[i, 2:ncol(GM04284_FMR1_allele_2)]))),collapse = "")
  i <- i + 1
}



# GM05538 HTT
i <- 1
while(i < (nrow(GM05538_HTT_allele_1)+1)) {
  Character_vector_GM05538_HTT_allele_1[i] <- paste(c(as.character(as.character(GM05538_HTT_allele_1[i, 2:ncol(GM05538_HTT_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM05538_HTT_allele_2)+1)) {
  Character_vector_GM05538_HTT_allele_2[i] <- paste(c(as.character(as.character(GM05538_HTT_allele_2[i, 2:ncol(GM05538_HTT_allele_2)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM05538_FMR1_allele_1)+1)) {
  Character_vector_GM05538_FMR1_allele_1[i] <- paste(c(as.character(as.character(GM05538_FMR1_allele_1[i, 2:ncol(GM05538_FMR1_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM05538_FMR1_allele_2)+1)) {
  Character_vector_GM05538_FMR1_allele_2[i] <- paste(c(as.character(as.character(GM05538_FMR1_allele_2[i, 2:ncol(GM05538_FMR1_allele_2)]))),collapse = "")
  i <- i + 1
}



# GM06891 HTT
i <- 1
while(i < (nrow(GM06891_HTT_allele_1)+1)) {
  Character_vector_GM06891_HTT_allele_1[i] <- paste(c(as.character(as.character(GM06891_HTT_allele_1[i, 2:ncol(GM06891_HTT_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM06891_HTT_allele_2)+1)) {
  Character_vector_GM06891_HTT_allele_2[i] <- paste(c(as.character(as.character(GM06891_HTT_allele_2[i, 2:ncol(GM06891_HTT_allele_2)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM06891_FMR1_allele_1)+1)) {
  Character_vector_GM06891_FMR1_allele_1[i] <- paste(c(as.character(as.character(GM06891_FMR1_allele_1[i, 2:ncol(GM06891_FMR1_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM06891_FMR1_allele_2)+1)) {
  Character_vector_GM06891_FMR1_allele_2[i] <- paste(c(as.character(as.character(GM06891_FMR1_allele_2[i, 2:ncol(GM06891_FMR1_allele_2)]))),collapse = "")
  i <- i + 1
}



# GM07541 HTT
i <- 1
while(i < (nrow(GM07541_HTT_allele_1)+1)) {
  Character_vector_GM07541_HTT_allele_1[i] <- paste(c(as.character(as.character(GM07541_HTT_allele_1[i, 2:ncol(GM07541_HTT_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM07541_HTT_allele_2)+1)) {
  Character_vector_GM07541_HTT_allele_2[i] <- paste(c(as.character(as.character(GM07541_HTT_allele_2[i, 2:ncol(GM07541_HTT_allele_2)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM07541_FMR1_allele_1)+1)) {
  Character_vector_GM07541_FMR1_allele_1[i] <- paste(c(as.character(as.character(GM07541_FMR1_allele_1[i, 2:ncol(GM07541_FMR1_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM07541_FMR1_allele_2)+1)) {
  Character_vector_GM07541_FMR1_allele_2[i] <- paste(c(as.character(as.character(GM07541_FMR1_allele_2[i, 2:ncol(GM07541_FMR1_allele_2)]))),collapse = "")
  i <- i + 1
}



# GM07861 HTT
i <- 1
while(i < (nrow(GM07861_HTT_allele_1)+1)) {
  Character_vector_GM07861_HTT_allele_1[i] <- paste(c(as.character(as.character(GM07861_HTT_allele_1[i, 2:ncol(GM07861_HTT_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM07861_HTT_allele_2)+1)) {
  Character_vector_GM07861_HTT_allele_2[i] <- paste(c(as.character(as.character(GM07861_HTT_allele_2[i, 2:ncol(GM07861_HTT_allele_2)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM07861_FMR1_allele_1)+1)) {
  Character_vector_GM07861_FMR1_allele_1[i] <- paste(c(as.character(as.character(GM07861_FMR1_allele_1[i, 2:ncol(GM07861_FMR1_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM07861_FMR1_allele_2)+1)) {
  Character_vector_GM07861_FMR1_allele_2[i] <- paste(c(as.character(as.character(GM07861_FMR1_allele_2[i, 2:ncol(GM07861_FMR1_allele_2)]))),collapse = "")
  i <- i + 1
}


# GM20239 HTT
i <- 1
while(i < (nrow(GM20239_HTT_allele_1)+1)) {
  Character_vector_GM20239_HTT_allele_1[i] <- paste(c(as.character(as.character(GM20239_HTT_allele_1[i, 2:ncol(GM20239_HTT_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM20239_HTT_allele_2)+1)) {
  Character_vector_GM20239_HTT_allele_2[i] <- paste(c(as.character(as.character(GM20239_HTT_allele_2[i, 2:ncol(GM20239_HTT_allele_2)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM20239_FMR1_allele_1)+1)) {
  Character_vector_GM20239_FMR1_allele_1[i] <- paste(c(as.character(as.character(GM20239_FMR1_allele_1[i, 2:ncol(GM20239_FMR1_allele_1)]))),collapse = "")
  i <- i + 1
}
i <- 1
while(i < (nrow(GM20239_FMR1_allele_2)+1)) {
  Character_vector_GM20239_FMR1_allele_2[i] <- paste(c(as.character(as.character(GM20239_FMR1_allele_2[i, 2:ncol(GM20239_FMR1_allele_2)]))),collapse = "")
  i <- i + 1
}








## Create the sequence logo
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/logos_plots/longshot_logos_plots")

# GM04284 HTT
jpeg("GM04284_HTT_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM04284_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM04284_HTT_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM04284_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM04284_FMR1_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM04284_FMR1_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM04284_FMR1_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM04284_FMR1_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

# GM05538 HTT
jpeg("GM05538_HTT_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM05538_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM05538_HTT_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM05538_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM05538_FMR1_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM05538_FMR1_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM05538_FMR1_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM05538_FMR1_allele_2, method = 'bits', seq_type = 'dna')
dev.off()


# GM06891 HTT
jpeg("GM06891_HTT_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM06891_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM06891_HTT_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM06891_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM06891_FMR1_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM06891_FMR1_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM06891_FMR1_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM06891_FMR1_allele_2, method = 'bits', seq_type = 'dna')
dev.off()


# GM07541 HTT
jpeg("GM07541_HTT_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM07541_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07541_HTT_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM07541_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07541_FMR1_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM07541_FMR1_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07541_FMR1_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM07541_FMR1_allele_2, method = 'bits', seq_type = 'dna')
dev.off()


# GM07861 HTT
jpeg("GM07861_HTT_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM07861_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07861_HTT_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM07861_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07861_FMR1_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM07861_FMR1_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07861_FMR1_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM07861_FMR1_allele_2, method = 'bits', seq_type = 'dna')
dev.off()


# GM20239 HTT
jpeg("GM20239_HTT_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM20239_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM20239_HTT_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM20239_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM20239_FMR1_allele_1_longshot.jpg")
ggseqlogo(Character_vector_GM20239_FMR1_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM20239_FMR1_allele_2_longshot.jpg")
ggseqlogo(Character_vector_GM20239_FMR1_allele_2, method = 'bits', seq_type = 'dna')
dev.off()






