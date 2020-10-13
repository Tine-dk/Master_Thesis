### This script creates logos plots for Illumina-Guided HP SNV DFs ###

#################################################################################################################################
## Set Working Directory ########################################################################################################
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/SNV_analyse/illumina_guided_hp/dataframes_and_genotypes")

## Load Libraries ###############################################################################################################
require(ggplot2)
require(ggseqlogo)

### Loading the table ###########################################################################################################
GM04284_HTT <- NA
GM05538_HTT <- read.table("G87-GM05538_HTT_bases_window40000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",")
colnames(GM05538_HTT)[1] <- c("read_id")

GM06891_HTT <- read.table("G87-GM06891_HTT_bases_window40000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",")
colnames(GM06891_HTT)[1] <- c("read_id")

GM07541_HTT <- read.table("G87-GM07541_HTT_bases_window40000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",")
colnames(GM07541_HTT)[1] <- c("read_id")

GM07861_HTT <- read.table("G87-GM07861_HTT_bases_window40000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",")
colnames(GM07861_HTT)[1] <- c("read_id")

GM20239_HTT <- read.table("G87-GM20239_HTT_bases_window40000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",")
colnames(GM20239_HTT)[1] <- c("read_id")

GM04284_FMR1 <- read.table("G87-GM04284_HTT_bases_window40000_AF_0_90_to_1_00_DP10.txt", header = TRUE, sep = ",")
colnames(GM04284_FMR1)[1] <- c("read_id")

GM05538_FMR1 <- read.table("G87-GM05538_HTT_bases_window40000_AF_0_90_to_1_00_DP10.txt", header = TRUE, sep = ",")
colnames(GM05538_FMR1)[1] <- c("read_id")

GM06891_FMR1 <- read.table("G87-GM06891_HTT_bases_window40000_AF_0_90_to_1_00_DP10.txt", header = TRUE, sep = ",")
colnames(GM06891_FMR1)[1] <- c("read_id")

GM07541_FMR1 <- read.table("G87-GM07541_FMR1_bases_window40000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",") #Female = heterozygous
colnames(GM07541_FMR1)[1] <- c("read_id")

GM07861_FMR1 <- read.table("G87-GM07861_FMR1_bases_window40000_AF_0_90_to_1_00_DP10.txt", header = TRUE, sep = ",") 
colnames(GM07861_FMR1)[1] <- c("read_id")

GM20239_FMR1 <- read.table("G87-GM20239_HTT_bases_window40000_AF_0_30_to_0_70_DP10.txt", header = TRUE, sep = ",") #Female = heterozygous
colnames(GM20239_FMR1)[1] <- c("read_id")

### Raw sequence logo containing all data (both alleles) ########################################################################
### Using ggseqlogo and character vectors

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
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/logos_plots/illumina_guided_logos_plots")
# GM04284 HTT
jpeg("GM04284_HTT_illumina_guided.jpg")
ggseqlogo(Character_vector_GM04284_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM05538 HTT
jpeg("GM05538_HTT_illumina_guided.jpg")
ggseqlogo(Character_vector_GM05538_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM06891 HTT
jpeg("GM06891_HTT_illumina_guided.jpg")
ggseqlogo(Character_vector_GM06891_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM07541 HTT
jpeg("GM07541_HTT_illumina_guided.jpg")
ggseqlogo(Character_vector_GM07541_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM07861 HTT
jpeg("GM07861_HTT_illumina_guided.jpg")
ggseqlogo(Character_vector_GM07861_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM20239 HTT
jpeg("GM20239_HTT_illumina_guided.jpg")
ggseqlogo(Character_vector_GM20239_HTT, method = 'bits', seq_type = 'dna')
dev.off()

# GM04284 FMR1
jpeg("GM04284_FMR1_illumina_guided.jpg")
ggseqlogo(Character_vector_GM04284_FMR1, method = 'bits', seq_type = 'dna')
dev.off()

# GM05538 FMR1
jpeg("GM05538_FMR1_illumina_guided.jpg")
ggseqlogo(Character_vector_GM05538_FMR1, method = 'bits', seq_type = 'dna')
dev.off()

# GM06891 FMR1
jpeg("GM06891_FMR1_illumina_guided.jpg")
ggseqlogo(Character_vector_GM06891_FMR1, method = 'bits', seq_type = 'dna')
dev.off()

# GM07541 FMR1
jpeg("GM07541_FMR1_illumina_guided.jpg")
ggseqlogo(Character_vector_GM07541_FMR1, method = 'bits', seq_type = 'dna')
dev.off()

# GM07861 FMR1
jpeg("GM07861_FMR1_illumina_guided.jpg")
ggseqlogo(Character_vector_GM07861_FMR1[2:10], method = 'bits', seq_type = 'dna')
dev.off()

# GM20239 FMR1
jpeg("GM20239_FMR1_illumina_guided.jpg")
ggseqlogo(Character_vector_GM20239_FMR1, method = 'bits', seq_type = 'dna')
dev.off()



### Making the consensus sequence for the two alleles ###############################################################################
## Subset the data frames according to haplotypes then make new logos plots
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/SNV_analyse/illumina_guided_hp/dataframes_and_genotypes/Output_haplotype_phasing_heterozygote")

## Load Read-Names for the alleles
GM04284_HTT_readnames_hap1 <- NA
GM04284_HTT_readnames_hap2 <- NA

GM05538_HTT_readnames_hap1_l <- read.table("G87-GM05538_HTT_reads_hap1_left.txt", header = T,  stringsAsFactors = F)
GM05538_HTT_readnames_hap1_rec_l <- read.table("G87-GM05538_HTT_reads_hap1_left_recovered.txt", header = T, sep = "")
GM05538_HTT_readnames_hap2_l <- read.table("G87-GM05538_HTT_reads_hap2_left.txt", header = T, sep = "")
GM05538_HTT_readnames_hap2_rec_l <- read.table("G87-GM05538_HTT_reads_hap2_left_recovered.txt", header = T, sep = "")
GM05538_HTT_readnames_hap1_r <- read.table("G87-GM05538_HTT_reads_hap1_left.txt", header = T,  stringsAsFactors = F)
GM05538_HTT_readnames_hap1_rec_r <- read.table("G87-GM05538_HTT_reads_hap1_left_recovered.txt", header = T, sep = "")
GM05538_HTT_readnames_hap2_r <- read.table("G87-GM05538_HTT_reads_hap2_left.txt", header = T, sep = "")
GM05538_HTT_readnames_hap2_rec_r <- read.table("G87-GM05538_HTT_reads_hap2_left_recovered.txt", header = T, sep = "")

GM06891_HTT_readnames_hap1 <- read.table("G87-GM06891_HTT_reads_hap1.txt", header = T, sep = "")
GM06891_HTT_readnames_hap2 <- read.table("G87-GM06891_HTT_reads_hap2.txt", header = T, sep = "")

GM07541_HTT_readnames_hap1 <- read.table("G87-GM07541_HTT_reads_hap1.txt", header = T, sep = "")
GM07541_HTT_readnames_hap1_rec <- read.table("G87-GM07541_HTT_reads_hap1_recovered.txt", header = T, sep = "")
GM07541_HTT_readnames_hap2 <- read.table("G87-GM07541_HTT_reads_hap2.txt", header = T, sep = "")

GM07861_HTT_readnames_hap1 <- read.table("G87-GM07861_HTT_reads_hap1.txt", header = T, sep = "")
GM07861_HTT_readnames_hap1_rec <- read.table("G87-GM07861_HTT_reads_hap1_recovered.txt", header = T, sep = "")
GM07861_HTT_readnames_hap2 <- read.table("G87-GM07861_HTT_reads_hap2.txt", header = T, sep = "")
GM07861_HTT_readnames_hap2_rec <- read.table("G87-GM07861_HTT_reads_hap2_recovered.txt", header = T, sep = "")

GM20239_HTT_readnames_hap1 <- read.table("G87-GM20239_HTT_reads_hap1.txt", header = T, sep = "")
GM20239_HTT_readnames_hap1_rec <- read.table("G87-GM20239_HTT_reads_hap1_recovered.txt", header = T, sep = "")
GM20239_HTT_readnames_hap2 <- read.table("G87-GM20239_HTT_reads_hap2.txt", header = T, sep = "")

GM04284_FMR1_readnames_hap1 <- NA
GM04284_FMR1_readnames_hap2 <- NA

GM05538_FMR1_readnames_hap1 <- NA
GM05538_FMR1_readnames_hap2 <- NA

GM06891_FMR1_readnames_hap1 <- NA
GM06891_FMR1_readnames_hap2 <- NA

GM07541_FMR1_readnames_hap1 <- read.table("G87-GM07541_FMR1_reads_hap1.txt", header = T, sep = "")
GM07541_FMR1_readnames_hap1_rec <- read.table("G87-GM07541_FMR1_reads_hap1_recovered.txt", header = T, sep = "")
GM07541_FMR1_readnames_hap2 <- read.table("G87-GM07541_FMR1_reads_hap2.txt", header = T, sep = "")

GM07861_FMR1_readnames_hap1 <- NA
GM07861_FMR1_readnames_hap2 <- NA

GM20239_FMR1_readnames_hap1 <- read.table("G87-GM20239_FMR1_reads_hap1.txt", header = T, sep = "")
GM20239_FMR1_readnames_hap1_rec <- read.table("G87-GM20239_FMR1_reads_hap1_recovered.txt", header = T, sep = "")
GM20239_FMR1_readnames_hap2 <- read.table("G87-GM20239_FMR1_reads_hap1.txt", header = T, sep = "")
GM20239_FMR1_readnames_hap2_rec <- read.table("G87-GM20239_FMR1_reads_hap1_recovered.txt", header = T, sep = "")


# Merge the recovered data with theit belonging haplotype
# (check the dimensions of the tables, remove the last column from the recovered as it is a tag)
GM05538_HTT_readnames_hap1_l <- rbind(GM05538_HTT_readnames_hap1_l, GM05538_HTT_readnames_hap1_rec_l[,1:dim(GM05538_HTT_readnames_hap1_l)[2]])
GM05538_HTT_readnames_hap2_l <- rbind(GM05538_HTT_readnames_hap2_l, GM05538_HTT_readnames_hap2_rec_l[,1:dim(GM05538_HTT_readnames_hap2_l)[2]])
GM05538_HTT_readnames_hap1_r <- rbind(GM05538_HTT_readnames_hap1_r, GM05538_HTT_readnames_hap1_rec_r[,1:dim(GM05538_HTT_readnames_hap1_r)[2]])
GM05538_HTT_readnames_hap2_r <- rbind(GM05538_HTT_readnames_hap2_r, GM05538_HTT_readnames_hap2_rec_r[,1:dim(GM05538_HTT_readnames_hap2_r)[2]])
GM05538_HTT_readnames_hap1 <- rbind(GM05538_HTT_readnames_hap1_r, GM05538_HTT_readnames_hap1_l)
GM05538_HTT_readnames_hap2 <- rbind(GM05538_HTT_readnames_hap2_r, GM05538_HTT_readnames_hap2_l)


GM07541_HTT_readnames_hap1 <- rbind(GM07541_HTT_readnames_hap1, GM07541_HTT_readnames_hap1_rec[, 1:dim(GM07541_HTT_readnames_hap1)[2]])

GM07861_HTT_readnames_hap1 <- rbind(GM07861_HTT_readnames_hap1, GM07861_HTT_readnames_hap1_rec[, 1:dim(GM07861_HTT_readnames_hap1)[2]])
GM07861_HTT_readnames_hap2 <- rbind(GM07861_HTT_readnames_hap2, GM07861_HTT_readnames_hap2_rec[, 1:dim(GM07861_HTT_readnames_hap2)[2]])

GM20239_HTT_readnames_hap1 <- rbind(GM20239_HTT_readnames_hap1, GM20239_HTT_readnames_hap1_rec[, 1:dim(GM20239_HTT_readnames_hap1)[2]])

GM07541_FMR1_readnames_hap1 <- rbind(GM07541_FMR1_readnames_hap1, GM07541_FMR1_readnames_hap1_rec[, 1:dim(GM07541_FMR1_readnames_hap1)[2]])

GM20239_FMR1_readnames_hap1 <- rbind(GM20239_FMR1_readnames_hap1, GM20239_FMR1_readnames_hap1_rec[, 1:dim(GM20239_FMR1_readnames_hap1)[2]])
GM20239_FMR1_readnames_hap2 <- rbind(GM20239_FMR1_readnames_hap2, GM20239_FMR1_readnames_hap2_rec[, 1:dim(GM20239_FMR1_readnames_hap2)[2]])


## Save the row numbers for the alleles as a vector
# (which row number is read-name from row 1 from hap1? etc. )

GM05538_HTT_allele_1_rows <- match(GM05538_HTT_readnames_hap1$read_id, GM05538_HTT$read_id, nomatch = NA)
GM05538_HTT_allele_2_rows <- match(GM05538_HTT_readnames_hap2$read_id, GM05538_HTT$read_id, nomatch = NA)

GM06891_HTT_allele_1_rows <- match(GM06891_HTT_readnames_hap1$read_id, GM06891_HTT$read_id, nomatch = NA)
GM06891_HTT_allele_2_rows <- match(GM06891_HTT_readnames_hap2$read_id, GM06891_HTT$read_id, nomatch = NA)

GM07541_HTT_allele_1_rows <- match(GM07541_HTT_readnames_hap1$read_id, GM07541_HTT$read_id, nomatch = NA)
GM07541_HTT_allele_2_rows <- match(GM07541_HTT_readnames_hap2$read_id, GM07541_HTT$read_id, nomatch = NA)

GM07541_FMR1_allele_1_rows <- match(GM07541_FMR1_readnames_hap1$read_id, GM07541_FMR1$read_id, nomatch = NA)
GM07541_FMR1_allele_2_rows <- match(GM07541_FMR1_readnames_hap2$read_id, GM07541_FMR1$read_id, nomatch = NA)

GM07861_HTT_allele_1_rows <- match(GM07861_HTT_readnames_hap1$read_id, GM07861_HTT$read_id, nomatch = NA)
GM07861_HTT_allele_2_rows <- match(GM07861_HTT_readnames_hap2$read_id, GM07861_HTT$read_id, nomatch = NA)

GM20239_HTT_allele_1_rows <- match(GM20239_HTT_readnames_hap1$read_id, GM20239_HTT$read_id, nomatch = NA)
GM20239_HTT_allele_2_rows <- match(GM20239_HTT_readnames_hap2$read_id, GM20239_HTT$read_id, nomatch = NA)

GM20239_FMR1_allele_1_rows <- match(GM20239_FMR1_readnames_hap1$read_id, GM20239_FMR1$read_id, nomatch = NA)
GM20239_FMR1_allele_2_rows <- match(GM20239_FMR1_readnames_hap2$read_id, GM20239_FMR1$read_id, nomatch = NA)


## Create allele character vectors for the samples
Character_vector_GM05538_HTT_allele_1 <- Character_vector_GM05538_HTT[GM05538_HTT_allele_1_rows]
Character_vector_GM05538_HTT_allele_2 <- Character_vector_GM05538_HTT[GM05538_HTT_allele_2_rows]

Character_vector_GM06891_HTT_allele_1 <- Character_vector_GM06891_HTT[GM06891_HTT_allele_1_rows]
Character_vector_GM06891_HTT_allele_2 <- Character_vector_GM06891_HTT[GM06891_HTT_allele_2_rows]

Character_vector_GM07541_HTT_allele_1 <- Character_vector_GM07541_HTT[GM07541_HTT_allele_1_rows]
Character_vector_GM07541_HTT_allele_2 <- Character_vector_GM07541_HTT[GM07541_HTT_allele_2_rows]

Character_vector_GM07541_FMR1_allele_1 <- Character_vector_GM07541_FMR1[GM07541_FMR1_allele_1_rows]
Character_vector_GM07541_FMR1_allele_2 <- Character_vector_GM07541_FMR1[GM07541_FMR1_allele_2_rows]

Character_vector_GM07861_HTT_allele_1 <- Character_vector_GM07861_HTT[GM07861_HTT_allele_1_rows]
Character_vector_GM07861_HTT_allele_2 <- Character_vector_GM07861_HTT[GM07861_HTT_allele_2_rows]

Character_vector_GM20239_HTT_allele_1 <- Character_vector_GM20239_HTT[GM20239_HTT_allele_1_rows]
Character_vector_GM20239_HTT_allele_2 <- Character_vector_GM20239_HTT[GM20239_HTT_allele_2_rows]

Character_vector_GM20239_FMR1_allele_1 <- Character_vector_GM20239_FMR1[GM20239_FMR1_allele_1_rows]
Character_vector_GM20239_FMR1_allele_2 <- Character_vector_GM20239_FMR1[GM20239_FMR1_allele_2_rows]


## Plot the logos plots
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/logos_plots/illumina_guided_logos_plots")
# GM05538 HTT 
jpeg("GM05538_HTT_allele_1.jpg")
ggseqlogo(Character_vector_GM05538_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM05538_HTT_allele_2.jpg")
ggseqlogo(Character_vector_GM05538_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

# GM06891 HTT 
jpeg("GM06891_HTT_allele_1.jpg")
ggseqlogo(Character_vector_GM06891_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM06891_HTT_allele_2.jpg")
ggseqlogo(Character_vector_GM06891_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

# GM07541 HTT 
jpeg("GM07541_HTT_allele_1.jpg")
ggseqlogo(Character_vector_GM07541_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07541_HTT_allele_2.jpg")
ggseqlogo(Character_vector_GM07541_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

# GM07541 FMR1 
jpeg("GM07541_FMR1_allele_1.jpg")
ggseqlogo(Character_vector_GM07541_FMR1_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07541_FMR1_allele_2.jpg")
ggseqlogo(Character_vector_GM07541_FMR1_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

# GM07861 HTT 
jpeg("GM07861_HTT_allele_1.jpg")
ggseqlogo(Character_vector_GM07861_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM07861_HTT_allele_2.jpg")
ggseqlogo(Character_vector_GM07861_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

# GM20239 HTT 
jpeg("GM20239_HTT_allele_1.jpg")
ggseqlogo(Character_vector_GM20239_HTT_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM20239_HTT_allele_2.jpg")
ggseqlogo(Character_vector_GM20239_HTT_allele_2, method = 'bits', seq_type = 'dna')
dev.off()

# GM20239 FMR1
jpeg("GM20239_FMR1_allele_1.jpg")
ggseqlogo(Character_vector_GM20239_FMR1_allele_1, method = 'bits', seq_type = 'dna')
dev.off()

jpeg("GM20239_FMR1_allele_2.jpg")
ggseqlogo(Character_vector_GM20239_FMR1_allele_2, method = 'bits', seq_type = 'dna')
dev.off()





