### Illumina-Guided HP and Longshot Comparison ####


# 1. Træk read-navne ud af Longshot BAM og skriv metode om det (ikke en del af dette skript)
# 2. Lav match imellem Longshot og Illumina HP 'read-names lister'.
# 3. Lav Venn Diagrammer




##################################


## Load tables with Illumina-Guided Read-names
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/SNV_analyse/illumina_guided_results_dataframes/Output_haplotype_phasing_heterozygote")

GM05538_HTT_Illumina_guided_HP1_left <- read.table("G87-GM05538_HTT_reads_hap1_left.txt", header = T)
GM05538_HTT_Illumina_guided_HP1_left_recovered <- read.table("G87-GM05538_HTT_reads_hap1_left_recovered.txt", header = T)
GM05538_HTT_Illumina_guided_HP2_left <- read.table("G87-GM05538_HTT_reads_hap2_left.txt", header = T)
GM05538_HTT_Illumina_guided_HP2_left_recovered <- read.table("G87-GM05538_HTT_reads_hap2_left_recovered.txt", header = T)

GM05538_HTT_Illumina_guided_HP1_right <- read.table("G87-GM05538_HTT_reads_hap1_right.txt", header = T)
GM05538_HTT_Illumina_guided_HP1_right_recovered <- read.table("G87-GM05538_HTT_reads_hap1_right_recovered.txt", header = T)
GM05538_HTT_Illumina_guided_HP2_right <- read.table("G87-GM05538_HTT_reads_hap2_right.txt", header = T)

GM06891_HTT_Illumina_guided_HP1 <- read.table("G87-GM06891_HTT_reads_hap1.txt", header = T)
GM06891_HTT_Illumina_guided_HP2 <- read.table("G87-GM06891_HTT_reads_hap2.txt", header = T)

GM07541_HTT_Illumina_guided_HP1 <- read.table("G87-GM07541_HTT_reads_hap1.txt", header = T)
GM07541_HTT_Illumina_guided_HP1_recovered <- read.table("G87-GM07541_HTT_reads_hap1_recovered.txt", header = T)
GM07541_HTT_Illumina_guided_HP2 <- read.table("G87-GM07541_HTT_reads_hap2.txt", header = T)

GM07541_FMR1_Illumina_guided_HP1 <- read.table("G87-GM07541_FMR1_reads_hap1.txt", header = T)
GM07541_FMR1_Illumina_guided_HP1_recovered <- read.table("G87-GM07541_FMR1_reads_hap1_recovered.txt", header = T)
GM07541_FMR1_Illumina_guided_HP2 <- read.table("G87-GM07541_FMR1_reads_hap2.txt", header = T)

GM07861_HTT_Illumina_guided_HP1 <- read.table("G87-GM07861_HTT_reads_hap1.txt", header = T)
GM07861_HTT_Illumina_guided_HP1_recovered <- read.table("G87-GM07861_HTT_reads_hap1_recovered.txt", header = T)
GM07861_HTT_Illumina_guided_HP2 <- read.table("G87-GM07861_HTT_reads_hap2.txt", header = T)
GM07861_HTT_Illumina_guided_HP2_recovered <- read.table("G87-GM07861_HTT_reads_hap2_recovered.txt", header = T)

GM20239_FMR1_Illumina_guided_HP1 <- read.table("G87-GM20239_FMR1_reads_hap1.txt", header = T)
GM20239_FMR1_Illumina_guided_HP1_recovered <- read.table("G87-GM20239_FMR1_reads_hap1_recovered.txt", header = T)
GM20239_FMR1_Illumina_guided_HP2 <- read.table("G87-GM20239_FMR1_reads_hap2.txt", header = T)
GM20239_FMR1_Illumina_guided_HP2_recovered <- read.table("G87-GM20239_FMR1_reads_hap2_recovered.txt", header = T)

GM20239_HTT_Illumina_guided_HP1 <- read.table("G87-GM20239_HTT_reads_hap1.txt", header = T)
GM20239_HTT_Illumina_guided_HP1_recovered <- read.table("G87-GM20239_HTT_reads_hap1_recovered.txt", header = T)
GM20239_HTT_Illumina_guided_HP2 <- read.table("G87-GM20239_HTT_reads_hap2.txt", header = T)



## Merge recovered Illumina-guided haplotyped tables with the "original table" if neccessary and remove unneccesary columns
# GM05538 HTT HP1 & HP2 left concatenating
GM05538_HTT_Illumina_guided_HP1_left <- as.data.frame(GM05538_HTT_Illumina_guided_HP1_left[,1])
GM05538_HTT_Illumina_guided_HP1_left_recovered <- as.data.frame(GM05538_HTT_Illumina_guided_HP1_left_recovered[,1])
colnames(GM05538_HTT_Illumina_guided_HP1_left) <- c("read_id")
colnames(GM05538_HTT_Illumina_guided_HP1_left_recovered) <- c("read_id")
GM05538_HTT_Illumina_guided_HP1_left_with_recovered <- rbind(GM05538_HTT_Illumina_guided_HP1_left, GM05538_HTT_Illumina_guided_HP1_left_recovered)

GM05538_HTT_Illumina_guided_HP2_left <- as.data.frame(GM05538_HTT_Illumina_guided_HP2_left[,1])
GM05538_HTT_Illumina_guided_HP2_left_recovered <- as.data.frame(GM05538_HTT_Illumina_guided_HP2_left_recovered[,1])
colnames(GM05538_HTT_Illumina_guided_HP2_left) <- c("read_id")
colnames(GM05538_HTT_Illumina_guided_HP2_left_recovered) <- c("read_id")
GM05538_HTT_Illumina_guided_HP2_left_with_recovered <- rbind(GM05538_HTT_Illumina_guided_HP2_left, GM05538_HTT_Illumina_guided_HP2_left_recovered)

# GM05538 HTT HP1 right concatenating and HP2 right cleanup
GM05538_HTT_Illumina_guided_HP1_right <- as.data.frame(GM05538_HTT_Illumina_guided_HP1_right[,1])
GM05538_HTT_Illumina_guided_HP1_right_recovered <- as.data.frame(GM05538_HTT_Illumina_guided_HP1_right_recovered[,1])
colnames(GM05538_HTT_Illumina_guided_HP1_right) <- c("read_id")
colnames(GM05538_HTT_Illumina_guided_HP1_right_recovered) <- c("read_id")
GM05538_HTT_Illumina_guided_HP1_right_with_recovered <- rbind(GM05538_HTT_Illumina_guided_HP1_right, GM05538_HTT_Illumina_guided_HP1_right_recovered)

GM05538_HTT_Illumina_guided_HP2_right <- as.data.frame(GM05538_HTT_Illumina_guided_HP2_right[,1])
colnames(GM05538_HTT_Illumina_guided_HP2_right) <- c("read_id")

# GM06891 HTT HP1 & HP2 cleanup
GM06891_HTT_Illumina_guided_HP1 <- as.data.frame(GM06891_HTT_Illumina_guided_HP1[,1])
GM06891_HTT_Illumina_guided_HP2 <- as.data.frame(GM06891_HTT_Illumina_guided_HP2[,1])
colnames(GM06891_HTT_Illumina_guided_HP1) <- c("read_id")
colnames(GM06891_HTT_Illumina_guided_HP2) <- c("read_id")

# GM07541 HTT HP1 concatenating and HP2 cleanup
GM07541_HTT_Illumina_guided_HP1 <- as.data.frame(GM07541_HTT_Illumina_guided_HP1[,1])
GM07541_HTT_Illumina_guided_HP1_recovered <- as.data.frame(GM07541_HTT_Illumina_guided_HP1_recovered[,1])
colnames(GM07541_HTT_Illumina_guided_HP1) <- c("read_id")
colnames(GM07541_HTT_Illumina_guided_HP1_recovered) <- c("read_id")
GM07541_HTT_Illumina_guided_HP1_with_recovered <- rbind(GM07541_HTT_Illumina_guided_HP1, GM07541_HTT_Illumina_guided_HP1_recovered)

GM07541_HTT_Illumina_guided_HP2 <- as.data.frame(GM07541_HTT_Illumina_guided_HP2[,1])
colnames(GM07541_HTT_Illumina_guided_HP2) <- c("read_id")

# GM07541 FMR1 HP1 concatenating and HP2 cleanup
GM07541_FMR1_Illumina_guided_HP1 <- as.data.frame(GM07541_FMR1_Illumina_guided_HP1[,1])
GM07541_FMR1_Illumina_guided_HP1_recovered <- as.data.frame(GM07541_FMR1_Illumina_guided_HP1_recovered[,1])
colnames(GM07541_FMR1_Illumina_guided_HP1) <- c("read_id")
colnames(GM07541_FMR1_Illumina_guided_HP1_recovered) <- c("read_id")
GM07541_FMR1_Illumina_guided_HP1_with_recovered <- rbind(GM07541_FMR1_Illumina_guided_HP1, GM07541_FMR1_Illumina_guided_HP1_recovered)

GM07541_FMR1_Illumina_guided_HP2 <- as.data.frame(GM07541_FMR1_Illumina_guided_HP2[,1])
colnames(GM07541_FMR1_Illumina_guided_HP2) <- c("read_id")

# GM07861 HTT HP1 & HP2 concatenating
GM07861_HTT_Illumina_guided_HP1 <- as.data.frame(GM07861_HTT_Illumina_guided_HP1[,1])
GM07861_HTT_Illumina_guided_HP1_recovered <- as.data.frame(GM07861_HTT_Illumina_guided_HP1_recovered[,1])
colnames(GM07861_HTT_Illumina_guided_HP1) <- c("read_id")
colnames(GM07861_HTT_Illumina_guided_HP1_recovered) <- c("read_id")
GM07861_HTT_Illumina_guided_HP1_with_recovered <- rbind(GM07861_HTT_Illumina_guided_HP1, GM07861_HTT_Illumina_guided_HP1_recovered)

GM07861_HTT_Illumina_guided_HP2 <- as.data.frame(GM07861_HTT_Illumina_guided_HP2[,1])
GM07861_HTT_Illumina_guided_HP2_recovered <- as.data.frame(GM07861_HTT_Illumina_guided_HP2_recovered[,1])
colnames(GM07861_HTT_Illumina_guided_HP2) <- c("read_id")
colnames(GM07861_HTT_Illumina_guided_HP2_recovered) <- c("read_id")
GM07861_HTT_Illumina_guided_HP2_with_recovered <- rbind(GM07861_HTT_Illumina_guided_HP2, GM07861_HTT_Illumina_guided_HP2_recovered)

# GM20239 FMR1 HP1 & HP2 concatenating and
GM20239_FMR1_Illumina_guided_HP1 <- as.data.frame(GM20239_FMR1_Illumina_guided_HP1[,1])
GM20239_FMR1_Illumina_guided_HP1_recovered <- as.data.frame(GM20239_FMR1_Illumina_guided_HP1_recovered[,1])
colnames(GM20239_FMR1_Illumina_guided_HP1) <- c("read_id")
colnames(GM20239_FMR1_Illumina_guided_HP1_recovered) <- c("read_id")
GM20239_FMR1_Illumina_guided_HP1_with_recovered <- rbind(GM20239_FMR1_Illumina_guided_HP1, GM20239_FMR1_Illumina_guided_HP1_recovered)

GM20239_FMR1_Illumina_guided_HP2 <- as.data.frame(GM20239_FMR1_Illumina_guided_HP2[,1])
GM20239_FMR1_Illumina_guided_HP2_recovered <- as.data.frame(GM20239_FMR1_Illumina_guided_HP2_recovered[,1])
colnames(GM20239_FMR1_Illumina_guided_HP2) <- c("read_id")
colnames(GM20239_FMR1_Illumina_guided_HP2_recovered) <- c("read_id")
GM20239_FMR1_Illumina_guided_HP2_with_recovered <- rbind(GM20239_FMR1_Illumina_guided_HP2, GM20239_FMR1_Illumina_guided_HP2_recovered)

# GM20239 HTT HP1 concatenating and HP2 cleanup
GM20239_HTT_Illumina_guided_HP1 <- as.data.frame(GM20239_HTT_Illumina_guided_HP1[,1])
GM20239_HTT_Illumina_guided_HP1_recovered <- as.data.frame(GM20239_HTT_Illumina_guided_HP1_recovered[,1])
colnames(GM20239_HTT_Illumina_guided_HP1) <- c("read_id")
colnames(GM20239_HTT_Illumina_guided_HP1_recovered) <- c("read_id")
GM20239_HTT_Illumina_guided_HP1_with_recovered <- rbind(GM20239_HTT_Illumina_guided_HP1, GM20239_HTT_Illumina_guided_HP1_recovered)

GM20239_HTT_Illumina_guided_HP2 <- as.data.frame(GM20239_HTT_Illumina_guided_HP2[,1])
colnames(GM20239_HTT_Illumina_guided_HP2) <- c("read_id")


#############################################

## Load tables with Longshot read-names
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/SNV_analyse/longshot_readnames_hp")

GM04284_HTT_longshot_readnames <- na.omit(read.table("GM04284_longshot_hp_chr4_readnames.txt", header = F))
GM04284_FMR1_longshot_readnames <- na.omit(read.table("GM04284_longshot_hp_chrX_readnames.txt", header = F))
GM05538_HTT_longshot_readnames <- na.omit(read.table("GM05538_longshot_hp_chr4_readnames.txt", header = F)) #edit NAs in data frame
GM05538_FMR1_longshot_readnames <- na.omit(read.table("GM05538_longshot_hp_chrX_readnames.txt", header = F)) #edit NAs in data frame
GM06891_HTT_longshot_readnames <- na.omit(read.table("GM06891_longshot_hp_chr4_readnames.txt", header = F)) #edit NAs in data frame
GM06891_FMR1_longshot_readnames <- na.omit(read.table("GM06891_longshot_hp_chrX_readnames.txt", header = F)) #edit NAs in data frame
GM07541_HTT_longshot_readnames <- na.omit(read.table("GM07541_longshot_hp_chr4_readnames.txt", header = F))
GM07541_FMR1_longshot_readnames <- na.omit(read.table("GM07541_longshot_hp_chrX_readnames.txt", header = F)) #edit NAs in data frame
GM07861_HTT_longshot_readnames <- na.omit(read.table("GM07861_longshot_hp_chr4_readnames.txt", header = F))
GM07861_FMR1_longshot_readnames <- na.omit(read.table("GM07861_longshot_hp_chrX_readnames.txt", header = F))
GM20239_HTT_longshot_readnames <- na.omit(read.table("GM20239_longshot_hp_chr4_readnames.txt", header = F)) #edit NAs in data frame
GM20239_FMR1_longshot_readnames <- na.omit(read.table("GM20239_longshot_hp_chrX_readnames.txt", header = F)) #edit NAs in data frame

## Subset Longshot data frames into haplotypes
GM04284_HTT_longshot_hp1 <- subset(GM04284_HTT_longshot_readnames[GM04284_HTT_longshot_readnames$V2 == "HP:i:1",])
GM04284_HTT_longshot_hp2 <- subset(GM04284_HTT_longshot_readnames[GM04284_HTT_longshot_readnames$V2 == "HP:i:2",])
GM04284_FMR1_longshot_hp1 <- subset(GM04284_FMR1_longshot_readnames[GM04284_FMR1_longshot_readnames$V2 == "HP:i:1",])
GM04284_FMR1_longshot_hp2 <- subset(GM04284_FMR1_longshot_readnames[GM04284_FMR1_longshot_readnames$V2 == "HP:i:2",])

GM05538_HTT_longshot_hp1 <- subset(GM05538_HTT_longshot_readnames[GM05538_HTT_longshot_readnames$V2 == "HP:i:1",])
GM05538_HTT_longshot_hp2 <- subset(GM05538_HTT_longshot_readnames[GM05538_HTT_longshot_readnames$V2 == "HP:i:2",])
GM05538_FMR1_longshot_hp1 <- subset(GM05538_FMR1_longshot_readnames[GM05538_FMR1_longshot_readnames$V2 == "HP:i:1",])
GM05538_FMR1_longshot_hp2 <- subset(GM05538_FMR1_longshot_readnames[GM05538_FMR1_longshot_readnames$V2 == "HP:i:2",]) 

GM06891_HTT_longshot_readnames <- na.omit(GM06891_HTT_longshot_readnames)
GM06891_HTT_longshot_hp1 <- subset(GM06891_HTT_longshot_readnames[GM06891_HTT_longshot_readnames$V2 == "HP:i:1",])
GM06891_HTT_longshot_hp2 <- subset(GM06891_HTT_longshot_readnames[GM06891_HTT_longshot_readnames$V2 == "HP:i:2",])
GM06891_FMR1_longshot_hp1 <- subset(GM06891_FMR1_longshot_readnames[GM06891_FMR1_longshot_readnames$V2 == "HP:i:1",])
GM06891_FMR1_longshot_hp2 <- subset(GM06891_FMR1_longshot_readnames[GM06891_FMR1_longshot_readnames$V2 == "HP:i:2",]) 

GM07541_HTT_longshot_hp1 <- subset(GM07541_HTT_longshot_readnames[GM07541_HTT_longshot_readnames$V2 == "HP:i:1",])
GM07541_HTT_longshot_hp2 <- subset(GM07541_HTT_longshot_readnames[GM07541_HTT_longshot_readnames$V2 == "HP:i:2",])
GM07541_FMR1_longshot_hp1 <- subset(GM07541_FMR1_longshot_readnames[GM07541_FMR1_longshot_readnames$V2 == "HP:i:1",])
GM07541_FMR1_longshot_hp2 <- subset(GM07541_FMR1_longshot_readnames[GM07541_FMR1_longshot_readnames$V2 == "HP:i:2",]) 

GM07861_HTT_longshot_hp1 <- subset(GM07861_HTT_longshot_readnames[GM07861_HTT_longshot_readnames$V2 == "HP:i:1",])
GM07861_HTT_longshot_hp2 <- subset(GM07861_HTT_longshot_readnames[GM07861_HTT_longshot_readnames$V2 == "HP:i:2",])
GM07861_FMR1_longshot_hp1 <- subset(GM07861_FMR1_longshot_readnames[GM07861_FMR1_longshot_readnames$V2 == "HP:i:1",])
GM07861_FMR1_longshot_hp2 <- subset(GM07861_FMR1_longshot_readnames[GM07861_FMR1_longshot_readnames$V2 == "HP:i:2",]) 

GM20239_HTT_longshot_hp1 <- subset(GM20239_HTT_longshot_readnames[GM20239_HTT_longshot_readnames$V2 == "HP:i:1",])
GM20239_HTT_longshot_hp2 <- subset(GM20239_HTT_longshot_readnames[GM20239_HTT_longshot_readnames$V2 == "HP:i:2",])
GM20239_FMR1_longshot_hp1 <- subset(GM20239_FMR1_longshot_readnames[GM20239_FMR1_longshot_readnames$V2 == "HP:i:1",])
GM20239_FMR1_longshot_hp2 <- subset(GM20239_FMR1_longshot_readnames[GM20239_FMR1_longshot_readnames$V2 == "HP:i:2",]) 


################################
## Create Venn diagrams
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/SNV_analyse/venndiagrams")
library(VennDiagram)
library(scales)
# The colors for the sets are set1=yellow and set2=blue

# GM05538 HTT Illumina left VS Longshot
GM05538_HTT_Illumina_guided_HP1_left_set <- paste(GM05538_HTT_Illumina_guided_HP1_left$read_id , sep="")
GM05538_HTT_Illumina_guided_HP2_left_set <- paste(GM05538_HTT_Illumina_guided_HP2_left$read_id , sep="")
GM05538_HTT_Illumina_guided_HP1_left_with_recovered_set <- paste(GM05538_HTT_Illumina_guided_HP1_left_with_recovered$read_id , sep="")
GM05538_HTT_Illumina_guided_HP2_left_with_recovered_set <- paste(GM05538_HTT_Illumina_guided_HP2_left_with_recovered$read_id , sep="")
GM05538_HTT_longshot_hp1_set <- paste(GM05538_HTT_longshot_hp1$V1 , sep="")
GM05538_HTT_longshot_hp2_set <- paste(GM05538_HTT_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_left_with_recovered_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_left_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 585), cat.fontface = "bold", cat.dist = c(0.125, 0.125)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_left_with_recovered_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_left_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(35, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP2_left_with_recovered_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP2_left_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 585), cat.fontface = "bold", cat.dist = c(0.125, 0.125)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP2_left_with_recovered_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP2_left_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(35, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_left_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_left_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 585), cat.fontface = "bold", cat.dist = c(0.125, 0.125)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_left_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_left_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(35, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP2_left_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP2_left_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 585), cat.fontface = "bold", cat.dist = c(0.125, 0.125)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP2_left_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP2_left_VS__Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(35, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)




# GM05538 HTT Illumina right VS Longshot With recovered reads
GM05538_HTT_Illumina_guided_HP1_right_with_recovered_set <- paste(GM05538_HTT_Illumina_guided_HP1_right_with_recovered$read_id , sep="")
GM05538_HTT_Illumina_guided_HP1_right_set <- paste(GM05538_HTT_Illumina_guided_HP1_right$read_id , sep="")
GM05538_HTT_Illumina_guided_HP2_right_set <- paste(GM05538_HTT_Illumina_guided_HP2_right$read_id , sep="")
GM05538_HTT_longshot_hp1_set <- paste(GM05538_HTT_longshot_hp1$V1 , sep="")
GM05538_HTT_longshot_hp2_set <- paste(GM05538_HTT_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_right_with_recovered_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_right_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)


venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_right_with_recovered_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_right_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_right_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_right_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)


venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_right_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_right_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP2_right_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP2_right_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP2_right_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP2_right_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)

  
  
  
# GM05538 HTT Illumina right VS Longshot Without recovered reads
GM05538_HTT_Illumina_guided_HP1_right_with_recovered_set <- paste(GM05538_HTT_Illumina_guided_HP1_right_with_recovered$read_id , sep="")
GM05538_HTT_Illumina_guided_HP1_right_set <- paste(GM05538_HTT_Illumina_guided_HP1_right$read_id , sep="")
GM05538_HTT_Illumina_guided_HP2_right_set <- paste(GM05538_HTT_Illumina_guided_HP2_right$read_id , sep="")
GM05538_HTT_longshot_hp1_set <- paste(GM05538_HTT_longshot_hp1$V1 , sep="")
GM05538_HTT_longshot_hp2_set <- paste(GM05538_HTT_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_right_with_recovered_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_right_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_right_with_recovered_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_right_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_right_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_right_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP1_right_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP1_right_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP2_right_set, GM05538_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP2_right_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 220), cat.fontface = "bold", cat.dist = c(0.105, 0.2)
)

venn.diagram(
  x = list(GM05538_HTT_Illumina_guided_HP2_right_set, GM05538_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM05538_HTT_Illumina_HP2_right_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)





# GM06891 HTT Illumina VS Longshot
GM06891_HTT_Illumina_guided_HP1_set <- paste(GM06891_HTT_Illumina_guided_HP1$read_id , sep="")
GM06891_HTT_Illumina_guided_HP2_set <- paste(GM06891_HTT_Illumina_guided_HP2$read_id , sep="")
GM06891_HTT_longshot_hp1_set <- paste(GM06891_HTT_longshot_hp1$V1 , sep="")
GM06891_HTT_longshot_hp2_set <- paste(GM06891_HTT_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM06891_HTT_Illumina_guided_HP1_set, GM06891_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM06891_HTT_Illumina_HP1_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)
  
venn.diagram(
  x = list(GM06891_HTT_Illumina_guided_HP1_set, GM06891_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM06891_HTT_Illumina_HP1_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM06891_HTT_Illumina_guided_HP2_set, GM06891_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM06891_HTT_Illumina_HP2_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM06891_HTT_Illumina_guided_HP2_set, GM06891_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM06891_HTT_Illumina_HP2_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)




# GM07541 HTT Illumina VS Longshot
GM07541_HTT_Illumina_guided_HP1_with_recovered_set <- paste(GM07541_HTT_Illumina_guided_HP1_with_recovered$read_id , sep="")
GM07541_HTT_Illumina_guided_HP1_set <- paste(GM07541_HTT_Illumina_guided_HP1$read_id , sep="")
GM07541_HTT_Illumina_guided_HP2_set <- paste(GM07541_HTT_Illumina_guided_HP2$read_id , sep="")
GM07541_HTT_longshot_hp1_set <- paste(GM07541_HTT_longshot_hp1$V1 , sep="")
GM07541_HTT_longshot_hp2_set <- paste(GM07541_HTT_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM07541_HTT_Illumina_guided_HP1_with_recovered_set, GM07541_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07541_HTT_Illumina_HP1_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_HTT_Illumina_guided_HP1_with_recovered_set, GM07541_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07541_HTT_Illumina_HP1_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_HTT_Illumina_guided_HP1_set, GM07541_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07541_HTT_Illumina_HP1_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_HTT_Illumina_guided_HP1_set, GM07541_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07541_HTT_Illumina_HP1_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_HTT_Illumina_guided_HP2_set, GM07541_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07541_HTT_Illumina_HP2_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_HTT_Illumina_guided_HP2_set, GM07541_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07541_HTT_Illumina_HP2_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)



# GM07541 FMR1 Illumina VS Longshot
GM07541_FMR1_Illumina_guided_HP1_with_recovered_set <- paste(GM07541_FMR1_Illumina_guided_HP1_with_recovered$read_id , sep="")
GM07541_FMR1_Illumina_guided_HP1_set <- paste(GM07541_FMR1_Illumina_guided_HP1$read_id , sep="")
GM07541_FMR1_Illumina_guided_HP2_set <- paste(GM07541_FMR1_Illumina_guided_HP2$read_id , sep="")
GM07541_FMR1_longshot_hp1_set <- paste(GM07541_FMR1_longshot_hp1$V1 , sep="")
GM07541_FMR1_longshot_hp2_set <- paste(GM07541_FMR1_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM07541_FMR1_Illumina_guided_HP1_with_recovered_set, GM07541_FMR1_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07541_FMR1_Illumina_HP1_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_FMR1_Illumina_guided_HP1_with_recovered_set, GM07541_FMR1_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07541_FMR1_Illumina_HP1_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_FMR1_Illumina_guided_HP1_set, GM07541_FMR1_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07541_FMR1_Illumina_HP1_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_FMR1_Illumina_guided_HP1_set, GM07541_FMR1_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07541_FMR1_Illumina_HP1_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_FMR1_Illumina_guided_HP2_set, GM07541_FMR1_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07541_FMR1_Illumina_HP2_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07541_FMR1_Illumina_guided_HP2_set, GM07541_FMR1_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07541_FMR1_Illumina_HP2_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)


# GM07861 HTT Illumina VS Longshot
GM07861_HTT_Illumina_guided_HP1_with_recovered_set <- paste(GM07861_HTT_Illumina_guided_HP1_with_recovered$read_id , sep="")
GM07861_HTT_Illumina_guided_HP2_with_recovered_set <- paste(GM07861_HTT_Illumina_guided_HP2_with_recovered$read_id , sep="")
GM07861_HTT_Illumina_guided_HP1_set <- paste(GM07861_HTT_Illumina_guided_HP1$read_id , sep="")
GM07861_HTT_Illumina_guided_HP2_set <- paste(GM07861_HTT_Illumina_guided_HP2$read_id , sep="")
GM07861_HTT_longshot_hp1_set <- paste(GM07861_HTT_longshot_hp1$V1 , sep="")
GM07861_HTT_longshot_hp2_set <- paste(GM07861_HTT_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM07861_HTT_Illumina_guided_HP1_with_recovered_set, GM07861_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07861_HTT_Illumina_HP1_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07861_HTT_Illumina_guided_HP1_with_recovered_set, GM07861_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07861_HTT_Illumina_HP1_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)


venn.diagram(
  x = list(GM07861_HTT_Illumina_guided_HP1_set, GM07861_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07861_HTT_Illumina_HP1_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07861_HTT_Illumina_guided_HP1_set, GM07861_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07861_HTT_Illumina_HP1_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07861_HTT_Illumina_guided_HP2_with_recovered_set, GM07861_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07861_HTT_Illumina_HP2_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07861_HTT_Illumina_guided_HP2_with_recovered_set, GM07861_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07861_HTT_Illumina_HP2_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)


venn.diagram(
  x = list(GM07861_HTT_Illumina_guided_HP2_set, GM07861_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM07861_HTT_Illumina_HP2_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM07861_HTT_Illumina_guided_HP2_set, GM07861_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM07861_HTT_Illumina_HP2_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)


# GM20239 HTT Illumina VS Longshot
GM20239_HTT_Illumina_guided_HP1_with_recovered_set <- paste(GM20239_HTT_Illumina_guided_HP1_with_recovered$read_id , sep="")
GM20239_HTT_Illumina_guided_HP1_set <- paste(GM20239_HTT_Illumina_guided_HP1$read_id , sep="")
GM20239_HTT_Illumina_guided_HP2_set <- paste(GM20239_HTT_Illumina_guided_HP2$read_id , sep="")
GM20239_HTT_longshot_hp1_set <- paste(GM20239_HTT_longshot_hp1$V1 , sep="")
GM20239_HTT_longshot_hp2_set <- paste(GM20239_HTT_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM20239_HTT_Illumina_guided_HP1_with_recovered_set, GM20239_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM20239_HTT_Illumina_HP1_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM20239_HTT_Illumina_guided_HP1_with_recovered_set, GM20239_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM20239_HTT_Illumina_HP1_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM20239_HTT_Illumina_guided_HP1_set, GM20239_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM20239_HTT_Illumina_HP1_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM20239_HTT_Illumina_guided_HP1_set, GM20239_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM20239_HTT_Illumina_HP1_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM20239_HTT_Illumina_guided_HP2_set, GM20239_HTT_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM20239_HTT_Illumina_HP2_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)

venn.diagram(
  x = list(GM20239_HTT_Illumina_guided_HP2_set, GM20239_HTT_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM20239_HTT_Illumina_HP2_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(45, 190), cat.fontface = "bold", cat.dist = c(0.145, 0.04)
)


# GM20239 FMR1 Illumina VS Longshot
GM20239_FMR1_Illumina_guided_HP1_with_recovered_set <- paste(GM20239_FMR1_Illumina_guided_HP1_with_recovered$read_id , sep="")
GM20239_FMR1_Illumina_guided_HP2_with_recovered_set <- paste(GM20239_FMR1_Illumina_guided_HP2_with_recovered$read_id , sep="")
GM20239_FMR1_Illumina_guided_HP1_set <- paste(GM20239_FMR1_Illumina_guided_HP1$read_id , sep="")
GM20239_FMR1_Illumina_guided_HP2_set <- paste(GM20239_FMR1_Illumina_guided_HP2$read_id , sep="")
GM20239_FMR1_longshot_hp1_set <- paste(GM20239_FMR1_longshot_hp1$V1 , sep="")
GM20239_FMR1_longshot_hp2_set <- paste(GM20239_FMR1_longshot_hp2$V1 , sep="")

venn.diagram(
  x = list(GM20239_FMR1_Illumina_guided_HP1_with_recovered_set, GM20239_FMR1_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM20239_FMR1_Illumina_HP1_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 575), cat.fontface = "bold", cat.dist = c(0.06, 0.125)
)

venn.diagram(
  x = list(GM20239_FMR1_Illumina_guided_HP1_with_recovered_set, GM20239_FMR1_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM20239_FMR1_Illumina_HP1_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(0, 565), cat.fontface = "bold", cat.dist = c(0.03, 0.125)
)

venn.diagram(
  x = list(GM20239_FMR1_Illumina_guided_HP1_set, GM20239_FMR1_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM20239_FMR1_Illumina_HP1_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(30, 575), cat.fontface = "bold", cat.dist = c(0.06, 0.125)
)

venn.diagram(
  x = list(GM20239_FMR1_Illumina_guided_HP1_set, GM20239_FMR1_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM20239_FMR1_Illumina_HP1_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(0, 565), cat.fontface = "bold", cat.dist = c(0.03, 0.125)
)



venn.diagram(
  x = list(GM20239_FMR1_Illumina_guided_HP2_with_recovered_set, GM20239_FMR1_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM20239_FMR1_Illumina_HP2_with_recovered_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(20, 575), cat.fontface = "bold", cat.dist = c(0.055, 0.15)
)

venn.diagram(
  x = list(GM20239_FMR1_Illumina_guided_HP2_with_recovered_set, GM20239_FMR1_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM20239_FMR1_Illumina_HP2_with_recovered_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(0, 565), cat.fontface = "bold", cat.dist = c(0.03, 0.125)
)


venn.diagram(
  x = list(GM20239_FMR1_Illumina_guided_HP2_set, GM20239_FMR1_longshot_hp1_set),
  category.names = c("" , ""),
  filename = 'GM20239_FMR1_Illumina_HP2_VS_Longshot_HP1.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(20, 575), cat.fontface = "bold", cat.dist = c(0.055, 0.15)
)

venn.diagram(
  x = list(GM20239_FMR1_Illumina_guided_HP2_set, GM20239_FMR1_longshot_hp2_set),
  category.names = c("" , ""),
  filename = 'GM20239_FMR1_Illumina_HP2_VS_Longshot_HP2.png',
  output=T, fill = c(alpha("#ffff00",0.3), alpha('#00ccff',0.3)), cex = 5, cat.cex = 1.75,
  cat.pos = c(0, 565), cat.fontface = "bold", cat.dist = c(0.03, 0.125)
)
























