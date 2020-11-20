# Libraries
library(tidyverse)
library(hrbrthemes)
library(ggplot2)

## Genome mean depth
#   counting only mapped (primary aligned) reads: samtools view -c -F 260 SAMPLE.bam
reference_length_hg19 <- 3101788170

# GM06891
GM06891_reads_mapped <- 13482122
GM06891_bases_generated <- 7380000000
GM06891_reads_generated <- 13980000
GM06891_mean_read_length <- GM06891_bases_generated/GM06891_reads_generated
GM06891_genome_depth <- (GM06891_reads_mapped*GM06891_mean_read_length)/reference_length_hg19 

# GM07541
GM07541_reads_mapped <- 11133665
GM07541_bases_generated <- 5.540.000.000+1.280.000.000
GM07541_reads_generated <- 9180000+2310000
GM07541_mean_read_length <- GM07541_bases_generated/GM07541_reads_generated
GM07541_genome_depth <- (GM07541_reads_mapped*GM07541_mean_read_length)/reference_length_hg19 

# GM07861
GM07861_reads_mapped <- 25260984
GM07861_bases_generated <- 17460000000
GM07861_reads_generated <- 26370000
GM07861_mean_read_length <- GM07861_bases_generated/GM07861_reads_generated
GM07861_genome_depth <- (GM07861_reads_mapped*GM07861_mean_read_length)/reference_length_hg19 

# GM20239
GM20239_reads_mapped <- 21326359
GM20239_bases_generated <- 15860000000
GM20239_reads_generated <- 21960000
GM20239_mean_read_length <- GM20239_bases_generated/GM20239_reads_generated
GM20239_genome_depth <- (GM20239_reads_mapped*GM20239_mean_read_length)/reference_length_hg19 






# Load target depth file
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/LAB/GM06891")
GM06891_targets <- read.table("targets_depth.txt", header = T)
GM06891_targets$fold_change <- c(GM06891_targets$Mean_Depth/GM06891_genome_depth)
GM06891_targets_subset <- GM06891_targets[,c(4,6)]
colnames(GM06891_targets_subset) <- c("x","y")


setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/LAB/GM07541")
GM07541_targets <- read.table("targets_depth.txt", header = T)
GM07541_targets$fold_change <- c(GM07541_targets$mean_depth/GM07541_genome_depth)
GM07541_targets_subset <- GM07541_targets[,c(4,6)]
colnames(GM07541_targets_subset) <- c("x","y")

setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/LAB/GM07861")
GM07861_targets <- read.table("fold-change.txt", header = T)
GM07861_targets$fold_change <- c(GM07861_targets$mean_depth/GM07861_genome_depth)
GM07861_targets_subset <- GM07861_targets[,c(4,6)]
colnames(GM07861_targets_subset) <- c("x","y")


setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/LAB/GM20239")
GM20239_targets <- read.table("GM20239_target_depth.txt", header = T)
GM20239_targets$fold_change <- c(GM20239_targets$mean_depth/GM20239_genome_depth)
GM20239_targets_subset <- GM20239_targets[,c(4,6)]
colnames(GM20239_targets_subset) <- c("x","y")





# plot

# GM06891
ggplot(GM06891_targets_subset, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + labs(x = "", y = "Fold-Change") +
  ggtitle("GM06891 Fast Basecalling") +
  theme(plot.title = element_text(hjust = 0.4))


# GM07541
ggplot(GM07541_targets_subset, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + labs(x = "", y = "Fold-Change") +
  ggtitle("GM07541 HAC Basecalling") +
  theme(plot.title = element_text(hjust = 0.4))


# GM07861
ggplot(GM07861_targets_subset, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + labs(x = "", y = "Fold-Change") +
  ggtitle("GM07861 HAC Basecalling") +
  theme(plot.title = element_text(hjust = 0.4))



# GM20239
ggplot(GM20239_targets_subset, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + labs(x = "", y = "Fold-Change") +
  ggtitle("GM20239 HAC Basecalling + Large Target") +
  theme(plot.title = element_text(hjust = 0.4))








