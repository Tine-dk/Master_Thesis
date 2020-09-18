library(dplyr)
library(proxy)
library(gower)

# Set working directory
setwd("C:/Users/Taine/OneDrive/Skrivebord/Nanopore/SNV_analyse/illumina_guided_results_dataframes")



############################################################################
### Helper functions ###
############################################################################
get_gene_information <- function(gene = "HTT"){
    # HTT gene information
    if (gene == "HTT"){
        chr <<- 4
        start_pos <<- 3076603
        end_pos <<- 3076774 # HTT repeat ends at 307660 but there are GCC repeats until 3076774 which is being ignored
    } else if (gene == "FMR1"){
        # FMR1
        chr <<- "X"
        start_pos <<- 146993569
        end_pos <<- 146993598
    } else {
        print("gene not supported")
    }
    # return(c(chr,start_pos,end_pos))
}
    
# Remove snps not in close proximity to gene
remove_distant_reads <- function(data, snp_left, snp_right){
    
    # Get reads with information on nearest snp on either side of gene
    if (nrow(snp_left)>0){
        reads_left = !is.na(data[,which(names(data)==snp_left$snp)])
    } else {
        reads_left = F
    }
    if (nrow(snp_right)>0){
        reads_right = !is.na(data[,which(names(data)==snp_right$snp)])
    } else {
        reads_right = F
    }
    
    # Keep only reads near or spanning gene
    data = data[reads_left|reads_right,]
    return(data)
}

# Keep only reads covering at least two snps
remove_uninformative_reads <- function(dat,min_snps){
    rows_to_keep=apply(dat, 1, function(row){
        # Check if at least "min_snps" snps exists in a read
        # note: snp_name is not NA and therefore counts +1 extra
        if (sum(!is.na(row))>=(min_snps+1)){
            return(TRUE)
        } else {
            # print(row)
            return(FALSE)
        }
    })
    rows_to_keep
}

# Check reads cover both snps
num_reads_cover_both_snps <- function(data2, snp_left, snp_right){
    
    # Get reads with information on nearest snp on either side of gene
    if (nrow(snp_left)>0 && nrow(snp_right)>0){
        return(sum(!is.na(data2[,snp_left$snp]) & !is.na(data2[,snp_right$snp])))
    } else if (nrow(snp_left)>0){
        return(sum(!is.na(data2[,snp_left$snp])))
    } else if (nrow(snp_right)>0){
        return(sum(!is.na(data2[,snp_right$snp])))
    } else {
        return(0)
    }

}

# Compute gower's distance
gower_distance <- function(x,y){
    # Check if no overlapping snps between two reads
    # if no overlapping snps, return 0.5 as their distance
    if (sum(!(is.na(x) | is.na(y)))==0){
        # print(x)
        # print(y)
        return(0.5)
    }
    # Compute and return gower's distance between two reads
    return(gower_dist(x,y))
}

# Hierarchical clustering of reads
hierarchical_clustering <- function(data, k=2){
    ### Hierarchical clustering ###
    dst = proxy::dist(data[,-1],method = gower_distance)
    # hc <- hclust(dst, method = "complete")
    hc <- hclust(dst, method = "ward.D2")
    cluster = cutree(hc, k=k)
    dend = as.dendrogram(hc)
    plot(dend)
    return(cluster)
}

# Construct haplotype from subset of reads
construct_haplotype <- function(h){
    apply(h[,-1], 2, function(col){
        d = sort(table(col), decreasing = T)
        # print(d)
        if (length(d)==0){
            return(NA)
        } else if (length(d)==1 || d[1]>d[2]){
            return(names(d[1]))
        } else if (d[1]==d[2]){
            return(paste0(names(d)[1], "/", names(d)[2]))
        }
    })
}

# Print some stats of the data
get_statistics <- function(data){
    # Count num snps for each read
    snps_per_read = apply(data,1,function(row){
        sum(!is.na(row))
    })
    print("Number of SNPs covered in each read")
    ## Print only snps_per_read count
    print(snps_per_read)
    ## Print read_id and snps_per_read count
    # print(data.frame(Read_ID=data$X, snps_per_read))
    
    # Count num reads covering every set of two snps
    # - if any is 0, then the reads must be split into separated groups 
    reads_covering_set_two_snps <- sapply(2:ncol(data)-1, function(i){
        sum(!is.na(data[,i]) & !is.na(data[,i+1]))
    })
    print("Number of reads covering every 2-set of snps")
    print(reads_covering_set_two_snps)
}

# Check if snp is found in haplotype
check_snp_in_hap <- function(read, hap, pos){
    if (!is.na(read[pos]) && read[pos] %in% strsplit(hap[pos], split = "/")[[1]]){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

############################################################################
## Main
############################################################################

input_files = list.files(pattern = "_bases_window40000_AF_0_30_to_0_70_DP20.txt")
#input_files = list.files(pattern = "_bases_")
for (file in input_files){
    
    # Define sample and gene from input file
    sample = strsplit(file, split = "_")[[1]][1]
    gene = strsplit(file, split = "_")[[1]][2]
    
    # read input
    data = read.table(file, header = T, sep = ",", na.strings = "NA")
    
    # Remove duplicate lines (that exist for some reason :S)
    data = distinct(data)
    # Get dimensionality of data
    dim(data)
    
    get_gene_information(gene)
    
    # Guess the most probable genotype at each snp
    bases = apply(data[,-1], 2, function(col){
        d = sort(table(col), decreasing = T)
        c(names(d)[1], names(d)[2])
    })
    
    # Get some stats on the data
    get_statistics(data)
    
    # Extract snp coordinates from SNP names
    snp_coords = data.frame(snp=names(data[,-1]), coord = as.integer(gsub(x = names(data[,-1]), pattern = "X", replacement = "")))
    # Get nearest snp to gene on either side
    snp_left = filter(snp_coords, coord <= start_pos) %>% slice_max(coord)
    snp_right = filter(snp_coords, coord >= end_pos) %>% slice_min(coord)
    # Remove reads not covering at least one of the closest snp on either side of the gene
    data = remove_distant_reads(data, snp_left, snp_right)
    
    # Remove "-" (deletions) as they are probably sequencing mistakes
    # data[data=="-"] = NA
    
    # read genotypes
    genotypes = read.table(gsub("bases", "Illumina_Dragen_genotypes", file), header = T)
    for(i in 1:nrow(genotypes)){
        genotypes_in_snp = c(genotypes$REF[i],genotypes$ALT[i])
        data[,i+1][!(data[,i+1] %in% genotypes_in_snp)] = NA
    }

    # Get some stats on the data
    get_statistics(data)
    
    # Remove reads covering only one snp
    rows_to_keep = remove_uninformative_reads(data, 2)
    data2 = data[rows_to_keep,]
    
    # Check if any reads exist covering both snps closest on either side to the gene, or only snps on one side exists
    reads_covering_both_snps = num_reads_cover_both_snps(data2, snp_left, snp_right)
    if (reads_covering_both_snps>0){
    
        # Cluster reads into two groups
        cluster = hierarchical_clustering(data2, 2)
        
        # Construct most likely haplotypes - possibly use logos plots
        hap1 = construct_haplotype(data2[cluster==1,])
        hap2 = construct_haplotype(data2[cluster==2,])
    
    } else {
        
        # Split data into two sets containing reads on either side of the gene
        data_left = data2[!is.na(data2[,snp_left$snp]),]
        data_right = data2[!is.na(data2[,snp_right$snp]),]
        
        # Cluster reads in each dataset into two groups
        cluster1 = hierarchical_clustering(data_left, 2)
        cluster2 = hierarchical_clustering(data_right, 2)
        
        # Construct haplotypes for cluster1
        hap1_left = construct_haplotype(data_left[cluster1==1,])
        hap2_left = construct_haplotype(data_left[cluster1==2,])
        # Construct haplotypes for cluster2
        hap1_right = construct_haplotype(data_right[cluster2==1,])
        hap2_right = construct_haplotype(data_right[cluster2==2,])
    }
    
    # Recover and cluster "uninformative" reads i.e. reads covering only one snp
    if (sum(!rows_to_keep)>0){
        uninformative_reads = data[!rows_to_keep,]
        uninformative_reads$haplotype = NA
        
        if (reads_covering_both_snps>0){
            for (i in 1:nrow(uninformative_reads)){
                # Get position of snp not being NA
                pos = which(!is.na(uninformative_reads[i,-1]))
                # If no snp found ignore read
                if (length(pos)==0){
                    next
                }
                
                # Check if snp matches both haplotypes, if so we cannot assign it to any cluster and we ignore it
                if (check_snp_in_hap(uninformative_reads[i,-1], hap1, pos) && check_snp_in_hap(uninformative_reads[i,-1], hap2, pos)){
                    next
                # Check if snp matches hap1, if so assign it to cluster1
                } else if (check_snp_in_hap(uninformative_reads[i,-1], hap1, pos)) {
                    uninformative_reads[i,"haplotype"] = "hap1"
                # Check if snp matches hap2, if so assign it to cluster1
                } else if (check_snp_in_hap(uninformative_reads[i,-1], hap2, pos)) {
                    uninformative_reads[i,"haplotype"] = "hap2"
                # If no matches we ignore the read
                } else {
                    next
                }
            }
        } else {
            for (i in 1:nrow(uninformative_reads)){
                # Get position of snp not being NA
                pos = which(!is.na(uninformative_reads[i,-1]))
                # If no snp found ignore read
                if (length(pos)==0){
                    next
                }
                
                # Check if snp matches both haplotypes, if so we cannot assign it to any cluster and we ignore it
                if (check_snp_in_hap(uninformative_reads[i,-1], hap1_left, pos) && check_snp_in_hap(uninformative_reads[i,-1], hap2_left, pos)){
                    next
                } 
                if (check_snp_in_hap(uninformative_reads[i,-1], hap1_right, pos) && check_snp_in_hap(uninformative_reads[i,-1], hap2_right, pos)){
                    next
                }
                
                # Check if snp matches hap1_left, if so assign it to hap1_left
                if (check_snp_in_hap(uninformative_reads[i,-1], hap1_left, pos)) {
                    uninformative_reads[i,"haplotype"] = "hap1_left"
                # Check if snp matches hap2_left, if so assign it to hap2_left
                } else if (check_snp_in_hap(uninformative_reads[i,-1], hap2_left, pos)) {
                    uninformative_reads[i,"haplotype"] = "hap2_left"
                # Check if snp matches hap1_right, if so assign it to hap1_right
                } else if (check_snp_in_hap(uninformative_reads[i,-1], hap1_right, pos)) {
                    uninformative_reads[i,"haplotype"] = "hap1_right"
                # Check if snp matches hap2_right, if so assign it to hap2_right
                } else if (check_snp_in_hap(uninformative_reads[i,-1], hap2_right, pos)) {
                    uninformative_reads[i,"haplotype"] = "hap1_right"
                }
            }
        }
    }
    
    # Write output files
    output_dir = "Output_haplotype_phasing_heterozygote/"
    dir.create(output_dir, showWarnings = F)
    
    output_files_for_haplotype <- function(data, hap, hap_name, cluster, clust_num){
        # Insert more sensible column names for the data table
        names(data) = c("read_id", as.character(snp_coords$coord))
        # Output reads related to haplotype 1
        write.table(data[cluster==clust_num,], file = paste0(output_dir, sample,"_",gene,"_reads_",hap_name,".txt"), quote = F, row.names = F)
        
        # Check if any reads in uninformative_reads
        if (sum(!rows_to_keep)>0){
            uninformative_reads_out = filter(uninformative_reads, haplotype == hap_name)
            # Check if any reads in subset of uninformative_reads
            if (nrow(uninformative_reads_out)>0){
                # Insert more sensible column names for the data table
                names(uninformative_reads_out) = c("read_id", as.character(snp_coords$coord), "haplotypes")
                # Write recovered reads to file
                write.table(uninformative_reads_out, file = paste0(output_dir, sample,"_",gene,"_reads_",hap_name,"_recovered.txt"), quote = F, row.names = F)
            }
        }
    }
    
    if (reads_covering_both_snps>0){
        # Output haplotypes
        hap1_out = paste(hap1, collapse = " ")
        hap2_out = paste(hap2, collapse = " ")
        write.table(paste0("Haplotype 1:\n",hap1_out,"\n\nHaplotype 2:\n", hap2_out), 
                    file = paste0(output_dir, sample,"_",gene,"_haplotypes.txt"), quote = F, row.names = F, col.names = F)
        
        # Write reads to files
        output_files_for_haplotype(data2, hap1, "hap1", cluster, 1)
        output_files_for_haplotype(data2, hap1, "hap2", cluster, 2)
        
    } else {
        # Output haplotypes
        hap1_left_out = paste(hap1_left, collapse = " ")
        hap2_left_out = paste(hap2_left, collapse = " ")
        hap1_right_out = paste(hap1_right, collapse = " ")
        hap2_right_out = paste(hap2_right, collapse = " ")
        write.table(paste0("hap1_left:\n",hap1_left_out,"\n\nhap2_left:\n", hap2_left_out, "\n\nhap1_right:\n",hap1_right_out,"\n\nhap2_right:\n", hap2_right_out), 
                    file = paste0(output_dir, sample,"_",gene,"_haplotypes.txt"), quote = F, row.names = F, col.names = F)
        
        # Write reads to files
        output_files_for_haplotype(data_left, hap1_left, "hap1_left", cluster1, 1)
        output_files_for_haplotype(data_left, hap2_left, "hap2_left", cluster1, 2)
        output_files_for_haplotype(data_right, hap1_right, "hap1_right", cluster2, 1)
        output_files_for_haplotype(data_right, hap2_right, "hap2_right", cluster2, 2)
    }
}
