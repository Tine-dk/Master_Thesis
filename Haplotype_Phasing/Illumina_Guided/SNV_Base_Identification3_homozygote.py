# Import libraries
import pysam
import pandas as pd
import sys
import argparse
path = "/work/sdularsen/tine/nanoporeData/"

# Add variables used in the code
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--sample', type = str, help = 'The sample to analyze')
parser.add_argument('--window', type = str, help = 'The window around the repeat')
parser.add_argument('--gene', type = str, help = 'The gene of interest')
parser.add_argument('--chr', type = str, help = 'The chromosome where the gene of interest is located')
parser.add_argument('--start', type = int, help = 'Start position on chromosome')
parser.add_argument('--end', type = int, help = 'End position on chromosome')
parser.add_argument('--AF_low', type = str, help = 'low end of AF interval')
parser.add_argument('--AF_high', type = str, help = 'high end of AF interval')
parser.add_argument('--DP', type = str, help = 'read depth from Illumina Data')
args = parser.parse_args()
sample = args.sample
gene = args.gene
window = args.window
chr = args.chr
start = args.start
end = args.end
AF_low = args.AF_low
AF_high = args.AF_high
DP = args.DP

# Load SNV positions into a set to be used in the data frame
snv_positions = set()
with open (path + "illumina_guided_hp/illumina_VCF_subset_filtered/"+sample+"_"+gene+"_Illumina_Dragen_Filtered_SNVs_window"+window+"_AF_0_"+AF_low+"_to_1_"+AF_high+"_DP"+DP+".vcf", "r") as f:
    for line in f:
        if not line.startswith("#"):
            snv_positions.add(int(line.split("\t")[1]))

# Load read names into a list to be used in the data frame
read_names = list()
with open (path + "illumina_guided_hp/ONT_BAM_subset_readnames/"+sample+"_"+gene+"_Read_Name_window"+window+"_AF_0_"+AF_low+"_to_1_"+AF_high+"_DP"+DP+".txt", "r") as f:
    for line in f:
        read_names.append(line.split("\t")[0])

# Create DataFrame object
snv_df = pd.DataFrame(pd.NA, columns = sorted(list(snv_positions)), index=read_names)

# Fill the dataframe by going through the SNV's and reads
bamfile = pysam.AlignmentFile(path + "illumina_guided_hp/ONT_BAM_subset/"+sample+"_"+gene+"_window"+window+"_AF_0_"+AF_low+"_to_1_"+AF_high+"_DP"+DP+".bam", "rb")
for pileupcolumn in bamfile.pileup(chr, start, end, min_base_quality=0):
    # python is 0-indexed, vcf-file is 1-indexed
    snv_pos = pileupcolumn.pos+1
    if snv_pos in snv_positions:
        for pileupread in pileupcolumn.pileups:
            read_name = pileupread.alignment.query_name
            ## checking for deletions
            if not pileupread.is_refskip:
                if pileupread.is_del:
                    snv_df.loc[read_name, snv_pos] = "-"
                else:
                    snv_df.loc[read_name, snv_pos] = pileupread.alignment.query_sequence[pileupread.query_position]

bamfile.close()

# Remove columns with only NAs
snv_df = snv_df.dropna(axis=1,how="all")
# Remove rows with only NAs
snv_df = snv_df.dropna(axis=0,how="all")


# Print data frame while keeping the NAs
print(snv_df, sep='\n')

# Save the data frame as a txt file
snv_df.to_csv(path + "illumina_guided_hp/results_dataframes/"+sample+"_"+gene+"_bases_window"+window+"_AF_0_"+AF_low+"_to_1_"+AF_high+"_DP"+DP+".txt", sep=',', mode='w',na_rep='NA')
