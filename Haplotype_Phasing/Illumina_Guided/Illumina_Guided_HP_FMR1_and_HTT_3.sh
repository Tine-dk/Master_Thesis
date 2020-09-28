### Illumina Guided HP ####################################
cd /work/sdularsen/tine/nanoporeData/

Window=40000
AF_low=30
AF_high=70
DP=10

# Go through the listed genes
for gene in "FMR1" "HTT"; do
	if [ $gene == "FMR1" ]; then
		# FMR1
		chromosome=X
		Start_Position=$(( 146993569 - $Window ))
		End_Position=$(( 146993598  + $Window ))
	elif [ $gene == "HTT" ]; then
		# HTT
		chromosome=4
		Start_Position=$(( 3076603 - $Window ))
		End_Position=$(( 3076774 + $Window ))
	else
		echo "Gene $gene not supported"
		break
	fi

	for file in Illumina_data/*vcf.gz; do
		vcf=`basename $file`
		SAMPLE=${vcf%_truseq-pcr-free-genome_*.vcf.gz}
		echo $SAMPLE $gene

		### SAMtools
		# Extracting a smaller BAM
		samtools view -bh alignedMinimap2/${SAMPLE}.bam "$chromosome:${Start_Position}-${End_Position}" > illumina_guided_hp/ONT_BAM_subset/${SAMPLE}_${gene}_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.bam
		samtools index illumina_guided_hp/ONT_BAM_subset/${SAMPLE}_${gene}_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.bam

		# Export the readnames into a txt file
		samtools view illumina_guided_hp/ONT_BAM_subset/${SAMPLE}_${gene}_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.bam | cut -f 1,3,4 > illumina_guided_hp/ONT_BAM_subset_readnames/${SAMPLE}_${gene}_Read_Name_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.txt



		### BCFtools
		# Creation of an index file
		INDEX_FILE=Illumina_data/${vcf}.tbi
		if [ ! -f "$INDEX_FILE" ]; then
			tabix -p vcf Illumina_data/${vcf}
		fi

		# Subset Illumina Dragen VCF file to ROI +/- the window
		bcftools filter -r $chromosome:${Start_Position}-${End_Position} Illumina_data/${vcf} > illumina_guided_hp/illumina_VCF_subset/${SAMPLE}_${gene}_Illumina_Dragen_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf

		# Filter VCF file
		bcftools filter -i "FORMAT/DP>=${DP} && FORMAT/AF>=0.${AF_low} && FORMAT/AF<=0.${AF_high} && STRLEN(REF) == 1 && STRLEN(ALT) == 1 && N_ALT==1" illumina_guided_hp/illumina_VCF_subset/${SAMPLE}_${gene}_Illumina_Dragen_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf \
		> illumina_guided_hp/illumina_VCF_subset_filtered/${SAMPLE}_${gene}_Illumina_Dragen_Filtered_SNVs_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf

		#Extract Illumina Genotypes from the filtered vcf file
		grep -v '^##' illumina_guided_hp/illumina_VCF_subset_filtered/${SAMPLE}_${gene}_Illumina_Dragen_Filtered_SNVs_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf | cut -f 1,2,4,5 | sed 's/#//g' > \
		illumina_guided_hp/illumina_VCF_genotypes/${SAMPLE}_${gene}_Illumina_Dragen_genotypes_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.txt

		# Antal linjer/SNVs tilbage efter filtrering
		snvs=$(bcftools filter -i "FORMAT/DP>=${DP} && FORMAT/AF>=0.${AF_low} && FORMAT/AF<=0.${AF_high} && STRLEN(REF) == 1 && STRLEN(ALT) == 1 && N_ALT==1" illumina_guided_hp/illumina_VCF_subset/${SAMPLE}_${gene}_Illumina_Dragen_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf | grep -v '^#' | wc -l)
		echo $snvs


		### SNV Bases
		if [ $snvs -gt 0 ]; then
			python /work/sdularsen/tine/nanoporeData/kommandoer/SNV_Base_Identification3.py \
				--sample $SAMPLE --window $Window --gene $gene \
				--chr $chromosome --start $Start_Position --end $End_Position \
				--AF_low $AF_low --AF_high $AF_high --DP $DP
		fi
	done
done
