### Illumina Guided HP ####################################
cd /work/sdularsen/tine/nanoporeData/

Window=40000
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

# Find VCF filer fra LongShot phasningen
	for file in longshot_ny/merged_vcf/*.vcf.gz; do
		vcf=`basename $file`
		SAMPLE=${vcf%_longshot_hp_*.vcf.gz}
		echo $SAMPLE $gene

		### SAMtools
		# Extracting a smaller BAM
		samtools view -bh alignedMinimap2/"G87-"${SAMPLE}.bam "$chromosome:${Start_Position}-${End_Position}" > longshot_ny2/ONT_BAM_subset/longshot_${SAMPLE}_${gene}_window${Window}_DP${DP}.bam
		samtools index longshot_ny2/ONT_BAM_subset/longshot_${SAMPLE}_${gene}_window${Window}_DP${DP}.bam

		# Export the readnames into a txt file
		samtools view longshot_ny2/ONT_BAM_subset/longshot_${SAMPLE}_${gene}_window${Window}_DP${DP}.bam | cut -f 1,3,4 > longshot_ny2/ONT_BAM_subset_readnames/longshot_${SAMPLE}_${gene}_Read_Name_window${Window}_DP${DP}.txt


		echo "vcf file name:" $vcf
		### BCFtools
		# Creation of an index file
		INDEX_FILE=longshot_ny/merged_vcf/${vcf}.tbi
		if [ ! -f "$INDEX_FILE" ]; then
			tabix -p vcf longshot_ny/merged_vcf/${vcf}
		fi

		# Subset Illumina Dragen VCF file to ROI +/- the window
		bcftools filter -r $chromosome:${Start_Position}-${End_Position} longshot_ny/merged_vcf/${vcf} > longshot_ny2/longshot_vcf_subset/${SAMPLE}_${gene}_longshot_window${Window}_DP${DP}.vcf

		# Filter VCF file
		bcftools filter -i "INFO/DP>=${DP} && STRLEN(REF) == 1 && STRLEN(ALT) == 1 && N_ALT==1" longshot_ny2/longshot_vcf_subset/${SAMPLE}_${gene}_longshot_window${Window}_DP${DP}.vcf \
		> longshot_ny2/longshot_vcf_subset_filtered/${SAMPLE}_${gene}_longshot_Filtered_SNVs_window${Window}_DP${DP}.vcf

		bcftools view -g het longshot_ny2/longshot_vcf_subset_filtered/${SAMPLE}_${gene}_longshot_Filtered_SNVs_window${Window}_DP${DP}.vcf > longshot_ny2/longshot_vcf_subset_filtered_het/${SAMPLE}_${gene}_longshot_Filtered_het_SNVs_window${Window}_DP${DP}.vcf

		#Extract Illumina Genotypes from the filtered vcf file
		grep -v '^##' longshot_ny2/longshot_vcf_subset_filtered_het/${SAMPLE}_${gene}_longshot_Filtered_het_SNVs_window${Window}_DP${DP}.vcf | cut -f 1,2,4,5 | sed 's/#//g' > \
		longshot_ny2/longshot_vcf_het_genotypes/${SAMPLE}_${gene}_longshot_het_genotypes_window${Window}_DP${DP}.txt

		# Antal linjer/SNVs tilbage efter filtrering
		#snvs=$(bcftools filter -i "INFO/DP>=${DP} && STRLEN(REF) == 1 && STRLEN(ALT) == 1 && N_ALT==1" longshot_ny2/longshot_vcf_subset/${SAMPLE}_${gene}_longshot_window${Window}_DP${DP}.vcf | grep -v '^#' | wc -l)
		#echo $snvs

		snvs=$(bcftools view -g het longshot_ny2/longshot_vcf_subset_filtered/${SAMPLE}_${gene}_longshot_Filtered_SNVs_window${Window}_DP${DP}.vcf | grep -v '^#' | wc -l)
		echo $snvs


		### SNV Bases
		# Skal chr, start og end initieres l ngere oppe?
		if [ $snvs -gt 0 ]; then
			python /work/sdularsen/tine/nanoporeData/kommandoer/Longshot_SNV_Base_Identification4_heterozygote.py \
				--sample $SAMPLE --window $Window --gene $gene \
				--chr $chromosome --start $Start_Position --end $End_Position \
				--DP $DP
		fi
	done
done
