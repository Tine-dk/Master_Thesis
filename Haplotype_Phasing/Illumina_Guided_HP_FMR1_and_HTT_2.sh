### Illumina Guided HP ####################################
cd /work/sdularsen/tine/nanoporeData/
# Skal vinduet være 0 eller 100.000 her? Hvis der står 0 kommer der vidst ikke nogle tabeller ud. Men hvis der står 100000 får vi så for mange reads med?
# altså, nede hvor read navne tages ud af den mindre bam skal vinduet da ikke være 0
# og så når der kigges på SNV positioner i Illumina VCF skal vinduet være 100.000
# Kan vi evt fjerne $window fra definitionen af start_position og end_position og tilføje den detalje nede ved Illumina VCF filen hvor vinduet bruges til at få nok SNV positioner med?

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

# Jeg forstår ikke dette. Hvor angives 'basename'?
	for file in Illumina_data/*vcf.gz; do
		vcf=`basename $file`
		SAMPLE=${vcf%_truseq-pcr-free-genome_*.vcf.gz}
		echo $SAMPLE $gene

		### SAMtools
		# Extracting a smaller BAM
		samtools view -bh alignedMinimap2/${SAMPLE}.bam "$chromosome:${Start_Position}-${End_Position}" > alignedMinimap2/${SAMPLE}_${gene}_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.bam
		samtools index alignedMinimap2/${SAMPLE}_${gene}_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.bam

		# Export the readnames into a txt file
		samtools view alignedMinimap2/${SAMPLE}_${gene}_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.bam | cut -f 1,3,4 > alignedMinimap2/${SAMPLE}_${gene}_Read_Name_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.txt



		### BCFtools
		# Creation of an index file
		INDEX_FILE=Illumina_data/${vcf}.tbi
		if [ ! -f "$INDEX_FILE" ]; then
			tabix -p vcf Illumina_data/${vcf}
		fi

		# Subset Illumina Dragen VCF file to ROI +/- the window
		bcftools filter -r $chromosome:${Start_Position}-${End_Position} Illumina_data/${vcf} > Illumina_data/${SAMPLE}_${gene}_Illumina_Dragen_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf

		# Filter VCF file
		bcftools filter -i "FORMAT/DP>=10 && FORMAT/AF>=0.${AF_low} && FORMAT/AF<=0.${AF_high} && STRLEN(REF) == 1 && STRLEN(ALT) == 1 && N_ALT==1" Illumina_data/${SAMPLE}_${gene}_Illumina_Dragen_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf \
		> Illumina_data/${SAMPLE}_${gene}_Illumina_Dragen_Filtered_SNVs_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf

		#Extract Illumina Genotypes from the filtered vcf file
		grep -v '^##' Illumina_data/${SAMPLE}_${gene}_Illumina_Dragen_Filtered_SNVs_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf | cut -f 1,2,4,5 | sed 's/#//g' > \
		Illumina_data/${SAMPLE}_${gene}_Illumina_Dragen_genotypes_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.txt

		# Antal linjer/SNVs tilbage efter filtrering
		snvs=$(bcftools filter -i "FORMAT/DP>=10 && FORMAT/AF>=0.${AF_low} && FORMAT/AF<=0.${AF_high} && STRLEN(REF) == 1 && STRLEN(ALT) == 1 && N_ALT==1" Illumina_data/${SAMPLE}_${gene}_Illumina_Dragen_window${Window}_AF_0_${AF_low}_to_0_${AF_high}_DP${DP}.vcf | grep -v '^#' | wc -l)
		echo $snvs


		### SNV Bases
		# Skal chr, start og end initieres længere oppe?
		if [ $snvs -gt 0 ]; then
			python /work/sdularsen/tine/nanoporeData/kommandoer/SNV_Base_Identification2.py \
				--sample $SAMPLE --window $Window --gene $gene \
				--chr $chromosome --start $Start_Position --end $End_Position \
				--AF_low $AF_low --AF_high $AF_high --DP $DP
		fi
	done
done
