#!/bin/bash

# input a merged BAM from Veritas Genetics / PGP
#   TODO: add steps to produce merged bam

DIR_PICARD=/usr/local/share/picard/lib
JAR_PICARD=/usr/local/share/picard/lib/picard.jar
JAR_GATK=/usr/share/gatk/lib/GenomeAnalysisTK.jar
JAVA=java
SAMTOOLS=samtools
BWA=/usr/bin/bwa

TMP_DIR=/home/bjp/tmp
OUT_DIR=/home/bjp/genome/PGP/bam/gatk

HG38_FASTA=Homo_sapiens_assembly38.fa
HG38_DICT=Homo_sapiens_assembly38.dict

HS38DH_DICT=hs38DH.dict
HS38DH_FASTA=hs38DH.fa

HG38_KNOWN_INDELS=Homo_sapiens_assembly38.known_indels.vcf
MILLS_INDELS=Mills_and_1000G_gold_standard.indels.hg38.vcf

DBSNP_HG38=Homo_sapiens_assembly38.dbsnp.vcf

HAPMAP38_VCF=hapmap_3.3.hg38.vcf
OMNI_VCF=1000G_omni2.5.hg38.vcf
G1000SNPS_VCF=1000G_phase1.snps.high_confidence.hg38.vcf

BAMBASE=$1

BAM="${BAMBASE}.bam"

BAM_RG="${BAMBASE}.rg.bam"
BAM_RG_UM="${BAMBASE}.u.bam"

BAM_MARKED="${BAMBASE}.MarkIlluminaAdapters.bam"
BAM_MARK_METRICS="${BAMBASE}.MarkIlluminAadapters_metrics.txt"

BAM_CLEAN="${BAMBASE}_mergebamalignment.bam"
BAM_SORTED="${BAMBASE}_sorted.bam"
BAM_NMUQ="${BAMBASE}_nmuq.bam"

BAM_INDELREALIGN="${BAMBASE}_indels_realigned.bam"
BAM_REALIGN_INTERVALS="${BAMBASE}_realignertargetcreator.intervals"

BAM_MARK_DUPLICATES="${BAMBASE}_markduplicates.bam"
BAM_DUPLICATE_METRICS="${BAMBASE}.MarkDuplicates_metrics.txt"

BAM_BQSR="${BAMBASE}_bqsr.bam"
BAM_BQSR_BAI="${BAMBASE}_bqsr.bai"
BAM_BQSR_TABLE="${BAMBASE}_recal_data.table"
BAM_BQSR_PLOTS="${BAMBASE}_recalibration_plots.pdf"
BAM_POST_BQSR_TABLE="${BAMBASE}_post_recal_data.table"

BAM_COUNT_LOCI="${BAMBASE}.countloci.txt"

RAW_VCF="${BAMBASE}_raw_variants.vcf"
RAW_GVCF="${BAMBASE}_raw.g.vcf"

VQSR_SNP_RECAL="${BAMBASE}_recalibrate_SNP.recal"
VQSR_SNP_RECAL_TRANCHES="${BAMBASE}_recalibrate_SNP.tranches"
VQSR_SNP_RECAL_PLOTS="${BAMBASE}_recalibrate_SNP_plots.R"
VQSR_RECAL_SNPS_RAW_INDELS="${BAMBASE}_recalibrated_snps_raw_indels.vcf"

VQSR_INDEL_RECAL="${BAMBASE}_recalibrate_INDEL.recal"
VQSR_INDEL_RECAL_TRANCHES="${BAMBASE}_recalibrate_INDEL.tranches"
VQSR_INDEL_RECAL_PLOTS="${BAMBASE}_recalibrate_INDEL_plots.R"

VQSR_RECAL_VCF="${BAMBASE}_recalibrated_variants.vcf"

if [[ ! -f "${BAMBASE}.bam" ]]; then
	echo "Cannot find input BAM file $BAM"
	exit 1
fi


echo 
echo "Using input Veritas Genetics PGP BAM file $BAM"
echo 
echo "============================================================="

if [[ ! -r "$BAM_BQSR" ]]; then
	if [[ ! -r "$BAM_CLEAN" ]]; then
		################################
		if [[ ! -r "$BAM_RG" ]]; then
			echo 
			echo "Add read groups to BAM (in $BAM out $BAM_RG)"
			echo
			echo $JAVA -Xmx8G -jar $JAR_PICARD AddOrReplaceReadGroups I=$BAM O=$BAM_RG RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=bjp TMP_DIR=$TMP_DIR
			time $JAVA -Xmx8G -jar $JAR_PICARD AddOrReplaceReadGroups I=$BAM O=$BAM_RG RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=bjp TMP_DIR=$TMP_DIR 
			echo "============================================================="
		fi
		################################


		################################
		if [[ ! -r "$BAM_RG_UM" ]]; then
			echo 
			echo "Sanitize and unmap BAM (in $BAM_RG out $BAM_RG_UM)"
			echo
			echo $JAVA -Xmx8G -jar $JAR_PICARD RevertSam INPUT=$BAM_RG OUTPUT=$BAM_RG_UM SANITIZE=true MAX_DISCARD_FRACTION=0.005 ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XT ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true TMP_DIR=$TMP_DIR OUTPUT_BY_READGROUP=$OUT_DIR
			time $JAVA -Xmx8G -jar $JAR_PICARD RevertSam INPUT=$BAM_RG OUTPUT=$BAM_RG_UM SANITIZE=true MAX_DISCARD_FRACTION=0.005 ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XT ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true TMP_DIR=$TMP_DIR OUTPUT_BY_READGROUP=$OUT_DIR 
			echo "============================================================="
		fi
		################################

		################################
		if [[ ! -r $BAM_MARKED ]]; then
			echo 
			echo "Mark Illumina adapters with Picard (in $BAM_RG_UM out $BAM_MARKED)"
			echo
			echo $JAVA -Xmx8G -jar $JAR_PICARD MarkIlluminaAdapters I=$BAM_RG_UM O=$BAM_MARKED M=$BAM_MARK_METRICS TMP_DIR=$TMP_DIR
			time $JAVA -Xmx8G -jar $JAR_PICARD MarkIlluminaAdapters I=$BAM_RG_UM O=$BAM_MARKED M=$BAM_MARK_METRICS TMP_DIR=$TMP_DIR
			echo "============================================================="
		fi
		################################


		################################
		if [[ ! -r "$HG38_DICT" ]]; then
			echo 
			echo "Create $HG38_DICT (in $HG38_FASTA out ${HG38_DICT})"
			echo
			echo $JAVA -jar $JAR_PICARD CreateSequenceDictionary R=$HG38_FASTA O=$HG38_DICT
			time $JAVA -jar $JAR_PICARD CreateSequenceDictionary R=$HG38_FASTA O=$HG38_DICT
			echo "============================================================="
		fi
		################################


		################################
		if [[ ! -r "${HG38_FASTA}.fai" ]]; then
			echo 
			echo "Index $HG38_FASTA with $SAMTOOLS (in $HG38_FASTA out ${HG38_FASTA}.fai)"
			echo
			echo $SAMTOOLS faidx $HG38_FASTA
			time $SAMTOOLS faidx $HG38_FASTA 
			echo "============================================================="
		fi
		################################


		################################
		if [[ ! -r "${HG38_FASTA}.bwt" ]]; then
			echo 
			echo "Index reference for BWA (in $HG38_FASTA out $HG38_FASTA*)"
			echo
			echo $BWA index $HG38_FASTA
			time $BWA index $HG38_FASTA
			echo "============================================================="
		fi
		################################

		################################
		if [[ ! -x "./bwa.kit/run-gen-ref" ]]; then
			echo 
			echo "Fetching bwa.kit"
			echo
			echo wget -O- https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.12_x64-linux.tar.bz2/download \| bzip2 -dc \| tar xf -
			wget -O- https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.12_x64-linux.tar.bz2/download | bzip2 -dc | tar xf -
			echo "============================================================="
		fi
		################################


		################################
		if [[ ! -r "$HS38DH_FASTA" ]]; then
			echo 
			echo "Fetching hs38DH.fa (GRCh38+ALT+DECOY)"
			echo
			echo ./bwa.kit/run-gen-ref hs38DH
			time ./bwa.kit/run-gen-ref hs38DH
			echo "============================================================="
		fi
		################################

		################################
		if [[ ! -r "hs38DH.fa.sa" ]]; then
			echo 
			echo "BWA indexing hs38DH.fa (GRCh38+ALT+DECOY)"
			echo
			echo $BWA index hs38DH.fa
			time $BWA index hs38DH.fa
			echo "============================================================="
		fi
		################################

		################################
		if [[ ! -r "$HS38DH_DICT" ]]; then
			echo 
			echo "Create $HS38DH_DICT (in $HS38DH_FASTA out ${HS38DH_DICT})"
			echo
			echo $JAVA -jar $JAR_PICARD CreateSequenceDictionary R=$HS38DH_FASTA O=$HS38DH_DICT
			time $JAVA -jar $JAR_PICARD CreateSequenceDictionary R=$HS38DH_FASTA O=$HS38DH_DICT
		fi
		################################

		################################
		if [[ ! -r "${HS38DH_FASTA}.fai" ]]; then
			echo 
			echo "Index $HS38DH_FASTA with $SAMTOOLS (in $HS38DH_FASTA out ${HS38DH_FASTA}.fai)"
			echo
			echo $SAMTOOLS faidx $HS38DH_FASTA
			time $SAMTOOLS faidx $HS38DH_FASTA 
			echo "============================================================="
		fi
		################################


		################################
		if [[ ! -r "$BAM_CLEAN" ]]; then
			echo 
			echo "Convert BAM to Fastq |"
			echo "Align reads using BWA |"
			echo "Merge BAM files"
			echo "    (in $BAM_MARKED out $BAM_CLEAN)"
			echo
			echo $JAVA -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx128m -jar $JAR_PICARD SamToFastq I=$BAM_MARKED FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=$TMP_DIR \| $BWA mem -M -t 7 -p $HS38DH_FASTA /dev/stdin \| $JAVA -Dsamjdk.buffer_size=131072 -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:+UseStringCache -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx5000m -jar $JAR_PICARD MergeBamAlignment R=$HS38DH_FASTA UNMAPPED_BAM=$BAM_RG_UM ALIGNED_BAM=/dev/stdin OUTPUT=$BAM_CLEAN CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA UNMAP_CONTAMINANT_READS=true TMP_DIR=$TMP_DIR SORT_ORDER=unsorted
			$JAVA -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx128m -jar $JAR_PICARD SamToFastq I=$BAM_MARKED FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=$TMP_DIR | $BWA mem -M -t 7 -p $HS38DH_FASTA /dev/stdin | $JAVA -Dsamjdk.buffer_size=131072 -Dsamjdk.use_async_io=true -Dsamjdk.compression_level=1 -XX:+UseStringCache -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx5000m -jar $JAR_PICARD MergeBamAlignment R=$HS38DH_FASTA UNMAPPED_BAM=$BAM_RG_UM ALIGNED_BAM=/dev/stdin OUTPUT=$BAM_CLEAN CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA UNMAP_CONTAMINANT_READS=true TMP_DIR=$TMP_DIR SORT_ORDER=unsorted
			echo "============================================================="
		fi
		################################

		################################
		if [[ -r "$BAM_CLEAN" ]]; then
			echo 
			echo "Safe to remove unmapped BAM and MarkIlluminaAdapters output"
			echo
			echo rm $BAM_RG $BAM_RG_UM $BAM_MARKED
			rm $BAM_RG $BAM_RG_UM $BAM_MARKED
			echo "============================================================="
		fi
		################################
	fi


	################################
	if [[ ! -r "$BAM_MARK_DUPLICATES" ]]; then
		echo 
		echo "Mark duplicate reads in BAM (in $BAM_CLEAN out $BAM_MARK_DUPLICATES)"
		echo
		echo $JAVA -Xmx14G -jar $JAR_PICARD MarkDuplicates INPUT=$BAM_CLEAN OUTPUT=$BAM_MARK_DUPLICATES METRICS_FILE=$BAM_DUPLICATE_METRICS OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$TMP_DIR ASSUME_SORT_ORDER=queryname
		time $JAVA -Xmx14G -jar $JAR_PICARD MarkDuplicates INPUT=$BAM_CLEAN OUTPUT=$BAM_MARK_DUPLICATES METRICS_FILE=$BAM_DUPLICATE_METRICS OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$TMP_DIR ASSUME_SORT_ORDER=queryname
		echo "============================================================="
	fi
	################################

	ls -lh $BAM_CLEAN
	ls -lh $BAM_MARK_DUPLICATES
	ls -lh $BAM

	################################
	if [[ -r "$BAM_MARK_DUPLICATES" ]]; then
		echo 
		echo "Safe to remove clean BAM "
		echo
		echo rm $BAM_CLEAN
		rm $BAM_CLEAN
		echo "============================================================="
	fi
	################################

	################################
	if [[ ! -r "$BAM_NMUQ" ]]; then
		echo 
		echo "Coordinate sort, fix NM/UQ tags and index BAM (in $BAM_MARK_DUPLICATES out $BAM_NMUQ)"
		echo
		echo $JAVA -Xmx8G -jar $JAR_PICARD SortSam INPUT=$BAM_MARK_DUPLICATES OUTPUT=/dev/stdout SORT_ORDER=coordinate TMP_DIR=$TMP_DIR \| $JAVA -Xmx8G -jar $JAR_PICARD SetNmAndUqTags INPUT=/dev/stdin OUTPUT=$BAM_NMUQ CREATE_INDEX=true R=$HS38DH_FASTA TMP_DIR=$TMP_DIR
		time $JAVA -Xmx8G -jar $JAR_PICARD SortSam INPUT=$BAM_MARK_DUPLICATES OUTPUT=/dev/stdout SORT_ORDER=coordinate TMP_DIR=$TMP_DIR | $JAVA -Xmx8G -jar $JAR_PICARD SetNmAndUqTags INPUT=/dev/stdin OUTPUT=$BAM_NMUQ CREATE_INDEX=true R=$HS38DH_FASTA TMP_DIR=$TMP_DIR
		echo "============================================================="
	fi
	################################

	################################
	if [[ -r "$BAM_NMUQ" ]]; then
		echo 
		echo "Safe to remove MarkDuplicates BAM "
		echo
		echo rm $BAM_MARK_DUPLICATES
		rm $BAM_MARK_DUPLICATES
		echo "============================================================="
	fi
	################################
fi


################################
if [[ ! -r "$BAM_REALIGN_INTERVALS" ]]; then
	echo 
	echo "Local realignment around indels - target selection (in $BAM_NMUQ out $BAM_REALIGN_INTERVALS)"
	echo
	echo $JAVA -Xmx8G -jar $JAR_GATK -T RealignerTargetCreator -R $HS38DH_FASTA -I $BAM_NMUQ -o $BAM_REALIGN_INTERVALS -known $HG38_KNOWN_INDELS -known $MILLS_INDELS
	time $JAVA -Xmx8G -jar $JAR_GATK -T RealignerTargetCreator -R $HS38DH_FASTA -I $BAM_NMUQ -o $BAM_REALIGN_INTERVALS -known $HG38_KNOWN_INDELS -known $MILLS_INDELS
	echo "============================================================="
fi
################################

################################
if [[ ! -r "$BAM_INDELREALIGN" ]]; then
	echo 
	echo "Local realignment around indels - alignment (in $BAM_NMUQ out $BAM_INDELREALIGN in $BAM_REALIGN_INTERVALS)"
	echo
	echo $JAVA -Xmx8G -jar $JAR_GATK -T IndelRealigner -R $HS38DH_FASTA -targetIntervals $BAM_REALIGN_INTERVALS -known $HG38_KNOWN_INDELS -known $MILLS_INDELS -I $BAM_NMUQ -o $BAM_INDELREALIGN
	time $JAVA -Xmx8G -jar $JAR_GATK -T IndelRealigner -R $HS38DH_FASTA -targetIntervals $BAM_REALIGN_INTERVALS -known $HG38_KNOWN_INDELS -known $MILLS_INDELS -I $BAM_NMUQ -o $BAM_INDELREALIGN
	echo "============================================================="
fi
################################

################################
if [[ ! -r "$BAM_BQSR" ]]; then
	if [[ ! -r "$BAM_BQSR_TABLE" ]]; then
		echo 
		echo "BQSR: Step 1 (in $BAM_INDELREALIGN out $BAM_BQSR_TABLE)"
		echo "BQSR: Using $DBSNP_HG38 and $MILLS_INDELS"
		echo
		echo $JAVA -Xmx8G -jar $JAR_GATK -T BaseRecalibrator -R $HS38DH_FASTA -I $BAM_INDELREALIGN -knownSites $DBSNP_HG38 -knownSites $MILLS_INDELS -o $BAM_BQSR_TABLE
		time $JAVA -Xmx8G -jar $JAR_GATK -T BaseRecalibrator -R $HS38DH_FASTA -I $BAM_INDELREALIGN -knownSites $DBSNP_HG38 -knownSites $MILLS_INDELS -o $BAM_BQSR_TABLE
	fi

	if [[ ! -r "$BAM_POST_BQSR_TABLE" ]]; then
		echo 
		echo "BQSR: Step 2 (in $BAM_INDELREALIGN in $BAM_BQSR_TABLE out $BAM_POST_BQSR_TABLE)"
		echo "BQSR: Using $DBSNP_HG38 and $MILLS_INDELS"
		echo
		echo $JAVA -Xmx8G -jar $JAR_GATK -T BaseRecalibrator -R $HS38DH_FASTA -I $BAM_INDELREALIGN -knownSites $DBSNP_HG38 -knownSites $MILLS_INDELS -BQSR $BAM_BQSR_TABLE -o $BAM_POST_BQSR_TABLE
		time $JAVA -Xmx8G -jar $JAR_GATK -T BaseRecalibrator -R $HS38DH_FASTA -I $BAM_INDELREALIGN -knownSites $DBSNP_HG38 -knownSites $MILLS_INDELS -BQSR $BAM_BQSR_TABLE -o $BAM_POST_BQSR_TABLE
	fi

	if [[ ! -r "$BAM_BQSR_PLOTS" ]]; then
		echo 
		echo "BQSR: Step 3 (in $BAM_BQSR_TABLE in $BAM_POST_BQSR_TABLE out $BAM_BQSR_PLOTS)"
		echo
		echo $JAVA -Xmx8G -jar $JAR_GATK -T AnalyzeCovariates -R $HS38DH_FASTA -before $BAM_BQSR_TABLE -after $BAM_POST_BQSR_TABLE -plots $BAM_BQSR_PLOTS 
		time $JAVA -Xmx8G -jar $JAR_GATK -T AnalyzeCovariates -R $HS38DH_FASTA -before $BAM_BQSR_TABLE -after $BAM_POST_BQSR_TABLE -plots $BAM_BQSR_PLOTS
	fi

	echo 
	echo "BQSR: Step 4 (in $BAM_INDELREALIGN in $BAM_BQSR_TABLE out $BAM_BQSR)"
	echo 
	echo $JAVA -Xmx8G -jar $JAR_GATK -T PrintReads -R $HS38DH_FASTA -I $BAM_INDELREALIGN -BQSR $BAM_BQSR_TABLE -o $BAM_BQSR
	time $JAVA -Xmx8G -jar $JAR_GATK -T PrintReads -R $HS38DH_FASTA -I $BAM_INDELREALIGN -BQSR $BAM_BQSR_TABLE -o $BAM_BQSR
	echo "============================================================="
fi
################################


################################
if [[ -r "$BAM_BQSR" ]]; then
	echo 
	echo "Safe to remove NMUQ, IndelRealign BAMs"
	echo
	echo rm $BAM_NMUQ $BAM_INDELREALIGN
	rm $BAM_NMUQ $BAM_INDELREALIGN
	echo "============================================================="
fi
################################

################################
if [[ ! -r "$BAM_BQSR_BAI" ]]; then
	echo 
	echo "Indexing $BAM_BQSR (in $BAM_BQSR out $BAM_BQSR_BAI)"
	echo 
	echo $JAVA -Xmx8G -jar $JAR_PICARD BuildBamIndex I=$BAM_BQSR 
	time $JAVA -Xmx8G -jar $JAR_PICARD BuildBamIndex I=$BAM_BQSR 
	echo "============================================================="
fi
################################

################################
if [[ ! -r "$RAW_VCF" ]]; then
	echo 
	echo "Calling variants (in $BAM_BQSR out $RAW_VCF)"
	echo 
	echo $JAVA -Xmx8g -jar $JAR_GATK -T HaplotypeCaller -R $HS38DH_FASTA -I $BAM_BQSR --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $RAW_VCF
	time $JAVA -Xmx8g -jar $JAR_GATK -T HaplotypeCaller -R $HS38DH_FASTA -I $BAM_BQSR --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $RAW_VCF
fi
################################

################################
if [[ ! -r "$RAW_GVCF" ]]; then
	echo 
	echo "Calling variants (gVCF mode) (in $BAM_BQSR out $RAW_GVCF)"
	echo 
	echo $JAVA -Xmx8g -jar $JAR_GATK -T HaplotypeCaller -R $HS38DH_FASTA -I $BAM_BQSR --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -ERC GVCF -o $RAW_GVCF -variant_index_type LINEAR -variant_index_parameter 128000
	time $JAVA -Xmx8g -jar $JAR_GATK -T HaplotypeCaller -R $HS38DH_FASTA -I $BAM_BQSR --minPruning 3 --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -ERC GVCF -o $RAW_GVCF -variant_index_type LINEAR -variant_index_parameter 128000
fi
################################



exit 0

################################
echo 
echo "Counting reads in $BAM_BQSR"
echo
echo $JAVA -Xmx8G -jar $JAR_GATK -T CountReads -R $HS38DH_FASTA -I $BAM_BQSR 
time $JAVA -Xmx8G -jar $JAR_GATK -T CountReads -R $HS38DH_FASTA -I $BAM_BQSR 
echo "============================================================="
################################

################################
echo 
echo "Counting loci in $BAM_BQSR"
echo
echo $JAVA -Xmx8G -jar $JAR_GATK -T CountLoci -R $HS38DH_FASTA -I $BAM_BQSR 
time $JAVA -Xmx8G -jar $JAR_GATK -T CountLoci -R $HS38DH_FASTA -I $BAM_BQSR  -o $BAM_COUNTLOCI
echo "============================================================="
################################







exit 0


###################################
###if [[ ! -s "$VQSR_RECAL_VCF" ]]; then
###	if [[ ! -s "$VQSR_SNP_RECAL_PLOTS" ]]; then
###		echo 
###		echo "VQSR - SNPs build model (in $RAW_VCF out $VQSR_SNP_RECAL out $VQSR_SNP_RECAL_TRANCHES out $VQSR_SNP_RECAL_PLOTS)"
###		echo 
###		echo $JAVA -Xmx8g -jar $JAR_GATK -T VariantRecalibrator -R $HS38DH_FASTA -input $RAW_VCF -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP38_VCF -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_VCF -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000SNPS_VCF -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP_HG38 -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $VQSR_SNP_RECAL -tranchesFile $VQSR_SNP_RECAL_TRANCHES -rscriptFile $VQSR_SNP_RECAL_PLOTS
###		time $JAVA -Xmx8g -jar $JAR_GATK -T VariantRecalibrator -R $HS38DH_FASTA -input $RAW_VCF -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP38_VCF -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_VCF -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000SNPS_VCF -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP_HG38 -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $VQSR_SNP_RECAL -tranchesFile $VQSR_SNP_RECAL_TRANCHES -rscriptFile $VQSR_SNP_RECAL_PLOTS
###	fi
###	exit 0
###	if [[ ! -s "$VQSR_RECAL_SNPS_RAW_INDELS" ]]; then
###		echo 
###		echo "VQSR - SNPs apply recalibration (in $RAW_VCF out $VQSR_RECAL_SNPS_RAW_INDELS)"
###		echo 
###		echo $JAVA -Xmx8g -jar $JAR_GATK -T ApplyRecalibration -R $HS38DH_FASTA -input $RAW_VCF -mode SNP --ts_filter_level 99.9 -recalFile $VQSR_SNP_RECAL -tranchesFile $VQSR_SNP_RECAL_TRANCHES -o $VQSR_RECAL_SNPS_RAW_INDELS
###		time $JAVA -Xmx8g -jar $JAR_GATK -T ApplyRecalibration -R $HS38DH_FASTA -input $RAW_VCF -mode SNP --ts_filter_level 99.9 -recalFile $VQSR_SNP_RECAL -tranchesFile $VQSR_SNP_RECAL_TRANCHES -o $VQSR_RECAL_SNPS_RAW_INDELS
###	fi
###	if [[ ! -s "$VQSR_INDEL_RECAL_PLOTS" ]]; then
###		echo 
###		echo "VQSR - INDELs build model (in $VQSR_RECAL_SNPS_RAW_INDELS out $VQSR_INDEL_RECAL out $VQSR_INDEL_RECAL_TRANCHES out $VQSR_INDEL_RECAL_PLOTS)"
###		echo 
###		echo $JAVA -Xmx8g -jar $JAR_GATK -T VariantRecalibrator -R $HS38DH_FASTA -input $VQSR_RECAL_SNPS_RAW_INDELS -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS_INDELS -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 99.0 --maxGaussians 4 -recalFile $VQSR_INDEL_RECAL -tranchesFile $VQSR_INDEL_RECAL_TRANCHES -rscriptfilt $VQSR_INDEL_RECAL_PLOTS
###		time $JAVA -Xmx8g -jar $JAR_GATK -T VariantRecalibrator -R $HS38DH_FASTA -input $VQSR_RECAL_SNPS_RAW_INDELS -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS_INDELS -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 99.0 --maxGaussians 4 -recalFile $VQSR_INDEL_RECAL -tranchesFile $VQSR_INDEL_RECAL_TRANCHES -rscriptfilt $VQSR_INDEL_RECAL_PLOTS
###	fi
###	echo 
###	echo "VQSR - INDELs apply recalibration (in $VQSR_RECAL_SNPS_RAW_INDELS out $VQSR_RECAL_VCF)"
###	echo 
###	echo $JAVA -Xmx8g -jar $JAR_GATK -T ApplyRecalibration -R $HS38DH_FASTA -input $VQSR_RECAL_SNPS_RAW_INDELS -mode INDEL --ts_filter_level 99.9 -recalFile $VQSR_INDEL_RECAL -tranchesFile $VQSR_INDEL_RECAL_TRANCHES -o $VQSR_RECAL_VCF
###	time $JAVA -Xmx8g -jar $JAR_GATK -T ApplyRecalibration -R $HS38DH_FASTA -input $VQSR_RECAL_SNPS_RAW_INDELS -mode INDEL --ts_filter_level 99.9 -recalFile $VQSR_INDEL_RECAL -tranchesFile $VQSR_INDEL_RECAL_TRANCHES -o $VQSR_RECAL_VCF
###echo "============================================================="
###fi
###################################
###

