GEM mapping
===========

# Generate genome index

gemtools index -i Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa -t 12

# Generate transcriptome index from genome annotation

gemtools t-index -a Taeniopygia_guttata.taeGut3.2.4.89.gtf -i Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.gem -t 12

# Start gemtools RNA-seq pipeline

gemtools rna-pipeline -i Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.gem -a Taeniopygia_guttata.taeGut3.2.4.89.gtf -f $FASTQ_R1 -q 33 -t 12

Read quantification
===================

featureCounts -p -t exon -g gene_id -a Taeniopygia_guttata.taeGut3.2.4.89.gtf -o featureCount *.bam -T30


Differential Regulation
=======================

See Differential_Regulation_edgeR.ipynb


GO enrichment
=============

See GO_enrichment_analysis.ipynb


Variant Calling
===============

STAR alignment
--------------

# Generate index

STAR --runMode genomeGenerate --genomeDir $INDEX_DIR --genomeFastaFiles $GENOME --sjdbGTFfile $GTF --runThreadN 12

# First pass mapping

STAR --genomeDir $INDEX_DIR --sjdbGTFfile $GTF --readFilesIn $FQ_R1 $FQ_R2 --outSAMtype None --runThreadN 12 --outFileNamePrefix $prefix_out --readFilesCommand zcat

# Second pass mapping

STAR --genomeDir $INDEX_DIR --sjdbGTFfile $GTF --readFilesIn $i ${i%_1.fastq.gz}_2.fastq.gz --runThreadN 12 --outFileNamePrefix $prefix_out --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --sjdbFileChrStartEnd ${FIRST_PASS_MAPPING_DIR}/*_SJ.out.tab

Bam post-processing
-------------------

# Adding read groups:

AddOrReplaceReadGroups I=$BAM O=${BAM%%.*}_rgid_sort.bam RGID=juncoOnlyLib RGLB=library RGPL=illumina RGPU=machine RGSM=${ID} SO=coordinate

# Mark duplicates:

MarkDuplicates I=${BAM%%.*}_rgid_sort.bam O=${BAM%%.*}_rgid_sort_DupMark.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${BAM%%.*}_output.metrics

# Split cigar

GenomeAnalysisTK -T SplitNCigarReads -R $FAS -I ${BAM%%.*}_rgid_sort_DupMark.bam -o ${BAM%.bam}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

# Merging bams:

MergeSamFiles O=$MERGED_BAM USE_THREADING=true VALIDATION_STRINGENCY=LENIENT I=$BAM1 I=$BAM2 I=...

SNP discovery
-------------

# First round of SNP discovery

GenomeAnalysisTK -T HaplotypeCaller -R $FAS -I all_bams_split_rg_merged.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $ROUND1_VCF -nct 5

# Extract indels for second round of SNP discovery

vcftools --vcf $ROUND1_VCF --keep-only-indels --recode --recode-INFO-all --out $INDELS_FIRST_ROUND

# Create intervals

GenomeAnalysisTK \
   -T RealignerTargetCreator \
   -R $FAS \
   -I $MERGED_BAM \
   -known $INDELS_FIRST_ROUND \
   -o forIndelRealigner.intervals
    
# Realignment around indels

GenomeAnalysisTK \
    -T IndelRealigner \
    -R $FAS \
    -I $MERGED_BAM \
    -known $INDELS_FIRST_ROUND \
    -targetIntervals forIndelRealigner.intervals \
    -o merged_split_realign.bam

# Base recalibration

GenomeAnalysisTK \
    -T BaseRecalibrator \
    -R $FAS \
    -I merged_split_realign.bam \
    -knownSites round1_raw_filtered.vcf \
    -o recal_data.table

# Exporting realigned reads

GenomeAnalysisTK \
    -T PrintReads \
    -R $FAS \
    -I merged_split_realign.bam \
    -BQSR recal_data.table \
    -o merged_split_realign_recal.bam

# Calling SNPs

GenomeAnalysisTK -T HaplotypeCaller -R $FAS -I merged_split_realign_recal.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o final.vcf -nct 5


SNP filtering
-------------

GenomeAnalysisTK -T VariantFiltration -R $FAS -V final.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o final_gatkFilt.vcf

vcftools --vcf final_gatkFilt.vcf --recode --out final_gatkFilt_vcfFilt --minDP 6 --minGQ 10

vcftools --vcf final_gatkFilt_vcfFilt.recode.vcf --remove-filtered-all --recode --stdout > final_gatkFilt_vcfFilt_onlyPASS.vcf


FST calculation
---------------

MAXMIS=0.25
MAC=1
MINAL=2

# Filter and calculate FSTs

vcftools --max-missing $MAXMIS --mac $MAC --min-alleles $MINAL --vcf $VCF --recode --stdout > $OUT_PREFIX.vcf

vcftools --vcf $OUT_PREFIX.vcf --weir-fst-pop $POP1 --weir-fst-pop $POP2 --stdout > ${OUT_PREFIX}_FST.tab

# Using in-house python script to localize SNPs

SNPs_annotator.py $GTF $OUT_PREFIX.vcf $MART_GO annotated_SNPs 0 ${OUT_PREFIX}_FST.tab
