#!/usr/bin/env nextflow


//############################################################################################################################
//
// Josh Campbell
// 1/20/2017
// Peforms alignment and preprocessing of paired-end RNA sequencing data specifically for variant calling.
// This pipline may not be suitable for expression estimation as it marks duplicates
// For all samples derived from the same individual, an indel co-cleaning step will be performed on all bams jointly
//
// GATK tutorials this pipeline is based on:
// Mapping: https://software.broadinstitute.org/gatk/guide/article?id=3891
// Marking Duplicates: http://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
// Realignment around indels: http://gatkforums.broadinstitute.org/gatk/discussion/2800/howto-perform-local-realignment-around-indels
// Base Quality Score Recalibration: http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr
// 
//############################################################################################################################


// Set up global variables for requried parameters:
inputFile = file(params.infile)
inputFileHeader = params.infile_header

// Set up global variables for parameters with preset defaults:
REF = file(params.ref_dir)
REF_FASTA = file(params.ref_fasta)
GATK = file(params.gatk_jar)
PICARD = file(params.picard_jar)
GOLD1 = file(params.gold_indels1)
GOLD2 = file(params.gold_indels2)
OUTDIR = file(params.output_dir)
DBSNP = file(params.dbsnp)
READ_LENGTH = params.read_length
ADAPTER = params.adapter
GENE_GTF = file(params.gene_gtf)
GENE_BED = file(params.gene_bed)
RSEM_REF = file(params.rsem_ref)
OVERHANG = READ_LENGTH - 1
STRANDED = params.stranded
PREFIX = params.prefix
CREATE_SE = params.create_SE_rscript

RSEM_FORWARD_PROB = 0.5
if(STRANDED == true) {
  RSEM_FORWARD_PROB = 0  
} 



logParams(params, "nextflow_parameters.txt")

VERSION = "1.0"

// Header log info
log.info ""
log.info "========================================="
log.info "GATK Best Practices for RNA-Seq Preprocessing v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="
log.info ""


//#############################################################################################################
//#############################################################################################################
//
// Main
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Send FASTQ files to two processes from input file: FastQC and FastqToSam
//
// ------------------------------------------------------------------------------------------------------------

Channel.from(inputFile)
  .splitCsv(sep: '\t', header: inputFileHeader)
  .into { readPairsFastQC; readPairsFastqToSTAR_1Pass; readPairsFastqToSTAR_2Pass }




// ------------------------------------------------------------------------------------------------------------
//
// Run STAR 2-pass to align reads to genome
//
// ------------------------------------------------------------------------------------------------------------

process runSTAR_1pass {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/STAR_1Pass/"
    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastqToSTAR_1Pass
    
    output:
    file(outfile_sj) into runSTAR_1PassOutput
    
    script:
    outfile_prefix = sampleID + "_" + libraryID + "_" + rgID + ".1pass."
    outfile_bam = outfile_prefix + "Aligned.out.bam"        
    outfile_sj = outfile_prefix + "SJ.out.tab"        
    
    """
    module load star/2.5.2b
    
	STAR --genomeDir ${REF}					\
     --readFilesIn ${fastqR1} ${fastqR2}	\
     --runThreadN 12						\
     --outFileNamePrefix ${outfile_prefix}	\
     --outSAMtype BAM Unsorted 				\
     --clip3pAdapterSeq ${ADAPTER}			\
	 --outFilterMultimapNmax 20 			\
	 --outFilterType BySJout				\
     --readFilesCommand zcat
     
    rm -vf ${outfile_bam}
    """
}

process runSTAR_GenomeGenerate {
    tag "Generating STAR genome reference with Splice Junctions"
    
    input:
    val sjdb_files from runSTAR_1PassOutput.flatten().toSortedList()
    
    output:
	file('Genome') into runSTAR_GenomeGenerateOutput
	
    script:
    
    """
    module load star/2.5.2b

	STAR --runMode genomeGenerate					\
      --genomeDir ./								\
      --genomeFastaFiles ${REF_FASTA}				\
      --sjdbFileChrStartEnd ${sjdb_files.join(' ')}	\
      --sjdbGTFfile ${GENE_GTF}						\
      --sjdbOverhang ${OVERHANG}					\
      --runThreadN 12					
    """
}

process runSTAR_2pass {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/STAR_2Pass/"
    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastqToSTAR_2Pass
    val genomeFile from runSTAR_GenomeGenerateOutput.first()
    
    output:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, file(outfile_bam) into runSTAR_2PassOutput
	set indivID, sampleID, file(outfile_tbam) into runSTAR_2PassOutput_For_RSEM
	
    script:
    outfile_prefix = sampleID + "_" + libraryID + "_" + rgID + ".2pass."
    outfile_bam = outfile_prefix + "Aligned.out.bam"        
    outfile_tbam = outfile_prefix + "Aligned.toTranscriptome.out.bam"
    genomeDir = genomeFile.getParent()
    
    """
    module load star/2.5.2b

    STAR --genomeDir ${genomeDir}					\
	  --readFilesIn ${fastqR1} ${fastqR2} 			\
	  --runThreadN 12                           	\
      --outFileNamePrefix ${outfile_prefix}			\
      --outSAMtype BAM Unsorted						\
	  --quantMode TranscriptomeSAM					\
	  --outFilterMultimapNmax 20 					\
	  --outFilterType BySJout						\
	  --outSAMunmapped Within						\
      --readFilesCommand zcat      
    """
}


process runRSEM {
    tag "${indivID}|${sampleID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/RSEM"
    
    input:
	set indivID, sampleID, tbam from runSTAR_2PassOutput_For_RSEM
    
    output:
    file outfile_plot into runRSEMOutput
    file genes_file into genesFileForSE
    file isoforms_file into isoformsFileForSE 

    script:
    outfile_plot_prefix = sampleID + "_RSEM"
    outfile_plot = sampleID + "_RSEM.pdf"
    genes_file = sampleID + ".genes.results"
    isoforms_file = sampleID + ".isoforms.results"
    
    """
    module load rsem/1.3.0

	rsem-calculate-expression 						\
    --calc-ci --estimate-rspd --no-bam-output --bam \
	--paired-end 									\
	--forward-prob ${RSEM_FORWARD_PROB}				\
	-p 12											\
	$tbam 											\
	${RSEM_REF}					 					\
	${sampleID}

	rsem-plot-model ${sampleID} ${outfile_plot}
	"""
}


// ------------------------------------------------------------------------------------------------------------
//
// Add read group information and sort
//
// ------------------------------------------------------------------------------------------------------------


process runAddReadGroupInfo {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/STAR_2Pass/"
    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, bam from runSTAR_2PassOutput
	
	output:
	set indivID, sampleID, file(outfile_bam) into runAddReadGroupInfoOutput, runAddReadGroupInfoOutput_For_RSeQC
	
    script:
    outfile_bam = sampleID + "_" + libraryID + "_" + rgID + ".rgid.bam"        
    
    """
    module load java/1.8.0_66

	java -Xmx5G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=tmp/ -jar ${PICARD} AddOrReplaceReadGroups \
		I=${bam}		 			\
		O=${outfile_bam}	 		\
		SO=coordinate 				\
		RGID=${rgID}				\
		RGLB=${libraryID}			\
		RGPL=${platform}	 		\
		RGPU=${platform_unit} 		\
		RGSM=${sampleID}			\
		RGDT=${run_date}			\
		RGCN=${center}				\
		RGPM=${platform_model}		\
		CREATE_INDEX=true	
	
    samtools index ${outfile_bam}			
    """
}
     



// ------------------------------------------------------------------------------------------------------------
//
// Combined libraries from the same Sample to send to MarkDuplicates
//
// ------------------------------------------------------------------------------------------------------------

runAddReadGroupInfoOutput_grouped_by_sample = runAddReadGroupInfoOutput.groupTuple(by: [0,1])


// ------------------------------------------------------------------------------------------------------------
//
// Run Picard MarkDuplicates
// Requires a lot of memory
// Need to set "ParallelGCThreads" otherwise it will "grab" extra available threads without asking (and potentially be terminated by SGE)
//
// ------------------------------------------------------------------------------------------------------------

process runMarkDuplicates {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/MarkDuplicates"
	
    input:
    set indivID, sampleID, aligned_bam_list from runAddReadGroupInfoOutput_grouped_by_sample
    
    output:
    set indivID, sampleID, file(outfile_bam) into runMarkDuplicatesOutput
    file(outfile_bam) into runMarkDuplicatesOutput_for_RSeQC
    file (outfile_metrics) into runMarkDuplicatesOutput_for_MultiQC

    script:
    outfile_bam = sampleID + ".dedup.bam"
    outfile_metrics = sampleID + "_duplicate_metrics.txt"	
	        
    """
    module load java/1.8.0_66
    
	java -Xmx25G -XX:ParallelGCThreads=5 -Djava.io.tmpdir=tmp/ -jar ${PICARD} MarkDuplicates \
		INPUT=${aligned_bam_list.join(" INPUT=")} \
		OUTPUT=${outfile_bam} \
		METRICS_FILE=${outfile_metrics} \
		CREATE_INDEX=true \
		TMP_DIR=tmp
	"""  
}



// ------------------------------------------------------------------------------------------------------------
//
// Split reads aligning to introns into sepearate reads for downstream analysis
//
// ------------------------------------------------------------------------------------------------------------

process runSplitNCigarReads {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/SplitNCigarReads"
	
    input:
    set indivID, sampleID, bam from runMarkDuplicatesOutput
    
    output:
    set indivID, sampleID, file(outfile_bam) into runSplitNCigarReadsOutput
    
    script:
    outfile_bam = sampleID + ".splitNreads.bam"
	        
    """
    module load java/1.8.0_66
	java -Xmx15g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T SplitNCigarReads 			\
		-R ${REF_FASTA}					\
		-I ${bam} 						\
		-o ${outfile_bam} 				\
		-rf ReassignOneMappingQuality 	\
		-RMQF 255 						\
		-RMQT 60 						\
		-U ALLOW_N_CIGAR_READS    
	"""  
}




// ------------------------------------------------------------------------------------------------------------
//
// Combine samples from the same Individual (e.g. tumor/normal pair) to send to runRealignerTargetCreator
//
// ------------------------------------------------------------------------------------------------------------

runSplitNCigarReadsOutput_grouped_by_sample = runSplitNCigarReadsOutput.groupTuple(by: [0])



// ------------------------------------------------------------------------------------------------------------
//
// Perform realignment around indels
// 1) Identify regions for realignement
// 2) Perform realignment
//
// ------------------------------------------------------------------------------------------------------------

process runRealignerTargetCreator {
    tag "${indivID}"
    publishDir "${OUTDIR}/${indivID}/Processing/RealignerTargetCreator/"
    
    input:
    set indivID, sampleID, dedup_bam_list from runSplitNCigarReadsOutput_grouped_by_sample
    
    output:
    set indivID, dedup_bam_list, file(target_file) into runRealignerTargetCreatorOutput
 	
    script:
    target_file = indivID + "_target_intervals.list"
	        
    """
    module load java/1.8.0_66

	java -Xmx15g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T RealignerTargetCreator \
		-R ${REF_FASTA} \
		-I ${dedup_bam_list.join(" -I ")} \
		-known ${GOLD1} \
		-known ${GOLD2} \
		-o ${target_file}
	"""  
}

process runIndelRealigner {
    tag "${indivID}"
	publishDir "${OUTDIR}/${indivID}/Processing/IndelRealigner/"
	    
    input:
    set indivID, dedup_bam_list, target_file from runRealignerTargetCreatorOutput
 	    
    output:
    set indivID, file('*.realign.bam') into runIndelRealignerOutput mode flatten
 
    script:
            
    """
    module load java/1.8.0_66

	java -Xmx25g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T IndelRealigner \
		-R ${REF_FASTA} \
		-I ${dedup_bam_list.join(" -I ")} \
		-targetIntervals ${target_file} \
		-known ${GOLD1} \
		-known ${GOLD2} \
		-nWayOut ".realign.bam"		
	"""  
}



// ------------------------------------------------------------------------------------------------------------
//
// Perform base quality score recalibration (BQSR) including
// 1) Generate a recalibration table
// 2) Generate a new table after applying recalibration
// 3) Compare differences between recalibration tables
// 4) Apply recalibration
//
// ------------------------------------------------------------------------------------------------------------

// First we need to recapture the SampleID from the filename
runIndelRealignerOutput_split = runIndelRealignerOutput.map { indivID, file -> tuple(indivID, file.baseName.replaceAll(".splitNreads.realign", ""), file) }

process runBaseRecalibrator {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/BaseRecalibrator/"
	    
    input:
    set indivID, sampleID, realign_bam from runIndelRealignerOutput_split
    
    output:
    set indivID, sampleID, realign_bam, file(recal_table) into runBaseRecalibratorOutput
    
    script:
    recal_table = sampleID + "_recal_table.txt" 
       
    """
    module load java/1.8.0_66
    
	java -XX:ParallelGCThreads=2 -Xmx25g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T BaseRecalibrator \
		-R ${REF_FASTA} \
		-I ${realign_bam} \
		-knownSites ${GOLD1} \
		-knownSites ${GOLD2} \
		-knownSites ${DBSNP} \
		-o ${recal_table}
	"""
}

process runPrintReads {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/"
	    
    input:
    set indivID, sampleID, realign_bam, recal_table from runBaseRecalibratorOutput 

    output:
    set indivID, sampleID, file(outfile_bam), file(outfile_bai) into runPrintReadsOutput
    set indivID, sampleID, realign_bam, recal_table into runPrintReadsOutput_for_PostRecal
            
    script:
    outfile_bam = sampleID + ".clean.bam"
    outfile_bai = sampleID + ".clean.bai"
           
    """
    module load java/1.8.0_66

	java -XX:ParallelGCThreads=1 -Xmx25g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T PrintReads \
		-R ${REF_FASTA} \
		-I ${realign_bam} \
		-BQSR ${recal_table} \
		-o ${outfile_bam}
    """
}    

process runBaseRecalibratorPostRecal {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/BaseRecalibratorPostRecal/"
	    
    input:
    set indivID, sampleID, realign_bam, recal_table from runPrintReadsOutput_for_PostRecal
    
    output:
    set indivID, sampleID, recal_table, file(post_recal_table) into runBaseRecalibratorPostRecalOutput_Analyze
        
    script:
    post_recal_table = sampleID + "_post_recal_table.txt" 
       
    """
    module load java/1.8.0_66
    
	java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T BaseRecalibrator \
		-R ${REF_FASTA} \
		-I ${realign_bam} \
		-knownSites ${GOLD1} \
		-knownSites ${GOLD2} \
		-knownSites ${DBSNP} \
		-BQSR ${recal_table} \
		-o ${post_recal_table}
	"""
}	

process runAnalyzeCovariates {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/AnalyzeCovariates/"
	    
    input:
    set indivID, sampleID, recal_table, post_recal_table from runBaseRecalibratorPostRecalOutput_Analyze

	output:
	set indivID, sampleID, recal_plots into runAnalyzeCovariatesOutput
	    
    script:
    recal_plots = sampleID + "_recal_plots.pdf" 

    """
    module load java/1.8.0_66
    
	java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T AnalyzeCovariates \
		-R ${REF_FASTA} \
		-before ${recal_table} \
		-after ${post_recal_table} \
		-plots ${recal_plots}
    """
}    





// ------------------------------------------------------------------------------------------------------------
//
// Call Variants for each sample separately
// 1) Run HaplotypeCaller
// 2) Filter variants
// 3) Combine into single VCF file
//
// ------------------------------------------------------------------------------------------------------------

process runHTC {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/"
	    
    input:
    set indivID, sampleID, bam, bai from runPrintReadsOutput

	output:
	set indivID, sampleID, file(outfile) into runHTC_Output
	    
    script:
    outfile = sampleID + ".raw.vcf" 

    """
    module load java/1.8.0_66
    
	java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T HaplotypeCaller 			\
		-R ${REF_FASTA} 			\
		-I ${bam} 					\
		-dontUseSoftClippedBases 	\
		-stand_call_conf 20.0 		\
		--dbsnp ${DBSNP} 			\
		-L ${GENE_BED}				\
		-o ${outfile}
    """
}    

process runFilterVariants {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/"
	    
    input:
	set indivID, sampleID, vcf from runHTC_Output

	output:
	set file(outfile) into runFilterVariantsOutput
	    
    script:
    outfile = sampleID + ".filter.vcf" 

    """
    module load java/1.8.0_66
    
	java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T VariantFiltration				\
		-R ${REF_FASTA} 					\
		-V ${vcf} 							\
		-window 35 -cluster 3 				\
		-filterName FS -filter "FS > 30.0" 	\
		-filterName QD -filter "QD < 2.0"	\
		-o ${outfile}
    """
}    

process runCombineVariants {
    tag "Combining VCFs"
	publishDir "${OUTDIR}/Summary/Variants"
	    
    input:
	val vcf from runFilterVariantsOutput.flatten().toSortedList()

	output:
	set file(outfile) into runCombineVariantsOutput
	    
    script:
    outfile = PREFIX + ".filter.vcf" 

    """
    module load java/1.8.0_66
    
	java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T CombineVariants		\
		-R ${REF_FASTA} 		\
		-V ${vcf.join(" -V ")} 	\
		-o ${outfile}
    """
}    





// ------------------------------------------------------------------------------------------------------------
//
// Perform QC:
// 1) Run FASTQC to assess read quality
// 2) MultiQC on STAR 1st pass and 2nd pass output
//
// ------------------------------------------------------------------------------------------------------------



process runFastQC {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastQC/"
	    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastQC

    output:
    set file("*.zip"), file("*.html") into FastQCOutput, FastQCOutputForSE
   
    script:

    """
    module load fastqc/0.11.3
    fastqc -t 1 -o . ${fastqR1} ${fastqR2}
    """
}


process runRSeQC {
    tag "${indivID}|${sampleID}"

    publishDir "${OUTDIR}/${indivID}/${sampleID}/RSeQC/"

    input:
    set indivID, sampleID, bam from runAddReadGroupInfoOutput_For_RSeQC

    output:
    file("${sampleID}*.bam_stat.txt") into rseqc_bam_stat_results
    file("${sampleID}*.inferred_experiment.txt") into rseqc_inferred_experiment_results
    file("${sampleID}*.read_distribution.txt") into rseqc_read_distribution_results
    file("${sampleID}*.tin.txt") into rseqc_tin_results
    file("${sampleID}*.junction_annotation.txt") into rseqc_junction_annotation_results    
    file("${sampleID}*") into rseqc_results

    script:
    outfile1 = sampleID + ".bam_stat.txt"
    outfile2 = sampleID + ".inferred_experiment.txt"
    outfile3 = sampleID + ".read_distribution.txt"
    outfile4 = sampleID + ".tin.txt"
    outfile5 = sampleID + ".junction_annotation.txt"
	
    """
    module load python/2.7.12
    module load rseqc/2.6.4
    
    bam_stat.py -i ${bam} > ${outfile1}
    geneBody_coverage.py -i ${bam} -r ${GENE_BED} -o ${sampleID}
    junction_annotation.py -i ${bam} -r ${GENE_BED} -o ${sampleID} 2> ${outfile5}
    junction_saturation.py -i ${bam} -r ${GENE_BED} -o ${sampleID} 
    tin.py -i ${bam} -r ${GENE_BED} > ${outfile4}
    inner_distance.py -i ${bam} -r ${GENE_BED} -o ${sampleID} 
    clipping_profile.py -i ${bam} -s "PE" -o ${sampleID}
    infer_experiment.py -i ${bam} -r ${GENE_BED} > ${outfile2}
    insertion_profile.py -s "PE" -i ${bam} -o ${sampleID} 
    deletion_profile.py -i ${bam} -l ${READ_LENGTH} -o ${sampleID}
    read_distribution.py -i ${bam} -r ${GENE_BED} > ${outfile3}    
    read_GC.py -i ${bam} -o ${sampleID}
    read_duplication.py -i ${bam} -o ${sampleID}
	read_NVC.py -i ${bam} -o ${sampleID}
	read_quality.py -i ${bam} -o ${sampleID}
    """
}


// ------------------------------------------------------------------------------------------------------------
//
// Plot results with multiqc
//
// ------------------------------------------------------------------------------------------------------------

process runMultiQCFastq {
    tag "Generating fastq level summary and QC plots"
	publishDir "${OUTDIR}/Summary/Fastq"
	    
    input:
    val fastqc_files from FastQCOutput.flatten().toSortedList()
    
    output:
    file("fastq_multiqc*") into runMultiQCFastqOutput
    	
    script:

    """
    module load python/2.7.12
    module load multiqc/0.9

    multiqc -n fastq_multiqc ${fastqc_files.join(' ')}
    """
}



process runMultiQCLibrary {
    tag "Generating library level summary and QC plots"
	publishDir "${OUTDIR}/Summary/Library"
	    
    input:
    val duplicate_files from runMarkDuplicatesOutput_for_MultiQC.flatten().toSortedList()

    output:
    file("library_multiqc*") into runMultiQCLibraryOutput
    file("${sampleID}.*multiqc_picard_dups.txt") into runMultiQCPicardOutputForSE
    file("${sampleID}.*multiqc_general_stats.txt") into runMultiQCGeneralStatsForSE
 	
    script:
    outfile1 = PREFIX + "multiqc_picard_dups.txt"
    outfile2 = PREFIX + "multiqc_general_stats.txt"
    """
    module load python/2.7.12
    module load multiqc/0.9

    multiqc -n library_multiqc ${duplicate_files.join(' ')}
    """
}


process runMultiQCSample {
    tag "Generating sample level summary and QC plots"
	publishDir "${OUTDIR}/Summary/Sample"
	    
    input:
    val rseqc_files from rseqc_results.flatten().toSortedList()

    output:
    file("sample_multiqc*") into runMultiQCSampleOutput
    	
    script:
    """
    module load python/2.7.12
    module load multiqc/0.9

    multiqc -n sample_multiqc ${rseqc_files.join(' ')}
    """
}


// ------------------------------------------------------------------------------------------------------------
//
// Combine results into SummarizedExperiment object
//
// ------------------------------------------------------------------------------------------------------------
process runCreateSE {
    tag "Combining results into SummarizedExperiment object"
	publishDir "${OUTDIR}/Summary/Expression"
	    
    input:
    val rseqc_bam_stat_files from rseqc_bam_stat_results.flatten().toSortedList()
    val fastqc_files from FastQCOutputForSE.flatten().toSortedList()
    val rseqc_inferred_experiment_files from rseqc_inferred_experiment_results.flatten().toSortedList()
    val rseqc_read_distribution_files from rseqc_read_distribution_results.flatten().toSortedList()
    val rseqc_junction_annotation_files from rseqc_junction_annotation_results.flatten().toSortedList()
    val mark_duplicates_files from runMultiQCPicardOutputForSE.flatten().toSortedList()
    val multiqc_files from runMultiQCGeneralStatsForSE.flatten().toSortedList()
    val genes_files from genesFileForSE.flatten().toSortedList()
    val isoforms_files from isoformsFileForSE.flatten().toSortedList()

    echo "${rseqc_bam_stat_files.join('\n')}" > rseqc_bam_stat.txt
    echo "${fastqc_files.join('\n')}" > fastqc_files.txt
    echo "${rseqc_inferred_experiment_files.join('\n')}" > rseqc_inferred_experiment.txt
    echo "${rseqc_read_distribution_files.join('\n')}" > rseqc_read_distribution.txt
    echo "${rseqc_junction_annotation_files.join('\n')}" > rseqc_junction_annotation.txt
    echo "${mark_duplicates_files.join('\n')}" > mark_duplicates.txt
    echo "${multiqc_files.join('\n')}" > multiqc_files.txt
    echo "${genes_files.join('\n')}" > genes_results_files.txt
    echo "${isoforms_files.join('\n')}" > isoforms_results_files.txt

    output:
    set gene_file, iso_files into runCreateSEOutput
    	
    script:
    gene_file = PREFIX + "_Gene_Expression.rds"
    iso_file = PREFIX + "_Isoform_Expression.rds"
    
    """
	module load R/3.3.2
	
    ${CREATE_SE} -a ${genes_results_files.txt} -b ${isoforms_results_files.txt} -c ${demo.txt} -d ${inputFile} -e ${fastqc_files.txt} -f ${mark_duplicates.txt} -g ${rseqc_bam_stat.txt} -z ${multiqc_files.txt} -i ${rseqc_inferred_experiment.txt} -x ${rseqc_junction_annotation.txt} -k ${rseqc_read_distribution.txt} -n ${GENE_GTF}
    """
}


workflow.onComplete {
  log.info ""
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}




//#############################################################################################################
//#############################################################################################################
//
// FUNCTIONS
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------

def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}
