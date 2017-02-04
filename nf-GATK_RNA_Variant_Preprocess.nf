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
OVERHANG = params.star_overhang
ADAPTER = params.adapter
GENE_GTF = file(params.gene_gtf)

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
     --sjdbGTFfile ${GENE_GTF}				\
     --readFilesCommand zcat
     
    rm -vf ${outfile_bam}
    """
}


process runSTAR_GenomeGenerate {
    tag "Generating STAR genome reference with Splice Junctions"
    
    input:
    val sjdb_files from runSTAR_1PassOutput.flatten().toList()
    
    output:
	file('Genome') into runSTAR_GenomeGenerateOutput
	
    script:
    
    """
    module load star/2.5.2b

	STAR --runMode genomeGenerate					\
      --genomeDir ./								\
      --genomeFastaFiles ${REF_FASTA}				\
      --sjdbFileChrStartEnd ${sjdb_files.join(' ')}	\
      --sjdbOverhang ${OVERHANG}					\
      --runThreadN 12							
    """
}


process runSTAR_2pass {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/STAR_2Pass/"
    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastqToSTAR_2Pass
    set genomeFile from runSTAR_GenomeGenerateOutput.first()
    
    output:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, file(outfile_bam) into runSTAR_2PassOutput
	
    script:
    outfile_prefix = sampleID + "_" + libraryID + "_" + rgID + ".2pass."
    outfile_bam = outfile_prefix + "Aligned.out.bam"        
    genomeDir = genomeFile.getParent()
    
    """
    module load star/2.5.2b

    STAR --genomeDir ${genomeDir}				\
	  --readFilesIn ${fastqR1} ${fastqR2} 		\
	  --runThreadN 12                           \
      --outFileNamePrefix ${outfile_prefix}		\
      --outSAMtype BAM Unsorted					\
      --clip3pAdapterSeq ${ADAPTER}				\
      --readFilesCommand zcat
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
	set indivID, sampleID, file(outfile_bam) into runAddReadGroupInfoOutput
	
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
		RGPM=${platform_model}		
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
    file(outfile_metrics) into runMarkDuplicatesOutput_for_MultiQC
    
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
runIndelRealignerOutput_split = runIndelRealignerOutput.map { indivID, file -> tuple(indivID, file.baseName.replaceAll(".dedup.realign", ""), file) }

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
    set file("*.zip"), file("*.html") into FastQCOutput
    	
    script:

    """
    module load fastqc/0.11.3
    fastqc -t 1 -o . ${fastqR1} ${fastqR2}
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
    val fastqc_files from FastQCOutput.flatten().toList()
    
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
    val duplicate_files from runMarkDuplicatesOutput_for_MultiQC.flatten().toList()

    output:
    file("library_multiqc*") into runMultiQCLibraryOutput
    	
    script:
    """
    module load python/2.7.12
    module load multiqc/0.9

    multiqc -n library_multiqc ${duplicate_files.join(' ')}
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
