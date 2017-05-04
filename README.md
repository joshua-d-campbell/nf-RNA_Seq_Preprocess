# nf-RNA_Seq_Preprocess

## 1. Description
This pipeline will preprocess RNA-seq data to quantify gene and isoform expression with STAR/RSEM and perform variant calling based on the GATK best practices.
```
nextflow nf-RNA_Seq_Preprocess -c nextflow.config -with-timeline -with-dag -with-trace
```

If the pipeline fails at any point and you fix the issue, the pipeline can be restarted with job avoidance using the command:
```
nextflow nf-RNA_Seq_Preprocess -c nextflow.config -with-timeline -with-dag -with-trace -resume
```

For more information about Nextflow commands, please refer to the following link:

https://www.nextflow.io


## 2. Documentation of tools and approaches used for preprocessing
#### Preprocessing:
- Processing for Variant Calling: https://software.broadinstitute.org/gatk/guide/article?id=3891
- STAR: https://github.com/alexdobin/STAR
- Marking Duplicates: http://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
- Realignment around indels: http://gatkforums.broadinstitute.org/gatk/discussion/2800/howto-perform-local-realignment-around-indels
- Base Quality Score Recalibration: http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr

#### Qualiy control metrics:
- FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Picard: https://broadinstitute.github.io/picard/
- RSeQC: http://rseqc.sourceforge.net/

#### Plotting:
- MultiQC: [multiqc.info](multiqc.info) (currently using v0.9)

## 3. Preparing reference files before running workflow
See here for downloading reference and vcf files:

https://software.broadinstitute.org/gatk/guide/article?id=1213

This workflow assumes that the FASTA reference file has been preprocessed and indices have been made according to:

http://gatkforums.broadinstitute.org/wdl/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk

## 4. Input Parameters
```
params.infile = "fastq_input.txt"
params.output_dir = "output_directory"
params.prefix = "prefix"
params.demofile = "demographics.txt"
params.read_length = 75
params.stranded = true
params.ref_fasta = "hg19.fa"
params.ref_dir = "STAR_reference_directory"
params.gatk_jar = "/share/pkg/gatk/3.5/install/GenomeAnalysisTK.jar"
params.picard_jar = "/share/pkg/picard/2.8.0/install/lib/picard.jar"
params.gold_indels1 = "1000G_phase1.indels.hg19.sites.vcf"
params.gold_indels2 = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
params.dbsnp = "dbsnp_138.hg19.vcf"
params.infile_header = true
params.gene_gtf = "Homo_sapiens.GRCh37.75.ucsc.base_random.gtf"
params.gene_bed = "Homo_sapiens.GRCh37.75.ucsc.base_random.bed"
params.rsem_ref = "Homo_sapiens.GRCh37.75.ucsc.base_random"
params.create_SE_Rscript = "createSEfromRSEM.R"
params.inferAncestry = "inferAncestry.R"
```
#### Parameter description
  - params.infile - file path to input file
  - params.ouutput_dir - path to output directory
  - params.prefix - prefix to 
  - params.demofile - file path to demographics file
  - params.read_length - 
  - params.stranded - true will cause
  - params.ref_fasta - 
#### Input file description
The input file needs to be tab delimited and contain 11 columns (Capitalized):
  1.  INDIVIDUAL_ID - The ID of the individual from which the sample was derived.
  2.  SAMPLE_ID - The ID of the sample. Not that more than one sample can come from the same individual (e.g. tumor/normal pair)
  3.  LIBRARY_ID - The ID of the DNA library. Multiple sequencing libraries can be prepared from the same sample.
  4.  RG_ID - Read group ID
  5.  PLATFORM_UNIT - Generally is the read group ID plus the library ID
  6.  PLATFORM - Sequencer (e.g. illumina)
  7.  PLATFORM_MODEL - Sequencer (e.g. HiSeq2500)
  8.  RUN_DATE - Date of sequencing run
  9.  CENTER - Location of sequencing run
  10. R1 - Full path to Fastq file 1
  11. R2 - Full path to Fastq file 2
  
For more information on how to properly form read group IDs, please see the following:

https://software.broadinstitute.org/gatk/guide/article?id=6472


#### Target and bait information for the hybrid capture protocol:
targets: The genomic location of gene/exon boundaries that are supposed to be captured (in interval list format)

baits: The genomic location of the actual baits/probes used for capture (in interval list format)

#### JAR files
gatk_jar: JAR of GATK (tested with v3.5)

picard_jar: JAR of Picard Tools (tested with v2.8.0)

#### Reference files:
ref: BWA Reference file in FASTA format

gold_indels1 = High quality indel vcf file for realignment around indels

gold_indels2 = High quality indel vcf file for realignment around indels

dbsnp = dbSNP vcf used in base quality score recalibration (BQSR)

#### Other workflow parameters
output_dir: Final output directory for linked files

infile_header: Whether or not the input file has a header


## 5. Config file
The config file "nextflow.config" is included which contains all of the input paramters. To run on a cluster, you may need to change the "executor" and ".clusterOptions" for each subtask to work on your own system. If you want to change the number of cpus or memory requirements for a subtask, you will need to change the code in the main script as these requirements are currenly hard coded in the actual Linux command. To adapt NextFlow workflows to your own cluster setup, see the following link: 

https://www.nextflow.io/docs/latest/executor.html

## 6. Program versions and dependencies
This pipeline has been successfully run with the following versions
  - star v2.5.2b
  - rsem v1.3.0
  - samtools (requires java v1.8)
  - Picard tools v2.8.0 (requires java v1.8)
  - GATK v3.5
  - FastQC v0.11.3
  - rseqc v2.6.4 (using python 2.7.12)
  - multiqc v0.9 (using python 2.7.12)

**Important note:** These programs are currently loaded using the "module load" command. However, this will vary from system to system depending on your local setup. Therefore you may need to delete these commands and make sure these programs are accessible in your path.

## 7. File cleanup
This workflow does not currently delete the intermediate bams produced during the various steps. There after workflow completion, the follow commands will delete all intermediate ".bam" files.

```
rm -rf work/*/*/*.out.bam
rm -rf work/*/*/*.dedup.bam
rm -rf work/*/*/*splitNreads.bam
rm -rf work/*/*/*realign.bam
rm -rf work/*/*/*clean.bam

```





