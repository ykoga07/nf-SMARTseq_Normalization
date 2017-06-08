# nf-SMARTseq_Normalization

## 1. Description
This pipeline will align and normalize SMART-seq data to quantify gene and isoform expression with STAR/RSEM and perform variant calling based on the GATK best practices.
```
nextflow nf-SMART_Seq_Normalization -c nextflow.config -with-timeline -with-dag -with-trace
```

If the pipeline fails at any point and you fix the issue, the pipeline can be restarted with job avoidance using the command:
```
nextflow nf-SMART_Seq_Normalization -c nextflow.config -with-timeline -with-dag -with-trace -resume
```

For more information about Nextflow commands, please refer to the following link:

https://www.nextflow.io


## 2. Documentation of tools and approaches used for preprocessing
#### Alignment:
- STAR: https://github.com/alexdobin/STAR

#### Qualiy control metrics:
- FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- RSeQC: http://rseqc.sourceforge.net/

#### Plotting:
- MultiQC: [multiqc.info](multiqc.info) (currently using v0.9)

#### Quantification:
- RSEM: https://deweylab.github.io/RSEM/

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
params.paired_end = false

params.ref_fasta = "hg19.fa"
params.ref_dir = "STAR_reference_directory"
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

#### Input file description
The input file needs to be tab delimited and contain 11 columns (Capitalized):
  1.  INDIVIDUAL_ID - The ID of the individual from which the sample was derived.
  2.  SAMPLE_ID - The ID of the sample. More than one sample can come from the same individual (e.g. tumor/normal pair)
  3.  LIBRARY_ID - The ID of the DNA library. Multiple sequencing libraries can be prepared from the same sample.
  4.  RG_ID - Read group ID
  5.  PLATFORM_UNIT - Generally is the read group ID plus the library ID
  6.  PLATFORM - Sequencer (e.g. illumina)
  7.  PLATFORM_MODEL - Sequencer (e.g. HiSeq2500)
  8.  RUN_DATE - Date of sequencing run
  9.  CENTER - Location of sequencing run
  10. R1 - Full path to Fastq file 1
  11. R2 - Full path to Fastq file 2, if single-end data, will be blank or NA

For more information on how to properly form read group IDs, please see the following:

https://software.broadinstitute.org/gatk/guide/article?id=6472

#### Main workflow parameters
infile: Input file

output_dir: Final output directory for linked files

prefix: Prefix to the output files

demofile: Demographics file, optional, first column must match SAMPLE_ID in the input file

read_length: Length of largest read

paired_end: true for paired-end sequencing data, false for single-end sequencing data

#### Reference/JAR files:
ref_fasta: BWA Reference file in FASTA format

ref_dir: Directory of the STAR Genome Directory

picard_jar: JAR of Picard Tools (tested with v2.8.0)

dbsnp: dbSNP vcf used in base quality score recalibration (BQSR)

infile_header: Whether or not the input file has a header

gene_gtf: Gene annotations in gtf format (used for STAR genome generation)

gene_bed: Reference gene model in bed format (used for HTC, RSeQC)

rsem_ref: Reference created by RSEM, used for quantifying gene/isoform expression 

#### R Scripts
create_SE_Rscript: Full file path to createSEfromRSEM.R script

## 5. Config file
The config file "nextflow.config" is included which contains all of the input paramters. To run on a cluster, you may need to change the "executor" and ".clusterOptions" for each subtask to work on your own system. If you want to change the number of cpus or memory requirements for a subtask, you will need to change the code in the main script as these requirements are currenly hard coded in the actual Linux command. To adapt NextFlow workflows to your own cluster setup, see the following link: 

https://www.nextflow.io/docs/latest/executor.html

## 6. Program versions and dependencies
This pipeline has been successfully run with the following versions
  - star v2.5.2b
  - rsem v1.3.0
  - samtools v1.4 (requires java v1.8)
  - FastQC v0.11.3
  - rseqc v2.6.4 (using python 2.7.12)
  - multiqc v0.9 (using python 2.7.12)

**Important note:** These programs are currently loaded using the "module load" command. However, this will vary from system to system depending on your local setup. Therefore you may need to delete these commands and make sure these programs are accessible in your path.

## 7. Making MultiQC plots across different runs
Due to the nature of the pipeline, there may not be enough storage space to run all of the samples at once. To Make MultiQC plots across different runs, use the following commands:

In this example, there are two batches, Batch 1 and Batch2. To make MultiQC plots from the outputs of both of the batches, we will concatenate the file paths that were used in the batches, stored inside each of the "multiqc_input_files.txt".  
```
module load python/2.7.12
module load multiqc/0.9
#For Fastq MultiQC
echo -e "Batch1/Output/QC/Fastq/fastq_multiqc_input_files.txt" >  all_fastq_multiqc_input_files.txt
echo -e "Batch2/Output/QC/Fastq/fastq_multiqc_input_files.txt" >>  all_fastq_multiqc_input_files.txt
multiqc -n all_fastq_multiqc --file-list all_fastq_multiqc_input_files.txt
```
```
#For library MultiQC
echo -e "Batch1/Output/QC/Library/library_multiqc_input_files.txt" >  all_library_multiqc_input_files.txt
echo -e "Batch2/Output/QC/Library/library_multiqc_input_files.txt" >>  all_library_multiqc_input_files.txt
multiqc -n all_library_multiqc --file-list all_library_multiqc_input_files.txt
```
```
#For sample MultiQC
echo -e "Batch1/Output/QC/Sample/sample_multiqc_input_files.txt" >  all_sample_multiqc_input_files.txt
echo -e "Batch2/Output/QC/Sample/sample_multiqc_input_files.txt" >>  all_sample_multiqc_input_files.txt
multiqc -n all_sample_multiqc --file-list all_sample_multiqc_input_files.txt
```


## 8. File cleanup
This workflow does not currently delete the intermediate bams produced during the various steps. There after workflow completion, the follow commands will delete all intermediate ".bam" files.

```
rm -rf work/*/*/*.out.bam
rm -rf work/*/*/*.dedup.bam
rm -rf work/*/*/*splitNreads.bam
rm -rf work/*/*/*realign.bam
rm -rf work/*/*/*clean.bam

```





