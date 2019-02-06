# Description of RNA-Seq data analysis

The analysis of alternative splicing events is a multi-stage process. This document shall be used as
a description and guideline to understand the deposited analysis code and to advise in running a
comparable analysis on different data.

For each sample, multiple steps of the analysis are run within a module. After completion of
this first stage for all samples, the result files need to be integrated in a joint analysis. 

## Stage 1 - Per-sample analysis

This first stage has several components and is based on the following list of input files
* RNA-Seq sample in fastq format
* reference genome in fasta format
* reference genome annotation in GTF format

The code for all analysis steps of this first stage is located in
[sample_processing_rna](../sample_processing_rna). Specifically,
[sample_processing_rna/process_samples.sh](../sample_processing_rna/process_samples.sh) orchestrates
the iteration over a list of samples, e.g. using HPC infrastructure (in our example with LSF as a
batch system), by calling
[sample_processing_rna/process_sample_one.sh](../sample_processing_rna/process_sample_one.sh), that
provides the individual steps of the module:
* preprocessing of input file
* collection of QC statistics
* alignment
* gene expression quantification
* generation of alternative splicing graphs

In the following, we will give a brief description of each of the steps.

### Preprocessing of input files

This step can be omitted, if the RNA-Seq files are already in FASTQ format. As the files downloaded
from GDC were in BAM format, this step converts them back into FASTQ to be compatible with our
pipelines. 

First, the input BAM file is sorted by read ID and subsequently split into pairs. If this strategy
is used, it is important that unaligned reads are part of the BAM file and both sequence and quality
information are present. Sorting was done using `samtools`, while conversion into FASTQ was done
using a custom python script
[sample_processing_rna/bam2fastq.py](../sample_processing_rna/bam2fastq.py). (Please note that a
specific version for single-end reads exists. However, our analysis excludes single-end RNA-Seq
files.)

### Collection of QC statistics

Quality measures if the FASTQ input are collected using the FastQC tool (version 0.11.6). Call and
parameters are collected in the script
[sample_processing_rna/run_fastqc_one.sh](../sample_processing_rna/run_fastqc_one.sh). We generate one
output file per input FASTQ and later aggregate the information for joint QC analysis.

Later in the pipeline, right after the following alignment step, the number of aligned reads is
counted (together with some statistics on the distribution of mismatches). For this step, the script
[sample_processing_rna/collect_quick_align_stats.py](../sample_processing_rna/collect_quick_align_stats.py)
is used.

### Alignment

For alignment, we use the STAR software with a set of parameters we had good experience with in our
evaluations and practical applications over the past years. The full alignment call is present in
[sample_processing_rna/process_sample_one.sh](../sample_processing_rna/process_sample_one.sh) in the
section ALIGNMENT. It is important for the further steps in our analysis, that the SAM attributes as
specified in the call (NH HI NM MD AS XS) are set - they are needed by the downstream analysis
steps. STAR is run in 2-pass mode. That is, in a first round of alignment novel junctions are
detected from the RNA-Seq data and subsequently integrated into a new temporary index. This
temporary (sample specific) index is then used in a second alignment round for sensitive junction
detection. From the alignment results, we keep the BAM file as well as the list of junctions for further
analysis.

As a preliminary step to this alignment, the reference genome and the reference annotation need 
to be integrated into a STAR index that is used for the alignment step above. All logic and
parameters to create the alignment index are present in
[sample_processing_rna/run_create_index.sh](../sample_processing_rna/run_create_index.sh). The
directory containing the alignment index will is represented by `$genome` in
`process_sample_one.sh`.

*Practical Notes:* To run the STAR alignment, we used 80GByte of memory and 6-8 parallel threads on
our compute servers. Lowering the number of parallel threads, will increase running time.

### Gene expression quantification

For gene expression counting, we use the custom script
[count_expression/count_expression.py](../count_expression/count_expression.py), which will count all
non-secondary alignments that overlap to at least one annotated exon position. All positions in the
annotation that are assigned to more than one gene model are masked out. Further, we create two
version of the alignment counts. One with only overlapping gene positions masked and a second one
with also all positions masked that are annotated as both exon and intron in the annotation.

### Generation of alternative splicing graphs

For each sample, we generate a splicing graph using the SplAdder tool. This will be called in the
SPLICING section
[sample_processing_rna/process_sample_one.sh](../sample_processing_rna/process_sample_one.sh). Based on
the alignment of the sample in BAM format as input, this script will generate an initial splicing
graph for each gene in the annotation based on all annotated transcripts. Using the RNA-Seq evidence
from the BAM file, each gene's splicing graph is then augmented with additional nodes and edges. 

Noteworthy parameters in this context are the following:
* `-c 2` sets the SplAdder confidence level to 2, affecting read filter criteria such as minimum number of
spliced alignments across junctions or the minimum anchor length on each side of a spliced
alignment.
* `-M single` runs SplAdder in `single`-mode generating a splicing graph from a single sample
* `--sparse bam y` generates a compressed hdf5 representation of the bam file that is later used for
more efficient access to the alignment information when quantifying the splicing graphs.
* `-P y` exclusively use primary alignments
* `-n 50` the read length of the sample, which is used to correctly adapt filter criteria of the
chosen filter level

### Post-processing

In a last post-processing step, all result files generated in the temporary working space of the
script are gathered, appropriately renamed, and moved to the final output location. Temporary files
will be removed.

## Stage 2 - Data integration across all samples

The second stage is concerned with integrating the data collected from all individual samples and
preparing them for a joint analysis.

This includes the following steps:
* integration of alternative splicing events
* collection of quantification data
* aggregation of quality control information

### Integration of alternative splicing events

The script [alt_splice/events/run_alt_splice.sh](events/run_alt_splice.sh) takes care of
integrating the single splicing graphs generated during stage 1, to form a joint splicing graph per
gene across all samples.

The above script will call the SplAdder executable, taking the following inputs:
* annotation: reference gene annotation in GTF format that was used also for stage 1
* alignments: a text file containing the absolute paths to all alignment files in BAM format (If in
  stage 1 sparse representations of the alignment files were created, these files can be used
  instead.)
* output directory: the directory that is used as the base directory for SplAdder analyses (If OUT
  would be this output directory, it is necessary that the individual splicing graphs that were
  generated in stage 1 are either copied or sym-linked to a subdirectory OUT/spladder. This step is
  necessary, as the computation of SplAdder can be made more independent this way.)

The SplAdder pipeline will then automatically iterate over all samples and form a joint splicing
graph per gene. Subsequently, the joint splicing graph is quantified per gene and sample and
alternative splicing events are extracted and quantified. In the following, we give a list of the
most important parameters used in that script:
* `-M merge graphs` runs SplAdder in integrative mode, generating a merged splicing graph across 
  samples from the single, per-sample graphs
* `-P y` the SplAdder routine allows for distributed computing across an HPC system, which is 
directly managed within the script. However, this is experimental and currently only supports a
subset of HPC batch systems (LSF, torque, SGE, ...). Switching this off will lead to linear
processing of all steps without dispatching to HPC.
* `-t exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons` the list of alternative splicing
  event_types to be extracted and quantified
* `-D 500` empirical threshold for the number of edges in a splicing graph such that the gene is  
  considered for event detection to prevent exceedingly long computation times


