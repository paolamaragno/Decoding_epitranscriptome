This repository is a workflow for decoding the whole epitranscriptome and determining whether it changes during the RNA life cycle.
It is implemented to analyse data obtained through Nanopore direct RNA sequencing of RNAs present within chromatin, nucleoplasm and cytoplasm of SUM159 cells, an in vitro model of triple negative breast cancer cells, and K562 lymphoblast cells. 

Two different tools - ELIGOS and m6Anet - were used to identify RNA marks by analysing Nanopore reads. 
The comparative analysis against an IVT baseline in which RNAs were devoid of any mark allowed characterising the whole epitranscriptome. In particular, 
this analysis focuses on all the marks whose removal impacted Nanopore signal and could be called by ELIGOS. 
ELIGOS hits were characterised by motif discovery analysis and associated with specific RNA marks using public databases. 

In this way, it was possible to decipher which modifications were present on the RNAs and where, as well as if the levels of the different RNA 
modifications and their localization inside the gene body change across the different cellular fractions.

To better understand the spatial and temporal evolution of RNA modifications, the nano-ID tool was used to classify the reads into two classes: 
those representing newly-synthetized nascent RNAs and those from pre-existing RNAs. 

Finally, the combinatorial co-occurrence of various modifications within the transcripts of the same genes was studied. 

# Repository content
* eligos-v2.1.0: folder containing the Dockerfile to assemble all the images required by the pipeline, a bash script to run the nextflow pipeline (nf-eligos.sh), the nextflow pipeline main script (nf-eligos.nf), the nextflow pipeline configuration file (nf-eligos.conf) and a sample.txt;
* m6Anet-v2.1.0: folder containing the Dockerfile to assemble all the images required by the pipeline, a bash script to run the nextflow pipeline (nf-m6anet.sh), the nextflow pipeline main script (nf-m6anet.nf), the nextflow pipeline configuration file (nf-m6anet.conf) and a sample.txt;
* nanoID: folder containing a bash script to run the nextflow pipeline (nanoID.sh), the nextflow pipeline main script (nanoID.nf), the nextflow pipeline configuration file (nanoID.conf) and a sample.txt;
* R_scripts: folder containing a set of R scripts for data pre and post-processing.

# Getting started
* [Nextflow](https://nf-co.re/docs/usage/installation)
* [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

# Usage
## Subsampling reads obtained through Nanopore dRNA sequencing 
To perform library-level and gene-level subsampling on the reads obtained through Nanopore dRNA sequencing of the RNAs from chromatin, nucleoplasm and cytoplasm fractions, create an environment with the installation of minimap2 (v2.26-r1175), samtools (v1.17) and seqtk (v1.4-r122).
Then, activate R environment and run the following command:
```
Rscript subsampling_reads.R [list of arguments]

--path_fastq_chr_ass                                     Path to the fastq file of chromatin fraction
--path_fastq_nucleo                                      Path to the fastq file of nucleoplasm fraction
--path_fastq_cyto                                        Path to the fastq file of cytoplasm fraction
--path_reference_genome                                  Path to the reference genome in fasta format (*)
--path_gtf_file                                          Path to the gtf file (*)
--threads                                                Number of threads for the mapping (*)
--num_reads                                              Minimum number of reads mapping on each gene (*)
--num_subsampling                                        Number of samplings to do (*)
--minimap2                                               Path to minimap2 (*)
--samtools                                               Path to samtools (*)
--seqtk                                                  Path to seqtk (*)
--cond1                                                  Condition of the first fastq file given in input
--cond2                                                  Condition of the second fastq file given in input
--cond3                                                  Condition of the third fastq file given in input

For all the parameters indicated with (*) the default values can be set directly inside the R script.

Example of execution with default parameters:
Rscript subsampling_reads.R path_fastq_chr_ass='/path/to/chromatin_associated_PASS_reads.fastq' path_fastq_nucleo='/path/to/nucleoplasmic_PASS_reads.fastq' path_fastq_cyto='/path/to/cytoplasmic_PASS_reads.fastq' cond1='chr' cond2='nucleo' cond3='cyto'
```

## Running ELIGOS
Open nf-eligos.conf configuration file and set the desired options.

```
nextflow -c nf-eligos.conf run nf-eligos.nf --samples="/path/to/samples.txt" --resultsDir="/path/to/resultsDir" -profile docker

Mandatory argument:
-profile                                                 Configuration profile to use. Available: docker, singularity

Other mandatory arguments which may be specified in the nf-eligos.conf file
--samples                                                Path to the tab-separated sample file including sample name, condition and path to fastq file
--baseline_condition                                     Condition to be considered as the baseline, must match one of the conditions in the samples file
--min_depth                                              Minimum number of reads
--max_depth                                              Maximum number of reads
--spliced_alignment_flag                                 Flag for splice-aware alignment, set to true for genome alignment and to false for transcriptome alignment
--min_mapq                                               Minimum mapping quality for filtering alignments
--resultsDir                                             Path to a folder where to store results
--reference_fasta                                        Path to the reference fasta file
--bed_file                                               Path to regions of interest bed file
--pval_thr                                               p-value threshold
--adjPval_thr                                            adjusted p-value threshold
--oddR_thr                                               Odds Ratio threshold
--esb_thr                                                Threshold on %Error of Specific Bases to be considered for de novo motifs discovery
--sb                                                     Selected basis for filtering modification of interest
--opt_args                                               Other optional arguments (e.g. "-bcf <file.bcf> -m <model.json>")
```

Compile the sample.txt file and nf-eligos.sh bash file and execute ELIGOS with the command
```
qsub nf-eligos.sh
```

## Running m6Anet
Open nf-m6anet.conf configuration file and set the desired options.

```
nextflow -c nf-m6anet.conf run nf-m6anet.nf --samples="/path/to/samples.txt" --resultsDir="/path/to/resultsDir" -profile docker

Mandatory argument:
-profile                                                 Configuration profile to use. Available: docker, singularity

Other mandatory arguments which may be specified in the nf-m6anet.conf file

--samples                                                Path to the tab-separated sample file including sample name, condition, path to fast5 folder and path to fastq file
--resultsDir                                             Path to a folder where to store results
--transcriptome_fasta                                    Path to the transcriptome fasta file
--gtf                                                    Path to genome annotation gtf file
--min_mapq                                               Minimum mapping quality
--prob_mod_thr                                           Probability modification threshold for calling a site as m6A+
--postprocessingScript                                   Path to Transcript_to_genome.R script
--bulkLevelScript                                        Path to Calculate_m6anet_bulk.R script
```

Compile the sample.txt file and nf-m6anet.sh bash file and execute m6Anet with the command
```
qsub nf-m6anet.sh
```

## Running nano-ID
Open nanoID.conf configuration file and set the desired options.

```
nextflow -c nanoID.conf run nanoID.nf --samples="/path/to/samples.txt" --resultsDir="/path/to/resultsDir" -profile docker

Mandatory argument:
-profile                                                 Configuration profile to use. Available: docker, singularity

Other mandatory arguments which may be specified in the nf-m6anet.conf file

--samples                                                Path to the tab-separated sample file including sample name, condition, path to workspace directory
--scripts                                                Path to the pipeline scripts folder
--nanoIDSupplementals                                    Path to nanoID supplemental files
--unlabeled_time                                         If the training is required, this must match the conditions to use as unlabeled
--fullylabeled_time                                      If the training is required, this must match the conditions to use as fully labeled
--nanoIDtrainedModel                                     If the training has been already done this should point to the resulting model
--resultsDir                                             Path to a folder where to store results
--fast5_slot                                             FAST5 slot containing the basecalled bases
--qvalue                                                 Q value for reads filtering
--alignmentqValue                                        Alignment quality score for reads filtering
--genome_fasta                                           Path to the reference fasta file
--genesBed                                               Path to genes bed file
--CCL                                                    Cell Cycle Length in minutes
--Nseeds                                                 Number of seeds for training
--Nsamplings                                             Number of subsamplings; linear between 0 and Max
--trainingSizes                                          Sizes to add to the training set
--moveSlot                                               Set to "true" to use the info contained in the Move slot
```

Compile the sample.txt file and nanoID.sh bash file and execute nano-ID with the command
```
qsub nanoID.sh
```

## Analysis on random sequences
To generate random sequences execute
```
qsub generate_random_sequences.sh
```

To perform the overlap between the random sequences and the databases of RNA marks execute
```
qsub overlap_with_RNA_marks.sh
```

To perform the overlap between the random sequences and the databases of effectors' binding sites execute
```
qsub overlap_with_RBPs_chr.sh
qsub overlap_with_RBPs_nucleo.sh
qsub overlap_with_RBPs_cyto.sh
```
