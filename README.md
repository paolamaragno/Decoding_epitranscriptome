This repository is a workflow for decoding the whole epitranscriptome and determining whether it changes during the RNA life cycle.
It is implemented to analyse data obtained through Nanopore direct RNA sequencing of RNA present within chromatin, nucleoplasm and cytoplasm of SUM159 cells,
an in vitro model of triple negative breast cancer. 

Two different tools - ELIGOS and m6Anet - were used to identify RNA marks by analysing Nanopore reads. 
The comparative analysis against an IVT baseline in which RNAs were devoid of any mark allowed characterising the whole epitranscriptome. In particular, 
this analysis focuses on all the marks whose removal impacted Nanopore signal and could be called by ELIGOS. 
ELIGOS hits were characterised by motif discovery analysis and associated with specific RNA marks using public databases. 

In this way, it was possible to decipher which modifications were present on the RNAs and where, as well as if the levels of the different RNA 
modifications and their localization inside the gene body change across the different cellular compartments.

To better understand the spatial and temporal evolution of RNA modifications, the nano-ID tool was used to classify the reads into two classes: 
those representing newly-synthetized nascent RNAs and those from pre-existing RNAs. 

Finally, the combinatorial co-occurrence of various modifications within the transcripts of the same genes was studied. 

# Repository content
* eligos-v2.1.0: folder containing the Dockerfile to assemble all the images required by the pipeline, a bash script to run the nextflow pipeline (nf-eligos.sh),
  the nextflow pipeline main script (nf-eligos.nf) the nextflow pipeline configuration file (nf-eligos.conf) and a sample.txt;
* m6Anet-v2.1.0: folder containing the Dockerfile to assemble all the images required by the pipeline, a bash script to run the nextflow pipeline (nf-m6anet.sh),
  the nextflow pipeline main script (nf-m6anet.nf) and the nextflow pipeline configuration file (nf-m6anet.conf) and a sample.txt;
* nanoID: folder containing a bash script to run the nextflow pipeline (nanoID.sh), the nextflow pipeline main script (nanoID.nf) and the nextflow
  pipeline configuration file (nanoID.conf) and a sample.txt.
* R_scripts: folder containing a set of R scripts for data pre and post-processing.

# Getting started
* [Nextflow](https://nf-co.re/docs/usage/installation)
* [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

# Usage
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




