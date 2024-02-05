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
