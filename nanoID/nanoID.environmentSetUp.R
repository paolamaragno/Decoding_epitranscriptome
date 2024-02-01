### This code is aimed at setting-up the R environment for nano-ID analyses ###

## Package to install specific versions
if(!("versions" %in% rownames(installed.packages())))install.packages("versions")
library("versions")

if(!("devtools" %in% rownames(installed.packages())))install.packages("devtools")
library("devtools")

## Bioconductor packages installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioconductorPackages <- c("rhdf5","GenomicRanges","GenomeInfoDb","SummarizedExperiment","GenomicAlignments","Rsamtools"
						 ,"Biobase","Biostrings","XVector","IRanges","S4Vectors","BiocGenerics","zlibbioc","BiocParallel","rtracklayer")

for(p in bioconductorPackages)
{
	if(!(p %in% rownames(installed.packages())))
	{
		print(p)
		BiocManager::install(p,suppressUpdate=TRUE,suppressAutoUpdate=TRUE)
	}
}

for(p in bioconductorPackages)
{
	library(p,verbose=FALSE,quietly=TRUE,character.only = TRUE)
}

## CRAN packages installation
cranPackages <- c("1.18-1","4.6-14","0.3.0","2.2.1","0.5.0","1.7-0","0.4-20","2.2-5.4","1.0.10","1.0.8","1.4.3","0.12.10","1.8.4","0.1-3","1.0-6","7.3-14","0.1.0","0.2.0","1.3.0","1.4","0.20-35","0.3.1","1.2-8","0.6-1","0.0.2","2.2.0","3.98-1.6","0.3-2","1.5","0.4.1","1.4","0.2-15","0.1","1.3-2","0.2.0","0.4.3","1.95-4.8","0.2")
names(cranPackages) <- c("dtw","randomForest","purrr","ggplot2","dplyr","e1071","proxy","sm","doParallel","iterators","foreach","Rcpp","plyr","base64enc","bitops","class","zeallot","gtable","tibble","jsonlite","lattice","rlang","Matrix","DBI","generics","R6","XML","whisker","magrittr","scales","tfruns","codetools","assertthat","colorspace","lazyeval","munsell","RCurl","vioplot")

for(p in seq_along(cranPackages))
{
	if(!(names(cranPackages)[[p]] %in% rownames(installed.packages())))
	{
		print(names(cranPackages)[[p]])
		# It is possible that the required version is not available
		# if so I install the available one
		try(install.versions(names(cranPackages)[[p]],cranPackages[[p]]))
		if(!(names(cranPackages)[[p]] %in% rownames(installed.packages())))
		{
			print("No version compatibility")
			install.packages(names(cranPackages)[[p]])
		}
	}
}

for(p in names(cranPackages))
{
	library(p,verbose=FALSE,quietly=TRUE,character.only = TRUE)
}

library(reticulate)
library(tensorflow)
library(keras)

# install.packages("largeList")
library("largeList")

# install.packages("data.table")
library("data.table")

library(parallel)
mcsapply <- function( X, FUN, ... ) do.call('cbind', mclapply( X, FUN, ... ))

# sessionInfo()

# print(paste0("In this section there are ",mc.cores," available cores."))

# R version 3.6.1 (2019-07-05)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)

# Matrix products: default
# BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.19.so

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
# [1] RCurl_1.98-1.2              munsell_0.5.0              
#  [3] lazyeval_0.2.2              colorspace_1.4-1           
#  [5] assertthat_0.2.1            codetools_0.2-16           
#  [7] tfruns_1.4                  scales_1.1.0               
#  [9] magrittr_1.5                whisker_0.4                
# [11] XML_3.98-1.20               R6_2.4.1                   
# [13] generics_0.0.2              DBI_1.0.0                  
# [15] Matrix_1.2-17               rlang_0.4.2                
# [17] lattice_0.20-38             jsonlite_1.6               
# [19] tibble_2.1.3                gtable_0.3.0               
# [21] zeallot_0.1.0               class_7.3-15               
# [23] bitops_1.0-6                base64enc_0.1-3            
# [25] plyr_1.8.4                  Rcpp_1.0.3                 
# [27] doParallel_1.0.10           iterators_1.0.8            
# [29] foreach_1.4.3               e1071_1.7-0                
# [31] dplyr_0.8.3                 ggplot2_3.2.1              
# [33] purrr_0.3.3                 randomForest_4.6-14        
# [35] vioplot_0.3.4               zoo_1.8-6                  
# [37] sm_2.2-5.4                  dtw_1.18-1                 
# [39] proxy_0.4-16                rtracklayer_1.46.0         
# [41] zlibbioc_1.32.0             GenomicAlignments_1.22.1   
# [43] Rsamtools_2.2.3             Biostrings_2.54.0          
# [45] XVector_0.26.0              SummarizedExperiment_1.16.1
# [47] DelayedArray_0.12.3         BiocParallel_1.20.1        
# [49] matrixStats_0.57.0          Biobase_2.46.0             
# [51] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
# [53] IRanges_2.20.2              S4Vectors_0.24.4           
# [55] BiocGenerics_0.32.0         rhdf5_2.30.1               
# [57] devtools_2.2.1              usethis_1.5.1              
# [59] versions_0.3               

# loaded via a namespace (and not attached):
#  [1] pkgload_1.0.2          BiocManager_1.30.10    GenomeInfoDbData_1.2.2
#  [4] remotes_2.1.0          sessioninfo_1.1.1      pillar_1.4.2          
#  [7] backports_1.1.5        glue_1.3.1             digest_0.6.23         
# [10] pkgconfig_2.0.3        processx_3.4.1         ellipsis_0.3.0        
# [13] withr_2.1.2            cli_1.1.0              crayon_1.3.4          
# [16] memoise_1.1.0          ps_1.3.0               fs_1.3.1              
# [19] pkgbuild_1.0.6         tools_3.6.1            prettyunits_1.0.2     
# [22] lifecycle_0.1.0        Rhdf5lib_1.8.0         callr_3.3.2           
# [25] compiler_3.6.1         grid_3.6.1             testthat_2.3.1        
# [28] rprojroot_1.3-2        desc_1.2.0             tidyselect_0.2.5           

## Expected session info
# R version 3.3.3 (2017-03-06)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Scientific Linux 7.3 (Nitrogen)

# locale:
#  [1] LC_CTYPE=C                 LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
#  [1] purrr_0.3.0                ggplot2_2.2.1             
#  [3] dplyr_0.5.0                keras_2.2.4               
#  [5] randomForest_4.6-14        e1071_1.7-0               
#  [7] dtw_1.18-1                 proxy_0.4-20              
#  [9] rhdf5_2.18.0               vioplot_0.2               
# [11] sm_2.2-5.4                 rtracklayer_1.34.2        
# [13] GenomicAlignments_1.10.1   Rsamtools_1.26.2          
# [15] SummarizedExperiment_1.4.0 Biobase_2.34.0            
# [17] GenomicRanges_1.26.4       GenomeInfoDb_1.10.3       
# [19] doParallel_1.0.10          iterators_1.0.8           
# [21] foreach_1.4.3              Biostrings_2.42.1         
# [23] XVector_0.14.1             IRanges_2.8.2             
# [25] S4Vectors_0.12.2           BiocGenerics_0.20.0       

# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.10       plyr_1.8.4         base64enc_0.1-3    bitops_1.0-6      
#  [5] class_7.3-14       tools_3.3.3        zlibbioc_1.20.0    zeallot_0.1.0     
#  [9] gtable_0.2.0       tibble_1.3.0       jsonlite_1.4       lattice_0.20-35   
# [13] rlang_0.3.1        Matrix_1.2-8       DBI_0.6-1          generics_0.0.2    
# [17] grid_3.3.3         reticulate_1.10    R6_2.2.0           XML_3.98-1.6      
# [21] BiocParallel_1.8.2 whisker_0.3-2      magrittr_1.5       scales_0.4.1      
# [25] tfruns_1.4         codetools_0.2-15   assertthat_0.1     colorspace_1.3-2  
# [29] tensorflow_1.10    lazyeval_0.2.0     munsell_0.4.3      RCurl_1.95-4.8 

