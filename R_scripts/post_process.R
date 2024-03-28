library('GenomicRanges')
library('GenomicFeatures')
library('ggplot2')
library('patchwork')
library('ggnewscale')
library('cowplot')
library('pheatmap')

# give in input, for each fraction, the 1,000 sets of DRACH- random hits (hits_without) and the 1,000 sets of DRACH+ random hits (hits_DRACH) 
# (both overlapped with the coordinates of RNA marks from the databases) and the percentage
# of ELIGOS DRACH-/DRACH+ hits annotated to the different RNA marks from the DB (ELIGOS_results_without/ELIGOS_results_DRACH).
# Path is the path to the directory in which saving the plot
post_process_mods <- function(hits_without, hits_DRACH, ELIGOS_results_without, ELIGOS_results_DRACH,name_pdf,t,st, path) {
  characterized_hits_without <- unlist(hits_without[[8]])
  hits_without[[8]] <- characterized_hits_without
  names(hits_without) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','Annotated hits')
  names(ELIGOS_results_without) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','Annotated hits')
  
  characterized_hits_DRACH <- unlist(hits_DRACH[[8]])
  hits_DRACH[[8]] <- characterized_hits_DRACH
  names(hits_DRACH) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','Annotated hits')
  names(ELIGOS_results_DRACH) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','Annotated hits')
  
  # compute the absolute difference between the percentage of ELIGOS DRACH- hits annotated to the different RNA marks/at least 
  # one RNA mark and the median percentage of DRACH- random hits annotated to the same RNA mod type divided by the interquantile 
  # distance computed on the 1,000 values of percentage of DRACH- random hits annotated to the same RNA mod type
  IQR_distance_without <- unlist(lapply(seq_along(1:8), function(x,n,i) {
    if (length(x[[i]]) != 1000) {print('error')}
    IQR <- summary(as.numeric(x[[i]]))[5] - summary(as.numeric(x[[i]]))[2]
    if (IQR==0) {
      0
    } else {
      round(abs((ELIGOS_results_without[names(ELIGOS_results_without) == n[i]] - summary(as.numeric(x[[i]]))[3]))/IQR,2)
    }
  }, x=hits_without[1:8], n=names(hits_without[1:8])))
  
  # the same for random DRACH+ hits
  IQR_distance_DRACH <- unlist(lapply(seq_along(1:8), function(x,n,i) {
    if (length(x[[i]]) != 1000) {print('error')}
    IQR <- summary(as.numeric(x[[i]]))[5] - summary(as.numeric(x[[i]]))[2]
    if (IQR == 0) {
      0
    } else {
      round(abs((ELIGOS_results_DRACH[names(ELIGOS_results_DRACH) == n[i]] - summary(as.numeric(x[[i]]))[3]))/IQR,2) 
    }
  }, x=hits_DRACH[1:8], n=names(hits_DRACH[1:8])))
  
  # compute how many of the 1,000 values of percentage of DRACH- random hits annotated to each RNA mod type/at least one RNA mod type
  # are higher than the percentage of ELIGOS DRACH- hits annotated to the same RNA mod type divided by 1,000
  p_value_without <- unlist(lapply(seq_along(1:8), function(x,n,i) {
    if (length(x[[i]][x[[i]]>ELIGOS_results_without[names(ELIGOS_results_without) == n[i]]])/1000 == 0) {
      paste0('< ',as.character(1e-3))
    } else{
      as.character(length(x[[i]][x[[i]]>ELIGOS_results_without[names(ELIGOS_results_without) == n[i]]])/1000)
    } }, x=hits_without[1:8], n=names(hits_without[1:8])))
  names(p_value_without) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','Annotated hits')
  
  # the same for random DRACH+ hits
  p_value_DRACH <- unlist(lapply(seq_along(1:8), function(x,n,i) {
    if (length(x[[i]][x[[i]]>ELIGOS_results_DRACH[names(ELIGOS_results_DRACH) == n[i]]])/1000 == 0) {
      paste0('< ',as.character(1e-3))
    } else{
      as.character(length(x[[i]][x[[i]]>ELIGOS_results_DRACH[names(ELIGOS_results_DRACH) == n[i]]])/1000)
    } }, x=hits_DRACH[1:8], n=names(hits_DRACH[1:8])))
  names(p_value_DRACH) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','Annotated hits')
  
  # collapse the two statistics
  stats_without <- unlist(lapply(seq_along(1:8), function(i) {
    paste0(IQR_distance_without[i],'\n',p_value_without[i])
  }))
  names(stats_without) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','Annotated hits')
  
  stats_DRACH <- unlist(lapply(seq_along(1:8), function(i) {
    paste0(IQR_distance_DRACH[i],'\n',p_value_DRACH[i])
  }))
  names(stats_DRACH) <- c('m6A','Y','m1A','m5C','m7G','A-I','Nm','Annotated hits')
  
  hits_non_m6A <- c(hits_without[!names(hits_without) %in% c('m6A','Annotated hits')],
                    hits_DRACH[!names(hits_DRACH) %in% c('m6A','Annotated hits')])
  hits_m6A <- c(hits_without[names(hits_without) %in% c('m6A','Annotated hits')],
                hits_DRACH[names(hits_DRACH) %in% c('m6A','Annotated hits')])
  
  results_non_m6A <- c(ELIGOS_results_without[!names(ELIGOS_results_without) %in% c('m6A','Annotated hits')],
                       ELIGOS_results_DRACH[!names(ELIGOS_results_DRACH) %in% c('m6A','Annotated hits')])
  results_m6A <- c(ELIGOS_results_without[names(ELIGOS_results_without) %in% c('m6A','Annotated hits')],
                   ELIGOS_results_DRACH[names(ELIGOS_results_DRACH) %in% c('m6A','Annotated hits')])
  
  stats_non_m6A <- c(stats_without[!names(stats_without) %in% c('m6A','Annotated hits')],
                     stats_DRACH[!names(stats_DRACH) %in% c('m6A','Annotated hits')])
  stats_m6A <- c(stats_without[names(stats_without) %in% c('m6A','Annotated hits')],
                 stats_DRACH[names(stats_DRACH) %in% c('m6A','Annotated hits')])
  
  # create a data frame reporting, for each mod type (not m6A), the 1,000 values of percentage of DRACH+/DRACH- random 
  # hits overlapping with that RNA mod type
  data_non_m6A <- data.frame(
    RNA_modification = factor(x = rep(rep(c('Y','m1A','m5C','m7G','A-I','Nm'), each=1000),2), 
                              levels=c('Y','m1A','m5C','m7G','A-I','Nm')),
    value = unlist(lapply(hits_non_m6A, function(x) {x})),
    Class = factor(c(rep('DRACH negative', 6000),rep('DRACH positive', 6000)), levels=c('DRACH positive','DRACH negative'))
  )
  
  # create a data frame reporting the 1,000 values of percentage of DRACH+/DRACH- random hits overlapping with m6A/at least one hit
  data_m6A <- data.frame(
    RNA_modification = factor(x = rep(rep(c('m6A','Annotated hits'), each=1000),2), 
                              levels=c('m6A','Annotated hits')),
    value = unlist(lapply(hits_m6A, function(x) {x})),
    Class = factor(c(rep('DRACH negative', 2000),rep('DRACH positive', 2000)), levels=c('DRACH positive','DRACH negative'))
  )
  
  # create a data frame reporting for each mod type (not m6A) the percentage of ELIGOS DRACH+/DRACH- 
  # hits overlapping with that RNA mod type and the relative statistics
  points_non_m6A <- data.frame(RNA_modification = factor(x = rep(c('Y','m1A','m5C','m7G','A-I','Nm'), 2), 
                                                         levels=c('Y','m1A','m5C','m7G','A-I','Nm')),
                               ELIGOS = results_non_m6A,
                               stats_without= c(stats_non_m6A[1:6],NA,NA,NA,NA,NA,NA),
                               stats_DRACH= c(NA,NA,NA,NA,NA,NA,stats_non_m6A[7:12]),
                               Class = factor(c(rep('DRACH negative', 6),rep('DRACH positive', 6)), levels=c('DRACH positive','DRACH negative')))
  
  # create a data frame reporting the percentage of ELIGOS DRACH+/DRACH- hits overlapping with 
  # m6A/at least one hit and the relative statistics
  points_m6A <- data.frame(RNA_modification = factor(x = rep(c('m6A','Annotated hits'), 2), 
                                                     levels=c('m6A','Annotated hits')),
                           ELIGOS = results_m6A,
                           stats_without= c(stats_m6A[1:2],NA,NA),
                           stats_DRACH= c(NA,NA,stats_m6A[3:4]),
                           Class = factor(c(rep('DRACH negative', 2),rep('DRACH positive', 2)), levels=c('DRACH positive','DRACH negative')))
  
  # generate a violin plot representing the distribution of the 1,000 values of percentage of DRACH+/DRACH- random hits
  # overlapping with each mod type (not m6A). A point is used to report the percentage of ELIGOS DRACH+/DRACH- 
  # hits overlapping with the same RNA mod type and the relative statistics
  p_non_m6A <- ggplot(data = data_non_m6A, mapping = aes(x=RNA_modification, y=value, fill = Class)) + 
    geom_violin(position=position_dodge(width = 0.9)) +
    geom_boxplot(position = position_dodge(width = 0.9),width=0.1, color='black') +
    scale_fill_manual(name="Random sequences",values=c('coral2','cornflowerblue')) +
    xlab('') + 
    ylab('')+
    new_scale_fill()+
    geom_point(data = points_non_m6A, aes(x=RNA_modification, y=ELIGOS, color=Class),size=12, shape=23,stroke = 2.5, position=position_dodge(width=0.9), fill='white')+ 
    scale_color_manual(name="ELIGOS values",values=c('coral2','cornflowerblue')) +
    geom_text(data = points_non_m6A,aes(x=RNA_modification,y=ELIGOS,label=stats_without), nudge_x = 0.43, nudge_y = 0.05, size=15) +
    geom_text(data = points_non_m6A,aes(x=RNA_modification,y=ELIGOS,label=stats_DRACH), nudge_x = 0, nudge_y = 0.05, size=15) +
    theme_classic() +  
    theme(legend.position = 'none', axis.text.x = element_text(size=40),
          axis.text.y = element_text(size=40), legend.title = element_text(size=50), legend.text =  element_text(size=40))
  
  # generate a violin plot representing the distribution of the 1,000 values of percentage of DRACH+/DRACH- random hits
  # overlapping with m6A and of the percentage of hits annotated at least with one RNA mod type. 
  # A point is used to report the percentage of ELIGOS DRACH+/DRACH- 
  # hits overlapping with m6A/at least one hit and the relative statistics
  p_m6A <- ggplot(data = data_m6A, mapping = aes(x=RNA_modification, y=value, fill = Class)) + 
    geom_violin(position=position_dodge(width = 0.9)) +
    geom_boxplot(position = position_dodge(width = 0.9),width=0.1, color='black') +
    scale_fill_manual(name="Random sequences",values=c('coral2','cornflowerblue')) +
    xlab('') + 
    ylab('')+
    new_scale_fill()+
    geom_point(data = points_m6A, aes(x=RNA_modification, y=ELIGOS, color=Class),size=12, shape=23, stroke = 2.5,position=position_dodge(width=0.9), fill='white')+ 
    scale_color_manual(name="ELIGOS values",values=c('coral2','cornflowerblue')) +
    geom_text(data = points_m6A,aes(x=RNA_modification,y=ELIGOS,label=stats_without), nudge_x = 0.4, nudge_y = 2.2, size=15) +
    geom_text(data = points_m6A,aes(x=RNA_modification,y=ELIGOS,label=stats_DRACH), nudge_x = -0.04, nudge_y = 2.2, size=15) +
    theme_classic() +  
    theme(legend.position = 'none', axis.text.x = element_text(size=40),
          axis.text.y = element_text(size=40), legend.title = element_text(size=50), legend.text =  element_text(size=40))
  
  # generate a unique plot
  p <- p_non_m6A / (p_m6A + plot_spacer()) + plot_annotation(title=t, subtitle=st) + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.justification = 'right',
                                                                                                                             plot.title = element_text(size=50,hjust = 0.5), plot.subtitle = element_text(size=45,hjust = 0.5))
  
  ggsave(paste0(path, name_pdf), plot=p,width = 40, height = 35)
  
  significant_RNA_mods_without_DRACH <- names(p_value_without)[p_value_without <= 0.05]
  significant_RNA_mods_with_DRACH <- names(p_value_DRACH)[p_value_DRACH <= 0.05]
  
  significant <- list(p_value_without, p_value_DRACH, significant_RNA_mods_without_DRACH[significant_RNA_mods_without_DRACH != 'Annotated hits'],
                      significant_RNA_mods_with_DRACH[significant_RNA_mods_with_DRACH != 'Annotated hits'])
  names(significant) <- c('p_value_without','p_value_DRACH','significant_RNA_mods_without_DRACH','significant_RNA_mods_with_DRACH')
  return(significant)
}


# create a heatmap reporting, for each RNA modification type, in which gene parts ELIGOS hits (considering DRACH+ and DRACH- hits together)
# annotated to that RNA modification type - both from the overlap with RNA marks and from the overlap with the effectors' binding sites 
# from the public databases - fall. 
# Limiting to the RNA mod significantly enriched in ELIGOS DRACH- (mods_significant_without_DRACH) and DRACH+ (mods_significant_with_DRACH) 
# hits with respect to the random data.
# This function is executed for each fraction passing the GRanges object with ELIGOS DRACH- (hits_non_DRACH) and DRACH+ hits 
# (hits_DRACH) that have been overlapped with RNA marks from the public databases.
# path_directory is the path to the directory containing the results of ELIGOS analysis. m6A is set to TRUE (plot
# a heatmap with the distribution of ELIGOS hits annotated to m6A) or FALSE (plot
# a heatmap with the distribution of ELIGOS hits annotated to all the RNA mod types except m6A)
generate_heatmap_gene <- function(hits_non_DRACH, hits_DRACH, fraction, path_directory, m6A, mods_significant_without_DRACH, mods_significant_with_DRACH) {
  
  if (!m6A) {
    # create a list reporting all ELIGOS DRACH- hits annotated to the significantly enriched RNA modifications
    hits_without_DRACH_with_mod <- unlist(lapply(as.list(mods_significant_without_DRACH), function(x) {
      unlist(as(unlist(lapply(seq_along(hits_non_DRACH), function(i) {
        if (x %in% unlist(strsplit(hits_non_DRACH[i]$mod_type, split=';'))) {
          hits_non_DRACH[i]
        }})), 'GRangesList'))
    }))
    names(hits_without_DRACH_with_mod) <- mods_significant_without_DRACH
    
    # create a list reporting all ELIGOS DRACH+ hits annotated to the significantly enriched RNA modifications
    hits_with_DRACH_with_mod <- unlist(lapply(as.list(mods_significant_with_DRACH), function(x) {
      unlist(as(unlist(lapply(seq_along(hits_DRACH), function(i) {
        if (x %in% unlist(strsplit(hits_DRACH[i]$mod_type, split=';'))) {
          hits_DRACH[i]
        }})), 'GRangesList'))
    }))
    names(hits_with_DRACH_with_mod) <- mods_significant_with_DRACH
    
    all_mods_significantly_enriched <- union(mods_significant_without_DRACH, mods_significant_with_DRACH)
    
    mod_sequence <- c('Y', 'm1A', 'm5C', 'm7G', 'A-I', 'Nm', 'm6Am')
    mod_sequence <- mod_sequence[mod_sequence %in% all_mods_significantly_enriched]
    
    # create a unique GRangesList reporting, for each RNA modification type significantly enriched in at least
    # DRACH- or DRACH+ hits, all the hits annotated to that RNA mod
    all_mods_significant <- lapply(as.list(mod_sequence), function (x) {
      if ((x %in% mods_significant_without_DRACH) & (x %in% mods_significant_with_DRACH)) {
        return(c(hits_without_DRACH_with_mod[[x]], hits_with_DRACH_with_mod[[x]]))
      } else if ((x %in% mods_significant_without_DRACH) & (!x %in% mods_significant_with_DRACH)) {
        return(hits_without_DRACH_with_mod[[x]])
      } else {
        return(hits_with_DRACH_with_mod[[x]])
      }
    })
    names(all_mods_significant) <- mod_sequence
    
    # initiate a matrix in which, for each RNA mod type (row) and gene part (column), the number of ELIGOS hits annotated
    # to that mod type and falling in that gene part is reported
    m <- matrix(0, ncol = 5, nrow = length(mod_sequence)) 
    colnames(m) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(m) <- mod_sequence
    
    # iterate over the RNA mod types
    for (i in 1:length(mod_sequence)) {
      mod <- mod_sequence[i]
      
      # identify the genes on which the hits annotated to this RNA mod type fall
      genes <- intersect(names(protein_coding_genes_5UTR_3UTR_introns_exons_stop),unique(all_mods_significant[[mod]]$gene_id))
      
      # iterate over the hits annotated to this RNA mod type and identify in which gene part they fall
      for (hit in 1:length(all_mods_significant[[mod]])) {
        gene <- all_mods_significant[[mod]][hit]$gene_id
        parts_gene <- protein_coding_genes_5UTR_3UTR_introns_exons_stop[gene][[1]]
        over <- findOverlaps(all_mods_significant[[mod]][hit], parts_gene)
        
        if (length(subjectHits(over)) != 0) {
          which_features <- unique(parts_gene[subjectHits(over)]$feature)    
          for (colname in colnames(m)) { 
            if (colname %in% which_features) {
              m[mod,colname] <- m[mod,colname] +1
            } 
          }
        }
      }
    }
    
    pdf(paste0(path_directory, '/mods_RBPs/mods_distribution_on_gene_', fraction, '.pdf'), height = 7, width = 5)
    pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, 
             treeheight_row = 0, fontsize = 10, color = colorRampPalette(c("white", "red"))(50), 
             legend_breaks = round(seq(0,max(m),length.out=13),0), legend_labels = as.character(round(seq(0,max(m),length.out=13),0)), fontsize_row = 10, fontsize_col = 10,display_numbers = m)
    dev.off()
  } else {
    # repeat the same for ELIGOS hits annotated with m6A
    all_mods <- 'm6A'
    
    hits_each_mod <- lapply(as.list(all_mods), function(x) {
      hits_non_DRACH_with_mod <- unlist(lapply(seq_along(hits_non_DRACH), function(i) {
        if (x %in% unlist(strsplit(hits_non_DRACH[i]$mod_type, split=';'))) {
          hits_non_DRACH[i]
        }
      }))
      
      hits_DRACH_with_mod <- unlist(lapply(seq_along(hits_DRACH), function(i) {
        if (x %in% unlist(strsplit(hits_DRACH[i]$mod_type, split=';'))) {
          hits_DRACH[i]
        }
      }))
      
      unlist(as(c(hits_non_DRACH_with_mod, hits_DRACH_with_mod), 'GRangesList'))
    })
    names(hits_each_mod) <- all_mods
    
    m <- matrix(0, ncol = 5, nrow = 1) 
    colnames(m) <- c('5UTR', 'coding exon', 'intron', 'stop codon', '3UTR')
    rownames(m) <- names(hits_each_mod)
    
    for (i in 1:length(hits_each_mod)) {
      mod <- names(hits_each_mod)[i]
      
      genes <- intersect(names(protein_coding_genes_5UTR_3UTR_introns_exons_stop),unique(hits_each_mod[[i]]$gene_id))
      
      for (hit in 1:length(hits_each_mod[[i]])) {
        gene <- hits_each_mod[[i]][hit]$gene_id
        parts_gene <- protein_coding_genes_5UTR_3UTR_introns_exons_stop[gene][[1]]
        over <- findOverlaps(hits_each_mod[[i]][hit], parts_gene)
        
        if (length(subjectHits(over)) != 0) {
          which_features <- unique(parts_gene[subjectHits(over)]$feature)  
          for (colname in colnames(m)) { 
            if (colname %in% which_features) {
              m[mod,colname] <- m[mod,colname] +1
            } 
          }
        }
      }
    }
    
    pdf(paste0(path_directory, '/mods_RBPs/m6A_distribution_on_gene_', fraction, '.pdf'), height = 1, width = 5)
    pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, legend = FALSE,
             treeheight_row = 0, main = fraction, fontsize = 10, color = colorRampPalette(c("white", "red"))(50), 
             fontsize_row = 10, fontsize_col = 10,display_numbers = m)
    dev.off()
  }
} 

# directory_hits is the path to the directory containing the results of ELIGOS analysis.
# mods_significant_without_DRACH reports, for each fraction, the RNA mods significantly enriched in ELIGOS DRACH- hits of that fraction
# with respect to the random data, mods_significant_with_DRACH those enriched in ELIGOS DRACH+ hits of each fraction
print_heatmap <- function(directory_hits, mods_significant_without_DRACH_chr, mods_significant_with_DRACH_chr,
                          mods_significant_without_DRACH_nucleo, mods_significant_with_DRACH_nucleo,
                          mods_significant_without_DRACH_cyto, mods_significant_with_DRACH_cyto) {
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  
  generate_heatmap_gene(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings, hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings, 'chromatin', directory_hits, FALSE, mods_significant_without_DRACH_chr, mods_significant_with_DRACH_chr)
  generate_heatmap_gene(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings, hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings, 'chromatin', directory_hits, TRUE, mods_significant_without_DRACH_chr, mods_significant_with_DRACH_chr)
  
  generate_heatmap_gene(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings, hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings, 'nucleoplasm', directory_hits, FALSE, mods_significant_without_DRACH_nucleo, mods_significant_with_DRACH_nucleo)
  generate_heatmap_gene(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings, hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings, 'nucleoplasm', directory_hits, TRUE, mods_significant_without_DRACH_nucleo, mods_significant_with_DRACH_nucleo)
  
  generate_heatmap_gene(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings, hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings, 'cytoplasm', directory_hits, FALSE, mods_significant_without_DRACH_cyto, mods_significant_with_DRACH_cyto)
  generate_heatmap_gene(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings, hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings, 'cytoplasm', directory_hits, TRUE, mods_significant_without_DRACH_cyto, mods_significant_with_DRACH_cyto)
}

# generate a heatmap reporting, for each fraction and RNA mod, the percentage of ELIGOS hits annotated to that RNA mod - after the overlap
# with the RNA marks reported in the public databases - and the p-value of the enrichment with respect to the random data
heatmap_enrichment <- function(directory_hits, percentage_chr, percentage_nucleo, percentage_cyto, 
                               p_value_chr, p_value_nucleo, p_value_cyto, DRACH) {
  
  df <- data.frame(RNA_mod = factor(rep(c('m6A', 'Y', 'm1A', 'm5C', 'm7G', 'A-I', 'Nm'),3), levels = rev(c('m6A', 'Y', 'm1A', 'm5C', 'm7G', 'A-I', 'Nm'))),
                   Fraction = factor(rep(c('chromatin', 'nucleoplasm', 'cytoplasm'), each =  7), levels=c('chromatin', 'nucleoplasm', 'cytoplasm')),
                   values = paste0(as.character(c(percentage_chr, percentage_nucleo, percentage_cyto)),'%'),
                   p_value = as.factor(c(as.numeric(gsub('< ', '', p_value_chr[names(p_value_chr) != 'Annotated hits'])), 
                                         as.numeric(gsub('< ', '', p_value_nucleo[names(p_value_nucleo) != 'Annotated hits'])), 
                                         as.numeric(gsub('< ', '', p_value_cyto[names(p_value_cyto) != 'Annotated hits'])))))
  
  p <- ggplot(df,aes(x=Fraction, y=RNA_mod, fill =p_value)) +
    geom_tile() +
    geom_text(aes(label = values), size=10)+
    theme_classic() +
    xlab('') + 
    ylab('') +
    scale_fill_manual(values= c("#FF8B00", "#FFAA00", "#FFB800", "#FFC600", "#FFD500", "#FFE300", "#FFF100", "#FFFF00", "#FFFF15", "#FFFF40", "#FFFF6A",
                                "#FFFF95", "#FFFFBF", "#FFFFEA")) +
    theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), legend.title = element_text(size=20), legend.text =  element_text(size=15))
  
  if (DRACH) {
    ggsave(paste0(directory_hits, '/RNA_mods_enrichment_DRACH.pdf'), plot=p,width = 13, height = 17)
  } else {
    ggsave(paste0(directory_hits, '/RNA_mods_enrichment.pdf'), plot=p,width = 13, height = 17)
  }
}

# identify the combinations of RNA modifications (limiting to the RNA mods significantly enriched in ELIGOS DRACH- 
# (mods_significant_without_DRACH) and DRACH+ (mods_significant_with_DRACH) hits with respect to the random data)
# co-occurring on the genes. Both the information about the modification type to which 
# each ELIGOS hit has been annotated during the overlap with the coordinates of RNA marks and the overlap with the effectors' 
# binding sites from RMBase3 and RMVar are considered.
# If condition = 'both' both ELIGOS DRACH- and DRACH+ hits are considered, if condition = 'DRACH' only ELIGOS DRACH+ hits are evaluated,
# otherwise only ELIGOS DRACH- hits
mods_combination_only_significant <- function(directory_hits,condition,mods_significant_without_DRACH_chr, mods_significant_with_DRACH_chr,
                                              mods_significant_without_DRACH_nucleo, mods_significant_with_DRACH_nucleo,
                                              mods_significant_without_DRACH_cyto, mods_significant_with_DRACH_cyto) {
  
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  
  if (condition == 'both') {
    # keep only ELIGOS DRACH- hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_chr]
      hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_nucleo]
      hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_cyto]
      hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_chr]
      hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_nucleo]
      hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_cyto]
      hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # for each fraction, concatenate ELIGOS DRACH+ and DRACH- hits
    hits_chr <- c(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings, hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings)
    hits_nucleo <- c(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings, hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings)
    hits_cyto<- c(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings, hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings)
    
  } else if (condition == 'DRACH') {
    # keep only ELIGOS DRACH+ hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_chr]
      hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_nucleo]
      hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_cyto]
      hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    hits_chr <- hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings
    hits_nucleo <- hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings
    hits_cyto<- hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings
  } else {
    # keep only ELIGOS DRACH- hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_chr]
      hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_nucleo]
      hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_cyto]
      hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    hits_chr <- hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings
    hits_nucleo <- hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings
    hits_cyto<- hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings
  }
  
  # initiate a matrix reporting, for each gene, the RNA modification types of the hits mapping on it
  all_genes_mods_chr <- matrix(ncol=1, nrow=length(unique(hits_chr$gene_id))) 
  rownames(all_genes_mods_chr) <- unique(hits_chr$gene_id)
  
  if (condition == 'both') {
    # iterate over the genes
    lapply(seq_along(rownames(all_genes_mods_chr)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      # identify with which RNA mod types are annotated ELIGOS hits mapping on that gene 
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          # if there are at least two specific RNA mod types they are added to the matrix
          all_genes_mods_chr[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_chr), y=hits_chr)
  } else {
    # iterate over the genes
    lapply(seq_along(rownames(all_genes_mods_chr)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      # identify with which RNA mod types are annotated ELIGOS DRACH+/DRACH- hits mapping on that gene 
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene aren't annotated with any specific modification type skip this gene
        NULL
      } else if (length(mods_gene) > 1) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene are annotated with multiple specific modification types they are reported in the matrix 
        all_genes_mods_chr[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene are annotated only with m6A they are reported in the matrix 
        all_genes_mods_chr[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_chr), y=hits_chr)
  }
  
  all_genes_mods_chr <- all_genes_mods_chr[!is.na(all_genes_mods_chr[,1]),1]
  
  # create a data frame reporting the frequency with which each combination of RNA modifications is found in the genes
  df_chr <- data.frame(table(all_genes_mods_chr))
  colnames(df_chr) <- c('Mods_combination', 'Frequency')
  
  all_genes_mods_nucleo <- matrix(ncol=1, nrow=length(unique(hits_nucleo$gene_id))) 
  rownames(all_genes_mods_nucleo) <- unique(hits_nucleo$gene_id)
  
  if (condition == 'both') {
    lapply(seq_along(rownames(all_genes_mods_nucleo)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          all_genes_mods_nucleo[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_nucleo), y=hits_nucleo)
  } else {
    lapply(seq_along(rownames(all_genes_mods_nucleo)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        NULL
      } else if (length(mods_gene) > 1) {
        all_genes_mods_nucleo[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        all_genes_mods_nucleo[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_nucleo), y=hits_nucleo)
  }
  
  all_genes_mods_nucleo <- all_genes_mods_nucleo[!is.na(all_genes_mods_nucleo[,1]),1]
  
  df_nucleo <- data.frame(table(all_genes_mods_nucleo))
  colnames(df_nucleo) <- c('Mods_combination', 'Frequency')
  
  all_genes_mods_cyto <- matrix(ncol=1, nrow=length(unique(hits_cyto$gene_id))) 
  rownames(all_genes_mods_cyto) <- unique(hits_cyto$gene_id)
  
  if (condition == 'both') {
    lapply(seq_along(rownames(all_genes_mods_cyto)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          all_genes_mods_cyto[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_cyto), y=hits_cyto)
  } else {
    lapply(seq_along(rownames(all_genes_mods_cyto)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        NULL
      } else if (length(mods_gene) > 1) {
        all_genes_mods_cyto[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        all_genes_mods_cyto[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_cyto), y=hits_cyto)
  }
  
  all_genes_mods_cyto <- all_genes_mods_cyto[!is.na(all_genes_mods_cyto[,1]),1]
  
  df_cyto <- data.frame(table(all_genes_mods_cyto))
  colnames(df_cyto) <- c('Mods_combination', 'Frequency')
  
  # identify all the possible combinations of RNA modifications found in at least one fraction
  all_mods_comb <- as.vector(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))
  all_mods_comb <- rep(all_mods_comb, each=3)
  Fractions <- factor(rep(c('chromatin', 'nucleoplasm', 'cytoplasm'), length(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))), levels = c('chromatin', 'nucleoplasm', 'cytoplasm'))
  number <- c()
  
  # for each fraction, report the frequency with which each combination of RNA modifications is found in the genes
  for (mod in as.vector(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))) {
    if (mod %in% df_chr$Mods_combination) {
      number_chr <- df_chr[df_chr$Mods_combination == mod,2]
    } else {
      number_chr <- 0
    }
    if (mod %in% df_nucleo$Mods_combination) {
      number_nucleo <- df_nucleo[df_nucleo$Mods_combination == mod,2]
    } else {
      number_nucleo <- 0
    }
    if (mod %in% df_cyto$Mods_combination) {
      number_cyto <- df_cyto[df_cyto$Mods_combination == mod,2]
    } else {
      number_cyto <- 0
    }
    number <- c(number, number_chr, number_nucleo, number_cyto)
  }
  
  df <- data.frame(Fractions, all_mods_comb, number)
  
  order <- lapply(unique(df$all_mods_comb), function(x) median(df$number[df$all_mods_comb == x]))
  names(order) <- unique(df$all_mods_comb)
  order <- sort(unlist(order), decreasing = TRUE)
  
  if (length(order) > 20) {
    df <- df[df$all_mods_comb %in% names(order)[1:20],]
    x_labels <- names(order)[1:20]
  } else {
    x_labels <- names(order)
  }
  
  # plot a barplot reporting, for the 10 combinations of RNA marks occurring in the genes of the different fractions with the highest median
  # frequency, the frequency of that combination in each fraction
  p <- ggplot(data = df, aes(fill=Fractions, y=number, x=all_mods_comb)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_classic() +
    theme(axis.text.x=element_text(size=45),axis.text.y=element_text(size=40), axis.title = element_text(size=45), legend.text=element_text(size=45), legend.title=element_text(size=45), legend.position = 'top') +
    xlab('')+
    ylab('Number of genes carrying a combination of modifications') +
    scale_x_discrete(limits=x_labels) 
  
  if (condition == 'both') {
    ggsave(paste0(directory_hits, 'mods_RBPs/mods_combinations.pdf'), plot = p, height = 20, width = 45)
  } else if (condition == 'DRACH') {
    ggsave(paste0(directory_hits, 'mods_RBPs/mods_combinations_with_DRACH.pdf'), plot = p, height = 20, width = 45)
  } else {
    ggsave(paste0(directory_hits, 'mods_RBPs/mods_combinations_without_DRACH.pdf'), plot = p, height = 20, width = 45)
  }
}

# identify which combinations of RNA marks have a median frequency at least equal to 4 in at least one of the analyses (either on nascent reads
# or on all the reads with the library-level subsampling threshold used for nascent reads).
# Generate a scatterplot reporting on the x axis the median frequency across the fractions of the selected combinations of RNA marks in nascent analysis
# and on the y axis the median frequency across the fractions of the selected combinations of RNA marks in the analysis on all the reads.
# If condition = 'both' both ELIGOS DRACH- and DRACH+ hits are considered, if condition = 'DRACH' only ELIGOS DRACH+ hits are evaluated,
# otherwise only ELIGOS DRACH- hits
comparison_mods_combination_nascent_total_only_significant <- function(directory_hits_nascent, directory_hits_total, condition,
                                                                       mods_significant_without_DRACH_chr_nascent, mods_significant_with_DRACH_chr_nascent,
                                                                       mods_significant_without_DRACH_nucleo_nascent, mods_significant_with_DRACH_nucleo_nascent,
                                                                       mods_significant_without_DRACH_cyto_nascent, mods_significant_with_DRACH_cyto_nascent,
                                                                       mods_significant_without_DRACH_chr_total, mods_significant_with_DRACH_chr_total,
                                                                       mods_significant_without_DRACH_nucleo_total, mods_significant_with_DRACH_nucleo_total,
                                                                       mods_significant_without_DRACH_cyto_total, mods_significant_with_DRACH_cyto_total) {
  
  load(paste0(directory_hits_nascent,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_nascent,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_nascent,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_nascent,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_nascent,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_nascent,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  
  if (condition == 'both') {
    # keep only ELIGOS DRACH- hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_chr_nascent]
      hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_nucleo_nascent]
      hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_cyto_nascent]
      hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_chr_nascent]
      hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_nucleo_nascent]
      hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_cyto_nascent]
      hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # for each fraction, concatenate ELIGOS DRACH+ and DRACH- hits
    hits_chr <- c(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings, hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings)
    hits_nucleo <- c(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings, hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings)
    hits_cyto<- c(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings, hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings)
    
  } else if (condition == 'DRACH') {
    # keep only ELIGOS DRACH+ hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_chr_nascent]
      hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_nucleo_nascent]
      hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_cyto_nascent]
      hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    hits_chr <- hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings
    hits_nucleo <- hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings
    hits_cyto<- hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings
    
  } else {
    # keep only ELIGOS DRACH- hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_chr_nascent]
      hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_nucleo_nascent]
      hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_cyto_nascent]
      hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    hits_chr <- hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings
    hits_nucleo <- hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings
    hits_cyto<- hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings
  }
  
  # initiate a matrix reporting, for each gene, the RNA modification types of the hits mapping on it
  all_genes_mods_chr <- matrix(ncol=1, nrow=length(unique(hits_chr$gene_id))) 
  rownames(all_genes_mods_chr) <- unique(hits_chr$gene_id)
  
  if (condition == 'both') {
    # iterate over the genes
    lapply(seq_along(rownames(all_genes_mods_chr)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      # identify with which RNA mod types are annotated ELIGOS hits mapping on that gene 
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          # if there are at least two specific RNA mod types they are added to the matrix
          all_genes_mods_chr[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_chr), y=hits_chr)
  } else {
    # iterate over the genes
    lapply(seq_along(rownames(all_genes_mods_chr)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      # identify with which RNA mod types are annotated ELIGOS DRACH+/DRACH- hits mapping on that gene 
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene aren't annotated with any specific modification type skip this gene
        NULL
      } else if (length(mods_gene) > 1) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene are annotated with multiple specific modification types they are reported in the matrix 
        all_genes_mods_chr[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene are annotated only with m6A they are reported in the matrix 
        all_genes_mods_chr[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_chr), y=hits_chr)
  }
  all_genes_mods_chr <- all_genes_mods_chr[!is.na(all_genes_mods_chr[,1]),1]
  
  # create a data frame reporting the frequency with which each combination of RNA modifications is found in the genes
  df_chr <- data.frame(table(all_genes_mods_chr))
  colnames(df_chr) <- c('Mods_combination', 'Frequency')
  
  all_genes_mods_nucleo <- matrix(ncol=1, nrow=length(unique(hits_nucleo$gene_id))) 
  rownames(all_genes_mods_nucleo) <- unique(hits_nucleo$gene_id)
  
  if (condition == 'both') {
    lapply(seq_along(rownames(all_genes_mods_nucleo)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          all_genes_mods_nucleo[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_nucleo), y=hits_nucleo)
  } else {
    lapply(seq_along(rownames(all_genes_mods_nucleo)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        NULL
      } else if (length(mods_gene) > 1) {
        all_genes_mods_nucleo[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        all_genes_mods_nucleo[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_nucleo), y=hits_nucleo)
  }
  all_genes_mods_nucleo <- all_genes_mods_nucleo[!is.na(all_genes_mods_nucleo[,1]),1]
  
  df_nucleo <- data.frame(table(all_genes_mods_nucleo))
  colnames(df_nucleo) <- c('Mods_combination', 'Frequency')
  
  all_genes_mods_cyto <- matrix(ncol=1, nrow=length(unique(hits_cyto$gene_id))) 
  rownames(all_genes_mods_cyto) <- unique(hits_cyto$gene_id)
  
  if (condition == 'both') {
    lapply(seq_along(rownames(all_genes_mods_cyto)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          all_genes_mods_cyto[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_cyto), y=hits_cyto)
  } else {
    lapply(seq_along(rownames(all_genes_mods_cyto)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        NULL
      } else if (length(mods_gene) > 1) {
        all_genes_mods_cyto[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        all_genes_mods_cyto[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_cyto), y=hits_cyto)
  }
  all_genes_mods_cyto <- all_genes_mods_cyto[!is.na(all_genes_mods_cyto[,1]),1]
  
  df_cyto <- data.frame(table(all_genes_mods_cyto))
  colnames(df_cyto) <- c('Mods_combination', 'Frequency')
  
  # identify all the possible combinations of RNA modifications found in at least one fraction
  all_mods_comb <- as.vector(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))
  all_mods_comb <- rep(all_mods_comb, each=3)
  Fractions <- factor(rep(c('chromatin', 'nucleoplasm', 'cytoplasm'), length(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))), levels = c('chromatin', 'nucleoplasm', 'cytoplasm'))
  number <- c()
  
  # for each fraction, report the frequency with which each combination of RNA modifications is found in the genes
  for (mod in as.vector(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))) {
    if (mod %in% df_chr$Mods_combination) {
      number_chr <- df_chr[df_chr$Mods_combination == mod,2]
    } else {
      number_chr <- 0
    }
    if (mod %in% df_nucleo$Mods_combination) {
      number_nucleo <- df_nucleo[df_nucleo$Mods_combination == mod,2]
    } else {
      number_nucleo <- 0
    }
    if (mod %in% df_cyto$Mods_combination) {
      number_cyto <- df_cyto[df_cyto$Mods_combination == mod,2]
    } else {
      number_cyto <- 0
    }
    number <- c(number, number_chr, number_nucleo, number_cyto)
  }
  
  df_nascent <- data.frame(Fractions, all_mods_comb, number)
  order_nascent <- lapply(unique(df_nascent$all_mods_comb), function(x) mean(df_nascent$number[df_nascent$all_mods_comb == x]))
  names(order_nascent) <- unique(df_nascent$all_mods_comb)
  order_nascent <- sort(unlist(order_nascent), decreasing = TRUE)
  # identify the combinations of RNA mods with a median frequency across the fractions at least equal to 4 (nascent reads)
  order_nascent <- order_nascent[order_nascent>=4]
  
  load(paste0(directory_hits_total,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_total,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_total,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_total,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_total,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits_total,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  
  if (condition == 'both') {
    # keep only ELIGOS DRACH- hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_chr_nascent]
      hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_nucleo_nascent]
      hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_cyto_nascent]
      hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_chr_nascent]
      hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_nucleo_nascent]
      hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_cyto_nascent]
      hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # for each fraction, concatenate ELIGOS DRACH+ and DRACH- hits
    hits_chr <- c(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings, hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings)
    hits_nucleo <- c(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings, hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings)
    hits_cyto<- c(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings, hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings)
    
  } else if (condition == 'DRACH') {
    # keep only ELIGOS DRACH+ hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_chr_nascent]
      hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_nucleo_nascent]
      hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH+ hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_with_DRACH_cyto_nascent]
      hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    hits_chr <- hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings
    hits_nucleo <- hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings
    hits_cyto<- hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings
    
  } else {
    # keep only ELIGOS DRACH- hits from chromatin associated RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_chr_nascent]
      hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from nucleoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_nucleo_nascent]
      hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    # keep only ELIGOS DRACH- hits from cytoplasmic RNAs annotated to significant RNA modifications
    lapply(seq_along(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings), function(i) {
      mods <- unlist(strsplit(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type, split =';'))
      mods_significant <- mods[mods %in% mods_significant_without_DRACH_cyto_nascent]
      hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[i]$mod_type <<- paste(mods_significant, collapse = ';')
    })
    
    hits_chr <- hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings
    hits_nucleo <- hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings
    hits_cyto<- hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings
  }
  
  # initiate a matrix reporting, for each gene, the RNA modification types of the hits mapping on it
  all_genes_mods_chr <- matrix(ncol=1, nrow=length(unique(hits_chr$gene_id))) 
  rownames(all_genes_mods_chr) <- unique(hits_chr$gene_id)
  
  if (condition == 'both') {
    # iterate over the genes
    lapply(seq_along(rownames(all_genes_mods_chr)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      # identify with which RNA mod types are annotated ELIGOS hits mapping on that gene 
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          # if there are at least two specific RNA mod types they are added to the matrix
          all_genes_mods_chr[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_chr), y=hits_chr)
  } else {
    # iterate over the genes
    lapply(seq_along(rownames(all_genes_mods_chr)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      # identify with which RNA mod types are annotated ELIGOS DRACH+/DRACH- hits mapping on that gene 
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene aren't annotated with any specific modification type skip this gene
        NULL
      } else if (length(mods_gene) > 1) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene are annotated with multiple specific modification types they are reported in the matrix 
        all_genes_mods_chr[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        # if ELIGOS DRACH+/DRACH- hits mapping on this gene are annotated only with m6A they are reported in the matrix 
        all_genes_mods_chr[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_chr), y=hits_chr)
  }
  all_genes_mods_chr <- all_genes_mods_chr[!is.na(all_genes_mods_chr[,1]),1]
  
  # create a data frame reporting the frequency with which each combination of RNA modifications is found in the genes
  df_chr <- data.frame(table(all_genes_mods_chr))
  colnames(df_chr) <- c('Mods_combination', 'Frequency')
  
  all_genes_mods_nucleo <- matrix(ncol=1, nrow=length(unique(hits_nucleo$gene_id))) 
  rownames(all_genes_mods_nucleo) <- unique(hits_nucleo$gene_id)
  
  if (condition == 'both') {
    lapply(seq_along(rownames(all_genes_mods_nucleo)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          all_genes_mods_nucleo[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_nucleo), y=hits_nucleo)
  } else {
    lapply(seq_along(rownames(all_genes_mods_nucleo)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        NULL
      } else if (length(mods_gene) > 1) {
        all_genes_mods_nucleo[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        all_genes_mods_nucleo[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_nucleo), y=hits_nucleo)
  }
  all_genes_mods_nucleo <- all_genes_mods_nucleo[!is.na(all_genes_mods_nucleo[,1]),1]
  
  df_nucleo <- data.frame(table(all_genes_mods_nucleo))
  colnames(df_nucleo) <- c('Mods_combination', 'Frequency')
  
  all_genes_mods_cyto <- matrix(ncol=1, nrow=length(unique(hits_cyto$gene_id))) 
  rownames(all_genes_mods_cyto) <- unique(hits_cyto$gene_id)
  
  if (condition == 'both') {
    lapply(seq_along(rownames(all_genes_mods_cyto)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      if (length(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]) > 1) {
        mods_gene <- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
        if (mods_gene != '') {
          all_genes_mods_cyto[x[i],] <<- mods_gene
        }
      }
    }, x=rownames(all_genes_mods_cyto), y=hits_cyto)
  } else {
    lapply(seq_along(rownames(all_genes_mods_cyto)), function(x,y,i) {
      genes <- y[y$gene_id == x[i]]
      mods_gene <- unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]
      if (length(mods_gene) == 0) {
        NULL
      } else if (length(mods_gene) > 1) {
        all_genes_mods_cyto[x[i],] <<- paste(sort(unique(unlist(strsplit(genes$mod_type, split=';')))[!unique(unlist(strsplit(genes$mod_type, split=';'))) %in% c('Unknown', 'Ambiguous')]), collapse='\n')
      } else if ((length(mods_gene) == 1) & (mods_gene == 'm6A')) {
        all_genes_mods_cyto[x[i],] <<- 'm6A'
      }
    }, x=rownames(all_genes_mods_cyto), y=hits_cyto)
  }
  all_genes_mods_cyto <- all_genes_mods_cyto[!is.na(all_genes_mods_cyto[,1]),1]
  
  df_cyto <- data.frame(table(all_genes_mods_cyto))
  colnames(df_cyto) <- c('Mods_combination', 'Frequency')
  
  # identify all the possible combinations of RNA modifications found in at least one fraction
  all_mods_comb <- as.vector(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))
  all_mods_comb <- rep(all_mods_comb, each=3)
  Fractions <- factor(rep(c('chromatin', 'nucleoplasm', 'cytoplasm'), length(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))), levels = c('chromatin', 'nucleoplasm', 'cytoplasm'))
  number <- c()
  
  # for each fraction, report the frequency with which each combination of RNA modifications is found in the genes
  for (mod in as.vector(unique(df_chr$Mods_combination, df_nucleo$Mods_combination, df_cyto$Mods_combination))) {
    if (mod %in% df_chr$Mods_combination) {
      number_chr <- df_chr[df_chr$Mods_combination == mod,2]
    } else {
      number_chr <- 0
    }
    if (mod %in% df_nucleo$Mods_combination) {
      number_nucleo <- df_nucleo[df_nucleo$Mods_combination == mod,2]
    } else {
      number_nucleo <- 0
    }
    if (mod %in% df_cyto$Mods_combination) {
      number_cyto <- df_cyto[df_cyto$Mods_combination == mod,2]
    } else {
      number_cyto <- 0
    }
    number <- c(number, number_chr, number_nucleo, number_cyto)
  }
  
  df_total <- data.frame(Fractions, all_mods_comb, number)
  order_total <- lapply(unique(df_total$all_mods_comb), function(x) mean(df_total$number[df_total$all_mods_comb == x]))
  names(order_total) <- unique(df_total$all_mods_comb)
  order_total <- sort(unlist(order_total), decreasing = TRUE)
  # identify the combinations of RNA mods with a median frequency across the fractions at least equal to 4 (all reads)
  order_total <- order_total[order_total >=4]
  
  all_comb <- unique(c(names(order_nascent), names(order_total)))
  df_tot <- matrix(ncol=2, nrow = length(all_comb))
  colnames(df_tot) <- c('Nascent', 'Total')
  rownames(df_tot) <- all_comb
  
  # if a combination of RNA marks doesn't have a frequency at least equal to 4 in the analysis on nascent reads/on total reads
  # set the corresponding median frequency to 0
  for (r in rownames(df_tot)) {
    if (r %in% names(order_nascent)) {
      df_tot[r,1] <- order_nascent[names(order_nascent) == r]
    } else {
      df_tot[r,1] <- 0
    }
    if (r %in% names(order_total)) {
      df_tot[r,2] <- order_total[names(order_total) == r]
    } else {
      df_tot[r,2] <- 0
    }
  }
  
  df_tot <- as.data.frame(df_tot)
  
  p <- ggplot(df_tot, aes(x=Nascent, y=Total)) +
    geom_point() +
    geom_text(aes(label=rownames(df_tot)),hjust = -0.1,vjust = -0.1, size =8) +
    #geom_text(aes(label=ifelse(Total==5,rownames(df_tot),'')),hjust = -0.1,vjust = +1.1, size=5) +
    geom_abline(data = data.frame(x=c(0,20), y=c(0,20)), aes(intercept=0, slope=1, colour ='grey'), show.legend = FALSE) +
    theme_classic() +
    #xlim(0,20) +
    #ylim(0,20) +
    labs(title='Co-occurring RNA modifications', subtitle = paste0('Total number of combinations of RNA modifications in nascent RNAs: ',
                                                                   as.character(nrow(df_nascent)), ';\nTotal number of combinations of RNA modifications in all the RNAs: ',
                                                                   as.character(nrow(df_total)))) +
    xlab('Nascent reads') +
    ylab('Nascent + pre-existing reads') +
    theme(text = element_text(size = 25), axis.text.x=element_text(size=20),axis.text.y=element_text(size=20)) +
    theme(plot.title = element_text(size=25)) +
    theme(plot.subtitle = element_text(size=20))
  
  if (condition == 'both') {
    ggsave(paste0(directory_hits_nascent,'/mods_RBPs/mods_frequency_nascent_total_DRACHneg_DRACHpos.pdf'), plot = p, height = 30, width = 30)
  } else if (condition == 'DRACH') {
    ggsave(paste0(directory_hits_nascent,'/mods_RBPs/mods_frequency_nascent_total_DRACHpos.pdf'), plot = p, height = 20, width = 20)
  } else {
    ggsave(paste0(directory_hits_nascent,'/mods_RBPs/mods_frequency_nascent_total_DRACHneg.pdf'), plot = p, height = 20, width = 20)
  }
}

load("/path/to/R_data/protein_coding_genes_5UTR_3UTR_introns_exons_stop.Rda")


####
# SUM159TP: nascent + pre-existent RNAs 

load('/path/to/folder_random_hits_cluster/chr_ass_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster/chr_ass_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_chr <- post_process_mods(chr_ass_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],chr_ass_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                  c(6.96,1.05,2.1,2.85,2.01,0.25,0.84,14.89), c(86.5,0.24,1.08,1.2,2.22,0.3,1.5,87.88),
                                  'violin_hit_per_mod_bothBD_chr.pdf','chromatin',
                                  paste(as.character(unique(lapply(chr_ass_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(chr_ass_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                  '/path/to/folder_random_hits_cluster/')

pvalues_without_DRACH_chr <- RNA_mods_chr[[1]]
pvalues_with_DRACH_chr <- RNA_mods_chr[[2]]
significant_RNA_mods_without_DRACH_chr <- RNA_mods_chr[[3]]
significant_RNA_mods_with_DRACH_chr <- RNA_mods_chr[[4]]

load('/path/to/folder_random_hits_cluster/nucleo_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster/nucleo_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_nucleo <- post_process_mods(nucleo_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],nucleo_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                     c(6.87,1.14,2.89,3.08,1.94,0.43,0.85,16.15), c(87.06,0.31,0.8,0.99,2.65,0.25,1.29,88.29),
                                     'violin_hit_per_mod_bothBD_nucleo.pdf','nucleoplasm',
                                     paste(as.character(unique(lapply(nucleo_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(nucleo_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                     '/path/to/folder_random_hits_cluster/')

pvalues_without_DRACH_nucleo <- RNA_mods_nucleo[[1]]
pvalues_with_DRACH_nucleo <- RNA_mods_nucleo[[2]]
significant_RNA_mods_without_DRACH_nucleo <- RNA_mods_nucleo[[3]]
significant_RNA_mods_with_DRACH_nucleo <- RNA_mods_nucleo[[4]]

load('/path/to/folder_random_hits_cluster/cyto_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster/cyto_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_cyto <- post_process_mods(cyto_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],cyto_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                   c(7.21,1,2.28,2.97,1.62,0.19,0.77,15.12), c(85.16,0.46,1.03,0.97,1.95,0.29,1.43,86.36),
                                   'violin_hit_per_mod_bothBD_cyto.pdf','cytoplasm',
                                   paste(as.character(unique(lapply(cyto_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(cyto_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                   '/path/to/folder_random_hits_cluster/')

pvalues_without_DRACH_cyto <- RNA_mods_cyto[[1]]
pvalues_with_DRACH_cyto <- RNA_mods_cyto[[2]]
significant_RNA_mods_without_DRACH_cyto <- RNA_mods_cyto[[3]]
significant_RNA_mods_with_DRACH_cyto <- RNA_mods_cyto[[4]]

print_heatmap(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/', 
              significant_RNA_mods_without_DRACH_chr, significant_RNA_mods_with_DRACH_chr, significant_RNA_mods_without_DRACH_nucleo, 
              significant_RNA_mods_with_DRACH_nucleo, significant_RNA_mods_without_DRACH_cyto, significant_RNA_mods_with_DRACH_cyto)

heatmap_enrichment(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/',
                   c(86.5,0.24,1.08,1.2,2.22,0.3,1.5), c(87.06,0.31,0.8,0.99,2.65,0.25,1.29), c(85.16,0.46,1.03,0.97,1.95,0.29,1.43),
                   pvalues_with_DRACH_chr, pvalues_with_DRACH_nucleo, pvalues_with_DRACH_cyto, TRUE)

heatmap_enrichment(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/',
                   c(6.96,1.05,2.1,2.85,2.01,0.25,0.84),c(6.87,1.14,2.89,3.08,1.94,0.43,0.85),c(7.21,1,2.28,2.97,1.62,0.19,0.77),
                   pvalues_without_DRACH_chr, pvalues_without_DRACH_nucleo, pvalues_without_DRACH_cyto, FALSE)


# both DRACH+ and DRACH- hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/',
                                  'both', significant_RNA_mods_without_DRACH_chr, significant_RNA_mods_with_DRACH_chr, significant_RNA_mods_without_DRACH_nucleo, 
                                  significant_RNA_mods_with_DRACH_nucleo, significant_RNA_mods_without_DRACH_cyto, significant_RNA_mods_with_DRACH_cyto)

# only DRACH+ hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/',
                                  'DRACH', significant_RNA_mods_without_DRACH_chr, significant_RNA_mods_with_DRACH_chr, significant_RNA_mods_without_DRACH_nucleo, 
                                  significant_RNA_mods_with_DRACH_nucleo, significant_RNA_mods_without_DRACH_cyto, significant_RNA_mods_with_DRACH_cyto)

# only DRACH- hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/',
                                  'non_DRACH', significant_RNA_mods_without_DRACH_chr, significant_RNA_mods_with_DRACH_chr, significant_RNA_mods_without_DRACH_nucleo, 
                                  significant_RNA_mods_with_DRACH_nucleo, significant_RNA_mods_without_DRACH_cyto, significant_RNA_mods_with_DRACH_cyto)

#################
# SUM159TP: nascent RNAs 

load('/path/to/folder_random_hits_cluster_nascent/chr_ass_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_nascent/chr_ass_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_chr_nascent <- post_process_mods(chr_ass_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],chr_ass_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                          c(9.13,0.96,1.98,3.15,1.62,0.41,1.62,7.55), c(71.34,0.25,1.14,1.89,3.41,0.38,3.54,74.62),
                                          'violin_hit_per_mod_bothBD_chr.pdf','chromatin',
                                          paste(as.character(unique(lapply(chr_ass_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(chr_ass_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                          '/path/to/folder_random_hits_cluster_nascent/')

pvalues_without_DRACH_chr_nascent <- RNA_mods_chr_nascent[[1]]
pvalues_with_DRACH_chr_nascent <- RNA_mods_chr[[2]]
significant_RNA_mods_without_DRACH_chr_nascent <- RNA_mods_chr_nascent[[3]]
significant_RNA_mods_with_DRACH_chr_nascent <- RNA_mods_chr_nascent[[4]]

load('/path/to/folder_random_hits_cluster_nascent/nucleo_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_nascent/nucleo_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_nucleo_nascent <- post_process_mods(nucleo_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],nucleo_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                             c(9.56,1.47,2.57,3.58,1.84,0.46,1.65,19.39), c(74.31,0.4,0.79,2.57,4.15,0.59,4.94,77.27),
                                             'violin_hit_per_mod_bothBD_nucleo.pdf','nucleoplasm',
                                             paste(as.character(unique(lapply(nucleo_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(nucleo_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                             '/path/to/folder_random_hits_cluster_nascent/')

pvalues_without_DRACH_nucleo_nascent <- RNA_mods_nucleo_nascent[[1]]
pvalues_with_DRACH_nucleo_nascent <- RNA_mods_nucleo_nascent[[2]]
significant_RNA_mods_without_DRACH_nucleo_nascent <- RNA_mods_nucleo_nascent[[3]]
significant_RNA_mods_with_DRACH_nucleo_nascent <- RNA_mods_nucleo_nascent[[4]]

load('/path/to/folder_random_hits_cluster_nascent/cyto_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_nascent/cyto_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_cyto_nascent <- post_process_mods(cyto_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],cyto_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                           c(7.72,1.08,2.17,3.86,1.83,0.14,1.56,16.87), c(71.07,0.68,1.02,2.2,2.2,1.18,3.89,73.6),
                                           'violin_hit_per_mod_bothBD_cyto.pdf','cytoplasm',
                                           paste(as.character(unique(lapply(cyto_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(cyto_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                           '/path/to/folder_random_hits_cluster_nascent/')

pvalues_without_DRACH_cyto_nascent <- RNA_mods_cyto_nascent[[1]]
pvalues_with_DRACH_cyto_nascent <- RNA_mods_cyto_nascent[[2]]
significant_RNA_mods_without_DRACH_cyto_nascent <- RNA_mods_cyto_nascent[[3]]
significant_RNA_mods_with_DRACH_cyto_nascent <- RNA_mods_cyto_nascent[[4]]

print_heatmap(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/',
              significant_RNA_mods_without_DRACH_chr_nascent, significant_RNA_mods_with_DRACH_chr_nascent, significant_RNA_mods_without_DRACH_nucleo_nascent, 
              significant_RNA_mods_with_DRACH_nucleo_nascent, significant_RNA_mods_without_DRACH_cyto_nascent, significant_RNA_mods_with_DRACH_cyto_nascent)

heatmap_enrichment(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/',
                   c(71.34,0.25,1.14,1.89,3.41,0.38,3.54), c(74.31,0.4,0.79,2.57,4.15,0.59,4.94), c(71.07,0.68,1.02,2.2,2.2,1.18,3.89),
                   pvalues_with_DRACH_chr_nascent, pvalues_with_DRACH_nucleo_nascent, pvalues_with_DRACH_cyto_nascent, TRUE)

heatmap_enrichment(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/',
                   c(9.13,0.96,1.98,3.15,1.62,0.41,1.62),c(9.56,1.47,2.57,3.58,1.84,0.46,1.65),c(7.72,1.08,2.17,3.86,1.83,0.14,1.56),
                   pvalues_without_DRACH_chr_nascent, pvalues_without_DRACH_nucleo_nascent, pvalues_without_DRACH_cyto_nascent, FALSE)

# both DRACH+ and DRACH- hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/',
                                  'both', significant_RNA_mods_without_DRACH_chr_nascent, significant_RNA_mods_with_DRACH_chr_nascent, significant_RNA_mods_without_DRACH_nucleo_nascent, 
                                  significant_RNA_mods_with_DRACH_nucleo_nascent, significant_RNA_mods_without_DRACH_cyto_nascent, significant_RNA_mods_with_DRACH_cyto_nascent)
# only DRACH+ hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/',
                                  'DRACH', significant_RNA_mods_without_DRACH_chr_nascent, significant_RNA_mods_with_DRACH_chr_nascent, significant_RNA_mods_without_DRACH_nucleo_nascent, 
                                  significant_RNA_mods_with_DRACH_nucleo_nascent, significant_RNA_mods_without_DRACH_cyto_nascent, significant_RNA_mods_with_DRACH_cyto_nascent)
# only DRACH- hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/',
                                  'non_DRACH', significant_RNA_mods_without_DRACH_chr_nascent, significant_RNA_mods_with_DRACH_chr_nascent, significant_RNA_mods_without_DRACH_nucleo_nascent, 
                                  significant_RNA_mods_with_DRACH_nucleo_nascent, significant_RNA_mods_without_DRACH_cyto_nascent, significant_RNA_mods_with_DRACH_cyto_nascent)

#################
# SUM159TP: nascent + pre-existent RNAs with library-level subsampling threshold used for nascent reads

load('/path/to/folder_random_hits_cluster_tot2/chr_ass_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_tot2/chr_ass_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_chr_tot2 <- post_process_mods(chr_ass_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],chr_ass_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                       c(9.18,1.41,2.81,3.09,2.16,0.19,1.22,18.84), c(84.04,0.3,1.37,1.37,3.65,0.3,2.58,86.17),
                                       'violin_hit_per_mod_bothBD_chr.pdf','chromatin',
                                       paste(as.character(unique(lapply(chr_ass_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(chr_ass_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                       '/path/to/folder_random_hits_cluster_tot2/')

pvalues_without_DRACH_chr_tot2 <- RNA_mods_chr_tot2[[1]]
pvalues_with_DRACH_chr_tot2 <- RNA_mods_chr_tot2[[2]]
significant_RNA_mods_without_DRACH_chr_tot2 <- RNA_mods_chr_tot2[[3]]
significant_RNA_mods_with_DRACH_chr_tot2 <- RNA_mods_chr_tot2[[4]]

load('/path/to/folder_random_hits_cluster_tot2/nucleo_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_tot2/nucleo_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_nucleo_tot2 <- post_process_mods(nucleo_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],nucleo_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                          c(6.82,1.86,3.23,3.23,2.23,0.37,0.62,17.37), c(87.73,0.72,0.54,1.08,3.25,0.36,2.35,88.81),
                                          'violin_hit_per_mod_bothBD_nucleo.pdf','nucleoplasm',
                                          paste(as.character(unique(lapply(nucleo_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(nucleo_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                          '/path/to/folder_random_hits_cluster_tot2/')

pvalues_without_DRACH_nucleo_tot2 <- RNA_mods_nucleo_tot2[[1]]
pvalues_with_DRACH_nucleo_tot2 <- RNA_mods_nucleo_tot2[[2]]
significant_RNA_mods_without_DRACH_nucleo_tot2 <- RNA_mods_nucleo_tot2[[3]]
significant_RNA_mods_with_DRACH_nucleo_tot2 <- RNA_mods_nucleo_tot2[[4]]

load('/path/to/folder_random_hits_cluster_tot2/cyto_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_tot2/cyto_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_cyto_tot2 <- post_process_mods(cyto_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],cyto_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                        c(6.9,1.72,2.78,4.12,1.82,0.38,0.77,17.15), c(84.63,0.62,0.78,1.55,3.11,1.24,2.95,86.49),
                                        'violin_hit_per_mod_bothBD_cyto.pdf','cytoplasm',
                                        paste(as.character(unique(lapply(cyto_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(cyto_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                        '/path/to/folder_random_hits_cluster_tot2/')

pvalues_without_DRACH_cyto_tot2 <- RNA_mods_cyto_tot2[[1]]
pvalues_with_DRACH_cyto_tot2 <- RNA_mods_cyto_tot2[[2]]
significant_RNA_mods_without_DRACH_cyto_tot2 <- RNA_mods_cyto_tot2[[3]]
significant_RNA_mods_with_DRACH_cyto_tot2 <- RNA_mods_cyto_tot2[[4]]

print_heatmap(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/',
              significant_RNA_mods_without_DRACH_chr_tot2, significant_RNA_mods_with_DRACH_chr_tot2, significant_RNA_mods_without_DRACH_nucleo_tot2, 
              significant_RNA_mods_with_DRACH_nucleo_tot2, significant_RNA_mods_without_DRACH_cyto_tot2, significant_RNA_mods_with_DRACH_cyto_tot2)

heatmap_enrichment(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/',
                   c(84.04,0.3,1.37,1.37,3.65,0.3,2.58), c(87.73,0.72,0.54,1.08,3.25,0.36,2.35), c(84.63,0.62,0.78,1.55,3.11,1.24,2.95),
                   pvalues_with_DRACH_chr_tot2, pvalues_with_DRACH_nucleo_tot2, pvalues_with_DRACH_cyto_tot2, TRUE)

heatmap_enrichment(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/',
                   c(9.18,1.41,2.81,3.09,2.16,0.19,1.22),c(6.82,1.86,3.23,3.23,2.23,0.37,0.62),c(6.9,1.72,2.78,4.12,1.82,0.38,0.77),
                   pvalues_without_DRACH_chr_tot2, pvalues_without_DRACH_nucleo_tot2, pvalues_without_DRACH_cyto_tot2, FALSE)

# both DRACH+ and DRACH- hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/',
                                  'both',significant_RNA_mods_without_DRACH_chr_tot2, significant_RNA_mods_with_DRACH_chr_tot2, significant_RNA_mods_without_DRACH_nucleo_tot2, 
                                  significant_RNA_mods_with_DRACH_nucleo_tot2, significant_RNA_mods_without_DRACH_cyto_tot2, significant_RNA_mods_with_DRACH_cyto_tot2)
# only DRACH+ hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/',
                                  'DRACH',significant_RNA_mods_without_DRACH_chr_tot2, significant_RNA_mods_with_DRACH_chr_tot2, significant_RNA_mods_without_DRACH_nucleo_tot2, 
                                  significant_RNA_mods_with_DRACH_nucleo_tot2, significant_RNA_mods_without_DRACH_cyto_tot2, significant_RNA_mods_with_DRACH_cyto_tot2)
# only DRACH- hits
mods_combination_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/',
                                  'non_DRACH',significant_RNA_mods_without_DRACH_chr_tot2, significant_RNA_mods_with_DRACH_chr_tot2, significant_RNA_mods_without_DRACH_nucleo_tot2, 
                                  significant_RNA_mods_with_DRACH_nucleo_tot2, significant_RNA_mods_without_DRACH_cyto_tot2, significant_RNA_mods_with_DRACH_cyto_tot2)


# both DRACH- and DRACH+ hits
comparison_mods_combination_nascent_total_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/', 
                                                           '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/', 'both',
                                                           significant_RNA_mods_without_DRACH_chr_nascent, significant_RNA_mods_with_DRACH_chr_nascent, significant_RNA_mods_without_DRACH_nucleo_nascent, 
                                                           significant_RNA_mods_with_DRACH_nucleo_nascent, significant_RNA_mods_without_DRACH_cyto_nascent, significant_RNA_mods_with_DRACH_cyto_nascent,
                                                           significant_RNA_mods_without_DRACH_chr_tot2, significant_RNA_mods_with_DRACH_chr_tot2, significant_RNA_mods_without_DRACH_nucleo_tot2, 
                                                           significant_RNA_mods_with_DRACH_nucleo_tot2, significant_RNA_mods_without_DRACH_cyto_tot2, significant_RNA_mods_with_DRACH_cyto_tot2)

# only DRACH+ hits
comparison_mods_combination_nascent_total_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/', 
                                                           '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/', 'DRACH',
                                                           significant_RNA_mods_without_DRACH_chr_nascent, significant_RNA_mods_with_DRACH_chr_nascent, significant_RNA_mods_without_DRACH_nucleo_nascent, 
                                                           significant_RNA_mods_with_DRACH_nucleo_nascent, significant_RNA_mods_without_DRACH_cyto_nascent, significant_RNA_mods_with_DRACH_cyto_nascent,
                                                           significant_RNA_mods_without_DRACH_chr_tot2, significant_RNA_mods_with_DRACH_chr_tot2, significant_RNA_mods_without_DRACH_nucleo_tot2, 
                                                           significant_RNA_mods_with_DRACH_nucleo_tot2, significant_RNA_mods_without_DRACH_cyto_tot2, significant_RNA_mods_with_DRACH_cyto_tot2)

# only DRACH- hits
comparison_mods_combination_nascent_total_only_significant('/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/', 
                                                           '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/', 'non_DRACH',
                                                           significant_RNA_mods_without_DRACH_chr_nascent, significant_RNA_mods_with_DRACH_chr_nascent, significant_RNA_mods_without_DRACH_nucleo_nascent, 
                                                           significant_RNA_mods_with_DRACH_nucleo_nascent, significant_RNA_mods_without_DRACH_cyto_nascent, significant_RNA_mods_with_DRACH_cyto_nascent,
                                                           significant_RNA_mods_without_DRACH_chr_tot2, significant_RNA_mods_with_DRACH_chr_tot2, significant_RNA_mods_without_DRACH_nucleo_tot2, 
                                                           significant_RNA_mods_with_DRACH_nucleo_tot2, significant_RNA_mods_without_DRACH_cyto_tot2, significant_RNA_mods_with_DRACH_cyto_tot2)

###################
# K562:nascent + pre-existent RNAs 

load('/path/to/folder_random_hits_cluster_K562/chr_ass_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_K562/chr_ass_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_chr <- post_process_mods(chr_ass_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],chr_ass_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                  c(10.45,1.79,2.84,4.48,2.24,0.15,0.9,20.9), c(83.74,0.69,1.38,2.08,2.42,0.69,3.81,85.47),
                                  'violin_hit_per_mod_bothBD_chr.pdf','chromatin',
                                  paste(as.character(unique(lapply(chr_ass_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(chr_ass_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                  '/path/to/folder_random_hits_cluster_K562/')

pvalues_without_DRACH_chr <- RNA_mods_chr[[1]]
pvalues_with_DRACH_chr <- RNA_mods_chr[[2]]
significant_RNA_mods_without_DRACH_chr <- RNA_mods_chr[[3]]
significant_RNA_mods_with_DRACH_chr <- RNA_mods_chr[[4]]

load('/path/to/folder_random_hits_cluster_K562/nucleo_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_K562/nucleo_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_nucleo <- post_process_mods(nucleo_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],nucleo_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                     c(7.43,1.18,2.87,5.24,1.52,0.34,1.18,17.74), c(76.63,1.03,1.37,4.12,1.37,1.03,2.06,78.35),
                                     'violin_hit_per_mod_bothBD_nucleo.pdf','nucleoplasm',
                                     paste(as.character(unique(lapply(nucleo_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(nucleo_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                     '/path/to/folder_random_hits_cluster_K562/')

pvalues_without_DRACH_nucleo <- RNA_mods_nucleo[[1]]
pvalues_with_DRACH_nucleo <- RNA_mods_nucleo[[2]]
significant_RNA_mods_without_DRACH_nucleo <- RNA_mods_nucleo[[3]]
significant_RNA_mods_with_DRACH_nucleo <- RNA_mods_nucleo[[4]]

load('/path/to/folder_random_hits_cluster_K562/cyto_mod_type_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_K562/cyto_mod_type_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each RNA mod type and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
RNA_mods_cyto <- post_process_mods(cyto_mod_type_without_DRACH[c(1,2,3,4,5,6,7,11)],cyto_mod_type_DRACH[c(1,2,3,4,5,6,7,11)],
                                   c(8.55,1.57,2.56,4.42,1.42,0.28,0.43,17.81), c(73.61,1.11,1.11,2.5,1.11,1.67,3.33,75.28),
                                   'violin_hit_per_mod_bothBD_cyto.pdf','cytoplasm',
                                   paste(as.character(unique(lapply(cyto_mod_type_DRACH[[12]], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(cyto_mod_type_without_DRACH[[12]], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                                   '/path/to/folder_random_hits_cluster_K562/')

pvalues_without_DRACH_cyto <- RNA_mods_cyto[[1]]
pvalues_with_DRACH_cyto <- RNA_mods_cyto[[2]]
significant_RNA_mods_without_DRACH_cyto <- RNA_mods_cyto[[3]]
significant_RNA_mods_with_DRACH_cyto <- RNA_mods_cyto[[4]]

print_heatmap(directory_hits = '/path/to/fractions_eligos_K562_all_reads/', 
              significant_RNA_mods_without_DRACH_chr, significant_RNA_mods_with_DRACH_chr, significant_RNA_mods_without_DRACH_nucleo, 
              significant_RNA_mods_with_DRACH_nucleo, significant_RNA_mods_without_DRACH_cyto, significant_RNA_mods_with_DRACH_cyto)

heatmap_enrichment(directory_hits = '/path/to/fractions_eligos_K562_all_reads/',
                   c(83.74,0.69,1.38,2.08,2.42,0.69,3.81),c(76.63,1.03,1.37,4.12,1.37,1.03,2.06),c(73.61,1.11,1.11,2.5,1.11,1.67,3.33),
                   pvalues_with_DRACH_chr, pvalues_with_DRACH_nucleo, pvalues_with_DRACH_cyto, TRUE)

heatmap_enrichment(directory_hits = '/path/to/fractions_eligos_K562_all_reads/',
                   c(10.45,1.79,2.84,4.48,2.24,0.15,0.9),c(7.43,1.18,2.87,5.24,1.52,0.34,1.18), c(8.55,1.57,2.56,4.42,1.42,0.28,0.43),
                   pvalues_without_DRACH_chr, pvalues_without_DRACH_nucleo, pvalues_without_DRACH_cyto, FALSE)


# both DRACH+ and DRACH- hits
mods_combination_only_significant('/path/to/fractions_eligos_K562_all_reads/',
                                  'both', significant_RNA_mods_without_DRACH_chr, significant_RNA_mods_with_DRACH_chr, significant_RNA_mods_without_DRACH_nucleo, 
                                  significant_RNA_mods_with_DRACH_nucleo, significant_RNA_mods_without_DRACH_cyto, significant_RNA_mods_with_DRACH_cyto)

# only DRACH+ hits
mods_combination_only_significant('/path/to/fractions_eligos_K562_all_reads/',
                                  'DRACH', significant_RNA_mods_without_DRACH_chr, significant_RNA_mods_with_DRACH_chr, significant_RNA_mods_without_DRACH_nucleo, 
                                  significant_RNA_mods_with_DRACH_nucleo, significant_RNA_mods_without_DRACH_cyto, significant_RNA_mods_with_DRACH_cyto)

# only DRACH- hits
mods_combination_only_significant('/path/to/fractions_eligos_K562_all_reads/',
                                  'non_DRACH', significant_RNA_mods_without_DRACH_chr, significant_RNA_mods_with_DRACH_chr, significant_RNA_mods_without_DRACH_nucleo, 
                                  significant_RNA_mods_with_DRACH_nucleo, significant_RNA_mods_without_DRACH_cyto, significant_RNA_mods_with_DRACH_cyto)

#################

# give in input, for each fraction, the 1,000 sets of DRACH- random hits (hits_effector_without_DRACH) and 
# the 1,000 sets of DRACH+ random hits (hits_effector_with_DRACH) (both overlapped
# with the coordinates of RNA marks+effectors' binding sites from the databases) and the percentage
# of ELIGOS DRACH-/DRACH+ hits overlapping with each category of effectors from the DB for the same fraction
# (ELIGOS_results_without_DRACH/ELIGOS_results_with_DRACH).
# Path is the path to the directory in which saving the plot
post_process_effectors <- function(hits_effector_without_DRACH, hits_effector_with_DRACH, ELIGOS_results_without_DRACH, ELIGOS_results_with_DRACH, name_pdf,t,st, path) {
  hits_effector <- c(hits_effector_without_DRACH[2],hits_effector_with_DRACH[2],
                     hits_effector_without_DRACH[1],hits_effector_with_DRACH[1],
                     hits_effector_without_DRACH[3],hits_effector_with_DRACH[3],
                     hits_effector_without_DRACH[4],hits_effector_with_DRACH[4])
  names(hits_effector) <- rep(c('% of hits overlapping with\nspecific effectors non m6A',
                                '% of hits overlapping with\nspecific effectors m6A',
                                '% of hits overlapping with\nspecific effectors',
                                '% of hits overlapping with\ntotal effectors'),each=2)
  
  names(ELIGOS_results_without_DRACH) <- c('% of hits overlapping with\nspecific effectors non m6A',
                                           '% of hits overlapping with\nspecific effectors m6A',
                                           '% of hits overlapping with\nspecific effectors',
                                           '% of hits overlapping with\ntotal effectors')
  names(ELIGOS_results_with_DRACH) <- c('% of hits overlapping with\nspecific effectors non m6A',
                                        '% of hits overlapping with\nspecific effectors m6A',
                                        '% of hits overlapping with\nspecific effectors',
                                        '% of hits overlapping with\ntotal effectors')
  
  # for each category of effectors, compute the absolute 
  # difference between the percentage of ELIGOS DRACH- hits overlapping with each category and the median percentage 
  # of DRACH- random hits overlapping with the same category divided by the interquantile distance computed on the 1,000 values of percentage of 
  # DRACH- random hits overlapping with the same category
  IQR_distance_without_DRACH <- unlist(lapply(seq_along(c(1,3,5,7)), function(x,n,i) {
    if (length(x[[i]]) != 1000) {print('error')}
    IQR <- summary(as.numeric(x[[i]]))[5] - summary(as.numeric(x[[i]]))[2]
    round(abs((ELIGOS_results_without_DRACH[names(ELIGOS_results_without_DRACH) == n[i]] - summary(as.numeric(x[[i]]))[3]))/IQR,2)
  },  x=hits_effector[c(1,3,5,7)], n=names(hits_effector[c(1,3,5,7)])))
  
  # the same for random DRACH+ hits
  IQR_distance_with_DRACH <- unlist(lapply(seq_along(c(2,4,6,8)), function(x,n,i) {
    if (length(x[[i]]) != 1000) {print('error')}
    IQR <- summary(as.numeric(x[[i]]))[5] - summary(as.numeric(x[[i]]))[2]
    round(abs((ELIGOS_results_with_DRACH[names(ELIGOS_results_with_DRACH) == n[i]] - summary(as.numeric(x[[i]]))[3]))/IQR,2)
  },  x=hits_effector[c(2,4,6,8)], n=names(hits_effector[c(2,4,6,8)])))
  
  # for each category of effectors, compute how many of the 1,000
  # values of percentage of DRACH- random hits overlapping with each category are higher than the percentage of ELIGOS DRACH- hits 
  # overlapping with the same category and then divide by 1,000
  p_value_without_DRACH <- unlist(lapply(seq_along(c(1,3,5,7)), function(x,n,i) {
    if (length(x[[i]][x[[i]]>ELIGOS_results_without_DRACH[names(ELIGOS_results_without_DRACH) == n[i]]])/1000 == 0) {
      paste0('< ',as.character(1e-3))
    } else{
      as.character(length(x[[i]][x[[i]]>ELIGOS_results_without_DRACH[names(ELIGOS_results_without_DRACH) == n[i]]])/1000)
    } },  x=hits_effector[c(1,3,5,7)], n=names(hits_effector[c(1,3,5,7)])))
  
  # the same for random DRACH+ hits
  p_value_with_DRACH <- unlist(lapply(seq_along(c(2,4,6,8)), function(x,n,i) {
    if (length(x[[i]][x[[i]]>ELIGOS_results_with_DRACH[names(ELIGOS_results_with_DRACH) == n[i]]])/1000 == 0) {
      paste0('< ',as.character(1e-3))
    } else{
      as.character(length(x[[i]][x[[i]]>ELIGOS_results_with_DRACH[names(ELIGOS_results_with_DRACH) == n[i]]])/1000)
    } },  x=hits_effector[c(2,4,6,8)], n=names(hits_effector[c(2,4,6,8)])))
  
  # collapse the two statistics
  stats_without_DRACH <- unlist(lapply(seq_along(1:4), function(i) {
    paste0(IQR_distance_without_DRACH[i],'\n',p_value_without_DRACH[i])
  }))
  
  stats_with_DRACH <- unlist(lapply(seq_along(1:4), function(i) {
    paste0(IQR_distance_with_DRACH[i],'\n',p_value_with_DRACH[i])
  }))
  
  # create a data frame reporting for each category of effectors the 1,000 values of percentage of DRACH-/DRACH+ random 
  # hits overlapping with that category
  df <- data.frame(
    EFFECTORS = factor(rep(c('% of hits overlapping with\nspecific effectors non m6A',
                             '% of hits overlapping with\nspecific effectors m6A',
                             '% of hits overlapping with\nspecific effectors',
                             '% of hits overlapping with\ntotal effectors'), each=2000), levels=c('% of hits overlapping with\nspecific effectors non m6A', 
                                                                                                  '% of hits overlapping with\nspecific effectors m6A',
                                                                                                  '% of hits overlapping with\nspecific effectors',
                                                                                                  '% of hits overlapping with\ntotal effectors')),
    Class = factor(rep(c(rep('DRACH negative',1000), rep('DRACH positive', 1000)),4),levels=c('DRACH positive','DRACH negative')),
    value = unlist(lapply(hits_effector, function(x) {x}))
  )
  
  # create a data frame reporting for each category of effectors the percentage of ELIGOS DRACH-/DRACH+
  # hits overlapping with that category and the relative statistics
  points <- data.frame(EFFECTORS = factor(rep(c('% of hits overlapping with\nspecific effectors non m6A',
                                                '% of hits overlapping with\nspecific effectors m6A',
                                                '% of hits overlapping with\nspecific effectors',
                                                '% of hits overlapping with\ntotal effectors'),each = 2), levels=c('% of hits overlapping with\nspecific effectors non m6A',
                                                                                                                   '% of hits overlapping with\nspecific effectors m6A',
                                                                                                                   '% of hits overlapping with\nspecific effectors',
                                                                                                                   '% of hits overlapping with\ntotal effectors')),
                       ELIGOS = c(ELIGOS_results_without_DRACH[1],ELIGOS_results_with_DRACH[1],
                                  ELIGOS_results_without_DRACH[2],ELIGOS_results_with_DRACH[2],
                                  ELIGOS_results_without_DRACH[3],ELIGOS_results_with_DRACH[3],
                                  ELIGOS_results_without_DRACH[4],ELIGOS_results_with_DRACH[4]),
                       stats_without= c(stats_without_DRACH[1],NA,stats_without_DRACH[2],NA,stats_without_DRACH[3],NA, stats_without_DRACH[4],NA),
                       stats_with= c(NA,stats_with_DRACH[1],NA,stats_with_DRACH[2],NA,stats_with_DRACH[3],NA,stats_with_DRACH[4]),
                       Class = factor(rep(c('DRACH negative', 'DRACH positive'),4), levels=c('DRACH positive','DRACH negative'))
  )
  
  # generate a violin plot representing the distribution of the 1,000 values of percentage of DRACH+/DRACH- random hits
  # overlapping with each category of effectors. A point is used to report the percentage of ELIGOS DRACH+/DRACH- 
  # hits overlapping with the same category and the relative statistics
  p <- ggplot(data = points, aes(x=EFFECTORS, y=ELIGOS, fill=Class)) + 
    geom_point(size=12, shape=23,stroke = 2.5, position=position_dodge(width=0.9))+ 
    scale_fill_manual(name="ELIGOS values",values=c('coral2','cornflowerblue')) +
    geom_text(data = points,aes(x=EFFECTORS,y=ELIGOS,label=stats_without), position = position_dodge(width=1.6), size=15) +
    geom_text(data = points,aes(x=EFFECTORS,y=ELIGOS,label=stats_with), position = position_dodge(width=0.3), size=15) +
    new_scale_fill()+
    geom_violin(data = df, mapping = aes(x=EFFECTORS, y=value, fill = Class),position=position_dodge(width = 0.9)) +
    geom_boxplot(data = df, mapping = aes(x=EFFECTORS, y=value, fill = Class),position = position_dodge(width = 0.9),width=0.1, color='black') +
    scale_fill_manual(name="Random sequences",values=c('coral2','cornflowerblue')) +
    xlab('') + 
    ylab('')+
    labs(title=t, subtitle = st) +
    theme_classic() +  
    theme(legend.position = 'bottom', axis.text.x = element_text(size=40),
          axis.text.y = element_text(size=35), plot.title = element_text(size=40,hjust = 0.5), plot.subtitle = element_text(size=35,hjust = 0.5), legend.title = element_text(size=45), legend.text =  element_text(size=40))
  
  ggsave(paste0(path, name_pdf), plot=p,width = 35, height = 20)
}


##########
# SUM159TP: nascent + pre-existent RNAs 

load('/path/to/folder_random_hits_cluster/hits_chr_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster/hits_chr_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_chr_mod_RBP_without_DRACH,hits_chr_mod_RBP_DRACH, c(18.85,63.72,65.13,88.35), c(27.08,88.72,89.51,97.64), 
                       'violin_hit_effectors_bothBD_chr.pdf', 'chromatin',
                       paste(as.character(unique(lapply(hits_chr_mod_RBP_DRACH[5:length(hits_chr_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_chr_mod_RBP_without_DRACH[5:length(hits_chr_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster/')

load('/path/to/folder_random_hits_cluster/hits_nucleo_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster/hits_nucleo_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_nucleo_mod_RBP_without_DRACH,hits_nucleo_mod_RBP_DRACH, c(18.91,64.51,65.8,88.39), c(27.08,87.15,87.87,97.18), 
                       'violin_hit_effectors_bothBD_nucleo.pdf', 'nucleoplasm',
                       paste(as.character(unique(lapply(hits_nucleo_mod_RBP_DRACH[5:length(hits_nucleo_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_nucleo_mod_RBP_without_DRACH[5:length(hits_nucleo_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster/')

load('/path/to/folder_random_hits_cluster/hits_cyto_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster/hits_cyto_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_cyto_mod_RBP_without_DRACH,hits_cyto_mod_RBP_DRACH, c(18.72,64.29,65.59,87.86), c(26.55,87.55,88.41,97.56), 
                       'violin_hit_effectors_bothBD_cyto.pdf', 'cytoplasm',
                       paste(as.character(unique(lapply(hits_cyto_mod_RBP_DRACH[5:length(hits_cyto_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_cyto_mod_RBP_without_DRACH[5:length(hits_cyto_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster/')

######
# SUM159TP: nascent RNAs 

load('/path/to/folder_random_hits_cluster_nascent/hits_chr_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_nascent/hits_chr_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_chr_mod_RBP_without_DRACH,hits_chr_mod_RBP_DRACH, c(25.67,70.72,72.29,92.35), c(34.67,85.39,86.53,97.42), 
                       'violin_hit_effectors_bothBD_chr.pdf', 'chromatin',
                       paste(as.character(unique(lapply(hits_chr_mod_RBP_DRACH[5:length(hits_chr_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_chr_mod_RBP_without_DRACH[5:length(hits_chr_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster_nascent/')

load('/path/to/folder_random_hits_cluster_nascent/hits_nucleo_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_nascent/hits_nucleo_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_nucleo_mod_RBP_without_DRACH,hits_nucleo_mod_RBP_DRACH, c(27.77,69.82,71.83,92.29), c(38.66,87.97,88.49,97.25), 
                       'violin_hit_effectors_bothBD_nucleo.pdf', 'nucleoplasm',
                       paste(as.character(unique(lapply(hits_nucleo_mod_RBP_DRACH[5:length(hits_nucleo_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_nucleo_mod_RBP_without_DRACH[5:length(hits_nucleo_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster_nascent/')

load('/path/to/folder_random_hits_cluster_nascent/hits_cyto_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_nascent/hits_cyto_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_cyto_mod_RBP_without_DRACH,hits_cyto_mod_RBP_DRACH, c(26.95,70.27,71.69,92.27), c(34.96,86.96,88.04,96.38), 
                       'violin_hit_effectors_bothBD_cyto.pdf', 'cytoplasm',
                       paste(as.character(unique(lapply(hits_cyto_mod_RBP_DRACH[5:length(hits_cyto_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_cyto_mod_RBP_without_DRACH[5:length(hits_cyto_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster_nascent/')

######
# SUM159TP: nascent + pre-existent RNAs with library-level subsampling threshold used for nascent reads

load('/path/to/folder_random_hits_cluster_tot2/hits_chr_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_tot2/hits_chr_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_chr_mod_RBP_without_DRACH,hits_chr_mod_RBP_DRACH, c(23.01,70.51,72,91.97), c(31.1,90.39,90.75,97.69), 
                       'violin_hit_effectors_bothBD_chr.pdf', 'chromatin',
                       paste(as.character(unique(lapply(hits_chr_mod_RBP_DRACH[5:length(hits_chr_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_chr_mod_RBP_without_DRACH[5:length(hits_chr_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster_tot2/')

load('/path/to/folder_random_hits_cluster_tot2/hits_nucleo_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_tot2/hits_nucleo_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_nucleo_mod_RBP_without_DRACH,hits_nucleo_mod_RBP_DRACH, c(24.41,69.15,71.11,91.99), c(37.1,89.85,90.52,97.84), 
                       'violin_hit_effectors_bothBD_nucleo.pdf', 'nucleoplasm',
                       paste(as.character(unique(lapply(hits_nucleo_mod_RBP_DRACH[5:length(hits_nucleo_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_nucleo_mod_RBP_without_DRACH[5:length(hits_nucleo_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster_tot2/')

load('/path/to/folder_random_hits_cluster_tot2/hits_cyto_mod_RBP_DRACH.Rda')
load('/path/to/folder_random_hits_cluster_tot2/hits_cyto_mod_RBP_without_DRACH.Rda')

# plot the distribution of the percentages of DRACH+/DRACH- random hits annotated to each category of effectors and the corresponding percentage 
# of ELIGOS DRACH+/DRACH- hits
post_process_effectors(hits_cyto_mod_RBP_without_DRACH,hits_cyto_mod_RBP_DRACH, c(24.17,67.74,69.57,90.3), c(32.6,91.39,91.89,98.14), 
                       'violin_hit_effectors_bothBD_cyto.pdf', 'cytoplasm',
                       paste(as.character(unique(lapply(hits_cyto_mod_RBP_DRACH[5:length(hits_cyto_mod_RBP_DRACH)], length))),'DRACH+ random sequences (10 nt) generated 1,000 times;', as.character(unique(lapply(hits_cyto_mod_RBP_without_DRACH[5:length(hits_cyto_mod_RBP_without_DRACH)], length))),'DRACH- random sequences (10 nt) generated 1,000 times'),
                       path = '/path/to/folder_random_hits_cluster_tot2/')

#######
