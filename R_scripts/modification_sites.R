library('GenomicRanges')
library('GenomicFeatures')
library('ggplot2')
library('Biostrings')
library('pheatmap')
library('xlsx')
library('readr')
library('dplyr')

## DATABASES OF RNA MARKS

# Download from RMBase3 https://rna.sysu.edu.cn/rmbase3/download.php (RNA modifications session) the bed files with the coordinates of
# m6A, m1A, m5C, m7G, Y, Nm, A-I and other RNA modifications selecting Mammal group, Homo sapiens hg38 genome.
# Save these bed files in a folder called /mod/, rename the file with the coordinates of other RNA modifications "other_RNA_mods_total.bed"
# and extract from it the coordinates of m5U and m6Am  
others <- read.table(file = '/path/to/mod/other_RNA_mods_total.bed')
m5U <- others[which(others$V7 == 'm5U'),]
write.table(x = m5U, file = '/path/to/mod/m5U.bed',row.names = FALSE, col.names = FALSE)
m6am <- others[which(others$V7 == 'm6Am'),]
write.table(x = m6am, file = '/path/to/mod/m6am.bed',row.names = FALSE, col.names = FALSE)
others <- others[-which(others$V7 == 'm6Am' | others$V7 == 'm5U'),]

# remove 'chr' from the annotation of the chromosome on which each RNA mark maps and 
# consider only the RNA marks on chromosomes 1-22, X and MT
others$V1 <- gsub('chr','',others$V1)
others <- others[which(others$V1 %in% c(as.character(seq(1,22)),'X','M')),]
write.table(others, '/path/to/mod/other_RNA_mods_total.bed',row.names = FALSE, col.names = FALSE)

# Download from RMVar https://rmvar.renlab.org/download.html the txt files with the coordinates of
# m6A, Nm, A-I, m1A, m5C, m5U, m6Am, m7G and Y relative to Human.
# These txt files were saved in the same folder /mod/

# For each RNA modification type the bed file and the txt file were renamed as "modification_type.bed" or "modification_type.txt"

files <- list.files(path = '/path/to/mod', pattern = 'txt', full.names = TRUE)

# iterate over the RNA modification types
for (f in files) {
  
  # read the coordinates of that RNA modification type identified by RMVar
  suppressWarnings(RMvar <- read_tsv(file = f,  col_types = cols())[,c(4,6,7,8,20,23)])
  if ('chrY' %in% unique(RMvar$chromosome)) {
    RMvar <- RMvar[-which(RMvar$chromosome =='chrY'),]
  }
  
  mod_grange_RMvar <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=RMvar$chromosome)),
                              ranges = IRanges(start = RMvar$position, end=RMvar$position),
                              strand = Rle(RMvar$strand),
                              mod_type =  RMvar$modification_type,
                              cell_line = RMvar$modification_sample_information,
                              DB='RMVar')
  
  # the sites that are reported multiple times with the same chromosome, coordinates, strand and same cell line 
  # (there are other metadata provided by RMVar that differ among these sites) are kept only once
  ind_unique <- which(isUnique(paste0(seqnames(mod_grange_RMvar), "_", start(mod_grange_RMvar), "_", strand(mod_grange_RMvar), "_", mod_grange_RMvar$cell_line)))
  dup <- unique(mod_grange_RMvar[-ind_unique])
  x_unique <- sort(c(mod_grange_RMvar[ind_unique], dup))
  mod_grange_RMvar <- x_unique
  
  print(paste0(unique(RMvar$modification_type), ': ', as.character(length(mod_grange_RMvar))))
  
  # read the coordinates of the same RNA modification type identified by RMBase3
  bed_file <- gsub('txt', 'bed', f)
  RMBase3 <- read.table(file = bed_file, fill = TRUE)
  RMBase3 <- RMBase3[which(RMBase3[,1] %in% RMvar$chromosome),]
  
  mod_grange_RMBase3 <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=RMBase3$V1)),
                                ranges = IRanges(start = RMBase3$V2, end=RMBase3$V3),
                                strand = Rle(RMBase3$V6),
                                mod_type = RMBase3$V7,
                                cell_line = RMBase3$V12,
                                DB='RMBase3')
  
  print(paste0(unique(RMvar$modification_type), ': ', as.character(length(mod_grange_RMBase3))))
  
  # compute how many RNA marks are reported in both databases (keeping the original coordinates)
  all_granges <- list(mod_grange_RMvar, mod_grange_RMBase3)
  names(all_granges) <- c('RMvar', 'RMBase3')
  order <- order(c(length(mod_grange_RMvar), length(mod_grange_RMBase3)))
  
  overlap <- findOverlaps(all_granges[[order[1]]],all_granges[[order[2]]], type = 'any')
  overlap_query <- all_granges[[order[1]]][unique(queryHits(overlap))]
  overlap_subject <- all_granges[[order[2]]][unique(subjectHits(overlap))]
  
  # print the percentage of marks of one DB also present in the other
  print(paste0(unique(RMvar$modification_type), ': ', as.character(round(length(overlap_query)/length(all_granges[[order[1]]])*100, 2)),'% in ', names(all_granges)[[order[1]]]))
  print(paste0(unique(RMvar$modification_type), ': ', as.character(round(length(overlap_subject)/length(all_granges[[order[2]]])*100, 2)),'% in ', names(all_granges)[[order[2]]]))
  
  # compute how many RNA marks are reported in both databases (resizing the original coordinates to 10 nucleotides)
  mod_grange_RMvar <- resize(mod_grange_RMvar, 10, fix='center')
  mod_grange_RMBase3 <- resize(mod_grange_RMBase3, 10, fix='center')
  
  all_granges <- list(mod_grange_RMvar, mod_grange_RMBase3)
  names(all_granges) <- c('RMvar', 'RMBase3')
  order <- order(c(length(mod_grange_RMvar), length(mod_grange_RMBase3)))
  
  overlap <- findOverlaps(all_granges[[order[1]]],all_granges[[order[2]]], type = 'any')
  overlap_query <- all_granges[[order[1]]][unique(queryHits(overlap))]
  overlap_subject <- all_granges[[order[2]]][unique(subjectHits(overlap))]
  
  # print the percentage of marks of one DB also present in the other
  print(paste0(unique(RMvar$modification_type), ': ', as.character(round(length(overlap_query)/length(all_granges[[order[1]]])*100, 2)),'% in ', names(all_granges)[[order[1]]]))
  print(paste0(unique(RMvar$modification_type), ': ', as.character(round(length(overlap_subject)/length(all_granges[[order[2]]])*100, 2)),'% in ', names(all_granges)[[order[2]]]))
}

# iterate over the RNA modification types  
for (f in files) {
  
  # read the coordinates of that RNA modification type identified by RMVar
  suppressWarnings(RMvar <- read_tsv(file=f,  col_types = cols())[,c(4,6,7,8,20,23)])
  if ('chrY' %in% unique(RMvar$chromosome)) {
    RMvar <- RMvar[-which(RMvar$chromosome =='chrY'),]
  }
  mod_grange_RMvar <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=RMvar$chromosome)),
                              ranges = IRanges(start = RMvar$position, end=RMvar$position),
                              strand = Rle(RMvar$strand),
                              mod_type =  RMvar$modification_type,
                              cell_line = RMvar$modification_sample_information,
                              DB='RMVar')
  
  # the sites that are reported multiple times with the same chromosome, coordinates, strand and same cell line 
  # (there are other metadata provided by RMVar that differ among these sites) are kept only once
  ind_unique <- which(isUnique(paste0(seqnames(mod_grange_RMvar), "_", start(mod_grange_RMvar), "_", strand(mod_grange_RMvar), "_", mod_grange_RMvar$cell_line)))
  dup <- unique(mod_grange_RMvar[-ind_unique])
  x_unique <- sort(c(mod_grange_RMvar[ind_unique], dup))
  mod_grange_RMvar <- x_unique
  
  # read the coordinates of the same RNA modification type identified by RMBase3
  bed_file <- gsub('txt', 'bed', f)
  RMBase3 <- read.table(file = bed_file, fill=TRUE)
  RMBase3 <- RMBase3[which(RMBase3$V1 %in% RMvar$chromosome),]
  
  if (f == '/path/to/mod/Nm.txt') {
    mod_grange_RMBase3 <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=RMBase3$V1)),
                                  ranges = IRanges(start = RMBase3$V2, end=RMBase3$V3),
                                  strand = Rle(RMBase3$V6),
                                  mod_type = 'Nm',
                                  cell_line = RMBase3$V12,
                                  DB='RMBase3')
  } else {
    mod_grange_RMBase3 <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=RMBase3$V1)),
                                  ranges = IRanges(start = RMBase3$V2, end=RMBase3$V3),
                                  strand = Rle(RMBase3$V6),
                                  mod_type = RMBase3$V7,
                                  cell_line = RMBase3$V12,
                                  DB='RMBase3')
  }
  
  # create a unique bed file with the coordinates of the RNA marks of that RNA modification types from both databases
  all_granges <- list(mod_grange_RMvar, mod_grange_RMBase3)
  names(all_granges) <- c('RMvar', 'RMBase3')
  order <- order(c(length(mod_grange_RMvar), length(mod_grange_RMBase3)))
  
  overlap <- findOverlaps(all_granges[[order[1]]],all_granges[[order[2]]], type = 'any')
  overlap_query <- all_granges[[order[1]]][unique(queryHits(overlap))]
  overlap_subject <- all_granges[[order[2]]][unique(subjectHits(overlap))]
  
  if (length(overlap) ==0) {
    print(unique(mod_grange_RMBase3$mod_type))
    total <- c(mod_grange_RMvar, mod_grange_RMBase3)
    print(length(total))
    total_table <- data.frame(chr = as.vector(seqnames(total)), start=start(total), end=end(total), strand=strand(total), mod_type= unique(mod_grange_RMBase3$mod_type), cell_line =total$cell_line, DB = total$DB)
    write.table(total_table, file = paste0('/path/to/mod/','/', unique(mod_grange_RMBase3$mod_type), '_total.bed'))
  } else {
    print(unique(mod_grange_RMBase3$mod_type))
    total <- c(all_granges[[order[1]]], all_granges[[order[2]]][-unique(subjectHits(overlap))])
    print(length(total))
    total_table <- data.frame(chr = as.vector(seqnames(total)), start=start(total), end=end(total), strand=strand(total), mod_type= unique(mod_grange_RMBase3$mod_type), cell_line =total$cell_line, DB = total$DB)
    write.table(total_table, file = paste0('/path/to/mod/','/', unique(mod_grange_RMBase3$mod_type), '_total.bed'))
  }
}
###############

gtf_file <- "/path/to/Homo_sapiens.GRCh38.104.gtf"
txdb <- makeTxDbFromGFF(gtf_file)

# function that produces a barplot, for each RNA modification type, reporting the number of RNA marks 
# identified in each cell line
num_cell_lines <- function(cell_lines_RMBase, cell_lines_RMvar, mod_type) {
  
  # substitute NA values with "Unknown"
  cell_lines_RMBase[is.na(cell_lines_RMBase)] <- 'Unknown'
  cell_lines_RMBase[cell_lines_RMBase == '^-'] <- 'Unknown'
  cell_lines_RMBase <- unlist(lapply(as.list(cell_lines_RMBase), function(x) {
    a <- unlist(strsplit(x, split=','))
    a<- gsub('^NA', 'Unknown', a)
    a<- gsub('^na', 'Unknown', a)
    a<- gsub('^-', 'Unknown', a)
    a
  }))
  
  cell_lines_RMvar[is.na(cell_lines_RMvar)] <- 'Unknown'
  cell_lines_RMvar[cell_lines_RMvar == '^-'] <- 'Unknown'
  cell_lines_RMvar <- unlist(lapply(as.list(cell_lines_RMvar), function(x) {
    a <- unlist(strsplit(x, split=';'))
    a<- gsub('^NA', 'Unknown', a)
    a<- gsub('^na', 'Unknown', a)
    a<- gsub('^-', 'Unknown', a)
    a
  }))
  
  # initiate a data frame reporting, for each cell line, the number of RNA marks of a specific RNA modification
  # type identified in that cell line
  counts <- data.frame(cell_line = unique(c(unique(cell_lines_RMBase), unique(cell_lines_RMvar))), num_mods = rep(0,length(unique(c(unique(cell_lines_RMBase), unique(cell_lines_RMvar))))))
  
  for (hit in cell_lines_RMBase) {
    a <- unlist(strsplit(hit, split=','))
    for (cell in a) {
      counts[counts$cell_line == cell,2] <- counts[counts$cell_line == cell,2] +1
    }
  }
  
  for (hit in cell_lines_RMvar) {
    a <- unlist(strsplit(hit, split=';'))
    for (cell in a) {
      counts[counts$cell_line == cell,2] <- counts[counts$cell_line == cell,2] +1
    }
  }
  
  freq <- counts$num_mods
  names(freq) <- counts$cell_line
  order <- sort(freq, decreasing = TRUE)
  counts <- counts[counts$cell_line %in% names(order)[1:20],]
  
  if (length(order) > 20) {
    x_labels <- names(order)[1:20]
  } else {
    x_labels <- names(order)
  }
  
  # generate a barplot reporting the number of RNA marks identified in each cell line only
  # considering the 20 cell lines with the highest number of RNA marks 
  p <- ggplot(data = counts, aes(x = cell_line, y = num_mods)) +
    geom_col(width = 0.1, position="dodge") +
    theme_classic() +
    xlab('Cell line')+
    ylab('Number of entries confirmed by a cell line') +
    theme(axis.text.x=element_text(size=18), axis.title = element_text(size=15), plot.title = element_text(size=20)) +
    labs(title=mod_type) +
    scale_x_discrete(limits=x_labels) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) 
  
  return(p)
}

mod <- list.files(path = '/path/to/mod', pattern = "total.bed", full.names = TRUE)
for (f in mod) {
  
  mod_file <- read.table(f)
  
  if (f == "/path/to/mod/other_RNA_mods_total.bed") {
    cell_line_RMBase <- mod_file[,12]
    p <- num_cell_lines(cell_line_RMBase, c(), 'Others')
  } else {
    cell_line_RMBase <- mod_file[mod_file$DB =='RMBase3',6]
    cell_line_RMvar <- mod_file[mod_file$DB =='RMvar',6]
    p <- num_cell_lines(cell_line_RMBase, cell_line_RMvar, unique(mod_file$mod_type))
  }
  
  if (f == "/path/to/mod/other_RNA_mods_total.bed") {
    mod <- 'Others'
  } else {
    mod <- unique(mod_file$mod_type)
  }
  
  ggsave(paste0('/path/to/mod/', mod, '_cell_lines_bothBD.pdf'), plot = p, height = 15, width = 40)
}

#################

# overlap between ELIGOS hits and RNA marks, of each RNA mod type, from the two databases
# directory_hits is the path to the directory containing the the results of ELIGOS analysis,
# create a new folder /mods_DB/ inside this directory to save the result of the overlap between
# ELIGOS hits of a fraction and the RNA marks databases.
# directory_mod is the path to the directory containing, for each RNA mod type, the bed file with
# the coordinates of all the RNA marks of that type from RMBase3+RMVar.
# m6A is TRUE if when analysing DRACH+ hits, otherwise FALSE
overlap_marks <- function(hits_chr_grange,hits_nucleo_grange,hits_cyto_grange,m6A,directory_hits,directory_mod) {
  
  # initiate a matrix that will report, for each fraction, the number (and percentage) of ELIGOS hits containing each RNA mod type
  # and the overall number (and percentage) of ELIGOS hits annotated at least to one RNA modification type. 
  # The matrix will also report the number of ELIGOS hits containing the same RNA mod type and shared by all the fractions
  number_hits_per_mod <- matrix(0,nrow=5, ncol=12)
  rownames(number_hits_per_mod) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm', 'Overlap','Number of annotated marks')
  if (m6A) {
    colnames(number_hits_per_mod) <- c('m6A','Y','m1A','m5C', 'm7G','A-I', 'Nm', 'm6Am', 'm5U', 'Others', 'DRACH+ hits annotated', 'Tot DRACH+ hits')
  } else {
    colnames(number_hits_per_mod) <- c('m6A','Y','m1A','m5C', 'm7G','A-I', 'Nm', 'm6Am', 'm5U', 'Others', 'DRACH- hits annotated', 'Tot DRACH- hits')
  }
  
  number_hits_per_mod[1,12] <- as.character(length(hits_chr_grange))
  number_hits_per_mod[2,12] <- as.character(length(hits_nucleo_grange))
  number_hits_per_mod[3,12] <- as.character(length(hits_cyto_grange))
  
  # initiate the vectors that, for each fraction, will report the index of ELIGOS hits overlapping with each RNA mod type
  query_matching_mark_chr <- c()
  query_matching_mark_nucleo <- c()
  query_matching_mark_cyto <- c()
  
  mod <- list.files(path = directory_mod, pattern = "total.bed", full.names = TRUE)
  
  for (f in mod) {
    
    mod_file <- read.table(f)
    
    if (f == "/path/to/mod/other_RNA_mods_total.bed") {
      mod_grange <- GRanges(seqnames = gsub(pattern = 'M', replacement = 'MT', x = gsub(pattern='chr', replacement = '', x=mod_file[,1])),
                            ranges = IRanges(start = mod_file[,2], end=mod_file[,3]),
                            strand = Rle(mod_file[,6]),
                            mod_type = mod_file[,7])
      number_hits_per_mod[5,'Others'] <- as.character(length(mod_grange))
    } else {
      mod_grange <- GRanges(seqnames = mod_file$chr,
                            ranges = IRanges(start = mod_file$start, end=mod_file$end),
                            strand = Rle(mod_file$strand),
                            mod_type =  mod_file$mod_type)
      number_hits_per_mod[5,unique(mod_grange$mod_type)] <- as.character(length(mod_grange))
    }
    
    # overlap between ELIGOS hits of a fraction and the coordinates of a specific RNA mark from the two databases
    suppressWarnings(overlap_chr <- findOverlaps(hits_chr_grange,mod_grange))
    hits_chr_with_mod <- c()
    if (length(overlap_chr) != 0) {
      # save inside the vector the index of ELIGOS hits overlapping with this RNA mod type
      query_matching_mark_chr <- c(query_matching_mark_chr, unique(queryHits(overlap_chr)))
      hits_chr_with_mod <- hits_chr_grange[unique(queryHits(overlap_chr))]
      # report in the matrix the number (and percentage) of hits that can be annotated with this RNA mod type
      if (f == "/path/to/mod/other_RNA_mods_total.bed") {
        number_hits_per_mod[1,'Others'] <- paste0(as.character(length(unique(queryHits(overlap_chr)))), ' - ', as.character(round(length(unique(queryHits(overlap_chr)))*100/as.numeric(number_hits_per_mod[1,12]),2)),'%')
      } else {
        number_hits_per_mod[1,unique(mod_grange$mod_type)] <- paste0(as.character(length(unique(queryHits(overlap_chr)))), ' - ', as.character(round(length(unique(queryHits(overlap_chr)))*100/as.numeric(number_hits_per_mod[1,12]),2)),'%')
      }
      # update the metadata field "mod_type" of ELIGOS hits overlapping with this modification type
      for (i in unique(queryHits(overlap_chr))) {
        if (hits_chr_grange[i]$mod_type == 'Unknown') {
          hits_chr_grange[i]$mod_type <- paste(unique(mod_grange[subjectHits(overlap_chr[queryHits(overlap_chr) == i])]$mod_type), collapse = ';')
        } else {
          hits_chr_grange[i]$mod_type <- paste(c(hits_chr_grange[i]$mod_type,unique(mod_grange[subjectHits(overlap_chr[queryHits(overlap_chr) == i])]$mod_type)), collapse = ';')
        }
      }
    } 
    
    suppressWarnings(overlap_nucleo <- findOverlaps(hits_nucleo_grange,mod_grange))   
    hits_nucleo_with_mod <- c()
    if (length(overlap_nucleo) != 0) {
      # save inside the vector the index of ELIGOS hits overlapping with this RNA mod type
      query_matching_mark_nucleo <- c(query_matching_mark_nucleo, unique(queryHits(overlap_nucleo)))
      hits_nucleo_with_mod <- hits_nucleo_grange[unique(queryHits(overlap_nucleo))]
      # report in the matrix the number (and percentage) of hits that can be annotated with this RNA mod type
      if (f == "/path/to/mod/other_RNA_mods_total.bed") {
        number_hits_per_mod[2,'Others'] <- paste0(as.character(length(unique(queryHits(overlap_nucleo)))), ' - ', as.character(round(length(unique(queryHits(overlap_nucleo)))*100/as.numeric(number_hits_per_mod[2,12]),2)),'%')
      } else {
        number_hits_per_mod[2,unique(mod_grange$mod_type)] <- paste0(as.character(length(unique(queryHits(overlap_nucleo)))), ' - ', as.character(round(length(unique(queryHits(overlap_nucleo)))*100/as.numeric(number_hits_per_mod[2,12]),2)),'%')
      }
      # update the metadata field "mod_type" of ELIGOS hits overlapping with this modification type
      for (i in unique(queryHits(overlap_nucleo))) {
        if (hits_nucleo_grange[i]$mod_type == 'Unknown') {
          hits_nucleo_grange[i]$mod_type <- paste(unique(mod_grange[subjectHits(overlap_nucleo[queryHits(overlap_nucleo) == i])]$mod_type), collapse = ';')
        } else {
          hits_nucleo_grange[i]$mod_type <- paste(c(hits_nucleo_grange[i]$mod_type,unique(mod_grange[subjectHits(overlap_nucleo[queryHits(overlap_nucleo) == i])]$mod_type)), collapse = ';')
        }
      }
    }
    
    suppressWarnings(overlap_cyto <- findOverlaps(hits_cyto_grange,mod_grange))
    hits_cyto_with_mod <- c()
    if (length(overlap_cyto) != 0) {
      # save inside the vector the index of ELIGOS hits overlapping with this RNA mod type
      query_matching_mark_cyto <- c(query_matching_mark_cyto, unique(queryHits(overlap_cyto)))
      hits_cyto_with_mod <- hits_cyto_grange[unique(queryHits(overlap_cyto))]
      # report in the matrix the number (and percentage) of hits that can be annotated with this RNA mod type
      if (f == "/path/to/mod/other_RNA_mods_total.bed") {
        number_hits_per_mod[3,'Others'] <- paste0(as.character(length(unique(queryHits(overlap_cyto)))), ' - ', as.character(round(length(unique(queryHits(overlap_cyto)))*100/as.numeric(number_hits_per_mod[3,12]),2)),'%')
      } else {
        number_hits_per_mod[3,unique(mod_grange$mod_type)] <- paste0(as.character(length(unique(queryHits(overlap_cyto)))), ' - ', as.character(round(length(unique(queryHits(overlap_cyto)))*100/as.numeric(number_hits_per_mod[3,12]),2)),'%')
      }
      # update the metadata field "mod_type" of ELIGOS hits overlapping with this modification type
      for (i in unique(queryHits(overlap_cyto))) {
        if (hits_cyto_grange[i]$mod_type == 'Unknown') {
          hits_cyto_grange[i]$mod_type <- paste(unique(mod_grange[subjectHits(overlap_cyto[queryHits(overlap_cyto) == i])]$mod_type), collapse = ';')
        } else {
          hits_cyto_grange[i]$mod_type <- paste(c(hits_cyto_grange[i]$mod_type,unique(mod_grange[subjectHits(overlap_cyto[queryHits(overlap_cyto) == i])]$mod_type)), collapse = ';')
        }
      }
    }
    
    # compute how many ELIGOS hits, overlapping with this RNA modification type, are shared by all the fractions and add
    # this information to the matrix
    if (length(hits_chr_with_mod) !=0 & length(hits_nucleo_with_mod) != 0 & length(hits_cyto_with_mod) != 0) {
      hits_with_mod <- list(hits_chr_with_mod, hits_nucleo_with_mod, hits_cyto_with_mod)
      order <- order(c(length(hits_chr_with_mod), length(hits_nucleo_with_mod), length(hits_cyto_with_mod)))
      
      overlap <- findOverlaps(hits_with_mod[[order[1]]],hits_with_mod[[order[2]]], type = 'any')
      overlap2 <- hits_with_mod[[order[1]]][unique(queryHits(overlap))]
      overlap_all <- findOverlaps(overlap2,hits_with_mod[[order[3]]], type = 'any')
      if (f == "/path/to/mod/other_RNA_mods_total.bed") {
        number_hits_per_mod[4,'Others'] <- length(overlap2[unique(queryHits(overlap_all))]) 
      } else {
        number_hits_per_mod[4,unique(mod_grange$mod_type)] <- length(overlap2[unique(queryHits(overlap_all))]) 
      }
    } 
  }
  
  # add to the matrix the number (and percentage) of ELIGOS hits, of each fraction, annotated at least with one RNA mod type
  number_hits_per_mod[1,11] <- paste0(as.character(length(unique(query_matching_mark_chr))), ' - ', as.character(round(length(unique(query_matching_mark_chr))*100/as.numeric(number_hits_per_mod[1,12]),2)), '%')
  number_hits_per_mod[2,11] <- paste0(as.character(length(unique(query_matching_mark_nucleo))), ' - ', as.character(round(length(unique(query_matching_mark_nucleo))*100/as.numeric(number_hits_per_mod[2,12]),2)), '%')
  number_hits_per_mod[3,11] <- paste0(as.character(length(unique(query_matching_mark_cyto))), ' - ', as.character(round(length(unique(query_matching_mark_cyto))*100/as.numeric(number_hits_per_mod[3,12]),2)), '%')
  
  # add to the matrix how many ELIGOS hits - overlapping at least with one RNA modification type - are shared by all the fractions
  hits_chr_annotated <- hits_chr_grange[unique(query_matching_mark_chr)]
  hits_nucleo_annotated <- hits_nucleo_grange[unique(query_matching_mark_nucleo)]
  hits_cyto_annotated <- hits_cyto_grange[unique(query_matching_mark_cyto)]
  
  hits_annotated <- list(hits_chr_annotated, hits_nucleo_annotated, hits_cyto_annotated)
  order <- order(c(length(hits_chr_annotated), length(hits_nucleo_annotated), length(hits_cyto_annotated)))
  
  overlap <- findOverlaps(hits_annotated[[order[1]]],hits_annotated[[order[2]]], type = 'any')
  overlap2 <- hits_annotated[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,hits_annotated[[order[3]]], type = 'any')
  number_hits_per_mod[4,11] <- length(overlap2[unique(queryHits(overlap_all))])
  
  if (m6A) {
    write.xlsx(x = data.frame(number_hits_per_mod),file = paste0(directory_hits, '/mods_DB/', 'number_hit_per_mod_DRACH_bothBD.xlsx'),col.names = TRUE, row.names=TRUE)
  } else {
    write.xlsx(x = data.frame(number_hits_per_mod),file = paste0(directory_hits, '/mods_DB/', 'number_hit_per_mod_non_DRACH_bothBD.xlsx'),col.names = TRUE, row.names=TRUE)
  }
  
  unique_mod_types <- function(hits) {
    hits$mod_type <- unlist(lapply(seq_along(hits), function(i,x) {
      paste(unique(unlist(strsplit(x[i]$mod_type, split =';'))), collapse = ';')
    }, x=hits))
    return(hits)
  }
  
  hits_chr_grange <- unique_mod_types(hits_chr_grange)
  hits_nucleo_grange <- unique_mod_types(hits_nucleo_grange)
  hits_cyto_grange <- unique_mod_types(hits_cyto_grange)
  
  l <- list(hits_chr_grange,hits_nucleo_grange,hits_cyto_grange)
  names(l) <- c('hits_chr_grange','hits_nucleo_grange','hits_cyto_grange')
  return(l)
}

ELIGOS_hits_DB_marks <- function(directory_mod, directory_hits) {
  
  # load the RData containing ELIGOS hits of each fraction, DRACH+ and DRACH- separately
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH.Rda'))
  
  # add a metadata column to each GRanges object in which the RNA modification types with which each hit overlaps will be reported (initiated 
  # with "Unknown")
  mcols(hits_eligos_chr_ass_confirmed_5_with_DRACH) <- cbind(mcols(hits_eligos_chr_ass_confirmed_5_with_DRACH), mod_type= 'Unknown')
  mcols(hits_eligos_nucleo_confirmed_5_with_DRACH) <- cbind(mcols(hits_eligos_nucleo_confirmed_5_with_DRACH), mod_type= 'Unknown')
  mcols(hits_eligos_cyto_confirmed_5_with_DRACH) <- cbind(mcols(hits_eligos_cyto_confirmed_5_with_DRACH), mod_type= 'Unknown')
  mcols(hits_eligos_chr_ass_confirmed_5_without_DRACH) <- cbind(mcols(hits_eligos_chr_ass_confirmed_5_without_DRACH), mod_type= 'Unknown')
  mcols(hits_eligos_nucleo_confirmed_5_without_DRACH) <- cbind(mcols(hits_eligos_nucleo_confirmed_5_without_DRACH), mod_type= 'Unknown')
  mcols(hits_eligos_cyto_confirmed_5_without_DRACH) <- cbind(mcols(hits_eligos_cyto_confirmed_5_without_DRACH), mod_type= 'Unknown')
  
  # perform the overlap between ELIGOS DRACH+ hits and the RNA marks, of each RNA mod type, from the two databases
  hits_DRACH <- overlap_marks(hits_eligos_chr_ass_confirmed_5_with_DRACH,hits_eligos_nucleo_confirmed_5_with_DRACH,hits_eligos_cyto_confirmed_5_with_DRACH,TRUE,directory_hits,directory_mod)
  hits_eligos_chr_ass_confirmed_5_with_DRACH <- hits_DRACH[[1]]
  hits_eligos_nucleo_confirmed_5_with_DRACH <- hits_DRACH[[2]]
  hits_eligos_cyto_confirmed_5_with_DRACH <- hits_DRACH[[3]]
  
  # perform the overlap between ELIGOS DRACH- hits and the RNA marks, of each RNA mod type, from the two databases
  hits_non_DRACH <- overlap_marks(hits_eligos_chr_ass_confirmed_5_without_DRACH,hits_eligos_nucleo_confirmed_5_without_DRACH,hits_eligos_cyto_confirmed_5_without_DRACH,FALSE,directory_hits,directory_mod)
  hits_eligos_chr_ass_confirmed_5_without_DRACH <- hits_non_DRACH[[1]]
  hits_eligos_nucleo_confirmed_5_without_DRACH <- hits_non_DRACH[[2]]
  hits_eligos_cyto_confirmed_5_without_DRACH <- hits_non_DRACH[[3]]
  
  save(hits_eligos_chr_ass_confirmed_5_with_DRACH, file=paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type.Rda'))
  save(hits_eligos_nucleo_confirmed_5_with_DRACH, file=paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type.Rda'))
  save(hits_eligos_cyto_confirmed_5_with_DRACH, file=paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type.Rda'))
  save(hits_eligos_chr_ass_confirmed_5_without_DRACH, file =paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type.Rda'))
  save(hits_eligos_nucleo_confirmed_5_without_DRACH, file=paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type.Rda'))
  save(hits_eligos_cyto_confirmed_5_without_DRACH, file=paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type.Rda'))
}

ELIGOS_hits_DB_marks(directory_mod = '/path/to/mod', 
                     directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/')

ELIGOS_hits_DB_marks(directory_mod = '/path/to/mod', 
                     directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/')

ELIGOS_hits_DB_marks(directory_mod = '/path/to/mod', 
                     directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/')

ELIGOS_hits_DB_marks(directory_mod = '/path/to/mod', 
                     directory_hits = '/path/to/fractions_eligos_STORM_K562/eligos_total_reads/')


## DATABASES OF EFFECTORS' BINDING SITES

# Download from RMBase3 https://rna.sysu.edu.cn/rmbase3/download.php (RNA Modifications related to RBPs session) the bed files with the 
# coordinates of the binding sites of m6A writers/readers/erasers, m5C writers/readers, Y writers, Nm writers, A-I writers and other effectors
# associated with other RNA modifications selecting Mammal group, Homo sapiens hg38 genome.
# These bed files were saved in a folder called /RBPs/ and renamed "modification_type_RBPs.bed" and "other_RNA_mods_RBPs.bed".

# concatenate the different bed files reporting the coordinates of the different categories of effectors (writer..)
# associated with the same RNA modification type
m6A_writers <- read.table('/path/to/RBPs/human.hg38.modrbp.m6A.writer.bed', fill =TRUE)
m6A_erasers <- read.table('/path/to/RBPs/human.hg38.modrbp.m6A.eraser.bed', fill =TRUE)
m6A_readers <- read.table('/path/to/RBPs/human.hg38.modrbp.m6A.reader.bed', fill =TRUE)
write.table(rbind(m6A_writers, m6A_erasers, m6A_readers), file = '/path/to/RBPs/m6A_RBPs.bed', row.names = FALSE, col.names = FALSE)

m5C_writers <- read.table('/path/to/RBPs/human.hg38.modrbp.m5C.writer.bed', fill =TRUE)
m5C_readers <- read.table('/path/to/RBPs/human.hg38.modrbp.m5C.reader.bed', fill =TRUE)
write.table(rbind(m5C_writers, m5C_readers), file = '/path/to/RBPs/m5C_RBPs.bed', row.names = FALSE, col.names = FALSE)

# From other_RNA_mods_RBPs.bed the coordinates of the binding sites associated with known RNA mod types are extracted and appended to 
# the already present files. Only the 4 columns of interest are kept
effectors_others_RMBase <- read.table('/path/to/RBPs/other_RNA_mods_RBPs.bed', fill=TRUE)[,c(2,3,10,17)]
m6A1 <- effectors_others_RMBase[which(effectors_others_RMBase$V10 == 'm6A'),]
m6A2 <- read.table('/path/to/RBPs/m6A_RBPs.bed', fill = TRUE)[,c(2,3,10,17)]
write.table(x = rbind(m6A1, m6A2), file = '/path/to/RBPs/m6A_RBPs.bed', row.names = FALSE, col.names = FALSE)

m5C1 <- effectors_others_RMBase[which(effectors_others_RMBase$V10 == 'm5C'),]
m5C2 <- read.table('/path/to/RBPs/m5C_RBPs.bed', fill = TRUE)[,c(2,3,10,17)]
write.table(x = rbind(m5C1, m5C2), file = '/path/to/RBPs/m5C_RBPs.bed', row.names = FALSE, col.names = FALSE)

AI1 <- effectors_others_RMBase[which(effectors_others_RMBase$V10 == 'A-I'),]
AI2 <- read.table('/path/to/RBPs/A-I_RBPs.bed', fill = TRUE)[,c(2,3,10,17)]
write.table(x = rbind(AI1, AI2), file = '/path/to/RBPs/A-I_RBPs.bed', row.names = FALSE, col.names = FALSE)

Y1 <- effectors_others_RMBase[which(effectors_others_RMBase$V10 == 'Y'),]
Y2 <- read.table('/path/to/RBPs/Y_RBPs.bed', fill = TRUE)[,c(2,3,10,17)]
write.table(x = rbind(Y1,Y2), file = '/path/to/RBPs/Y_RBPs.bed', row.names = FALSE, col.names = FALSE)

m1A <- effectors_others_RMBase[which(effectors_others_RMBase$V10 == 'm1A'),]
write.table(x = m1A, file = '/path/to/RBPs/m1A_RBPs.bed',  row.names = FALSE, col.names = FALSE)

m6Am <- effectors_others_RMBase[which(effectors_others_RMBase$V10 == 'm6Am'),]
write.table(x = m6Am, file = '/path/to/RBPs/m6Am_RBPs.bed', row.names = FALSE, col.names = FALSE)

m7G <- effectors_others_RMBase[which(effectors_others_RMBase$V10 == 'm7G'),]
write.table(x = m7G, file = '/path/to/RBPs/m7G_RBPs.bed', row.names = FALSE, col.names = FALSE)

m5U <- effectors_others_RMBase[which(effectors_others_RMBase$V10 == 'm5U'),]
write.table(x = m5U, file = '/path/to/RBPs/m5U_RBPs.bed', row.names = FALSE, col.names = FALSE)

Nm1 <- effectors_others_RMBase[which(effectors_others_RMBase$V10 %in% c('Um', 'Am', 'Cm', 'Gm')),]
Nm2 <- read.table('/path/to/RBPs/Nm_RBPs.bed', fill = TRUE)[,c(2,3,10,17)]
write.table(x = rbind(Nm1,Nm2), file = '/path/to/RBPs/Nm_RBPs.bed', row.names = FALSE, col.names = FALSE)

effectors_others_RMBase <- effectors_others_RMBase[-which(effectors_others_RMBase$V10 %in% c('m6A', 'm5C', 'A-I', 'Y', 'm1A', 'm6Am', 'm7G', 'm5U', 'Um', 'Am', 'Cm', 'Gm')),]
write.table(effectors_others_RMBase, '/path/to/RBPs/other_RNA_mods_RBPs.bed', row.names = FALSE, col.names = FALSE)

# from each bed file with the coordinates of the binding sites associated with each RNA mod type extract the 
# unique rows (effector, coordinates, mod, type_effector)
RBPs <- list.files(path = '/path/to/RBPs', pattern = 'bed', full.names = TRUE)
for (i in 1:length(RBPs)) {
  f <- read.table(RBPs[i])
  if (RBPs[i] == '/path/to/RBPs/Nm_RBPs.bed') {
    f[,3] <- 'Nm'
  }
  f_unique <- f %>% distinct()
  write.table(f_unique, file = RBPs[i],row.names = FALSE, col.names = FALSE)
}

# compute the overall number of binding sites (the same effector may have the same binding site reported multiple times 
# because it is associated with different RNA mods), the number of binding sites for each RNA mod and the overall number of 
# effectors associated with RNA modification types reported in RMBase3
enzymes_RMBase <- c()
tot_sites <- c()
for (f in RBPs) {
  enzymes_RMBase <- c(enzymes_RMBase, read.table(f)[,1])
  tot_sites <- rbind(tot_sites, unique(read.table(f)[,c(1,2)]))
  print(f)
  print(nrow(unique(read.table(f)[,c(1,2,3)])))
}

nrow(tot_sites %>% distinct())
unique_RMBase3_binding_sites <- unique((tot_sites %>% distinct())[,2])
enzymes_RMBase <- unique(enzymes_RMBase)

# Download from RMVar https://rmvar.renlab.org/download.html the txt file with the coordinates of the binding sites of
# effectors associated with RNA modifications. This file is called RMVar_Human_RBP_info.txt and was saved in the same folder /RBPs/
RBPs_RMvar <- read_tsv('/path/to/RBPs/RMVar_Human_RBP_info.txt',  col_types = cols())
# the columns of interest and the unique binding sites are extracted 
RBPs_RMvar <- RBPs_RMvar[,c(2,6)]
RBPs_RMvar <- RBPs_RMvar %>% distinct()
# compute the overall number of binding sites and of effectors associated with RNA modification types reported in RMVar
nrow(RBPs_RMvar)
enzymes_RMvar <- unique(RBPs_RMvar$RBP)

# create a data frame reporting the name of the effector, the function of the effector (writer, eraser, reader, other) and 
# the modification type associated to it (on the base of the information reported in RMBase3)
enzyme_type_mod <- data.frame()
for (f in RBPs) {
  enzyme_type_mod <- rbind(enzyme_type_mod, read.table(f)[,c(1,4,3)])
}

length(unique(enzyme_type_mod[,3]))
# 33 different modifications reported in RMBase3
length(unique(enzyme_type_mod[,1]))
# 203 different enzymes reported in RMBase3
length(enzymes_RMvar[enzymes_RMvar %in% unique(enzyme_type_mod[,1])])
# 197 enzymes in common between RMBase3 and RMvar

enzyme_type_mod_unique <- enzyme_type_mod %>% distinct()
colnames(enzyme_type_mod_unique) <- c('Name_effector', 'Function', 'Mod_type')
nrow(enzyme_type_mod_unique[,c(1,2)] %>% distinct())
# 203 -> the same enzyme is always annotated to the same function 

# identify which effectors are associated only with one RNA mod type (specific) and
# which with more than one (ambiguous).
# in each RMBase file there are no duplicates of the same binding site (the binding site of each effector is not duplicated); 
# but the same binding site of the same effector can be associated with different modifications (contained in different RMBase3 files)
specific_effectors <- names(table(enzyme_type_mod_unique$Name_effector)[table(enzyme_type_mod_unique$Name_effector)==1])
# 31 specific effectors 
save(specific_effectors, file = '/path/to/RBPs/specific_effectors.RDa')

specific_effectors_only_m6A <- enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% specific_effectors,][enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% specific_effectors,]$Mod_type =='m6A',]$Name_effector
# 22 specific effectors associated with m6A
save(specific_effectors_only_m6A, file = '/path/to/RBPs/specific_effectors_only_m6A.RDa')

specific_effectors_non_m6A <- specific_effectors[!specific_effectors %in% specific_effectors_only_m6A]
# 9 specific effectors associated with a RNA mod different from m6A
save(specific_effectors_non_m6A, file = '/path/to/RBPs/specific_effectors_non_m6A.RDa')

# identify which effectors are reported only in RMVar
enzymes_only_RMvar <-setdiff(unique(RBPs_RMvar$RBP),unique(enzyme_type_mod_unique$Name_effector))

# add to the data frame the names of the effectors only present in RMVar, in this case the relative function and modification type are set to "Unknown"
enzyme_type_mod_unique <- rbind(enzyme_type_mod_unique, data.frame(Name_effector=enzymes_only_RMvar,Function=rep('Unknown', length(enzymes_only_RMvar)),Mod_type=rep('Unknown', length(enzymes_only_RMvar))))

length(unique(enzyme_type_mod_unique$Name_effector))
# 207 unique effectors annotated in at least one of the two databases

ambiguous_effectors <- unique(setdiff(unique(enzyme_type_mod_unique$Name_effector), specific_effectors))
# 176 ambiguous effectors 
save(ambiguous_effectors, file = '/path/to/RBPs/ambiguous_effectors.RDa')

save(enzyme_type_mod_unique, file = '/path/to/RBPs/enzyme_type_mod_unique.RDa')

# distribution of the dimension of the binding sites (narrow peak) from all the bed files from RMBase3
dim_RMBase <- unlist(lapply(as.list(unique_RMBase3_binding_sites), function(x) {
  coord <- strsplit(x, split=':')[[1]][2]
  start <- strsplit(coord, split = '-')[[1]][1]
  end <- strsplit(coord, split = '-')[[1]][2]
  as.numeric(end)-as.numeric(start)
}))

summary(dim_RMBase)

# saturation at 200 nucleotides
dim_RMBase[dim_RMBase>200] <- 200

# plot the distribution of the dimension of the binding sites (narrow peak) from RMBase3
p <- ggplot(as.data.frame(dim_RMBase), aes(x=dim_RMBase)) +
  geom_histogram(bins = 100) +
  theme_classic() +
  labs(title='Distribution of binding site dimension - RMBase3') +
  xlab('Binding site dimension')+
  ylab('Number of effectors')

ggsave('/path/to/RBPs/binding_site_dimension_RMBase.jpeg', plot = p, height = 7, width = 7)

# distribution of the dimension of the binding sites (narrow peak) from RMVar
dim_RMVar <- unlist(lapply(as.list(unique(RBPs_RMvar$binding_region)), function(x) {
  coord <- strsplit(x, split=':')[[1]][2]
  start <- strsplit(coord, split = '\\..')[[1]][1]
  end <- strsplit(coord, split = '\\..')[[1]][2]
  as.numeric(end)-as.numeric(start)
}))

summary(dim_RMVar)

# saturation at 200 nucleotides
dim_RMVar[dim_RMVar>200] <- 200

# plot the distribution of the dimension of the binding sites (narrow peak) from RMVar
p <- ggplot(as.data.frame(dim_RMVar), aes(x=dim_RMVar)) +
  geom_histogram(bins = 100) +
  theme_classic() +
  labs(title='Distribution of binding site dimension - RMvar') +
  xlab('Binding site dimension')+
  ylab('Number of effectors')

ggsave('/path/to/RBPs/binding_site_dimension_RMvar.jpeg', plot = p, height = 7, width = 7)

# resize the binding sites reported by RMVar and RMBase3 to a width at most of 50 nucleotides to compute how many binding 
# sites are confirmed by both
grange_all_sites_RMvar <-  GRanges(seqnames = gsub('M', 'MT', gsub('chr', '',unlist(lapply(as.list(RBPs_RMvar$binding_region), function(x) strsplit(x, split=':')[[1]][1])))),
                                   ranges = IRanges(start = as.numeric(unlist(lapply(as.list(RBPs_RMvar$binding_region), 
                                                                                     function(x) strsplit(strsplit(x, split=':')[[1]][2], split='\\..')[[1]][1]))), 
                                                    end=as.numeric(unlist(lapply(as.list(RBPs_RMvar$binding_region), function(x) strsplit(strsplit(x, split=':')[[1]][2], split='\\..')[[1]][2])))),
                                   RBP_name = RBPs_RMvar$RBP)

grange_all_sites_RMvar_min50 <- grange_all_sites_RMvar[width(grange_all_sites_RMvar) <=50]
grange_all_sites_RMvar_mag50 <- grange_all_sites_RMvar[width(grange_all_sites_RMvar) >50]
grange_all_sites_RMvar_mag50 <- resize(grange_all_sites_RMvar_mag50, 50, fix = 'center')
grange_all_sites_RMvar <- c(grange_all_sites_RMvar_min50,grange_all_sites_RMvar_mag50)

grange_all_sites_RMBase3 <- grange <- GRanges(seqnames = gsub('M', 'MT', gsub('chr', '', unlist(lapply(as.list((tot_sites %>% distinct())[,2]), function(x) strsplit(x, split=':')[[1]][1])))),
                                              ranges = IRanges(start = as.numeric(unlist(lapply(as.list((tot_sites %>% distinct())[,2]), function(x) strsplit(strsplit(x, split=':')[[1]][2], split='-')[[1]][1]))), 
                                                               end=as.numeric(unlist(lapply(as.list((tot_sites %>% distinct())[,2]), function(x) strsplit(strsplit(x, split=':')[[1]][2], split='-')[[1]][2])))),
                                              RBP_name = (tot_sites %>% distinct())[,1])

grange_min50 <- grange_all_sites_RMBase3[width(grange_all_sites_RMBase3) <=50]
grange_mag50 <- grange_all_sites_RMBase3[width(grange_all_sites_RMBase3) >50]
grange_mag50 <- resize(grange_mag50, 50, fix = 'center')
grange_all_sites_RMBase3 <- c(grange_min50,grange_mag50)

# compute how many binding sites are confirmed by both databases
tot_confirmed <- 0
for (enzyme in unique(grange_all_sites_RMBase3$RBP_name)) {
  enzyme_RMBase <- grange_all_sites_RMBase3[grange_all_sites_RMBase3$RBP_name == enzyme]
  enzyme_RMvar <- grange_all_sites_RMvar[grange_all_sites_RMvar$RBP_name == enzyme]
  tot_confirmed <- tot_confirmed + length(which(countOverlaps(enzyme_RMBase,enzyme_RMvar) != 0))
}
print(tot_confirmed)

# resize each binding site from RMBase3 and RMVar to its midpoint to merge the binding sites from 
# the two databases
grange_all_sites_RMBase3 <- GRanges()
for (f in RBPs) {
  t <- read.table(f)
  grange <- GRanges(seqnames = gsub('M', 'MT', gsub('chr', '', unlist(lapply(as.list(t[,2]), function(x) strsplit(x, split=':')[[1]][1])))),
                    ranges = IRanges(start = as.numeric(unlist(lapply(as.list(t[,2]), function(x) strsplit(strsplit(x, split=':')[[1]][2], split='-')[[1]][1]))), end=as.numeric(unlist(lapply(as.list(t[,2]), function(x) strsplit(strsplit(x, split=':')[[1]][2], split='-')[[1]][2])))),
                    RBP_name = t[,1],
                    mod_type = t[,3],
                    type_enzyme = t[,4])
  
  grange <- resize(grange, 1, fix = 'center')
  grange_all_sites_RMBase3 <- c(grange_all_sites_RMBase3, grange)
}

# report each effector annotated in RMBase3 only once 
# indicating in the "mod_type" field all the modification types associated with it
grange_all_sites_RMBase_unique <- GRanges()
for (enzyme in unique(grange_all_sites_RMBase3$RBP_name)) {
  grange_enzyme_RMBase <- grange_all_sites_RMBase3[grange_all_sites_RMBase3$RBP_name==enzyme]
  ind_unique <- which(isUnique(paste0(seqnames(grange_enzyme_RMBase), "_", start(grange_enzyme_RMBase), "_", end(grange_enzyme_RMBase))))
  grange_enzyme_RMBase_unique <- grange_enzyme_RMBase[ind_unique]
  if (length(ind_unique) != length(grange_enzyme_RMBase)) {
    dup <- grange_enzyme_RMBase[-ind_unique]
    for (r in 1:length(ranges(unique(dup)))) {
      duplicates <- dup[ranges(dup) == ranges(unique(dup))[r]] 
      s <- unique(duplicates) 
      s$mod_type <- paste(unique(duplicates$mod_type), collapse = ';')
      grange_enzyme_RMBase_unique <- c(grange_enzyme_RMBase_unique, s)
    }}
  grange_all_sites_RMBase_unique <- c(grange_all_sites_RMBase_unique, grange_enzyme_RMBase_unique)}

# check if all the effectors reported by RMBase3 have each binding site reported only once
for (enzyme in unique(grange_all_sites_RMBase_unique$RBP_name)) {
  grange_enzyme_RMBase <- grange_all_sites_RMBase_unique[grange_all_sites_RMBase_unique$RBP_name==enzyme]
  ind_unique <- which(isUnique(paste0(seqnames(grange_enzyme_RMBase), "_", start(grange_enzyme_RMBase), "_", end(grange_enzyme_RMBase))))
  if (length(ind_unique) != length(grange_enzyme_RMBase)) {
    print(enzyme)
    print(grange_enzyme_RMBase[-ind_unique])
  }}

which(grange_all_sites_RMBase_unique$RBP_name=='PTBP1' & start(grange_all_sites_RMBase_unique) == 17309625 & seqnames(grange_all_sites_RMBase_unique) ==11)
grange_all_sites_RMBase_unique <- grange_all_sites_RMBase_unique[-1716511]
which(grange_all_sites_RMBase_unique$RBP_name=='PTBP1' & start(grange_all_sites_RMBase_unique) == 17309625 & seqnames(grange_all_sites_RMBase_unique) ==19)
grange_all_sites_RMBase_unique <- grange_all_sites_RMBase_unique[-1716511]
save(grange_all_sites_RMBase_unique, file ='/path/to/RBPs/grange_all_sites_RMBase_unique.Rda')

# resize each binding site from RMVar to its midpoint
grange_all_sites_RMvar <- resize(grange_all_sites_RMvar, 1, fix='center')

# add the "mod_type" field to grange_all_sites_RMvar in which, for each effector present both in RMVar and RMBase3, the 
# modification types associated with that effector are reported (using the information from RMBase3)
mcols(grange_all_sites_RMvar) <- cbind(mcols(grange_all_sites_RMvar), mod_type='Unknown', type_enzyme='Unknown')
all_enzymes <- unique(union(grange_all_sites_RMBase_unique$RBP_name, grange_all_sites_RMvar$RBP_name))

# create a unique GRanges reporting all the binding sites from both databases: the overlap is done on the binding sites of 1 nt 
all_sites <- GRanges()
for (enzyme in all_enzymes) {
  if (enzyme %in% unique(grange_all_sites_RMBase_unique$RBP_name) & enzyme %in% unique(grange_all_sites_RMvar$RBP_name)) {
    grange_enzyme_RMBase <- grange_all_sites_RMBase_unique[grange_all_sites_RMBase_unique$RBP_name==enzyme]
    grange_enzyme_RMvar <- grange_all_sites_RMvar[grange_all_sites_RMvar$RBP_name==enzyme]
    # the same enzyme from RMVar may have different binding sites of different dimension that, when reduced to 1 nt, are identical
    grange_enzyme_RMvar <- unique(grange_enzyme_RMvar)
    if (enzyme %in% specific_effectors) {
      grange_enzyme_RMvar$mod_type <- enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector==enzyme,3]
      grange_enzyme_RMvar$type_enzyme <- enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector==enzyme,2]
    }
    
    over <- findOverlaps(grange_enzyme_RMBase,grange_enzyme_RMvar)
    if (length(over) != 0) {
      all_sites_new <- c(grange_enzyme_RMBase,grange_enzyme_RMvar[-unique(subjectHits(over))]) 
    } else {
      all_sites_new <- c(grange_enzyme_RMBase,grange_enzyme_RMvar)
    }
    all_sites <- c(all_sites,all_sites_new)
  } else if ((enzyme %in% unique(grange_all_sites_RMBase_unique$RBP_name)) & !(enzyme %in% unique(grange_all_sites_RMvar$RBP_name))) {
    grange_enzyme_RMBase <- grange_all_sites_RMBase_unique[grange_all_sites_RMBase_unique$RBP_name==enzyme]
    all_sites <- c(all_sites,grange_enzyme_RMBase)
  } else if (!(enzyme %in% unique(grange_all_sites_RMBase_unique$RBP_name)) & (enzyme %in% unique(grange_all_sites_RMvar$RBP_name))) {
    grange_enzyme_RMvar <- grange_all_sites_RMvar[grange_all_sites_RMvar$RBP_name==enzyme]
    grange_enzyme_RMvar <- unique(grange_enzyme_RMvar)
    all_sites <- c(all_sites,grange_enzyme_RMvar)
  }
}

length(unique(all_sites$RBP_name))

# check if there aren't effector-binding site associations that are duplicated after the joining of RMBase3 and RMvar databases
for (en in unique(all_sites$RBP_name)) {
  gr <- all_sites[all_sites$RBP_name==en]
  if (length(gr[-which(isUnique(paste0(seqnames(gr), "_", start(gr), "_", end(gr))))]) !=0) {
    print(en)
  }
}

# resize the resulting binding sites to 50 nucleotides
all_sites_RMBase_RMvar <- resize(all_sites, 50, fix = 'center')
save(all_sites_RMBase_RMvar, file ='/path/to/RBPs/all_sites_RMBase_RMvar.Rda')

# compute how many binding sites, from both databases, are associated with specific effectors associated with m6A
length(all_sites_RMBase_RMvar[all_sites_RMBase_RMvar$RBP_name %in% specific_effectors_only_m6A])

# overlap (ignoring the strand) between ELIGOS hits and the binding sites of the effectors associated with RNA mod types 
# from the two databases
overlap_binding_sites <- function(hits, n, m, p, q, r, m6A) {
  
  over <-findOverlaps(hits,all_sites_RMBase_RMvar, ignore.strand=TRUE, type='any')
  hits_overlapping_specific_effector_m6A <- c()
  hits_overlapping_specific_effector_non_m6A <- c()
  # report the index of ELIGOS hits overlapping with specific effectors associated with m6A or specific effectors not associated
  # with m6A
  for (i in 1:length(over)) {
    if (all_sites_RMBase_RMvar[subjectHits(over[i])]$RBP_name %in% specific_effectors_only_m6A) {
      hits_overlapping_specific_effector_m6A <- c(hits_overlapping_specific_effector_m6A, queryHits(over[i]))
    }
    if (all_sites_RMBase_RMvar[subjectHits(over[i])]$RBP_name %in% specific_effectors_non_m6A) {
      hits_overlapping_specific_effector_non_m6A <- c(hits_overlapping_specific_effector_non_m6A, queryHits(over[i]))
    }
  }
  if (m6A) {
    number_hits_per_mod_RBP_DRACH[n,m] <<-paste0(as.character(length(unique(hits_overlapping_specific_effector_non_m6A))), ' - ', as.character(round(length(unique(hits_overlapping_specific_effector_non_m6A))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_DRACH[n,p] <<- paste0(as.character(length(unique(hits_overlapping_specific_effector_m6A))),' - ', as.character(round(length(unique(hits_overlapping_specific_effector_m6A))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_DRACH[n,q] <<-paste0(as.character(length(unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)))),' - ', as.character(round(length(unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_DRACH[n,r] <<- paste0(as.character(length(unique(queryHits(over)))), ' - ', as.character(round(length(unique(queryHits(over)))/length(hits)*100,2)),'%')
  } else {
    number_hits_per_mod_RBP_non_DRACH[n,m] <<-paste0(as.character(length(unique(hits_overlapping_specific_effector_non_m6A))), ' - ',as.character(round(length(unique(hits_overlapping_specific_effector_non_m6A))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_non_DRACH[n,p] <<- paste0(as.character(length(unique(hits_overlapping_specific_effector_m6A))),' - ', as.character(round(length(unique(hits_overlapping_specific_effector_m6A))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_non_DRACH[n,q] <<-paste0(as.character(length(unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)))),' - ',as.character(round(length(unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)))/length(hits)*100,2)),'%')
    number_hits_per_mod_RBP_non_DRACH[n,r] <<- paste0(as.character(length(unique(queryHits(over)))), ' - ',as.character(round(length(unique(queryHits(over)))/length(hits)*100,2)),'%')
  }
  
  # update the field "mod_type" of ELIGOS hits that didn't overlap with any RNA mark from the databases
  for (i in 1:length(hits)) {
    if (hits[i]$mod_type == 'Unknown') {
      if (i %in% queryHits(over)) {
        binding_sites <- unique(all_sites_RMBase_RMvar[subjectHits(over[queryHits(over)==i])])
        RBP_specific <- unique(binding_sites$RBP_name[which(binding_sites$RBP_name %in% specific_effectors)])
        RBP_ambiguous <- unique(binding_sites$RBP_name[which(!binding_sites$RBP_name %in% specific_effectors)])
        if (length(RBP_specific) == 0) {
          # if the hit overlaps only with one or more ambiguous effectors, the names of these effectors are added in the metadata separated by “;” 
          # and the modification type of the hit is annotated as “Ambiguous”
          hits[i]$effector <- paste(RBP_ambiguous, collapse = ';')
          hits[i]$mod_type <- 'Ambiguous'
          hits[i]$type_effector <- 'Ambiguous'
        } else if (length(RBP_specific) == length(unique(binding_sites$RBP_name))) {
          # if the hit overlaps only with one or more specific effectors (associated either with m6A or with another RNA mark), 
          # the names of these effectors are added in the metadata separated by “;” and the modification type of the hit is annotated with the 
          # names of the RNA marks related to these specific effectors separated by “;”
          hits[i]$effector <- paste(RBP_specific, collapse = ';')
          hits[i]$mod_type <- paste(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,3]), collapse=';')
          hits[i]$type_effector <- paste(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,2]), collapse=';')
        } else {
          # if the hit overlaps both with one or more specific effectors and with one or more ambiguous effectors, the names of these effectors 
          # are added in the metadata separated by “;” and the modification type of the hit is annotated with the names of the RNA marks related 
          # to the specific effectors with which the hit overlaps separated by “;” also including the modification type “Ambiguous”
          hits[i]$effector <- paste(c(RBP_specific,RBP_ambiguous), collapse = ';')
          hits[i]$mod_type <- paste(c(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,3]), 'Ambiguous'), collapse=';')
          hits[i]$type_effector <- paste(c(unique(enzyme_type_mod_unique[enzyme_type_mod_unique$Name_effector %in% RBP_specific,2]), 'Ambiguous'), collapse=';')
        }
      }
    }
  }
  
  hits$mod_type <- unlist(lapply(seq_along(hits), function(i,x) {
    paste(unique(unlist(strsplit(x[i]$mod_type, split =';'))), collapse = ';')
  }, x=hits))
  
  l <- list(unique(hits_overlapping_specific_effector_non_m6A), unique(hits_overlapping_specific_effector_m6A),unique(c(hits_overlapping_specific_effector_non_m6A,hits_overlapping_specific_effector_m6A)),unique(queryHits(over)), hits)
  names(l) <- c('specific non m6A', 'specific m6A','specific', 'total', 'hits')
  return(l)
}

# directory_hits is the path to the directory containing the results of ELIGOS analysis.
# create a new folder /mods_RBPs/ inside this directory to save the result of the overlap between
# ELIGOS hits of a fraction and the effectors' binding sites databases
ELIGOS_hits_DB_effectors <- function(directory_hits) {
  
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type.Rda'))
  
  # initiate a matrix to report, for each fraction, the percentage of ELIGOS DRACH- hit overlapping with each category of effectors and
  # the percentage of ELIGOS DRACH- hit overlapping at least with one category. For each category of effectors, the number of ELIGOS DRACH- hits 
  # overlapping with that category and shared by all the fractions is computed
  number_hits_per_mod_RBP_non_DRACH <<- matrix(0,nrow=5, ncol=5)
  rownames(number_hits_per_mod_RBP_non_DRACH) <<- c('Chromatin', 'Nucleoplasm', 'Cytoplasm', 'Overlap', 'Number of RBPs annotated')
  colnames(number_hits_per_mod_RBP_non_DRACH) <<- c('% of hits overlapping with\nspecific effectors non m6A', 
                                                    '% of hits overlapping with\nspecific effectors m6A',
                                                    '% of hits overlapping with\nspecific effectors',
                                                    '% of hits overlapping with\ntotal effectors',
                                                    'Tot number of DRACH- hits')
  
  # the same for ELIGOS DRACH+ hits
  number_hits_per_mod_RBP_DRACH <<- matrix(0,nrow=5, ncol=5)
  rownames(number_hits_per_mod_RBP_DRACH) <<- c('Chromatin', 'Nucleoplasm', 'Cytoplasm','Overlap','Number of RBPs annotated')
  colnames(number_hits_per_mod_RBP_DRACH) <<- c('% of hits overlapping with\nspecific effectors non m6A', 
                                                '% of hits overlapping with\nspecific effectors m6A',
                                                '% of hits overlapping with\nspecific effectors',
                                                '% of hits overlapping with\ntotal effectors',
                                                'Tot number of DRACH+ hits')
  
  number_hits_per_mod_RBP_non_DRACH[1,5] <<- length(hits_eligos_chr_ass_confirmed_5_without_DRACH)
  number_hits_per_mod_RBP_non_DRACH[2,5] <<- length(hits_eligos_nucleo_confirmed_5_without_DRACH)
  number_hits_per_mod_RBP_non_DRACH[3,5] <<- length(hits_eligos_cyto_confirmed_5_without_DRACH)
  number_hits_per_mod_RBP_non_DRACH[5,1] <<-length(specific_effectors_non_m6A)
  number_hits_per_mod_RBP_non_DRACH[5,2] <<-length(specific_effectors_only_m6A)
  number_hits_per_mod_RBP_non_DRACH[5,3] <<-length(specific_effectors)
  number_hits_per_mod_RBP_non_DRACH[5,4] <<-length(unique(enzyme_type_mod_unique[,1]))
  number_hits_per_mod_RBP_DRACH[1,5] <<- length(hits_eligos_chr_ass_confirmed_5_with_DRACH)
  number_hits_per_mod_RBP_DRACH[2,5] <<-length(hits_eligos_nucleo_confirmed_5_with_DRACH)
  number_hits_per_mod_RBP_DRACH[3,5] <<-length(hits_eligos_cyto_confirmed_5_with_DRACH)
  number_hits_per_mod_RBP_DRACH[5,1] <<-length(specific_effectors_non_m6A)
  number_hits_per_mod_RBP_DRACH[5,2] <<-length(specific_effectors_only_m6A)
  number_hits_per_mod_RBP_DRACH[5,3] <<-length(specific_effectors)
  number_hits_per_mod_RBP_DRACH[5,4] <<-length(unique(enzyme_type_mod_unique[,1]))
  
  # add a metadata column to each GRanges object in which the names of the effectors with which each hit overlaps will 
  # be reported as well as the function of these effectors (both initiated with "Unknown")
  mcols(hits_eligos_chr_ass_confirmed_5_without_DRACH) <- cbind(mcols(hits_eligos_chr_ass_confirmed_5_without_DRACH), effector = 'Unknown', type_effector = 'Unknown')
  mcols(hits_eligos_nucleo_confirmed_5_without_DRACH) <- cbind(mcols(hits_eligos_nucleo_confirmed_5_without_DRACH), effector = 'Unknown', type_effector = 'Unknown')
  mcols(hits_eligos_cyto_confirmed_5_without_DRACH) <- cbind(mcols(hits_eligos_cyto_confirmed_5_without_DRACH),effector = 'Unknown', type_effector = 'Unknown')
  mcols(hits_eligos_chr_ass_confirmed_5_with_DRACH) <- cbind(mcols(hits_eligos_chr_ass_confirmed_5_with_DRACH), effector = 'Unknown', type_effector = 'Unknown')
  mcols(hits_eligos_nucleo_confirmed_5_with_DRACH) <- cbind(mcols(hits_eligos_nucleo_confirmed_5_with_DRACH), effector = 'Unknown', type_effector = 'Unknown')
  mcols(hits_eligos_cyto_confirmed_5_with_DRACH) <- cbind(mcols(hits_eligos_cyto_confirmed_5_with_DRACH), effector = 'Unknown', type_effector = 'Unknown')
  
  results_chr <- overlap_binding_sites(hits_eligos_chr_ass_confirmed_5_without_DRACH,1,1,2,3,4, FALSE)
  chr_specific_non_m6A <- hits_eligos_chr_ass_confirmed_5_without_DRACH[results_chr[[1]]]
  chr_specific_m6A <- hits_eligos_chr_ass_confirmed_5_without_DRACH[results_chr[[2]]]
  chr_specific <- hits_eligos_chr_ass_confirmed_5_without_DRACH[results_chr[[3]]]
  chr_total <- hits_eligos_chr_ass_confirmed_5_without_DRACH[results_chr[[4]]]
  hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings <- results_chr[[5]]
  
  results_nucleo <- overlap_binding_sites(hits_eligos_nucleo_confirmed_5_without_DRACH,2,1,2,3,4, FALSE)
  nucleo_specific_non_m6A <- hits_eligos_nucleo_confirmed_5_without_DRACH[results_nucleo[[1]]]
  nucleo_specific_m6A <- hits_eligos_nucleo_confirmed_5_without_DRACH[results_nucleo[[2]]]
  nucleo_specific <- hits_eligos_nucleo_confirmed_5_without_DRACH[results_nucleo[[3]]]
  nucleo_total <- hits_eligos_nucleo_confirmed_5_without_DRACH[results_nucleo[[4]]]
  hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings <- results_nucleo[[5]]
  
  results_cyto <- overlap_binding_sites(hits_eligos_cyto_confirmed_5_without_DRACH,3,1,2,3,4, FALSE)
  cyto_specific_non_m6A <- hits_eligos_cyto_confirmed_5_without_DRACH[results_cyto[[1]]]
  cyto_specific_m6A <- hits_eligos_cyto_confirmed_5_without_DRACH[results_cyto[[2]]]
  cyto_specific <- hits_eligos_cyto_confirmed_5_without_DRACH[results_cyto[[3]]]
  cyto_total <- hits_eligos_cyto_confirmed_5_without_DRACH[results_cyto[[4]]]
  hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings <- results_cyto[[5]]
  
  # compute, for each category of effectors, the number of ELIGOS DRACH- hit overlapping with that category and 
  # shared by all the fractions 
  all_hits_specific_non_m6A <- list(chr_specific_non_m6A,nucleo_specific_non_m6A,cyto_specific_non_m6A)
  order <- order(c(length(chr_specific_non_m6A), length(nucleo_specific_non_m6A), length(cyto_specific_non_m6A)))
  overlap <- findOverlaps(all_hits_specific_non_m6A[[order[1]]],all_hits_specific_non_m6A[[order[2]]], type = 'any')
  overlap2 <- all_hits_specific_non_m6A[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits_specific_non_m6A[[order[3]]], type = 'any')
  number_hits_per_mod_RBP_non_DRACH[4,1] <- length(overlap2[unique(queryHits(overlap_all))])
  
  all_hits_specific_m6A <- list(chr_specific_m6A,nucleo_specific_m6A,cyto_specific_m6A)
  order <- order(c(length(chr_specific_m6A), length(nucleo_specific_m6A), length(cyto_specific_m6A)))
  overlap <- findOverlaps(all_hits_specific_m6A[[order[1]]],all_hits_specific_m6A[[order[2]]], type = 'any')
  overlap2 <- all_hits_specific_m6A[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits_specific_m6A[[order[3]]], type = 'any')
  number_hits_per_mod_RBP_non_DRACH[4,2] <- length(overlap2[unique(queryHits(overlap_all))])
  
  all_hits_specific <- list(chr_specific,nucleo_specific,cyto_specific)
  order <- order(c(length(chr_specific), length(nucleo_specific), length(cyto_specific)))
  overlap <- findOverlaps(all_hits_specific[[order[1]]],all_hits_specific[[order[2]]], type = 'any')
  overlap2 <- all_hits_specific[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits_specific[[order[3]]], type = 'any')
  number_hits_per_mod_RBP_non_DRACH[4,3] <- length(overlap2[unique(queryHits(overlap_all))])
  
  all_hits_total <- list(chr_total,nucleo_total,cyto_total)
  order <- order(c(length(chr_total), length(nucleo_total), length(cyto_total)))
  overlap <- findOverlaps(all_hits_total[[order[1]]],all_hits_total[[order[2]]], type = 'any')
  overlap2 <- all_hits_total[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits_total[[order[3]]], type = 'any')
  number_hits_per_mod_RBP_non_DRACH[4,4] <- length(overlap2[unique(queryHits(overlap_all))])
  
  results_chr <- overlap_binding_sites(hits_eligos_chr_ass_confirmed_5_with_DRACH,1,1,2,3,4, TRUE)
  chr_specific_non_m6A <- hits_eligos_chr_ass_confirmed_5_with_DRACH[results_chr[[1]]]
  chr_specific_m6A <- hits_eligos_chr_ass_confirmed_5_with_DRACH[results_chr[[2]]]
  chr_specific <- hits_eligos_chr_ass_confirmed_5_with_DRACH[results_chr[[3]]]
  chr_total <- hits_eligos_chr_ass_confirmed_5_with_DRACH[results_chr[[4]]]
  hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings <- results_chr[[5]]
  
  results_nucleo <- overlap_binding_sites(hits_eligos_nucleo_confirmed_5_with_DRACH,2,1,2,3,4, TRUE)
  nucleo_specific_non_m6A <- hits_eligos_nucleo_confirmed_5_with_DRACH[results_nucleo[[1]]]
  nucleo_specific_m6A <- hits_eligos_nucleo_confirmed_5_with_DRACH[results_nucleo[[2]]]
  nucleo_specific <- hits_eligos_nucleo_confirmed_5_with_DRACH[results_nucleo[[3]]]
  nucleo_total <- hits_eligos_nucleo_confirmed_5_with_DRACH[results_nucleo[[4]]]
  hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings <- results_nucleo[[5]]
  
  results_cyto <- overlap_binding_sites(hits_eligos_cyto_confirmed_5_with_DRACH,3,1,2,3,4, TRUE)
  cyto_specific_non_m6A <- hits_eligos_cyto_confirmed_5_with_DRACH[results_cyto[[1]]]
  cyto_specific_m6A <- hits_eligos_cyto_confirmed_5_with_DRACH[results_cyto[[2]]]
  cyto_specific <- hits_eligos_cyto_confirmed_5_with_DRACH[results_cyto[[3]]]
  cyto_total <- hits_eligos_cyto_confirmed_5_with_DRACH[results_cyto[[4]]]
  hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings <- results_cyto[[5]]
  
  # compute, for each category of effectors, the number of ELIGOS DRACH+ hits overlapping with that category 
  # and shared by all the fractions 
  all_hits_specific_non_m6A <- list(chr_specific_non_m6A,nucleo_specific_non_m6A,cyto_specific_non_m6A)
  order <- order(c(length(chr_specific_non_m6A), length(nucleo_specific_non_m6A), length(cyto_specific_non_m6A)))
  overlap <- findOverlaps(all_hits_specific_non_m6A[[order[1]]],all_hits_specific_non_m6A[[order[2]]], type = 'any')
  overlap2 <- all_hits_specific_non_m6A[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits_specific_non_m6A[[order[3]]], type = 'any')
  number_hits_per_mod_RBP_DRACH[4,1] <- length(overlap2[unique(queryHits(overlap_all))])
  
  all_hits_specific_m6A <- list(chr_specific_m6A,nucleo_specific_m6A,cyto_specific_m6A)
  order <- order(c(length(chr_specific_m6A), length(nucleo_specific_m6A), length(cyto_specific_m6A)))
  overlap <- findOverlaps(all_hits_specific_m6A[[order[1]]],all_hits_specific_m6A[[order[2]]], type = 'any')
  overlap2 <- all_hits_specific_m6A[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits_specific_m6A[[order[3]]], type = 'any')
  number_hits_per_mod_RBP_DRACH[4,2] <- length(overlap2[unique(queryHits(overlap_all))])
  
  all_hits_specific <- list(chr_specific,nucleo_specific,cyto_specific)
  order <- order(c(length(chr_specific), length(nucleo_specific), length(cyto_specific)))
  overlap <- findOverlaps(all_hits_specific[[order[1]]],all_hits_specific[[order[2]]], type = 'any')
  overlap2 <- all_hits_specific[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits_specific[[order[3]]], type = 'any')
  number_hits_per_mod_RBP_DRACH[4,3] <- length(overlap2[unique(queryHits(overlap_all))])
  
  all_hits_total <- list(chr_total,nucleo_total,cyto_total)
  order <- order(c(length(chr_total), length(nucleo_total), length(cyto_total)))
  overlap <- findOverlaps(all_hits_total[[order[1]]],all_hits_total[[order[2]]], type = 'any')
  overlap2 <- all_hits_total[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,all_hits_total[[order[3]]], type = 'any')
  number_hits_per_mod_RBP_DRACH[4,4] <- length(overlap2[unique(queryHits(overlap_all))])
  
  write.xlsx(x = data.frame(number_hits_per_mod_RBP_non_DRACH),file = paste0(directory_hits, '/mods_RBPs/', 'hits_overlapping_RBPs_non_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
  write.xlsx(x = data.frame(number_hits_per_mod_RBP_DRACH),file = paste0(directory_hits, '/mods_RBPs/', 'hits_overlapping_RBPs_DRACH.xlsx'),col.names = TRUE, row.names = TRUE)
  
  save(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings, file =paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  save(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings, file=paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  save(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings, file=paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  save(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings, file=paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  save(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings, file=paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  save(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings, file=paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
}

ELIGOS_hits_DB_effectors(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/')

ELIGOS_hits_DB_effectors(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/')

ELIGOS_hits_DB_effectors(directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/')

ELIGOS_hits_DB_effectors(directory_hits = '/path/to/fractions_eligos_STORM_K562/eligos_total_reads/')

#############
# combine the information from the overlap with the RNA marks and the overlap with the binding sites of the
# effectors from the databases.
# directory_hits is the path to the directory containing the results of ELIGOS analysis
final_summary <- function(all_mods, directory_hits) {
  
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
  
  # initiate a matrix reporting, for each RNA mod type and each fraction, the number (and percentage) of ELIGOS DRACH- hits annotated 
  # to that RNA mod. Report also the number of ELIGOS DRACH- hits annotated at least with one RNA mod and 
  # the number of ELIGOS DRACH- hits annotated to the same RNA mod and shared by all the fractions
  number_hits_non_DRACH <- matrix(0,nrow=4, ncol=length(all_mods)+2)
  rownames(number_hits_non_DRACH) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm', 'Overlap')
  colnames(number_hits_non_DRACH) <- c(all_mods, 'Tot number of DRACH- annotated hits', 'Tot number of DRACH- hits')
  
  number_hits_non_DRACH[1,'Tot number of DRACH- hits'] <- length(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings)
  number_hits_non_DRACH[2,'Tot number of DRACH- hits'] <- length(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings)
  number_hits_non_DRACH[3,'Tot number of DRACH- hits'] <- length(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings)
  
  # initiate a matrix reporting, for each RNA mod type and each fraction, the number (and percentage) of ELIGOS DRACH+ hits annotated 
  # to that RNA mod. Report also the number of ELIGOS DRACH+ hits annotated at least with one RNA mod and 
  # the number of ELIGOS DRACH+ hits annotated to the same RNA mod and shared by all the fractions
  number_hits_DRACH <- matrix(0,nrow=4, ncol=length(all_mods)+2)
  rownames(number_hits_DRACH) <- c('Chromatin', 'Nucleoplasm', 'Cytoplasm', 'Overlap')
  colnames(number_hits_DRACH) <- c(all_mods, 'Tot number of DRACH+ annotated hits', 'Tot number of DRACH+ hits')
  
  number_hits_DRACH[1,'Tot number of DRACH+ hits'] <- length(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings)
  number_hits_DRACH[2,'Tot number of DRACH+ hits'] <- length(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings)
  number_hits_DRACH[3,'Tot number of DRACH+ hits'] <- length(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings)
  
  count_mods_non_DRACH <- function(hits_chr, hits_nucleo, hits_cyto) {
    mod_types_chr <- unlist(strsplit(hits_chr$mod_type, split=';'))
    mod_types_nucleo <- unlist(strsplit(hits_nucleo$mod_type, split=';'))
    mod_types_cyto <- unlist(strsplit(hits_cyto$mod_type, split=';'))
    # iterate over the modification types
    for (mod in colnames(number_hits_non_DRACH)[1:length(all_mods)]) {
      # for each fraction, compute the number and percentage of ELIGOS DRACH- hits annotated to that RNA mod type
      number_hits_non_DRACH[1,mod] <<- paste0(as.character(length(mod_types_chr[mod_types_chr==mod])),' - ', as.character(round(length(mod_types_chr[mod_types_chr==mod])*100/as.numeric(number_hits_non_DRACH[1,'Tot number of DRACH- hits']),2)),'%')
      number_hits_non_DRACH[2,mod] <<- paste0(as.character(length(mod_types_nucleo[mod_types_nucleo==mod])),' - ', as.character(round(length(mod_types_nucleo[mod_types_nucleo==mod])*100/as.numeric(number_hits_non_DRACH[2,'Tot number of DRACH- hits']),2)),'%')
      number_hits_non_DRACH[3,mod] <<- paste0(as.character(length(mod_types_cyto[mod_types_cyto==mod])),' - ', as.character(round(length(mod_types_cyto[mod_types_cyto==mod])*100/as.numeric(number_hits_non_DRACH[3,'Tot number of DRACH- hits']),2)),'%')
      
      # compute how many ELIGOS DRACH- hits annotated to the same RNA mod type are shared by all the fractions
      hits_mod_chr <- hits_chr[unlist(lapply(seq_along(hits_chr),function(i) {mod %in% unlist(strsplit(hits_chr[i]$mod_type, split=';'))}))] 
      hits_mod_nucleo <- hits_nucleo[unlist(lapply(seq_along(hits_nucleo),function(i) {mod %in% unlist(strsplit(hits_nucleo[i]$mod_type, split=';'))}))] 
      hits_mod_cyto <-  hits_cyto[unlist(lapply(seq_along(hits_cyto),function(i) {mod %in% unlist(strsplit(hits_cyto[i]$mod_type, split=';'))}))] 
      hits_mod <- list(hits_mod_chr, hits_mod_nucleo, hits_mod_cyto)
      order <- order(c(length(hits_mod_chr), length(hits_mod_cyto), length(hits_mod_cyto)))
      
      overlap <- findOverlaps(hits_mod[[order[1]]],hits_mod[[order[2]]], type = 'any')
      overlap2 <- hits_mod[[order[1]]][unique(queryHits(overlap))]
      overlap_all <- findOverlaps(overlap2,hits_mod[[order[3]]], type = 'any')
      number_hits_non_DRACH[4,mod] <<- length(overlap2[unique(queryHits(overlap_all))])
    }
  }
  
  count_mods_DRACH <- function(hits_chr, hits_nucleo, hits_cyto) {
    mod_types_chr <- unlist(strsplit(hits_chr$mod_type, split=';'))
    mod_types_nucleo <- unlist(strsplit(hits_nucleo$mod_type, split=';'))
    mod_types_cyto <- unlist(strsplit(hits_cyto$mod_type, split=';'))
    # iterate over the modification types
    for (mod in colnames(number_hits_DRACH)[1:length(all_mods)]) {
      # for each fraction, compute the number and percentage of ELIGOS DRACH+ hits annotated to that RNA mod type
      number_hits_DRACH[1,mod] <<- paste0(as.character(length(mod_types_chr[mod_types_chr==mod])),' - ', as.character(round(length(mod_types_chr[mod_types_chr==mod])*100/as.numeric(number_hits_DRACH[1,'Tot number of DRACH+ hits']),2)),'%')
      number_hits_DRACH[2,mod] <<- paste0(as.character(length(mod_types_nucleo[mod_types_nucleo==mod])),' - ', as.character(round(length(mod_types_nucleo[mod_types_nucleo==mod])*100/as.numeric(number_hits_DRACH[2,'Tot number of DRACH+ hits']),2)),'%')
      number_hits_DRACH[3,mod] <<- paste0(as.character(length(mod_types_cyto[mod_types_cyto==mod])),' - ', as.character(round(length(mod_types_cyto[mod_types_cyto==mod])*100/as.numeric(number_hits_DRACH[3,'Tot number of DRACH+ hits']),2)),'%')
      
      # compute how many ELIGOS DRACH+ hits annotated to the same RNA mod type are shared by all the fractions
      hits_mod_chr <- hits_chr[unlist(lapply(seq_along(hits_chr),function(i) {mod %in% unlist(strsplit(hits_chr[i]$mod_type, split=';'))}))] 
      hits_mod_nucleo <- hits_nucleo[unlist(lapply(seq_along(hits_nucleo),function(i) {mod %in% unlist(strsplit(hits_nucleo[i]$mod_type, split=';'))}))] 
      hits_mod_cyto <-  hits_cyto[unlist(lapply(seq_along(hits_cyto),function(i) {mod %in% unlist(strsplit(hits_cyto[i]$mod_type, split=';'))}))] 
      hits_mod <- list(hits_mod_chr, hits_mod_nucleo, hits_mod_cyto)
      order <- order(c(length(hits_mod_chr), length(hits_mod_cyto), length(hits_mod_cyto)))
      
      overlap <- findOverlaps(hits_mod[[order[1]]],hits_mod[[order[2]]], type = 'any')
      overlap2 <- hits_mod[[order[1]]][unique(queryHits(overlap))]
      overlap_all <- findOverlaps(overlap2,hits_mod[[order[3]]], type = 'any')
      number_hits_DRACH[4,mod] <<- length(overlap2[unique(queryHits(overlap_all))])
    }
  }
  
  count_mods_non_DRACH(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings,hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings,hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings)
  count_mods_DRACH(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings,hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings,hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings)
  
  # remove ELIGOS DRACH- hits that don't overlap with any RNA mark/effectors' binding site (Unknown) and those that overlap only with one or more ambiguous
  # effectors' binding sites (Ambiguous)
  hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings_marks_specific <- hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings[!hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings$mod_type %in% c('Unknown', 'Ambiguous')]
  hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings_marks_specific <- hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings[!hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings$mod_type %in% c('Unknown', 'Ambiguous')]
  hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings_marks_specific <- hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings[!hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings$mod_type %in% c('Unknown', 'Ambiguous')]
  
  # compute how many ELIGOS DRACH- hits are annotated at least with one RNA mark or overlap with the binding site of at least one specific effector
  number_hits_non_DRACH[1,'Tot number of DRACH- annotated hits'] <- paste0(as.character(length(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings_marks_specific)), ' - ', as.character(round(length(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings_marks_specific)*100/as.numeric(number_hits_non_DRACH[1,'Tot number of DRACH- hits']),2)),'%')
  number_hits_non_DRACH[2,'Tot number of DRACH- annotated hits'] <- paste0(as.character(length(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings_marks_specific)), ' - ', as.character(round(length(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings_marks_specific)*100/as.numeric(number_hits_non_DRACH[2,'Tot number of DRACH- hits']),2)),'%')
  number_hits_non_DRACH[3,'Tot number of DRACH- annotated hits'] <- paste0(as.character(length(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings_marks_specific)), ' - ', as.character(round(length(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings_marks_specific)*100/as.numeric(number_hits_non_DRACH[3,'Tot number of DRACH- hits']),2)),'%')
  
  # compute how many ELIGOS DRACH- hits are annotated at least with one RNA mark or overlap with the binding site of at least one specific effector
  # and are shared by all the fractions
  hits <- list(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings_marks_specific, hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings_marks_specific, hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings_marks_specific)
  order <- order(c(length(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings_marks_specific), length(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings_marks_specific), length(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings_marks_specific)))
  
  overlap <- findOverlaps(hits[[order[1]]],hits[[order[2]]], type = 'any')
  overlap2 <- hits[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,hits[[order[3]]], type = 'any')
  number_hits_non_DRACH[4,'Tot number of DRACH- annotated hits'] <- length(overlap2[unique(queryHits(overlap_all))])
  
  # remove ELIGOS DRACH+ hits that don't overlap with any RNA mark/effectors' binding site (Unknown) and those that overalp only with one or more ambiguous
  # effectors' binding sites (Ambiguous)
  hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings_marks_specific <- hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings[!hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings$mod_type  %in% c('Unknown', 'Ambiguous')]
  hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings_marks_specific <- hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings[!hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings$mod_type %in% c('Unknown', 'Ambiguous')]
  hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings_marks_specific <- hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings[!hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings$mod_type %in% c('Unknown', 'Ambiguous')]
  
  # compute how many ELIGOS DRACH+ hits are annotated at least with one RNA mark or overlap with the binding site of at least one specific effector
  number_hits_DRACH[1,'Tot number of DRACH+ annotated hits'] <- paste0(as.character(length(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings_marks_specific)), ' - ', as.character(round(length(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings_marks_specific)*100/as.numeric(number_hits_DRACH[1,'Tot number of DRACH+ hits']),2)),'%')
  number_hits_DRACH[2,'Tot number of DRACH+ annotated hits'] <- paste0(as.character(length(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings_marks_specific)), ' - ', as.character(round(length(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings_marks_specific)*100/as.numeric(number_hits_DRACH[2,'Tot number of DRACH+ hits']),2)),'%')
  number_hits_DRACH[3,'Tot number of DRACH+ annotated hits'] <- paste0(as.character(length(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings_marks_specific)), ' - ', as.character(round(length(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings_marks_specific)*100/as.numeric(number_hits_DRACH[3,'Tot number of DRACH+ hits']),2)),'%')
  
  # compute how many ELIGOS DRACH+ hits are annotated at least with one RNA mark or overlap with the binding site of at least one specific effector
  # and are shared by all the fractions
  hits <- list(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings_marks_specific, hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings_marks_specific, hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings_marks_specific)
  order <- order(c(length(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings_marks_specific), length(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings_marks_specific), length(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings_marks_specific)))
  
  overlap <- findOverlaps(hits[[order[1]]],hits[[order[2]]], type = 'any')
  overlap2 <- hits[[order[1]]][unique(queryHits(overlap))]
  overlap_all <- findOverlaps(overlap2,hits[[order[3]]], type = 'any')
  number_hits_DRACH[4,'Tot number of DRACH+ annotated hits'] <- length(overlap2[unique(queryHits(overlap_all))])
  
  write.xlsx(x = data.frame(number_hits_non_DRACH),file = paste0(directory_hits, '/mods_RBPs/', 'counts_mod_RBPs_non_DRACH.xlsx'),col.names = TRUE, row.names=TRUE)
  write.xlsx(x = data.frame(number_hits_DRACH),file = paste0(directory_hits, '/mods_RBPs/', 'counts_mod_RBPs_DRACH.xlsx'),col.names = TRUE, row.names=TRUE)
}

# all the reads
directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/'

load(paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))

# identify which are the RNA mod types with which ELIGOS hits are annotated after the overlap with the RNA marks and 
# the effectors from the databases
all_mods <- c(unique(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings$mod_type),
              unique(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings$mod_type),
              unique(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings$mod_type))
all_mods <- unique(all_mods)
all_mods <- unique(unlist(strsplit(all_mods, split=';')))
all_mods <- c('m6A', 'Y', 'm1A','m5C', 'm7G', 'A-I','Nm', 'm6Am', 'Ambiguous')

final_summary(all_mods, directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_min05_min05_mag1/')

# nascent reads
directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/'

load(paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))

all_mods <- c(unique(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings$mod_type),
              unique(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings$mod_type),
              unique(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings$mod_type))
all_mods <- unique(all_mods)
all_mods <- unique(unlist(strsplit(all_mods, split=';')))
all_mods <- c('m6A', 'Y', 'm1A','m5C', 'm7G', 'A-I','Nm', 'Ambiguous')

final_summary(all_mods,directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_nascent_min05_min05_mag1/')

# all the reads with the same library-level subsampling threshold used for nascent reads
directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/'

load(paste0(directory_hits,'/without_DRACH/hits_eligos_chr_ass_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/without_DRACH/hits_eligos_nucleo_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/without_DRACH/hits_eligos_cyto_confirmed_5_without_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_chr_ass_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_nucleo_confirmed_5_with_DRACH_mod_type_RBP.Rda'))
load(paste0(directory_hits,'/hits_ELIGOS/hits_eligos_cyto_confirmed_5_with_DRACH_mod_type_RBP.Rda'))

all_mods <- c(unique(hits_eligos_chr_ass_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_nucleo_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_cyto_confirmed_5_without_DRACH_with_bindings$mod_type),
              unique(hits_eligos_chr_ass_confirmed_5_with_DRACH_with_bindings$mod_type),
              unique(hits_eligos_nucleo_confirmed_5_with_DRACH_with_bindings$mod_type),
              unique(hits_eligos_cyto_confirmed_5_with_DRACH_with_bindings$mod_type))
all_mods <- unique(all_mods)
all_mods <- unique(unlist(strsplit(all_mods, split=';')))
all_mods <- c('m6A', 'Y', 'm1A','m5C', 'm7G', 'A-I','Nm', 'm6Am', 'Ambiguous')

final_summary(all_mods,directory_hits = '/path/to/fractions_eligos_4sU_library_gene_subsampling_total_510645_min05_min05_mag1/')

