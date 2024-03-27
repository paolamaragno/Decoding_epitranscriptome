# create a folder /stech/ in which saving the results of m6Anet on each of the 5 sets of reads coming from the 5 gene-level subsamplings
# for each cellular fraction. In particular save the tsv file with the coordinates of the sites analysed by m6Anet on the 
# transcriptomic coordinates

m6Anet_K562 <- list.files(path = '/path/to/fractions_m6anet_K562_all_reads_prob0.75/stech/', pattern = 'tsv', full.names = TRUE)
m6Anet_SUM159 <- list.files(path = '/path/to/fractions_m6anet_4sU_library_gene_subsampling_prob0.75/stech/', pattern = 'tsv', full.names = TRUE)

# function to retrieve the stoichiometry of the sites analysed by m6Anet
m6Anet_output <- function(dir) {
  
  stech <- c()
  
  for (file in dir) {
    t <- read.table(file, header=TRUE, sep = ',')
    
    stech <- c(stech, t$mod_ratio)
  }
  
 return(stech)
}

stech_k562 <- m6Anet_output(m6Anet_K562)
stech_SUM159 <- m6Anet_output(m6Anet_SUM159)

pdf('/path/to/stechiometry.pdf')
par(mfrow = c(1,2))

boxplot(stech_k562, main = ("Stoichiometry m6Anet\nanalysed sites - K562"))
boxplot(stech_SUM159, main = ("Stoichiometry m6Anet\nanalysed sites - SUM159"))

dev.off()
