bootstrapCancerGenomes <-
function(genomes) {
  genome.col.sums <- apply(genomes, 2, sum)
  norm.genomes <- genomes / matrix(genome.col.sums,
                                   ncol = ncol(genomes),
                                   nrow = nrow(genomes),
                                   byrow = TRUE)
  bootstrapGenomes <- sapply(1:length(genome.col.sums), (function(i){
    stats::rmultinom(1,
                     size = genome.col.sums[i],
                     prob = norm.genomes[,i])
  }))
  #
  return(bootstrapGenomes)
}
