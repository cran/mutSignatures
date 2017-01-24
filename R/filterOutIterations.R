filterOutIterations <-
function (wall, hall, cnt.errors, cnt.reconstructed, params) {
  #
  num.processes.toextract <-  params$num.processes.toextract
  remove.last.percent <- params$remove.last.percent
  #
  tot.iterations <- ncol(wall)/num.processes.toextract
  tot.rm.iterations <- round(remove.last.percent * tot.iterations)
  if (tot.rm.iterations > 0) {
    closeness.mutCounts <- matrix(0, nrow = tot.iterations, ncol = 1)
    for (i in 1:tot.iterations) {
      closeness.mutCounts[i,] <- norm(cnt.errors[[i]], 'F')
    }
    indexClosenessGenomes <- order(closeness.mutCounts, decreasing = TRUE)
    removeIterations <- indexClosenessGenomes[1:tot.rm.iterations]
    
    removeIterationSets <- matrix(0, nrow = (num.processes.toextract * tot.rm.iterations), ncol = 1)
    for (i in 1: tot.rm.iterations){
      iStart <- num.processes.toextract * ( removeIterations[i] - 1) + 1
      iEnd <- num.processes.toextract * removeIterations[i]
      removeIterationSets[(num.processes.toextract * ( i - 1) + 1):(num.processes.toextract * i) ,] <- iStart:iEnd
    }
    wall <- wall[,-removeIterationSets]
    hall <- hall[-removeIterationSets, ]
    cnt.errors <- cnt.errors[-removeIterations]
    cnt.reconstructed <- cnt.reconstructed[-removeIterations]
  }
  #
  res.list <- list()
  res.list$Wall <- wall
  res.list$Hall <- hall
  res.list$mutCounts.errors <- cnt.errors
  res.list$mutCounts.reconstructed <- cnt.reconstructed
  return(res.list)
}
