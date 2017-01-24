extractSignatures <-
function(mutCountMatrix, params, verbose = TRUE) {
  #
  tot.iterations <- params$tot.iterations
  num.processes.toextract <- params$num.processes.toextract
  eps <- params$eps
  #
  total.mutationTypes <- nrow(mutCountMatrix)
  total.samples <- ncol(mutCountMatrix)
  #
  Wall <- matrix(0, nrow = total.mutationTypes, ncol = num.processes.toextract * tot.iterations)
  Hall <- matrix(0, nrow = num.processes.toextract * tot.iterations, ncol = total.samples)
  mutCountMatrix.errors <- lapply(1: tot.iterations, (function(i){
    matrix(0, nrow = total.mutationTypes, ncol = total.samples)
  }))
  mutCountMatrix.reconstructed <- lapply(1: tot.iterations, (function(i){
    matrix(0, nrow = total.mutationTypes, ncol = total.samples)
  }))
  #
  processCount <- 1
  #
  for (j in 1:tot.iterations) {
    #
    # message(paste("Iter#...", leadZeros(j, tot.iterations), ": ", sep = ""), appendLF = FALSE)
    message(">", appendLF = FALSE)
    #
    # generate bootstrapped mut counts (genomes)
    bstrpd.result <- bootstrapCancerGenomes(mutCountMatrix)  
    for(iii in 1:ncol(bstrpd.result)){
      bstrpd.result[bstrpd.result[,iii]<eps,iii] <- eps
    }
    # solving NMF for these mut counts
    nmf.results <- do.nmf(v = bstrpd.result, r = num.processes.toextract, params = params)
    tmp.w <- nmf.results$w
    tmp.h <- nmf.results$h
    #
    for (jj in 1: num.processes.toextract) {
      tmp.tot <- sum(tmp.w[,jj])
      tmp.w[,jj] <- tmp.w[,jj] / tmp.tot
      tmp.h[jj,] <- tmp.h[jj,] * tmp.tot
    }
    #
    mutCountMatrix.errors[[j]] <- bstrpd.result - (tmp.w %*% tmp.h)
    mutCountMatrix.reconstructed[[j]] <- tmp.w %*% tmp.h
    #
    Wall[,processCount : (processCount + num.processes.toextract - 1)] <- tmp.w
    Hall[processCount : (processCount + num.processes.toextract - 1),] <- tmp.h
    processCount <- processCount + num.processes.toextract
    #
  }
  result.list <- list()
  result.list$Wall <- Wall
  result.list$Hall <- Hall
  result.list$mutCounts.errors <- mutCountMatrix.errors
  result.list$mutCounts.reconstructed <- mutCountMatrix.reconstructed
  message("", appendLF = TRUE)
  return(result.list)
}
