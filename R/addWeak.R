addWeak <-
function(mutationTypesToAddSet,
                    processes_I,
                    processesStd_I,
                    Wall_I,
                    genomeErrors_I,
                    genomesReconstructed_I)
{
  # check if there are indeed low-freq mutations to add back in
  # if not, skip the whole process and return inputs
  #
  if(length(mutationTypesToAddSet) > 0 & mutationTypesToAddSet[1] > 0) {
    totalMutTypes <- nrow(Wall_I) + length(mutationTypesToAddSet)
    processes <- matrix(0, nrow = totalMutTypes, ncol = ncol(processes_I))
    processesStd <- matrix(0, nrow = totalMutTypes, ncol = ncol(processesStd_I))
    Wall <- matrix(0, nrow = totalMutTypes, ncol = ncol(Wall_I))
    genomeErrors <- lapply(1: length(genomeErrors_I), (function(i){
      matrix(0, nrow = totalMutTypes, ncol = ncol(genomeErrors_I[[1]]))
    }))
    genomesReconstructed <- lapply(1: length(genomesReconstructed_I), (function(i){
      matrix(0, nrow = totalMutTypes, ncol = ncol(genomesReconstructed_I[[1]]))
    }))
    #
    origArrayIndex <- 1
    for (i in 1 : totalMutTypes) {
      #
      if (! (i %in% mutationTypesToAddSet)) {
        #
        processes[i,] <- processes_I[origArrayIndex, ]
        processesStd[i, ] <- processesStd_I[origArrayIndex, ]
        Wall[i, ] <- Wall_I[origArrayIndex, ]
        for(j in 1: length(genomeErrors_I)) {
          genomeErrors[[j]][i,] <- genomeErrors_I[[j]][origArrayIndex,]
        }
        for(j in 1: length(genomesReconstructed_I)) {
          genomesReconstructed[[j]][i,] <- genomesReconstructed_I[[j]][origArrayIndex, ]
        }
        origArrayIndex <- origArrayIndex + 1;
      }
    }
  } else {
    processes <- processes_I
    processesStd <- processesStd_I
    Wall <- Wall_I
    genomeErrors <- genomeErrors_I
    genomesReconstructed <- genomesReconstructed_I
  }
  #
  # Wrap up
  weakAdded.list <- list()
  weakAdded.list$processes <- processes
  weakAdded.list$processesStd <- processesStd
  weakAdded.list$Wall <- Wall
  weakAdded.list$mutCountErrors <- genomeErrors
  weakAdded.list$mutCountReconstructed <- genomesReconstructed
  return(weakAdded.list)
}
