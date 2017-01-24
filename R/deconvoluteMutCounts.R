deconvoluteMutCounts <-
function(input.mutCounts, params) {
  # params
  tot.iterations <- params$tot.iterations
  num.processes.toextract <- params$num.processes.toextract
  process.distance <- params$process.distance
  remove.weak.muttypes <- params$remove.weak.muttypes
  tot.cores <- params$tot.cores
  tot.Replicates <- params$tot.Replicates
  remove.last.percent <- params$remove.last.percent
  #
  # Check and prepare tables/matrices
  #
  colnames(input.mutCounts) <- NULL
  rownames(input.mutCounts) <- NULL
  input.mutCounts <- as.matrix(input.mutCounts)
  #
  bckgrnd.removed.mutCounts <- removeWeak(input.mutCounts, params)
  bckgrnd.removed.mutset <- bckgrnd.removed.mutCounts$removed.mutset
  bckgrnd.removed.mutCounts <- bckgrnd.removed.mutCounts$output.mutCounts
  #
  # take care of signatures to be trimmed
  #
  total.mutationTypes <- nrow(bckgrnd.removed.mutCounts)
  total.samples <- ncol(bckgrnd.removed.mutCounts)
  #
  # start parallel processing (eventually)
  #
  max.cores <- parallel::detectCores()
  max.cores <- max.cores - 1
  max.cores <- ifelse(max.cores < 1, 1, max.cores)
  use.cores <- ifelse(1 <= tot.cores & tot.cores <= max.cores, tot.cores, max.cores)
  cl <- parallel::makeCluster(use.cores, outfile = "")
  message(paste("Extracting", num.processes.toextract, "mutational signatures X", tot.iterations, "iterations on", use.cores, "cores"))
  message(paste("Total cycles [iterat. * cores]:", (tot.iterations * use.cores), " {This should be a number >= 1000}"))
  #
  doParallel::registerDoParallel(cl)
  stuffToExp <- c("do.nmf", "leadZeros", "extractSignatures", "bootstrapCancerGenomes")
  parallel::clusterExport(cl, stuffToExp)
  #
  muCounts.check <- tryCatch(foreach::foreach(j = (1:use.cores), .verbose = TRUE, .packages = "stats") %dopar%
  { #message(paste("Core", j, "up and running... Iteration #:", params$tot.iterations))
    extractSignatures(mutCountMatrix = bckgrnd.removed.mutCounts, params = params) },
  error = (function(e){ print(e) }),
  finally = (function(f){parallel::stopCluster(cl)}))
  #
  # Build look-up tables
  nu.Wall <- matrix(0, nrow = total.mutationTypes, ncol = num.processes.toextract * tot.iterations* length(muCounts.check))
  nu.Hall <- matrix(0, nrow = num.processes.toextract * tot.iterations * length(muCounts.check), ncol = total.samples)
  #
  nu.mutCounts.errors <- lapply(1:(tot.iterations * length(muCounts.check)), (function(iii){
    matrix(0, nrow = total.mutationTypes, ncol = total.samples)
  }))
  nu.mutCounts.reconstructed <- lapply(1:(tot.iterations * length(muCounts.check)), (function(iii){
    matrix(0, nrow = total.mutationTypes, ncol = total.samples)
  }))
  #
  #
  stepAll <- num.processes.toextract * tot.iterations
  stp <- 1
  tmp.elems <- seq(from = 1, to = ncol(nu.Wall), by = stepAll)
  for (startAll in tmp.elems) {
    endAll <- startAll + stepAll - 1
    nu.Wall[,startAll:endAll] <-  muCounts.check[[stp]]$Wall
    nu.Hall[startAll:endAll,] <-  muCounts.check[[stp]]$Hall
    stp <- stp + 1
  }
  #
  stp <- 1
  stepAll <- tot.iterations
  tmp.elems <- seq(from = 1, to = length(nu.mutCounts.errors), by = stepAll)
  for (startAll in tmp.elems) {
    endAll <- startAll + stepAll - 1
    ## add an extra for loop to process list.elements one-by-one
    tmp.mutcnt.err <- muCounts.check[[stp]]$mutCounts.errors
    for (ijk in 1:length(tmp.mutcnt.err)) {
      nu.mutCounts.errors[[(startAll + ijk - 1)]] <-  tmp.mutcnt.err[[ijk]]
    }
    tmp.mutcnt.rec <- muCounts.check[[stp]]$mutCounts.reconstructed
    for (ijk in 1:length(tmp.mutcnt.rec)) {
      nu.mutCounts.reconstructed[[(startAll + ijk - 1)]] <-  tmp.mutcnt.rec[[ijk]]
    }
  }
  #
  tmp.mutcnt.rec <- tmp.mutcnt.err <- muCounts.check <- NULL
  #
  fltr.mutCounts.data <- filterOutIterations(wall = nu.Wall,
                                             hall = nu.Hall,
                                             cnt.errors = nu.mutCounts.errors,
                                             cnt.reconstructed = nu.mutCounts.reconstructed,
                                             params)
  # 
  stability.check <- evaluateStability(wall = fltr.mutCounts.data$Wall,
                                       hall = fltr.mutCounts.data$Hall,
                                       params)
  #
  final.mutCounts.data <- addWeak(mutationTypesToAddSet = bckgrnd.removed.mutset,
                                  processes_I = stability.check$centroids,
                                  processesStd_I = stability.check$centroidStd,
                                  Wall_I = fltr.mutCounts.data$Wall,
                                  genomeErrors_I = fltr.mutCounts.data$mutCounts.errors,
                                  genomesReconstructed_I = fltr.mutCounts.data$mutCounts.reconstructed)
  #
  deconvoluted.results <- list()
  deconvoluted.results$Wall <- final.mutCounts.data$Wall
  deconvoluted.results$Hall <- fltr.mutCounts.data$Hall
  deconvoluted.results$mutCountErrors <- final.mutCounts.data$mutCountErrors
  deconvoluted.results$mutCountReconstructed <- final.mutCounts.data$mutCountReconstructed
  deconvoluted.results$idx <- stability.check$idx
  deconvoluted.results$idxS <- stability.check$idxS
  deconvoluted.results$processes <- final.mutCounts.data$processes
  deconvoluted.results$processesStd <-final.mutCounts.data$processesStd
  deconvoluted.results$exposure <- stability.check$exposure
  deconvoluted.results$exposureStd <-stability.check$exposureStd
  deconvoluted.results$processStab <- stability.check$processStab
  deconvoluted.results$processStabAvg <- stability.check$processStabAvg
  deconvoluted.results$clusterCompactness <- stability.check$clusterCompactness
  return(deconvoluted.results)
}
