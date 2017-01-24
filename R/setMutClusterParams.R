setMutClusterParams <-
function(num.processes.toextract = 2, # number of de-novo signatures to extract
                                tot.iterations = 4,
                                tot.cores = 5,
                                remove.weak.muttypes = 0.01,
                                remove.last.percent = 0.07,
                                process.distance = "cosine",
                                tot.Replicates = 100,
                                eps = 2.2204e-16,
                                stopconv = 10000,
                                niter = 1000000
)
{
  #
  # Step-by-step Parameter Validation and preparation
  paramList <- list()
  #
  if (!((is.numeric(num.processes.toextract[1]) & num.processes.toextract[1] > 0) ))
    stop("Provide a reasonable number of signatures/processes to extract")
  #
  paramList$num.processes.toextract <- round(num.processes.toextract[1]) 
  #
  if (!(is.numeric(tot.iterations[1]) & tot.iterations[1] > 0))
    stop("Provide a reasonable number of iterations to run (Bootstrapping)")
  paramList$tot.iterations <- round(tot.iterations[1])
  if (paramList$tot.iterations  < 2 ) {
    paramList$tot.iterations <- 2
    message("At least 2 iterations have to be run on each core. 'tot.iterations' defaulted to 2")
  }
  #
  if (!(is.numeric(tot.cores[1]) & tot.cores[1] > 0))
    stop("Provide a reasonable number of CPU cores to use for the analysis")
  paramList$tot.cores <- round(tot.cores[1])
  #
  if (!(is.numeric(remove.weak.muttypes[1]) & remove.weak.muttypes[1] >= 0 & remove.weak.muttypes[1] < 1))
    stop("Provide a reasonable (0.00-0.99) number of low-occurring mutation types to remove from the input before starting the analysis")
  paramList$remove.weak.muttypes <- remove.weak.muttypes[1]
  #
  if (!(is.numeric(remove.last.percent[1]) & remove.last.percent[1] >= 0 & remove.last.percent[1] < 1))
    stop("Provide a reasonable (0.00-0.99) number for filtering out poor iterations")
  paramList$remove.last.percent <- remove.last.percent[1]
  #
  allowed.dist.methods <- c("Braun-Blanquet", "Chi-squared", "correlation", "cosine", "Cramer", "Dice",
                            "eDice", "eJaccard", "Fager", "Faith", "Gower", "Hamman", "Jaccard",
                            "Kulczynski1", "Kulczynski2", "Michael", "Mountford", "Mozley", "Ochiai",
                            "Pearson", "Phi", "Phi-squared", "Russel", "simple matching", "Simpson",
                            "Stiles", "Tanimoto", "Tschuprow", "Yule", "Yule2", "Bhjattacharyya",
                            "Bray", "Canberra", "Chord", "divergence", "Euclidean", "fJaccard",
                            "Geodesic", "Hellinger", "Kullback", "Levenshtein", "Mahalanobis",
                            "Manhattan", "Minkowski", "Podani", "Soergel","supremum", "Wave", "Whittaker")
  if (!(is.character(process.distance[1]) & process.distance[1] %in% allowed.dist.methods))
    stop("Unknown method for calculating distances. For options, run: <<summary(proxy::pr_DB)>>")
  paramList$process.distance <- process.distance[1]
  #
  if (!(is.numeric(tot.Replicates[1]) & tot.Replicates[1] > 99))
    stop("Provide a reasonable number of replicates for stability evaluation of the results")
  paramList$tot.Replicates <- round(tot.Replicates[1])
  #
  if (!(is.numeric(eps[1]) & eps[1] > 0 & eps[1] < 0.0001))
    stop("Provide a reasonably small number (0 < n < 0.0001) for data overflow prevention")
  paramList$eps <- eps[1]
  #
  if (!(is.numeric(stopconv[1]) & stopconv[1] > 1000))
    stop("Provide a reasonable large number: number of 'conn-matrix-stable' iterations before stopping NMF")
  paramList$stopconv <- round(stopconv[1])
  #
  if (!(is.numeric(niter[1]) & niter[1] > 200000))
    stop("Provide a reasonable large number: total NMF iterations")
  paramList$niter <- round(niter[1])
  #
  return(paramList)
}
