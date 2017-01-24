evaluateStability <-
function(wall,
                              hall,
                              params)
{
  # Defining function constants
  BIG_NUMBER <- 100
  CONVERG_ITER <- 10
  CONVERG_CUTOFF <- 0.005 # cosine distance
  TOTAL_INIT_CONDITIONS <- 5
  #
  num.processes.toextract <- params$num.processes.toextract
  tot.Replicates <- params$tot.Replicates
  process.distance <- params$process.distance
  #
  # Clustering mutational processes using custom clustering procedure
  minClusterDist <- BIG_NUMBER
  totalIter <- ncol(wall) / num.processes.toextract
  idx = matrix(0, nrow = nrow(hall), ncol = 1)
  clusterCompactness <- matrix(0, nrow = num.processes.toextract, ncol = totalIter)
  iStartDataSet = seq(1 , ncol(wall), by = num.processes.toextract)
  iStartingDataSet = iStartDataSet[sample(1:totalIter)]
  #
  for (iInitData in 1 : min(c(TOTAL_INIT_CONDITIONS, totalIter))){
    iStartingData <- iStartingDataSet[iInitData]
    iEnd <- iStartingData + num.processes.toextract - 1
    centroids <- wall[, iStartingData:iEnd]
    #
    centroidsTest <- sapply(1:ncol(centroids), (function(kk){
      runif(nrow(centroids))
    }))
    countIRep <- 0
    #
    for (iRep in 1 : tot.Replicates) {
      tmp.tab <- t(cbind(centroids, wall))
      tmp.pdist <- as.vector(proxy::dist(tmp.tab, process.distance))   
      tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
      allDist <- pracma::squareform(tmp.pdist)
      centroidDist = t(allDist[ 1:ncol(centroids), (ncol(centroids)+1): ncol(allDist)])
      #
      jRange <- sort(1:num.processes.toextract)
      for (jIndex in 1 : num.processes.toextract) {
        j <-  jRange[jIndex]
        for (i in seq(1,ncol(wall), by = num.processes.toextract)){
          iRange = i: (i + num.processes.toextract - 1)
          #
          tmp.min <- min(centroidDist[iRange, j] , na.rm = TRUE)
          Ind <- which(centroidDist[iRange, j] ==  tmp.min)[1]
          centroidDist[iRange[Ind], ] <- BIG_NUMBER
          idx[iRange[Ind], 1] <- j 
        }
      }
      maxDistToNewCentroids <- 0
      #
      for (i in 1 : ncol(centroids)){
        tmp.dset <- wall[,as.vector(idx == i)]
        centroids[, i] <- apply(tmp.dset, 1, mean)
        tmp.dset <- t(cbind(centroids[, i], centroidsTest[, i]))
        tmp.pdist <- as.vector(proxy::dist(tmp.dset, process.distance))   
        tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
        maxDistToNewCentroids <- max(maxDistToNewCentroids, tmp.pdist, na.rm = TRUE)
      }
      if (maxDistToNewCentroids < CONVERG_CUTOFF){
        countIRep <-countIRep + 1
      } else {
        countIRep <- 0
        centroidsTest <- centroids
      }
      if (countIRep == CONVERG_ITER){
        break
      }
    }
    for (i in 1 : ncol(centroids)){
      tmp.tab <- t(cbind(centroids[,i], wall[,as.vector(idx == i)]))
      tmp.pdist <- as.vector(proxy::dist(tmp.tab, process.distance))   
      tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
      clusterDist <- pracma::squareform(tmp.pdist)
      clusterCompactness[i,] = clusterDist[1, 2:ncol(clusterDist)]
    }
    #
    # mean returns column-wise means in matlab
    # all comparisons have to be true to return true
    dist.test <- apply(clusterCompactness, 2, (function(clm) { mean(clm, na.rm=TRUE) }))
    if (sum(minClusterDist > dist.test) == length(dist.test))  {
      centroidsFinal <- centroids
      idxFinal <- idx
      clusterCompactnessFinal <- clusterCompactness
    }
  }
  #
  centroids <- t(centroidsFinal)
  idx <- idxFinal
  clusterCompactness <- clusterCompactnessFinal
  #
  centDist <- apply(clusterCompactness, 1, (function(tmprw){mean(tmprw, na.rm = TRUE)}))
  centDistInd <- order(centDist, decreasing = FALSE)
  clusterCompactness <- clusterCompactness[centDistInd,]
  centroids <-centroids[centDistInd, ]
  idxNew <- idx
  #
  for (i in 1 : num.processes.toextract){
    idxNew[as.vector(idx == centDistInd[i]),1] <- i
  }
  idx <- idxNew
  #
  if (num.processes.toextract > 1) {
    processStab <- silhouetteMLB(data = t(wall), fac = idx, process.distance)
    processStabAvg <- matrix(0, nrow = 1, ncol = num.processes.toextract)
    for (i in 1 : num.processes.toextract) {
      processStabAvg[1,i] = mean(processStab[idx==i])
    }
  } else {
    #
    tmp.tab <- t(cbind(t(centroids), wall))
    tmp.pdist <- as.vector(proxy::dist(tmp.tab, process.distance))
    tmp.pdist[tmp.pdist == 1 | tmp.pdist == 0] <- NA
    allDist <- pracma::squareform(tmp.pdist)
    #processStab = 1 - allDist( 1:size(centroids', 2), (size(centroids', 2)+1): size(allDist, 2) )';
    processStab <- 1 - t(allDist[1:nrow(centroids), (nrow(centroids) + 1): ncol(allDist)])
    processStabAvg <- apply(processStab, 2, (function(clmn) {mean(clmn, na.rm = TRUE)}))
  }
  #
  centroidStd <- matrix(0, nrow = nrow(centroids), ncol = ncol(centroids))
  for (i in 1 : num.processes.toextract) {
    centroidStd[i,] <- apply(wall[, idx == i], 1, (function(rw){sd(rw, na.rm = TRUE)}))
  }
  #
  centroids <- t(centroids)
  centroidStd <- t(centroidStd)
  idxS <- matrix(0, nrow = length(idx), ncol = 1)
  #
  for (i in seq(1, ncol(wall), by = num.processes.toextract)){
    iEnd <- i + num.processes.toextract - 1
    idxG <- idx[i : iEnd]
    #
    for (j in 1 : num.processes.toextract){
      idxS[(i + j - 1), ] = which(idxG == j)
    }
  }
  #
  exposure <- matrix(0, nrow = max(idxS), ncol(hall))
  exposureStd <- matrix(0, nrow = max(idxS), ncol(hall))
  #
  for (i in 1 : max(idxS)) {
    exposure[i, ] <- apply(hall[idx==i,], 2, (function(cl){mean(cl, na.rm = TRUE)}))
    exposureStd[i, ] <- apply(hall[idx==i,], 2, (function(cl){sd(cl, na.rm = TRUE)}))
  }
  result.list <- list()
  result.list$centroids <- centroids
  result.list$centroidStd <- centroidStd
  result.list$exposure <- exposure
  result.list$exposureStd <- exposureStd
  result.list$idx <- idx
  result.list$idxS <-idxS
  result.list$processStab <- processStab
  result.list$processStabAvg <- processStabAvg
  result.list$clusterCompactness <- clusterCompactness
  return(result.list)
}
