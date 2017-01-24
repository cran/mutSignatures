do.nmf <-
function(v, r, params, verbose = TRUE) {
  rownames(v) <- NULL
  v <- as.matrix(v)
  #test for neg values in v
  stopconv <- params$stopconv
  niter <- params$niter
  eps <- params$eps
  #
  if (min(v) < 0 )
    stop("Matrix entries caan not be negative")
  if (min(apply(v, 1, sum) ) == 0 )
    stop("Entries can not all be equal to 0")
  #
  tmp.n <- nrow(v)
  tmp.m <- ncol(v)
  #
  cons <- matrix(0, nrow = tmp.m, ncol = tmp.m)
  consold <- cons
  tmp.inc <- 0
  tmp.j <- 0
  #
  #initialize random w and h
  tmp.w <- sapply(1:r, (function(z){
    runif(tmp.n)
  }))
  tmp.h <- sapply(1:tmp.m, (function(z){
    runif(r)
  }))
  #
  for (i in 1:niter) {
    #divergence-reducing NMF iterations
    x1 <- sapply(1:tmp.m, (function(i2){
      apply(tmp.w, 2, sum)
    }))
    #
    tmp.h <- tmp.h * ( t(tmp.w) %*% (v/(tmp.w %*% tmp.h)) ) / x1
    x2 <- do.call(rbind, lapply(1:tmp.n, (function(z){
      apply(tmp.h, 1, sum)
    })))
    tmp.w <- tmp.w * ((v / (tmp.w %*% tmp.h)) %*% t(tmp.h)) / x2
    #
    if (i %% 5000 == 0 & verbose == TRUE){
      message(".", appendLF = FALSE)
    }
    # test convergence every 10 iterations
    #
    if (i %% 10 == 0){
      tmp.j <- tmp.j + 1
      #
      # adjust small values to avoid overflow
      tmp.h[tmp.h<eps] <- eps
      tmp.w[tmp.w<eps] <- eps
      #
      # build connectivity matrix
      y <- apply(tmp.h, 2, max)
      index <- apply(tmp.h, 2, (function(dt){
        which.max(dt)[1]
      }))
      #
      mat1 = t(sapply(1:tmp.m, (function(ii){
        index
      })))
      mat2 = sapply(1:tmp.m, (function(ii){
        index
      }))
      cons <- mat1 == mat2
      if (sum(cons != consold) == 0 ) {
        tmp.inc <- tmp.inc + 1
      } else {
        tmp.inc <- 0
      }
      if (tmp.inc > stopconv){
        break()
      }
      consold <- cons
    }
  }
  my.list <- list()
  my.list$w <- tmp.w
  my.list$h <- tmp.h
  return(my.list)
}
