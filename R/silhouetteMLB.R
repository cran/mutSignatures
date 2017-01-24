silhouetteMLB <-
function(data, fac, method = "cosine", plot = TRUE){
  #
  # Damiano Fantini f(x)
  #
  if (nrow(data) != length(fac))
    stop ("Bad input!")
  dist.matrix <- as.matrix(proxy::dist(x = data, method = method))
  sil.check <- cluster::silhouette(x = as.numeric(as.factor(fac)), dist = dist.matrix)
  if (plot == TRUE) {
    tmp <- lapply(unique(sil.check[,1]), (function(clid){
      part.out <-  sil.check[sil.check[,1] == clid,]
      part.out[order(part.out[,3], decreasing = TRUE),]
    }))
    tmp <- do.call(rbind, tmp)
    barplot(tmp[nrow(tmp):1,3],
            col = as.factor(tmp[nrow(tmp):1,1]),
            horiz = TRUE,
            xlab = "Silhouette Value",
            ylab = "Samples by Cluster",
            main = "Silhouette Plot",
            border = as.factor(tmp[nrow(tmp):1,1]))
  }
  return(as.vector(sil.check[,3]))
}
