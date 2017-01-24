removeWeak <-
function(input.mutCounts, params){
  #
  remove.weak.muttypes <- params$remove.weak.muttypes
  #
  sum.counts <- apply(input.mutCounts, 1, sum)
  sum.counts.idx <- order(sum.counts, decreasing = FALSE)
  sorted.sum.counts <- sum.counts[sum.counts.idx]
  #
  tot.mut.counts <- sum(input.mutCounts)
  #
  tot.muttypes.toremove <- sum((sapply(1:length(sorted.sum.counts), (function(i){
    sum(sorted.sum.counts[1:i])
  })) / tot.mut.counts) < remove.weak.muttypes)
  #
  return.list <- list()
  #
  if (tot.muttypes.toremove > 0) {
    removed.mutset <- sum.counts.idx[c(1:tot.muttypes.toremove)]
    input.mutCounts <- input.mutCounts[-removed.mutset,]
    #
    return.list$removed.mutset <- removed.mutset
  } else {
    return.list$removed.mutset <- (-1)
  }
  #
  return.list$output.mutCounts <- input.mutCounts
  return(return.list)
}
