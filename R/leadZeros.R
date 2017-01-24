leadZeros <-
function(n, m) {
  max.zeros <- nchar(as.character(round(m)))
  tmp.digits <- nchar(as.character(round(n)))
  zeros.toAdd <- max.zeros - tmp.digits
  returnVect <- c(rep("0", zeros.toAdd), as.character(round(n)))
  return(paste(returnVect, sep = "", collapse = ""))
}
