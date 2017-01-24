plotSignatureExposures <-
function(mutSign.analysis.result)
{
  mut.by.sample <- apply(mutSign.analysis.result$exposures, 2, sum)
  reord.id <- order(mut.by.sample, decreasing = TRUE)
  re.exposures <- mutSign.analysis.result$exposures[,reord.id]
  #
  sig.id <- order(apply(re.exposures, 1, sum), decreasing = TRUE)
  #
  # Workaround to avoid problems at the ccheck step
  plt.sample <- plt.signat <- plt.mutcnt <- NULL
  #
  procExpo.bySample <- data.frame(do.call(rbind, lapply(1:ncol(re.exposures), (function(i){
    t(sapply(sig.id, (function(j){
      c(plt.sample = i, # totMut-ordered sample ID
        plt.signat = j, # mutation signature (process) #
        plt.mutcnt = re.exposures[j,i] # mut count explained by signature j
      )
    })))
  }))))
  procExpo.bySample$plt.signat <- factor(procExpo.bySample$plt.signat, levels = rev(sig.id))
  bp <- ggplot(data=procExpo.bySample, aes(x=plt.sample, y=plt.mutcnt, fill=plt.signat)) +
    geom_bar(stat="identity")
  bp <- bp + theme_minimal() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  return(bp)
}
