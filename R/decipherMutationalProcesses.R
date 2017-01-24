decipherMutationalProcesses <-
function(input, params)
{
  options(warn = -1)
  # First, validate the input
  # Input is a list containing the following elements:
  # 1. cancerType, 2. mutCounts, 3. sampleNames, 4. mutTypes
  requiredElements <- c("cancerType", "mutCounts", "sampleNames", "mutTypes")
  if (! ("input" %in% ls() ) |
      ! (is.list(input)) |
      sum(requiredElements %in% names(input)) != length(requiredElements))
    stop("Malformed Input / Input does not include all the required fields!")
  input <- input
  #
  # Next, check which analysis will be performed: de-novo vs. custom signatures
  # First, a quick extra check!
  #
  if (is.numeric(params$num.processes.toextract)) {
    params$analyticApproach <- "denovo"
    # de novo assignment
  } else {
    # Something is wrong!
    stop("An error occurred!")
  }
  #
  # Validate system parallelization capacity
  if (parallel::detectCores()  < 4)
    stop("Weak computational environment. At least 4 cores are required for running this framework...")
  #
  # Run the analysis (1 line of code.. lol!)
  if (params$analyticApproach == "denovo") {
    deconvData <- deconvoluteMutCounts(input.mutCounts = input$mutCounts, params = params) 
  } else {
    stop("An error occurred!")
  }
  #
  # Capture the output and prepare for return
  mutProcesses <- list()
  mutProcesses$input <- input
  mutProcesses$params <- params
  mutProcesses$allProcesses <- deconvData$Wall
  mutProcesses$allExposures <- deconvData$Hall
  mutProcesses$idx <- deconvData$idx
  mutProcesses$processes <- deconvData$processes
  mutProcesses$exposures <- deconvData$exposure
  mutProcesses$mutCountErrors <- deconvData$mutCountErrors
  mutProcesses$mutCountReconstructed <- deconvData$mutCountReconstructed
  mutProcesses$processStab <- deconvData$processStab
  mutProcesses$processStabAvg <- deconvData$processStabAvg
  #
  # Done! (Hopefully) success! Return results
  options(warn = 0)
  return(mutProcesses)
  #
}
