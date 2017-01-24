setMutCountObject <-
function(mutCountMatrix, mutationTypes = NULL, sampleNames = NULL, datasetName = NULL) {
  #
  # Initialize collector
  input <- list()
  #
  # Data set name
  if (is.character(datasetName)) {
    input$cancerType <- datasetName[1]
  } else {
    input$cancerType <- "Custom_Dataset"
  }
  #
  # mutCountMatrix
  if (is.matrix(mutCountMatrix) | is.data.frame(mutCountMatrix) ) { 
    input$mutCounts <- mutCountMatrix
  } else {
    stop("Bad input")
  }
  #
  # Colnames // sample names
  if (is.null(sampleNames) & !is.null(colnames(mutCountMatrix))) {
    input$sampleNames <- colnames(mutCountMatrix)  
  } else if (is.character(sampleNames) & length(sampleNames) == ncol(mutCountMatrix)) {
    input$sampleNames <- sampleNames
  } else {
    input$sampleNames <- paste("sample", 10001:(10000 + ncol(mutCountMatrix)), sep = ".")
  }
  #
  # MutTypes // rownames
  if (is.null(mutationTypes) & !is.null(rownames(mutCountMatrix))) {
    input$mutTypes <- rownames(mutCountMatrix)  
  } else if (!is.null(mutationTypes) & is.character(mutationTypes) & length(mutationTypes) == nrow(mutCountMatrix)) {
    input$mutTypes <- mutationTypes
  } else {
    stop("Please, make sure to specify the mutation types")
  }
  #
  # Final warning if non-standard mutTypes
  if (sum(regexpr("^(A|C|T|G)\\[(A|C|T|G)>(A|C|T|G)\\](A|C|T|G)$",  input$mutTypes) < 0) > 0) {
    message("Alert!!! Please, note that one or more mutationTypes are in a non-standard format...")
  }
  return(input)
}
