## ----setup, include=FALSE------------------------------------------------
set.seed(1234567)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, fig.align = "center", results = "asis")

my.t0 <- Sys.time()

library(dplyr)
library(reshape2)
library(kableExtra)
library(ggplot2)
library(gridExtra)
library(mutSignatures)

# load data
data("mutSigData")

# Extra cols
head.muty <- c("T[C>A]G", "T[C>T]A", "G[C>T]A", "C[C>T]T", "T[C>G]T", "T[C>G]C", 
               "T[C>G]A", "T[C>T]T", "T[C>G]T", "T[C>T]G", "T[C>A]G", "T[C>T]A", 
               "A[T>C]T", "A[T>G]G", "C[T>C]T", "T[C>T]A", "T[C>T]A", "G[T>A]C")

## ----eval = FALSE--------------------------------------------------------
#  install.packages("mutSignatures")

## ----eval = FALSE--------------------------------------------------------
#  devtools::install_github("dami82/mutSignatures", force = TRUE,
#                           build_opts = NULL, build_vignettes = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  # Required libs
#  library(dplyr)
#  library(reshape2)
#  library(kableExtra)
#  library(ggplot2)
#  library(gridExtra)
#  library(BSgenome.Hsapiens.UCSC.hg19)
#  
#  # Load mutSignatures
#  library(mutSignatures)
#  
#  # prep hg19
#  hg19 <- BSgenome.Hsapiens.UCSC.hg19
#  
#  # load data
#  data("mutSigData")

## ------------------------------------------------------------------------
# Import data (VCF-like format)
x <- mutSigData$inputC

# Filter non SNV
x <- filterSNV(dataSet = x,  seq_colNames = c("REF", "ALT"))

# Visualize head
head(x) %>% kable() %>% kable_styling(bootstrap_options = "striped")

## ----eval=FALSE----------------------------------------------------------
#  # Attach context
#  x <- attachContext(mutData = x,
#                     chr_colName = "CHROM",
#                     start_colName = "POS",
#                     end_colName = "POS",
#                     nucl_contextN = 3,
#                     BSGenomeDb = hg19)

## ----include=FALSE, eval=TRUE, echo=FALSE--------------------------------
x <- mutSignatures::mutSigData$inputC.ctx

## ------------------------------------------------------------------------
# Visualize head
head(x) %>% kable() %>% kable_styling(bootstrap_options = "striped")
# Remove mismatches
x <- removeMismatchMut(mutData = x,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")    

## ----eval=FALSE----------------------------------------------------------
#  # Compute mutType
#  x <- attachMutType(mutData = x,                      # as above
#                     ref_colName = "REF",              # column name for ref base
#                     var_colName = "ALT",              # column name for mut base
#                     context_colName = "context")

## ----include=FALSE, eval=TRUE, echo=FALSE--------------------------------
x <- x[1:18, ]
x$mutType <- head.muty

## ------------------------------------------------------------------------
# Visualize head
head(x) %>% kable() %>% kable_styling(bootstrap_options = "striped")

## ----eval=FALSE----------------------------------------------------------
#  # Count
#  blca.counts <- countMutTypes(mutTable = x,
#                               mutType_colName = "mutType",
#                               sample_colName = "SAMPLEID")

## ----include=FALSE, eval=TRUE, echo=FALSE--------------------------------
# Count
blca.counts <- mutSignatures::as.mutation.counts(mutSigData$blcaMUTS)

## ----results='markup'----------------------------------------------------
# Mutation Counts
print(blca.counts)

## ------------------------------------------------------------------------
# how many signatures should we extract? 
num.sign <- 4

# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
blca.params <- 
  setMutClusterParams(
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 25,               # bootstrapping: usually 500-1000
    num_parallelCores = 1)                # total num of cores to use (parallelization)

# Extract new signatures - may take a while
blca.analysis <- 
  decipherMutationalProcesses(input = blca.counts,
                              params = blca.params)

## ------------------------------------------------------------------------
# Retrieve signatures (results)
blca.sig <- blca.analysis$Results$signatures

# Retrieve exposures (results)
blca.exp <- blca.analysis$Results$exposures

# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
msigPlot(blca.sig, signature = 1, ylim = c(0, 0.10))
# Plot exposures (ggplot2 object, you can customize as any other ggplot2 object)
msigPlot(blca.exp, main = "BLCA samples") + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))
# Export Signatures as data.frame
xprt <- mutSignatures::as.data.frame(blca.sig) 
head(xprt) %>% kable() %>% kable_styling(bootstrap_options = "striped")
# Get signatures from data (imported as data.frame) 
# and then convert it to mutSignatures object
cosmixSigs <- mutSigData$blcaSIGS %>% 
  dplyr::select(starts_with("COSMIC")) %>% 
  as.mutation.signatures()

blcaKnwnSigs <- mutSigData$blcaSIGS %>% 
  dplyr::select(starts_with("BLCA")) %>% 
  as.mutation.signatures()

# Compare de-novo signatures with selected COSMIC signatures
msig1 <- matchSignatures(mutSign = blca.sig, reference = cosmixSigs, 
                         threshold = 0.45, plot = TRUE) 
msig2 <- matchSignatures(mutSign = blca.sig, reference = blcaKnwnSigs, 
                         threshold = 0.45, plot = TRUE)

## ----fig.height=5, fig.width=12, fig.align='center'----------------------
# Visualize match
# signature 1 is similar to COSMIC ; 
# signatures 2 and 3 are similar to COSMIC
# Here, we should probably extract only 2 mutational signatures
hm1 <- msig1$plot + ggtitle("Match to COSMIC signs.")
hm2 <- msig2$plot + ggtitle("Match to known BLCA signs.")

# Show
grid.arrange(hm1, hm2, ncol = 2)

## ------------------------------------------------------------------------
# Retrieve a mutation.counts data.frame
x <- mutSigData$blcaMUTS

# Visualize header
x[1:10, 1:5] %>% kable() %>% kable_styling(bootstrap_options = "striped")
# Convert it
xx <- as.mutation.counts(x)

## ----results='markup'----------------------------------------------------
# Print to screen
print(xx)

## ------------------------------------------------------------------------
# Obtain 4 COSMIC signatures
cosmx <- mutSigData$blcaSIGS %>% dplyr::select(starts_with("COSMIC"))
cosmx <- as.mutation.signatures(cosmx)

# Obtain 4 BLCA signatures
blcmx <- mutSigData$blcaSIGS %>% dplyr::select(starts_with("BLCA"))
blcmx <- as.mutation.signatures(blcmx)

## ----results='markup'----------------------------------------------------
# Visualize cosmx
print(cosmx)
# Visualize cosmx
print(blcmx)

## ------------------------------------------------------------------------
# Run analysis
blca.expo1 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = cosmx)

blca.expo2 <- resolveMutSignatures(mutCountData = xx, 
                                   signFreqData = blcmx)

## ----fig.width=12, fig.height=5------------------------------------------
# Retrieve exposures (results)
blca.exp.1x <- blca.expo1$results$count.result
blca.exp.2x <- blca.expo2$results$count.result

# Plot exposures
bp1 <- msigPlot(blca.exp.1x, main = "BLCA | COSMIC sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))

bp2 <- msigPlot(blca.exp.2x, main = "BLCA | pre-blca sigs.") + 
  scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#fb9a99", "#1f78b4"))

# Visualize
grid.arrange(bp1, bp2, ncol = 2)

# Compare sigs

# Export Exposures as data.frame
xprt <- as.data.frame(blca.exp.1x, transpose = TRUE)
head(xprt) %>% round() %>% kable() %>% 
  kable_styling(bootstrap_options = "striped")

## ----include=FALSE-------------------------------------------------------
my.t1 <- Sys.time()
tdif <- difftime(time1 = my.t1, time2 = my.t0, units = "mins") %>% 
  as.numeric() %>% round(digits = 3)

## ----results='markup'----------------------------------------------------
sessionInfo()

