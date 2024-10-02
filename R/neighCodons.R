# tocID <- "R/neighCodons.R"
#
# Purpose: Define a function that returns all 9 neighbouring codons for an
#          input codons.- i.e. codons that can be derived from the input
#          with a single point mutation.
#          This script is source()'d from .util.R
#
# Version: 1.1
# Date:    2024-10-01
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.1    2024 - standalone and sourced from .util.R
#   1.0    Written in 2023 as part of sampleSolutionGeneticCode.R
#
# Input:   a codon
#
# Output:  a vector of its nine neighbouring codons
#
# ToDo:
#
# Notes:
#
# ==============================================================================


neighCodons <- function(inCodon) {

  # validate input
  if (! grepl("^[ACGT]{3}$", inCodon)) {
    stop(sprintf("%s%s%s\n",
                 "Input codon \"",
                 inCodon,
                 "\" does not contain exactly three nucleotides from {ACGT}."))
  }

  NUC <- c("A", "C", "G", "T")

  outCodons <- character(9)
  for (i in 1:3) {
    x <- NUC[! NUC == substring(inCodon, i, i)]
    for (j in 1:3) {
      s <- inCodon
      substring(s, i, i) <- x[j]
      outCodons[((i-1) * 3) + j] <- s
    }
  }
  return(outCodons)
}



# === USAGE, TESTS, AND VALIDATION =============================================

if (FALSE) {

  inC <- "AAA"
  (nC <- neighCodons(inC))
  length(nC) == 9                              # TRUE if nC has 9 codons.
  length(unique(nC)) == 9                              # TRUE if nC has 9 codons.


  sum(xGC == "*") == 3                         # TRUE if xGC has 3 "*".
  all(sort(names(xGC)) == sort(GCdf$Codon))    # TRUE if all codons are used.
  all(sort(table(xGC)) == sort(table(GCdf$A))) # TRUE if symbol counts match.


}


# [END]
