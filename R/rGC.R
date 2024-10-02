# tocID <- "R/rGC.R"
#
# Purpose: Define a function that creates a random genetic code with the same
#          redundancy as the standard genetic code, and three stop codons.
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
# Precondition: GCdf in the global name space
#
# Output:  a random genetic code
#
# ToDo:
#
# Notes:
#
# ==============================================================================


rGC <- function() {
  # Returns a random Genetic Code with three stop codons, and the same
  # redundancy as the universal code, but randomly replaced amino
  # acids.

  A <- GCdf$A                       # Fetch the single-letter symbols.
  A <- A[A != "*"]                  # Remove the "*" symbols.
  tGC <- table(A)                   # table() returns named counts.
  names(tGC) <- sample(names(tGC))  # Shuffle the names.
  GC <- rep(names(tGC), tGC)        # Expand tGC into a vector with rep().
  GC <- c(GC, "*", "*", "*")        # Add the three "*" symbols.
  GC <- sample(GC)                  # Shuffle the vector.
  names(GC) <- GCdf$Codon           # Give it codon names.
  return(GC)                        # Voilá - a random code
}



# === USAGE, TESTS, AND VALIDATION =============================================

if (FALSE) {

  xGC <- rGC()

  printGCtable(myGC = xGC)                     # Inspect visually
  length(xGC) == 64                            # TRUE if xGC has 64 symbols.
  sum(xGC == "*") == 3                         # TRUE if xGC has 3 "*".
  all(sort(names(xGC)) == sort(GCdf$Codon))    # TRUE if all codons are used.
  all(sort(table(xGC)) == sort(table(GCdf$A))) # TRUE if symbol counts match.


}


# [END]
