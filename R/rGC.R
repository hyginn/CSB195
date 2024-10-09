# tocID <- "R/rGC.R"
#
# Purpose: Define a function that creates a random genetic code with the same
#          redundancy as the standard genetic code, and three stop codons.
#          This script is source()'d from .util.R
#
# Version: 1.2.1
# Date:    2024-10-08
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.2.1  Turn restoration of RNG off by default to prevent cycling
#          the same RNG state in a looped call to the function and
#          RNG dependent generation of the seed.
#   1.2    support setting a seed
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


rGC <- function(seed, noSideEffects = FALSE) {
  # Returns a random Genetic Code with three stop codons, and the same
  # redundancy as the universal code, but randomly replaced amino
  # acids.

  # If noSideEffects is TRUE, the original state of the RNG is saved, and
  # restored on exit. Use this for applications in which the state of the
  # RNG should not be modified by rGC(). However,
  # code that produces the seed from a RNG dependent function like sample()
  # will always produce the same output in a loop if the RNG is not modified.
  if (noSideEffects) {
    oSeed <- .Random.seed
    on.exit(assign(".Random.seed", oSeed, envir = globalenv()))
  }
  if (! missing(seed)) {
    set.seed(seed)
  }

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

  rGC(112358)
    #    AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG
    #    "L" "Y" "M" "F" "T" "E" "A" "N" "A" "C" "R" "A" "*" "D" "N"
    #    ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC
    #    "W" "K" "M" "G" "E" "K" "D" "D" "N" "S" "D" "H" "W" "N" "V"
    #    CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA
    #    "Q" "A" "H" "I" "H" "E" "M" "A" "T" "E" "K" "R" "R" "P" "*"
    #    GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT
    #    "S" "*" "S" "F" "A" "E" "S" "I" "Y" "G" "V" "P" "M" "L" "E"
    #    TTA TTC TTG TTT
    #    "M" "M" "K" "H"
}


# [END]
