# neighCodons.R
#
# Purpose: Input: a DNA codon.
#          Output: its nine neighbouring codons in a vector
#
# Version: 1.0
# Date:    2023-10-11
# Author:  boris.steipe@utoronto.ca
#
#  THIS SCRIPT CAN BE SOURCE()'D WITHOUT SIDE EFFECT TO LOAD THE FUNCTION.
#
# ToDo:
#    - Write test
#
# Notes:
#    - We fail on any input not matching "^[ACGT]{3}$" - we do NOT make
#      assumptions to try and fix input.
#
# ==============================================================================

neighCodons <- function(inCodon) {

  if (! grepl("^[ACGT]{3}$", inCodon)) {
    stop(sprintf("%s%s%s\n",
                 "Input codon \"",
                 inCodon,
                 "\" does not contain exactly three nucleotides from {ACGT}."))
  }

  NUC <- c("A", "C", "G", "T")

  outCodons <- character()
  for (i in 1:3) {
    for (x in NUC[! NUC == substring(inCodon, i, i)]) {
      s <- inCodon
      substring(s, i, i) <- x
      outCodons <- c(outCodons, s)
    }
  }
  return(outCodons)
}



# ====  TESTS  =================================================================
#

if (FALSE) {

  all(neighCodons("TCG") == c("ACG", "CCG", "GCG",  # pos 1
                              "TAG", "TGG", "TTG",  # pos 2
                              "TCA", "TCC", "TCT")) # pos 3

  all(neighCodons("GGG") == c("AGG", "CGG", "TGG",
                              "GAG", "GCG", "GTG",
                              "GGA", "GGC", "GGT"))
  neighCodons()
  neighCodons("xx")
  neighCodons("yyy")
  neighCodons("zzzz")

}

# [END]
