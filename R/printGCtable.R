# tocID <- "R/printGCtable.R"
#
# Purpose: Define the function printGCtable() to print any genetic code in
#          either the standard or any alternative layout.
#
# Version: 2.0
# Date:    2023-09 - 2024-10
# Author:  boris.steipe@utoronto.ca; CSB195 2023 with ChatGPT-4; ChatGPT-4o
#
#  THIS SCRIPT CAN BE source()'D WITHOUT SIDE EFFECT TO LOAD THE FUNCTION.
#
# Versions:
#      2.0  2024 code refactoring. Rewrite of processing logic to work from
#           a 3D data cube instead of the ad hoc tables that were used
#           previously. Support amino acid names as output format.
#      1.2  Make standard code layout the default. Accept a code vector
#           as input.
#      1.1  Change to work with both one-letter and three-letter codes
#      1.0  Basic running code
#
# ToDo:
#    - Comment the code
#    - Refactor for coding style (variable names!)
#    - Refactor to use alternative column names
#    - Return the reordered data in a format that is suitable for further
#      processing / printing / analysis
#    - Write tests
#
# Notes:
#   The v 1.0 code was developed in class, in a ChatGPT dialogue. The dialogue
#   is here:
#     https://chat.openai.com/share/53e3b709-7067-4842-b61e-5a7010181f28
#
#   The function does not change any external objects (no side effects).
#   That is good practice.
#
#   The function validates user input. We did not ask ChatGPT to do that.
#   Doing so is good practice.
#
#   ChatGPT likes to use pothole_style for names. We prefer camelCase. But
#   when we implicitly suggested printGCtable() for our function in
#   the example function signature we gave it, it used that name
#   going forward. That is good dialogue.
#
#   ChatGPT added an empty line after each fourth row of the table.
#   That structures the output and makes it more readable. We did not ask
#   for that, it came up with that on its own. Doing so is a good idea.
#
#   The variable names are not very good. For example: don't use "c",
#   it is easy to misread for the inbuilt function "c()".
#
# ==============================================================================


printGCtable <- function(myGC,
                         nuc = c("T", "C", "A", "G"),
                         order = c(3, 2, 1),
                         dat = GCdf,                  # loaded from .util.R
                         format = "Aaa") {
  # Parameters:
  #   myGC:      a genetic code vector with 64 symbols name()'d with the
  #                64 codons. If missing, the standard code (NCBI code "1")
  #                is used and taken from GCdf.
  #   nuc:       a permutation of {A, C, G, T} in the order they should appear
  #                in the table. Default: T C A G.
  #   order:     a permutation of {1, 2, 3} corresponding to which nucleotide
  #                should vary in the table's blocks, columns, and rows,
  #                respectively. Default: 3 2 1.
  #   dat:      a data frame with the genetic code: required columns are
  #                "Codon" and the column that is named for the
  #                encoded amino acid symbol - either "Aaa" (default) or "A".
  #                Default is GCdf provided by .util.R
  #   format:    passed to A2aa - either "A", "Aaa", or "Name".
  #
  #   Valid parameter combinations reflect three major use cases:
  #    (a) to print a genetic code vector in the standard table, use:
  #        printGCtable(myGC).
  #    (b) to print the universal code in a different table, use something
  #        like  printGCtable(order = 3:1).

  # Check input consistency
  if (! all(sort(nuc) == c("A", "C", "G", "T")) ) {
    stop("\"nuc\" is not a permutation of A C G T.")
  }
  if (! all(sort(order) == 1:3)) {
    stop("\"order\" is not a permutation of 1 2 3.")
  }
  if (! any(format == colnames(dat))) {
    stop(sprintf("format \"%s\" does not appear in the dataframe \"dat\"",
                 format))
  }

  # Load the standard genetic code if myGC is missing
  if (missing(myGC)) {
    myGC <- dat[ , format]
    names(myGC) <- dat$Codon
  } else {
    # Check myGC's consistency
    if (length(myGC) != 64 ) {
      stop("\"myGC\" does not contain 64 elements.")
    }
    if (! all(sort(names(myGC)) == sort(dat$Codon)) ) {
      stop("\"myGC\" is not named with the 64 codons.")
    }
    for (i in 1:length(myGC)) {
      dat[names(myGC)[i], format] <- myGC[i]
    }
    # Overwrite column "format" with the symbols in the
    # input vector
  }

  # create a 3D matrix with the codons in the desired order.

  myCodons <- character(64)
  dim(myCodons) <- c(4,4,4)

  codon <- character(3)
  for (i in 1:4) {
    iPos <- order[1]
    codon[iPos] <- nuc[i]
    for (j in 1:4) {
      jPos <- order[2]
      codon[jPos] <- nuc[j]
      for (k in 1:4) {
        kPos <- order[3]
        codon[kPos] <- nuc[k]
        myCodons[i, j, k] <- paste0(codon, collapse = "")
      }
    }
  }

  for (iBlock in 1:4) {
    for (iRow in 1:4) {
      for (iCell in 1:4) {
        thisCodon <- myCodons[iRow, iCell, iBlock]
        thisAA <- myGC[thisCodon]
        thisAA <- A2Aaa(thisAA, out = format)
        if (format == "Name") {
          thisAA <- sprintf("%-13s", thisAA)
        }
        cat(sprintf("%s %s  ", thisCodon, thisAA))
      }
      cat("\n")
    }
    cat("\n\n")
  }
}


# ====  TESTS AND USAGE ========================================================
#

if (FALSE) {

  # This should produce the standard way to print the code:
  # ... all defaults
  printGCtable()

  # ... three-letter code (should be the same)
  printGCtable(nuc = c("T", "C", "A", "G"), order = c(3, 2, 1), dat = GCdf)

  # ... one-letter code
  printGCtable(nuc = c("T", "C", "A", "G"),
               order = c(3, 2, 1),
               dat = GCdf,
               format = "A")

  # ... amino acid names
  printGCtable(format = "Name")


  # ... Three letter code, reversed, using defaults:
  printGCtable(nuc = c("G", "A", "C", "T"),
               order = c(1, 2, 3))

  # print from Biostrings supplied NCBI code "1"
  printGCtable(Biostrings::getGeneticCode("1"))

  # ... a shuffled code with all-upper case three-letter symbols, in
  #     standard order
  myGC <- toupper(sample(GCdf$Aaa))
  names(myGC) <- GCdf$Codon
  printGCtable(myGC = myGC)


}


# [END]
