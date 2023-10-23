# tocID <- "aaSim.R"
#
# Purpose: Define a function that computes a pairwise similarity score for
#          amino acids.
#
# Version: 1.0
# Date:    2023-10
# Author:  boris.steipe@utoronto.ca; CSB195 2023 Class; ChatGPT-4
#
# Versions:
#   1.0    Define aaSim() efficiently as a closure.
#   0.5    Split this into two portions:
#            - aaSim.R  defines the actual functions for general use.
#            - aminoAcidSimilarity.R  contains the development process.
#   0.4    Minor change: Replace all usages of one-letter / three-letter
#          "code" with "symbol" to conform to IUPAC conventions.
#   0.3    Scale dimensions of the aaFeatureSpace according to the variance
#          of the principal components. Previously they all had the same
#          mean of 0 and sd() of 0.229 ... but that gives inappropriately high
#          weights to the higher-order PC's. The old feature space is in
#          "data/aaFeatureSpace.1.0.RData". The new (improved) feature
#          space is in "data/aaFeatureSpace.2.0.RData".
#   0.2    Implemented much of the pseudocode in co-development
#          with ChatGPT 4. This happened in two parts.
#            Part 1 is here - It covers the initial prompt, and ends up
#            with a cleaned and validated dataset of amino acid
#            property indices.
#            https://chat.openai.com/share/0cbd454e-2fa3-40e8-89d5-0e3f88a427b7
#
#            Part 2 is here - It analyzes the data, transforms it with
#            PCA, defines a distance function, and then leads into
#            Exploratory Data Analysis. In the end, we discovered some non-
#            obvious things about biology.
#            https://chat.openai.com/share/330a27ba-c40c-4d40-b8c3-ad96439c0764
#
#   0.1.1  Pseudocode format update. Use structuring keywords. Separate out
#          function definition from function use.
#   0.1    First workflow pseudocode co-developed with ChatGPT-4
#          https://chat.openai.com/share/44ebd17d-6e74-4bb0-8ed2-5d15cff78828
#
# Input:   a vector with two amino acid one-letter symbols
# Output:  a scalar distance in a feature space between the amino acids
#
# ToDo:   Code cleanup: make sure the file is safe-to source and defines
#         the aaSim() function upon source()'ing.
#
# Notes:  This script will contain BOTH the code for the function development
#         and the actual function definition. Thus all relevant information
#         is sanely kept in one single file.
#
# ==============================================================================


aaSimConstructor <- function() {

  # Value: a function that computes pairwise amino acid distances.

  # Note:
  # This function returns another function as its output. The returned function
  # is a so-called "closure" - a combination of data in its environment, and
  # code instructions. By doing this, the auxiliary data for our function
  # is created only once, when the function is defined, not every time it is
  # called. Yet, since the objects are defined locally, we do not risk
  # overwriting objects in our workspace when we define the function.
  #
  # Note:
  #    The resulting function takes as its input two amino acid
  #    one-letter symbols and returns the Euclidian distance between
  #    the two vectors in the feature space defined in SPACEFILE.
  #    a1, a1: two letters from "ACDEFGHIKLMNPQRSTVWY*"
  #    Value:  distance

  # === Parameters ==============================

  # Feature space:
  # This file is based on selected AAINDEX indices, curated by CSB195-2023
  # to remove indices that were considered to be _derived_ from the genetic
  # code itself. This makes the resulting feature space useable to evaluate
  # the robustness of the genetic code against point mutations. The indices
  # were reduced to 14 Principal Components after scaling them, and the
  # principal components were scaled by multiplying them with their
  # contribution to the total variance.
  SPACEFILE <- "data/aaFeatureSpace.2.0.Rds"

  # Stop codon distance:
  # The distance of an amino acid to a stop codon is STOPDIST times the
  # maximum distance in the distance matrix.
  STOPDIST <- 1.5

  # Recreate a 14-dimensional feature space of amino acids
  # (see aminoAcidSimilarity.R for details). Normalize rownames to single
  # letter symbols.
  AASPACE <- readRDS(SPACEFILE)
  rownames(AASPACE) <- A2Aaa(rownames(AASPACE), out = "A")

  # Compute a 21 x 21 matrix of Euclidian distances between any pair of vectors
  # in AASPACE.
  AADMAT <- matrix(numeric(21 * 21), nrow = 21)
  # Fill the first 20 x 20 values with amino acid pair distances
  for (i in 1:20) {
    for (j in 1:20) {
      AADMAT[i, j] <- sqrt(sum((AASPACE[i, ] - AASPACE[j, ])^2))
    }
  }

  # Define distance to stop codons: distance of stop codon "*" to
  # any other codon as STOPDIST times the maximum distance in the
  # distance matrix.

  stopDist <- STOPDIST * max(AADMAT)
  AADMAT[ 21, ] <- stopDist
  AADMAT[ , 21] <- stopDist
  AADMAT[21,21] <- 0

  rownames(AADMAT) <- c(rownames(AASPACE), "*")
  colnames(AADMAT) <- c(rownames(AASPACE), "*")

  # Define a function to return the pairwise distance between two points
  # in the distance space
  f <- function(a1, a2) {
    return(AADMAT[a1, a2])
  }

  return(f)  # return the function
}

aaSim <- aaSimConstructor()   # Define aaSim()
rm(aaSimConstructor)        # No need to keep the function in the workspace




# === USAGE, TESTS, AND VALIDATION =============================================

if (FALSE) {

  aaSim("Q", "Q")   # same
  aaSim("Q", "F")   # different
  aaSim("F", "Q")   # same as above
  aaSim("F", "I")   # quite similar
  aaSim("Q", "*")   # quite different
  aaSim("*", "G")   # same as above

  # Sanity check: Find the smallest and largets pairwise distance
  # recreate the distance matrix
  aa <- unlist(strsplit("ACDEFGHIKLMNPQRSTVWY", ""))
  x <- matrix(numeric(400), nrow = 20)
  for (i in 1:20) {
    for (j in 1:20) {
      x[i, j] <- aaSim(aa[i], aa[j])
    }
  }
  # Print minima and maxima (most and least similar)
  dMin <- min(x[x > 0])  # ignore 0, these are just the identical
  dMax <- max(x[x > 0])  # amino acids
  for (i in 1:19) {
    for (j in i:20) {
      if (x[i, j] == dMin) {
        cat(sprintf("min: (%s, %s) = %.3f\n", aa[i], aa[j], x[i, j]))
      } else if (x[i, j] == dMax) {
        cat(sprintf("max: (%s, %s) = %.3f\n", aa[i], aa[j], x[i, j]))
      }
    }
  }
  # This may vary with the distance matrix details, but will probably turn out
  # to make (I, V) the minimum distance pair, (I, P) the maximum distance pair.


}  # end if (FALSE) ...


# [END]
