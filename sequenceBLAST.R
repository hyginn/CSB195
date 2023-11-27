# tocID <- "./sequenceBLAST.R"
#
#
# Purpose: Explore homologues of a sequence (PPP2R1A).
#
#
# Version: 1.0
# Date:    2023-11
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.0  First lecture version.
#
# ToDo:
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                  Line
#TOC> ------------------------------------------------------
#TOC>   1        Preparation: packages                    55
#TOC>   1.1        Preparation: Helper functions          78
#TOC>   2        Introduction: BLAST                      83
#TOC>   3        Coverage of the Tree of Life            123
#TOC>   4        Conservation patterns                   135
#TOC> 
#TOC> ==========================================================================


################################################################################
#                                                                              #
#                            S O U R C E   S A F E                             #
#                                                                              #
#     This script is  "source safe".  source()'ing the  code will  define      #
#     (refresh) the global parameters, and define the functions. It won't      #
#     actually run code. Execute the statements in the if (FALSE) { ... }      #
#     blocks to run the interactive parts of this script.                      #
#                                                                              #
################################################################################



# NOTE: You need to read and explore this script, execute the functions but also
# try what happens if you change things.

# OPEN A NEW R SCRIPT into which you write your code experiments. Then you can
# easily reproduce what you did, and the experiments are automatically
# documented for your report.



# =    1  Preparation: packages  ===============================================

# Many things we do in computational biology require tools of the BioConductor
# project. These are loaded from their own installation framework: the
# "Bioconductor manager" ...
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Package information:
#  library(help = Biostrings)       # basic information
#  browseVignettes("Biostrings")    # available vignettes
#  data(package = "Biostrings")     # available datasets


# Once the BiocManager is installed, we can use it to install packages like
# we would install from CRAN. Biostrings is one of the foundational packages
# that has functions to handle the basic objects of biological sequences.
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}



# ==   1.1  Preparation: Helper functions  =====================================

# ...


# =    2  Introduction: BLAST  =================================================

# In this script we explore the evolutionary descent  of the PP2A
# (Protein Phosphatase 2) subunit A.

# Ask ChatGPT to explain the principles evolutionary relationships:

if (FALSE) {

  t2c("I am studying PPP2R1A. Please explain the concept of \"homology\", and in particular the difference betwwen \"orthologues\" and \"paralogues\". If sequences are othologous (or paralogous), what difference does it make for their sequence conservation and conservation of function? And what are \"isoforms\"?")


  t2c("Following up on that: what do I need to know about \"neo-functionalization\" and \"sub-functionalization\"? Also, do you know whether either or both apply to human PPP2R1A?")


  t2c("What is the best tool that I can use to find protein sequences that are related to PPP2R1A to explore these issues? How should I use it? Can you walk me through an example?")


  t2c("Once I find sequences related to PPP2R1A, how would I estimate what portion of the \"Tree of Life\" shares the PPP2R1A protein? What would be the most likely sequence(s) present in the the \"Last Common Ancestor\"?")

  # Now perform an actual BLAST search. Make sure you are running the
  # BLASTp program (protein sequence) against the "landmark" database. Use
  # default parameters otherwise. You don't need to enter the actual sequence
  # into the search field, the RefSeq ID (NP_055040) suffices. Once you
  # have a result, define the RID in your script. Then you can ask ChatGPT-4
  # (or Bing) to help you interpret the complex but very comprehensive
  # and informative results.

  myRID <- "P8911YEM013"

  t2c(sprintf("I have run a BLAST search but find the results very complex. Here is the link to the result page: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=%s Can you be so kind and have a look and walk me through the highlights?",
              myRID))


  # Conclusion?

  }



# =    3  Coverage of the Tree of Life  ========================================


if (FALSE) {
  # ...


  # Conclusion?
}



# =    4  Conservation patterns  ===============================================


if (FALSE) {
  # ...

  # Conclusion?
}





# [END]
