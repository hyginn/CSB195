# tocID <- "courseScripts/07-geneticCodeOptimization.R"
#
# Purpose: Optimizing the Genetic Code
#
# Version: 1.0
# Date:    2024-10-09
# Author:  boris.steipe@utoronto.ca; ChatGPT-4 and 4o
#
# Versions:
#   1.0    In tutorial - Tutorial session 05: new script taking the
#          insights from 05-geneticCodeExperiments.R further to arrive at
#          better-than-natural codes.
#
# Preconditions:
#    AADAT and GCdf are defined via .util.R, functions aaSim(), neighCodons()
#    and rGC() are available.
#
# ToDo:
#
# Notes:
#
# ==============================================================================
#                                                                              #
#   THIS SCRIPT IS NOT "SAFE TO SOURCE". SOURCING IT WILL HAVE SIDE EFFECTS.   #
#                                                                              #
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                    Line
#TOC> ----------------------------------------
#TOC>   1        INITIALIZATIONS            42
#TOC>   1.1        Parameters               44
#TOC>   1.2        Packages                 47
#TOC>   2        RATIONALE                  50
#TOC>   3        PROCESS                    71
#TOC> 
#TOC> ==========================================================================


# =    1  INITIALIZATIONS  =====================================================

# ==   1.1  Parameters  ========================================================


# ==   1.2  Packages  ==========================================================


# =    2  RATIONALE  ===========================================================

# We have previously constructed code that evaluates the "quality" of a non-
# natural genetic code. We have also seen that the standard genetic code is
# much better than the randomized codes. But how could nature have found such a
# highly performant code, if not by random chance? And can we do even better?
# can we do even better? Let's try our hand at optimizing genetic codes.

# There are many approaches to optimization ... Let's discuss.
gAIinit()

# Task: Navigate to your AI tutor of choice and paste the initialization prompt.
# Then continue:

# Prompt:
#
#    I would like to optimize the genetic code under some objective function.
#    What options for such a multidimensional optimization problem do I have
#    in principle?


# =    3  PROCESS  =============================================================


# Prompt:
#    This is code I have for my objective function:

qGC <- function(GC) {

  sumDist <- 0                      # Initialize the objective variable

  for (codonX in names(GC)) {       # For each codon in the code
    aaX <- GC[codonX]               # ... get the encoded amino acid

    for (codonY in neighCodons(codonX)) { # For all nine neighbors of the codon
      aaY <- GC[codonY]             # ... get the encoded amino acid
      dist <- aaSim(aaX, aaY)       # ... compute distance in feature space
      sumDist <- sumDist + dist     # ... add to the sum of distances
    }
  }
  return(sumDist)
}

#    Also, I want to optimize according to the following two constraints:
#      (a) the number of stop codons "*" should not change;
#      (b) the level of redundancy in the code should not change
#    Given that, which of the optimization options is the most straightforward
#    to implement, and how would you go about writing the code?
# --



# Note: The most important part of the task is to ensure that the produced
# code is logically sound. We would also like to visualize the progress of
# the optimization, and, of course, to be able to see the resulting code.


# [END]
