# tocID <- "./sequenceExplorations.R"
#
#
# Purpose: Explore a protein sequence to understand the effects of
#            evolution that shaped it.
#
#
# Version: 1.0
# Date:    2023-11
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.0  Tutorial version
#   0.5  First draft of contents. Used material from BIN-SEQA-Composition and
#          BIN-ALI-Dotplot.
#
# ToDo:
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                                   Line
#TOC> -----------------------------------------------------------------------
#TOC>   1        Preparation: packages                                     70
#TOC>   2        Introduction                                              93
#TOC>   2.01       Terminology:                                            97
#TOC>   2.02       Sources of Information                                 124
#TOC>   3        Download and inspect a PPP2R1A sequence                  158
#TOC>   3.01       Download and save a PPP2R1A sequence                   160
#TOC>   3.02       Read a sequence from file                              165
#TOC>   3.03       Analyse a protein sequence                             194
#TOC>   3.04       Barplot, and side-by-side barplot                      232
#TOC>   3.05       Plotting ratios                                        267
#TOC>   3.06       Plotting log ratios                                    285
#TOC>   3.07       Sort by frequency                                      301
#TOC>   3.08       Color by amino acid type                               319
#TOC>   3.09       Dotplot: Sequence comparison                           348
#TOC>   3.10       Customizing the dotplot                                410
#TOC>   3.11       Finding repeats in sequences                           473
#TOC>   4        Find and inspect a PPP2A and PPP2R1A structure           480
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



# =    2  Introduction  ========================================================

# In this script we explore the PP2A (Protein Phosphatase 2) subunit A.

# ==   2.01  Terminology:  =====================================================

# Phosphatase - an enzyme with a catalytic activity:
# ===========---------------------------------------
# ... the activity of a phosphatase is the removal of phosphate groups.

# Kinase - an enzyme with a catalytic activity:
# ======---------------------------------------
# ... the activity of a kinase is the addition of phosphate groups.

# ... together, phosphatases and kinases build reversible molecular on- and off
# switches. Bbiology needs to make these switches specific. This is generally
# achieved by recruiting catalytic domains or proteins to their intended site of
# action with adapter molecules. The adapter molecules can be constructed from
# simple, variable components, and the very sensitive catalytic activity is not
# jeopardized by the high variability of the adapter domains. It seems intuitive
# that such adapter proteins would be constructed from a scaffold of
# interchangeable building blocks.

# Ask ChatGPT to explain this principle:
if (FALSE) {
t2c("Please explain how phosphatases and kinases build adaptable and evolvable signalling networks for a non-specialist in which phosphate is used as a generic - indeed abstract - signalling device. Let's use the metaphor of the cell as a human city, individual people have particular skills and roles for the community, and the addition and removal of phosphate groups is like handing out an official hat, cap, beret or helmet, and retrieving it again, and through this attribute a person shifts from their private existence into their official role - and back again. With this metaphor, please explain how such a generic system well suited for both global and detailed control of large-scale processes, such as responding to growth signals, and passing through replication checkpoints; much better than delegating control to individuals' specific, bespoke decisions would be. Emphasize the crucial distinction between generic activity of kinases and phosphatases, and how they are directed through \"clerks\", i.e. the reaction-specific adapter molecules. What do these need to \"know\" so that everything ends up happening at the right time?")
}

# ...


# ==   2.02  Sources of Information  ===========================================

# In order to analyse protein sequence information, we need to access and download the sequence. There are several large, public databases that provide gene and protein information.

# The NCBI - Biomolecular data:
# ========---------------------
# Ask ChatGPT to explain:
if (FALSE) {
t2c("Can you tell me why the NCBI is important for computational biology, and perhaps walk me through the steps I need to find and download the sequence of the human PP2A (Protein Phosphatase 2) subunit A?")
}



# UniProt - Protein Sequences at the EBI:
# =======--------------------------------
# Ask ChatGPT to explain:
if (FALSE) {
t2c("Can you tell me why UniProt is important for computational biology, and perhaps walk me through the steps I need to find and download the sequence of the human PP2A (Protein Phosphatase 2) subunit A?")
}





# Ensembl - EBI Genome data:
# =======-------------------
# Ask ChatGPT to explain:
if (FALSE) {
t2c("Can you tell me about the relationship of Ensembl and UniProt? Why is Ensembl important for computational biology? And perhaps walk me through the steps I need to find and download the sequence of the human PP2A (Protein Phosphatase 2) subunit A gene? Also: what is \"Biomart\"?")
}




# =    3  Download and inspect a PPP2R1A sequence  =============================

# ==   3.01  Download and save a PPP2R1A sequence  =============================

# Instructions were presented by ChatGPT - the actual sequence is in the ./data
# directory.

# ==   3.02  Read a sequence from file  ========================================

# Generally, your sequence should be "FASTA" formatted.
if (FALSE) {
t2c("Can you explain why everybody uses FASTA format rather than a \"proper\" data grammar? If I need annotations for individual residues, how could I manage that? If you can suggest alternate or additional formats, and define their structure, that would be great.")

# When you look at the alternatives, you should be able to identify that GFF is
# a suitable format to store annotations for our generic sequence analysis use
# case. Perhaps you would like to try ...
t2c("My professor says that GFF is a suitable format to store annotations to specific sequences that we are analysing in our computational biology course. Could you imagine why?")
}

# Ok. Let's remember that for later. But now let's look at the PP2R1A sequence
# ...
if (FALSE) {

  # Read a FASTA formatted sequence.
  PPP2R1A <- Biostrings::readAAStringSet("./data/PPP2R1A.fa")

  # Turn it into a vector of one-letter symbols:
  (vPP2Aa <- unlist(strsplit(Biostrings::toString(PPP2R1A), "" )))
# Read a FASTA formatted sequence.
PPP2R1A <- Biostrings::readAAStringSet("./data/PPP2R1A.fa")

# Turn it into a vector of one-letter symbols:
(vPP2Aa <- unlist(strsplit(Biostrings::toString(PPP2R1A), "" )))

}

# ==   3.03  Analyse a protein sequence  =======================================

if (FALSE) {

(countsObs <- table(vPP2Aa))  # Use table to get the counts of each amino acid

# Is that a lot? Not very much? Is it what we would expect? Actually,
# what would we expect?
t2c(sprintf("My human PPA2 subunit A has the following amino acid composition:\n```\n %s\n%s\n```\n Is that about what I should expect?", paste(names(countsObs), collapse = "  "), paste(countsObs, collapse = " ")))

# Probably a somewhat generic response. In order to learn more, we need a
# reference data set. Several of the aaindex tables have such information:

aaindex[[459]]$D

aaindex[[459]]$I


}

# Note that the index data is percent, not counts, and that the names are
# three-letter symbols, not one letter. We need to convert this. First we change
# our PP2AA counts to percent and we call the counts "fObs" ("observed
# frequencies") ...

if (FALSE) {

  (obsData <- 100 * (countsObs / sum(countsObs)))

  # Then we convert the aaindex #459 names into one-letter symbols;
  # we call that "fRef" ("reference frequencies").

  refData <- aaindex[[459]]$I
  names(refData) <- A2Aaa(names(refData)) # Here: 3-letter to 1-letter symbols
  refData

}

# ==   3.04  Barplot, and side-by-side barplot  ================================

barplot(obsData, col = "#CCCCCC", cex.names = 0.7)
abline(h = 100/20, col="#BB0000")

barplot(refData, col = "#BB0000", cex.names = 0.7)
abline(h = 100/20, col="#555555")

# Ok: first problem - the values in obsData are in alphabetical order. But the
# values in refData are in alphabetical order of amino acid name: alanine,
# arginine, asparagine, aspartic acid ... A, R, N, D, E ... you will see this
# order a lot - one of the old biochemistry tropes in the field. So we need to
# re-order one of the vectors to match the other. That's easy though:
refData
(refData <- refData[names(obsData)])

barplot(refData, col = "#BB0000", cex.names = 0.7)
abline(h = 100/20, col="#555555")

# To compare the values, we want to see them in a barplot, side-by-side ...
barplot(rbind(obsData, refData),
        ylim = c(0, 12),
        beside = TRUE,
        col = c("#CCCCCC", "#BB0000"),
        cex.names = 0.7)
abline(h = 100/20, col="#00000044")

# ... and add a legend
legend (x = 1, y = 12,
        legend = c("PPP2R1A", "aaIndex #459"),
        fill = c("#CCCCCC", "#BB0000"),
        cex = 0.7,
        bty = "n")


# ==   3.05  Plotting ratios  ==================================================

# To better compare the values, we'll calculate ratios between
# obsData and refData

barplot(obsData / refData,
        col = "#CCCCCC",
        ylab = "Sequence / Average",
        ylim = c(0, 2.5),
        cex.names = 0.7)
abline(h = 1, col="#BB0000")
abline(h = c(1/2, 2), lty = 2, col="#BB000055")

# ... but  ratios are not very good here, since the difference in height on the
# plot now depends on the order we compare in: ratios of 1/2 and 2 (i.e. the
# dotted lines) are exactly the same "fold-difference" !


# ==   3.06  Plotting log ratios  ==============================================

# A better way to display this is to plot log(ratios).

barplot(log(obsData / refData),
        col = "#CCCCCC",
        ylab = "log(Sequence / Average)",
        ylim = log(c(1/3, 3)),
        cex.names = 0.7)
abline(h = log(1), col="#BB0000")
abline(h = log(c(1/2, 2)), lty = 2, col="#BB000055")

# Note how the two-fold difference lines are now the same distance from the
# line of equal ratio.


# ==   3.07  Sort by frequency  ================================================

barplot(sort(log(obsData / refData), decreasing = TRUE),
        ylim = log(c(1/3, 3)),
        col = "#CCCCCC",
        ylab = "log(Sequence / Average)",
        cex.names = 0.7)
abline(h = log(1), col="#BB0000")
abline(h = log(c(1/2, 2)), lty = 2, col="#BB000055")

yTxt <- log(0.9)
arrows(4, yTxt, 0, yTxt, length = 0.07)
text(5.5, yTxt, "Enriched", cex = 0.7)
yTxt <- log(1.1)
arrows(20, yTxt, 24, yTxt, length = 0.07)
text(19.5, yTxt, "Depleted", pos = 2, cex = 0.7)


# ==   3.08  Color by amino acid type  =========================================

# Color the bars by amino acid type. Use AACOLS , defined in the .utilities.R
# script, or define your own.

barplot(rep(1, 20), names.arg = names(AACOLS), col = AACOLS, cex.names = 0.5)

lR <- sort(log(obsData / refData), decreasing = TRUE)
barplot(lR,
        ylim = log(c(1/3, 3)),
        col = AACOLS[names(lR)],
        ylab = "log(Sequence / Average)",
        cex.names = 0.7)
abline(h = log(1), col="#00000055")
abline(h = log(c(1/2, 2)), lty = 2, col="#00000033")

yTxt <- log(0.9)
arrows(4, yTxt, 0, yTxt, length = 0.07)
text(5.5, yTxt, "Enriched", cex = 0.7)
yTxt <- log(1.1)
arrows(20, yTxt, 24, yTxt, length = 0.07)
text(19.5, yTxt, "Depleted", pos = 2, cex = 0.7)


# Task:
#   Interpret this plot. (Can you?) Which types of amino acids are enriched?
#   Depleted?


# ==   3.09  Dotplot: Sequence comparison  =====================================

if (FALSE) {
  t2c("What is a \"dotplot\" and what is it useful for?")

  # Our brief exploration of the structure suggested that there are many
  # repetitions in PPPP2R1A. Do they show up in the sequence? A dotplot might
  # tell...

  # Creating a dotplot to compare two sequences is actually quite simple. We
  # build a matrix that has as many rows as one sequence, as many columns as
  # another. (Obviously, if we compare a sequence with itself to find repeats,
  # this will be a square matrix.) Then we go through every cell of the matrix
  # and enter a pairwise similarity score for the amino acid pair whose position
  # corresponds to the row and column index. Finally we visualize the matrix in
  # a plot.

  # Our sequence was vPP2Aa
  str(vPP2Aa)

  # How do we get the pairwise similarity scores? Well - aaSim() was written
  # exactly for that, right? Consider:
  vPP2Aa[21]
  vPP2Aa[34]
  aaSim(vPP2Aa[21], vPP2Aa[34])

  # So the dotplot is simple:

  # First we build an empty matrix that will hold all pairscores ...
  dotMat <- matrix(numeric(length(vPP2Aa) * length(vPP2Aa)),
                   nrow = length(vPP2Aa))

  # ... then we loop over the sequences and store the scores in the matrix.
  #
  for (i in 1:length(vPP2Aa)) {
    for (j in 1:length(vPP2Aa)) {
      dotMat[i, j] <- aaSim(vPP2Aa[i], vPP2Aa[j])
    }
  }

  # Even though this is a large matrix, this does not take much time ...
  # Let's have a look at a small block of the values:

  dotMat[1:6, 1:6]

  # Note the 0 values along the diagonal.

  # To plot this, we use the image() function. Here, with default parameters.

  image(dotMat)

  # Be patient, this may takes a few moments to render: more than 34,000 values.

  # Nice.
  # What do you expect?
  # What do you see?
  # DotPlot


}


# ==   3.10  Customizing the dotplot  ==========================================

if (FALSE) {
  # Let's magnify this a bit by looking at only the first 200 amino acids ...
  image(dotMat[1:200, 1:200])

  # ... and, according to our normal writing convention, we would like the
  # diagonal to run from top-left to bottom-right since we write from left to
  # right and from top to bottom...
  image(dotMat[1:200, 1:200], ylim = 1.0:0.0)

  # ... and we would like the range of the x- and y- axis to correspond to the
  # sequence position ...
  image(x = 1:200, y = 1:200,  dotMat[1:200, 1:200], ylim=c(200,1))

  # ... and labels! Axis labels would be nice ...
  image(x = 1:200, y = 1:200,  dotMat[1:200, 1:200], ylim=c(200,1),
        xlab = "PPP2R1A", ylab = "PPP2R1A" )

  # ... and why don't we have axis-numbers on all four sides? Go, make that
  # right too ...
  len <- 200
  image(x = 1:len, y = 1:len,  dotMat[1:len, 1:len], ylim=c(len,1),
        xlab = "PPP2R1A", ylab = "PPP2R1A", axes = FALSE)
  box()
  axis(1, at = c(1, seq(10, len, by=10)))
  axis(2, at = c(1, seq(10, len, by=10)))
  axis(3, at = c(1, seq(10, len, by=10)))
  axis(4, at = c(1, seq(10, len, by=10)))

  # ... you get the idea, we can infinitely customize our plot. However a good
  # way to do this is to develop a particular view for, say, a report or
  # publication in a script and then put it into a function. I have put a
  # function into the utilities file and called it dotPlot2(). Why not dotPlot()
  # ... that's because there already is a dotplot function in the seqinr
  # package:

  seqinr::dotPlot(vPP2Aa, vPP2Aa) # seqinr
  dotPlot2(vPP2Aa, vPP2Aa, xlab = "PPP2R1A", ylab = "PPP2R1A")

  # Let's see if we can enhance the contrast between distributed noise and the
  # actual alignment of conserved residues. We can filter the dot matrix with a
  # pattern that enhances diagonally repeated values. Every value in the matrix
  # will be replaced by a weighted average of its neighborhood. Here is  a
  # diagonal-filter:

  myFilter <- matrix(numeric(25), nrow = 5)
  myFilter[1, ] <- c( 1, 0, 0, 0, 0)
  myFilter[2, ] <- c( 0, 1, 0, 0, 0)
  myFilter[3, ] <- c( 0, 0, 1, 0, 0)
  myFilter[4, ] <- c( 0, 0, 0, 1, 0)
  myFilter[5, ] <- c( 0, 0, 0, 0, 1)

  # I have added the option to read such filters (or others that you could
  # define on your own) as a parameter of the function.
  dotPlot2(vPP2Aa, vPP2Aa, xlab = "PPP2R1A", ylab = "PPP2R1A", f = myFilter)

  # I think the result finally shows up some segments of repetitions.
  # But is there not a better way?

}


# ==   3.11  Finding repeats in sequences  =====================================

t2c("I think there may be sequence repeats in the PP2A subunit A sequence. I am working in R - can you figure out how to identify sequence repeats? The sequence is in a Biostrings AA stringset object called \"PPP2R1A\", or alternatively in a single letter symbol vector called vPP2Aa. Thanks!")




# =    4  Find and inspect a PPP2A and PPP2R1A structure  ======================



# What do we learn




# [END]
