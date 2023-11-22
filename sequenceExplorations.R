# tocID <- "./sequenceExplorations.R"
#
#
# Purpose: Explore a protein sequence to understand the effects of
#            evolution that shaped it.
#
#
# Version: 1.3
# Date:    2023-11
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.3  Explore 1B3U via ChimeraX scripting (cf. RPR-ChimeraX_remote.R)
#   1.2  Add RADAR results
#   1.1  Updated the "RADAR" prompt
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
#TOC>   1        Preparation: packages                                     81
#TOC>   1.01       Preparation: Helper functions                          122
#TOC>   2        Introduction                                             139
#TOC>   2.01       Terminology:                                           143
#TOC>   2.02       Sources of Information                                 170
#TOC>   3        Download and inspect a PPP2R1A sequence                  204
#TOC>   3.01       Download and save a PPP2R1A sequence                   206
#TOC>   3.02       Read a sequence from file                              211
#TOC>   3.03       Analyse a protein sequence                             235
#TOC>   3.04       Barplot, and side-by-side barplot                      273
#TOC>   3.05       Plotting ratios                                        311
#TOC>   3.06       Plotting log ratios                                    332
#TOC>   3.07       Sort by frequency                                      351
#TOC>   3.08       Colour by amino acid type                              372
#TOC>   3.09       Dotplot: Sequence comparison                           404
#TOC>   3.10       Customizing the dotplot                                466
#TOC>   3.11       Finding repeats in sequences                           529
#TOC>   4        Find and inspect a PPP2A and PPP2R1A structure           583
#TOC>   4.01       The trouble with sequence numbers                      637
#TOC>   4.02       Mapping RADAR repeats to structure                     704
#TOC>   5        Sequence alignment / Structure superposition             743
#TOC>   5.01       A multiple structure superposition                     752
#TOC>   5.01.1         Defining secondary structure boundaries            767
#TOC>   5.02       Split the chain                                        894
#TOC>   5.03       Superpositions                                         923
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


# We need the httr:: package for remote control of ChimeraX
if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
# Package information:
#  library(help = httr)       # basic information
#  browseVignettes("httr")    # available vignettes


# Additional computations with protein structures can be done with the bio3d
# package:
if (! requireNamespace("bio3d", quietly = TRUE)) {
  install.packages("bio3d")
}
# Package information:
#  library(help = bio3d)       # basic information
#  browseVignettes("bio3d")    # available vignettes


# ==   1.01  Preparation: Helper functions  ====================================

# ... Here is a helper function to highlight specific elements of a ChimeraX
#  scene.
#
highlightCX <- function(sel, hl = "#aa5555", bg = "#aaaab5") {
  # Color everything in the scene with the bg colour, then colour sel with the
  # hl colour.
  CX("select")
  CX(sprintf("color sel %s", bg))
  CX(sprintf("select %s", sel))
  CX(sprintf("color sel %s", hl))
  CX("~select")
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

if (FALSE) {

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

}

# ==   3.05  Plotting ratios  ==================================================

if (FALSE) {

}
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

if (FALSE) {
}

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

if (FALSE) {
}

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


# ==   3.08  Colour by amino acid type  ========================================

if (FALSE) {
}

# Colour the bars by amino acid type. Use AACOLS , defined in the .utilities.R
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

if (FALSE) {
}

t2c("The protein structure looks like there may be sequence repeats in the PP2A subunit A sequence, of one or two dozen residues in length or so. Could you kindly explain the nature and significance of sequence repeats? \n\nI am working in R - can you figure out how to identify such sequence repeats? The sequence is in a Biostrings AA stringset object called \"PPP2R1A\", or alternatively in a single letter symbol vector called vPP2Aa. I have heard of a utility called \"RADAR\". Is that available for R - or is there an online server? Thanks!")

# You should be given a general overview, and hints how to access the RADAR
# server at the EBI, among the "tools" they are making available on their main
# Web page. Run a RADAR analysis. Look at the results, what do you think?

# A RADAR analysis of the PPP2R1A sequence finds 10 repeats.
# https://www.ebi.ac.uk/jdispatcher/pfa/radar

# ---------------------------------------------------------------------------
# No. of Repeats|Total Score|Length  |Diagonal| BW-From|   BW-To|   Level
#             10|     380.33|      35|      36|     162|     196|       1
# ---------------------------------------------------------------------------
#   47 -  80 (40.26/25.68)        	T.RSELLPFLTDTIYDE.DEVLLALAEQLGTFTTLV
#  162 - 196 (48.04/32.14)        	V.KAELRQYFRNLCSDDTPMVRRAAASKLGEFAKVL
#  201 - 235 (43.08/28.02)        	V.KSEIIPMFSNLASDEQDSVRLLAVEACVNIAQLL
#  240 - 274 (41.80/26.96)        	L.EALVMPTLRQAAEDKSWRVRYMVADKFTELQKAV
#  278 - 313 (39.73/25.24)        	ItKTDLVPAFQNLMKDCEAEVRAAASHKVKEFCENL
#  322 - 356 (37.87/23.70)        	I.MSQILPCIKELVSDANQHVKSALASVIMGLSPIL
#  362 - 395 (42.71/27.72)        	I.EHLLPLFLAQL.KDECPEVRLNIISNLDCVNEVI
#  400 - 422 (20.98/ 9.69)        	L.SQSLLPAIVELAEDAKWRVRLA............
#  518 - 551 (41.12/26.40)        	T.KHMLPTVLR.MAGDPVANVRFNVAKSLQKIGPIL
#  556 - 584 (24.74/12.81)        	L.QSEVKPILEKLTQDQDVDVKYFAQEALT......
# ---------------------------------------------------------------------------

# Let's put these results into a list, so we have them available. I don't know
# what the number of structure repeats will be, so I won't use a data frame
# for this case.

Annot <- list()
Annot$RADAR$start <- c( 47, 162, 201, 240, 278, 322, 362, 400, 518, 556)
Annot$RADAR$end   <- c( 80, 196, 235, 274, 313, 356, 395, 422, 551, 584)
Annot$RADAR$seq   <- c( "T.RSELLPFLTDTIYDE.DEVLLALAEQLGTFTTLV",
                        "V.KAELRQYFRNLCSDDTPMVRRAAASKLGEFAKVL",
                        "V.KSEIIPMFSNLASDEQDSVRLLAVEACVNIAQLL",
                        "L.EALVMPTLRQAAEDKSWRVRYMVADKFTELQKAV",
                        "ItKTDLVPAFQNLMKDCEAEVRAAASHKVKEFCENL",
                        "I.MSQILPCIKELVSDANQHVKSALASVIMGLSPIL",
                        "I.EHLLPLFLAQL.KDECPEVRLNIISNLDCVNEVI",
                        "L.SQSLLPAIVELAEDAKWRVRLA............",
                        "T.KHMLPTVLR.MAGDPVANVRFNVAKSLQKIGPIL",
                        "L.QSEVKPILEKLTQDQDVDVKYFAQEALT......")

# Here are the repeat lengths:
Annot$RADAR$end - Annot$RADAR$start




# =    4  Find and inspect a PPP2A and PPP2R1A structure  ======================

if (FALSE) {

  t2c("I am analysing the human PPP2R1A protein, the A subunit of PP2A. How can I find a 3D structure of the protein at the PDB? If I find several structures, how can I compare them to choose? If I find no structures, what can I do then?")

  # Among the many structures we find, the file 2IAE contains the complex, and
  # the file 1B3U. Let's explore:

  # One of the cool features of ChimeraX is that it can be driven by Python
  # code, both within a running session and through Python scripts. What I find
  # even cooler though is that ChimeraX can be driven from any programming
  # language via its remote control function that can listen to commands sent
  # from any other application. The interface that is used here is the standard
  # REST (method) - the GET and POST verbs that ubiquitously underly the
  # communication of clients and servers on the Web.

  # In order to establish the communication between this script and ChimeraX,
  # all we need to do is:
  #  - open ChimeraX;
  #  - tell it to listen on a specific "port";
  #  - send commands to that port via httr::

  # I have abstracted all of that away into a function CX() which is
  # loaded from the .utils.R script on startup.

  # (1) Open a fresh, new session in a recently updated version of ChimeraX
  # (2) Execute this command to place the string into the clipboard:
  t2c("remotecontrol rest start port 61803")
  # (3) PASTE this into the ChimeraX console and hit return. Then you can
  #     send commands to ChimeraX directly from your RStudio session.

  CX("open 1B3U")

  # The structure file contains two copies of the same protein, we delete
  # the second one:

  CX("select /B")          # First selecting, then deleting, allows you a sanity
  CX("delete atoms sel")   # check: is this really what you wanted to delete?


  CX("cofr /A")            # cofr: center-of-rotation
  CX("camera sbs")
  CX("lighting soft")
  CX("color sequential #1 & protein target abc palette powderblue:orchid:white")

  CX("label defaultheight 2")
  CX("label /A:1, 588")

  # Spend some time with this ...

}


# ==   4.01  The trouble with sequence numbers  ================================

if (FALSE) {

  # The first thing we need to consider is whether the "sequence numbers" of the
  # 3D file match the sequential numbering of the sequence. This actually really
  # important, because if they do not, we are making wrong assumptions about
  # what is what. And there are many reasons why they might be different.
  # Actually, the sequence numbers are more likely to be different than
  # identical.

  CX("sequence chain /A")

  # Unfortunately, there does not seem to be a good way to return selected
  # sequences from ChimeraX, but there is an R package - bio3d - that we can use
  # to extract the sequence from the PDB file. The process that I use below is a
  # bit "involved". But this could all be done quite easily by hand, its just a
  # bit tedious. And, TBH, it does not happen often that we have to consider as
  # many repeats as we do in our PPP2R1A sequence.

  pdb1B3U <- bio3d::read.pdb("1B3U")
  idx <- bio3d::atom.select(pdb1B3U, chain = "A", "calpha")
  pdbSeq <- bio3d::aa321(pdb1B3U$atom$resid[ idx$atom ])
  pdbSeq <- Biostrings::AAString(paste0(pdbSeq, collapse = ""))

  # Now we need to perform pairwise alignments and extract the ranges of
  # the alignment. This is best done with a function:

  getAlignmentRange <- function(fragment, sequence) {
    # Helper function to find the position of a query fragment
    # in a sequence. Returns a named vector with two integers for start and
    # end of the aligned range.
    alignment <- Biostrings::pairwiseAlignment(fragment,
                                               sequence,
                                               type="local")
    # The alignment returns an IRanges object that needs special Biostrings::
    # accessor functions - subject(), start(), end() - for processing.
    alignedSeq <- Biostrings::subject(alignment)
    range <- alignedSeq@range
    range <- c(Biostrings::start(range), Biostrings::end(range))
    names(range) <- c("start", "end")
    return(range)
  }

  # Check the RADAR reported ranges against the aligned positions in the PDB
  # sequence:

  for (i in 1:length(Annot$RADAR$seq)) {
    fragment <- toupper(gsub("[^A-Z,a-z]", "", Annot$RADAR$seq[i]))
    range <- getAlignmentRange(fragment, pdbSeq)
    cat(sprintf("%d RADAR %d:%d aligns to PDB: %d:%d\n", i,
                Annot$RADAR$start[i], Annot$RADAR$end[i],
                range["start"], range["end"]))
    Annot$PDB$start[i] <- range["start"]
    Annot$PDB$end[i]   <- range["end"]
  }

  # Notice that these indices are not the same? In this case, the difference
  # is only one residue - but it may be much larger. NEVER asssume that
  # the same sequence numbers label the same residues, unless you are working
  # from the same source. There is no authoritative reference for sequence
  # numbers! You always have to compare things based on the actual sequence
  # itself.

}


# ==   4.02  Mapping RADAR repeats to structure  ===============================

if (FALSE) {
  # Let us inspect where the RADAR segments map on the structure:
  CX("select /A")
  CX("color sequential sel palette #333339:#9999AA")  # Gradient colour
  CX("lighting flat")
  CX("lighting shadows true intensity 0.5")
  CX("graphics silhouettes true")

  # Define some colours: I am using a mauve-to-orange gradient, but "striping"
  # the colours with a cool light-gray to emphasize the separate fragments:
  myCols <- colorRampPalette(c("#662946", "#B1554D"))(length(Annot$PDB$start))
  myCols[ as.logical((0:9)%%2) ] <- "#DDDDE9"

  # Show each fragment.
  for (i in 1:length(Annot$PDB$start)) {
    CX(sprintf("select /A:%d-%d", Annot$PDB$start[i], Annot$PDB$end[i]))
    CX(sprintf("color sel %s", myCols[i]))
    CX("~select")
  }

  # This is very informative:
  #   - A canonical repeat seems to consist of two helices connected by a turn.
  #   - RADAR had identified 10 such repeats, but the structure contains
  #     five more. Repeat # 8 has only been partially identified.
  #   - Although we must assume that the repeats were at some point derived
  #     from duplications of sequence (i.e. at that point in time the sequences
  #     were identical), their sequences have drifted apart and their ancestral
  #     relationship has become nearly unrecognizable over time.
  #   - This is not at all the case when we consider structure - structure is
  #     much more highly conserved than sequence.
  #   - Still, at this point we can't be sure whether a repeated _unit_
  #     should be considered to comprise 2, 4, 6, or even more helices.

}



# =    5  Sequence alignment / Structure superposition  ========================

if (FALSE) {

  t2c("Can you explain how results from sequence alignment and structure superposition give complementary but distinct insights on protein evolution and function? In particular, how they represent information about the past and the present of a protein?")

}


# ==   5.01  A multiple structure superposition  ===============================

if (FALSE) {

  # We can perform a multiple superposition of fragments, but we need to define
  # the coordinates separately. I am doing that with four residues of overlap,
  # relative to the RADAR annotations, since we don't know precisely where we
  # should place the fragment boundaries - this will come out of the
  # superposition. And obviously, I will do that for all fifteen fragments in
  # the structure.

}



# ===   5.01.1  Defining secondary structure boundaries

if (FALSE) {

  # I started to do this by hand, but it got tedious, so I automated this. The
  # ChimeraX dssp command annotates every range and residue with its secondary
  # structure. We just need to take the helix summary information, get the
  # boundaries of secondary structure and define which ranges we want.

  CX("select #1")
  dsspReport <- CX("dssp sel minHelixLen 4 report true")
  dsspReport <- unlist(strsplit(dsspReport, '\\n'))
  sel <- grepl("^$", dsspReport)     # select empty lines
  dsspReport <- dsspReport[! sel]    # remove empty lines
  head(dsspReport, 70)               # inspect

  # The information we need is between the "Helix Summary" tag and the
  # "Ladder Summary" tag.
  iStart <- grep("Helix Summary",  dsspReport)
  iEnd   <- grep("Ladder Summary", dsspReport)
  dsspReport <- dsspReport[(iStart + 1):(iEnd - 1)]
  head(dsspReport)

  # We need to parse out two integers from each line. Let's just remove
  # everything else (except for the hyphen):

  dsspReport <- gsub("[^0-9-]", "", dsspReport)

  x <- strsplit(dsspReport, "-")  # This is one of the cases where we DON'T
  # unlist() a strsplit() result.

  for (i in 1:(length(x))) {      # Enter this into our Annotation list
    Annot$dssp$start[i] <- as.integer(x[[i]][1])
    Annot$dssp$end[i]   <- as.integer(x[[i]][2])
  }

  # ... and we can have a look at those annotations:

  CX("select /A")
  CX("color sequential sel palette #333339:#9999AA")  # Gradient colour

  # Same colour scheme as above:
  myCols <- colorRampPalette(c("#662946", "#B1554D"))(length(Annot$dssp$start))
  myCols[ as.logical((0:9)%%2) ] <- "#DDDDE9"

  # Show each fragment.
  for (i in 1:length(Annot$dssp$start)) {
    CX(sprintf("select /A:%d-%d", Annot$dssp$start[i], Annot$dssp$end[i]))
    CX(sprintf("color sel %s", myCols[i]))
    CX(sprintf("label /A:%d", Annot$dssp$start[i]))
    CX("~select")
  }

  # Now, this is seriously interesting. What appeared to be a boring bundle
  # of helices actually has fascinating internal structure. The helices are
  # kinked, there are breaks and irregularities, and it appears that
  # towards the end, in repeats 11 and 12 one helix each got lost.

  # Using the dssp annotations, we can make quick work of defining the ranges
  # (with some manual adjustments, see below). I define them so the ranges are
  # not overlapping (ChimeraX doesn't allow them to overlap, even though in this
  # case it would make sense - but I could work around that by selecting
  # segments and writing the coordinates to separate PDB files.):

  Annot$MDL$start[ 1] <-   1
  Annot$MDL$start[ 2] <-  43
  Annot$MDL$start[ 3] <-  82
  Annot$MDL$start[ 4] <- 120
  Annot$MDL$start[ 5] <- 159
  Annot$MDL$start[ 6] <- 197
  Annot$MDL$start[ 7] <- 237    # - 2
  Annot$MDL$start[ 8] <- 275
  Annot$MDL$start[ 9] <- 318    # - 2
  Annot$MDL$start[10] <- 357
  Annot$MDL$start[11] <- 396
  Annot$MDL$start[12] <- 436    # - 4
  Annot$MDL$start[13] <- 474
  Annot$MDL$start[14] <- 513
  Annot$MDL$start[15] <- 552

  Annot$MDL$end[ 1] <-  42
  Annot$MDL$end[ 2] <-  81
  Annot$MDL$end[ 3] <- 119
  Annot$MDL$end[ 4] <- 158
  Annot$MDL$end[ 5] <- 196
  Annot$MDL$end[ 6] <- 236
  Annot$MDL$end[ 7] <- 274
  Annot$MDL$end[ 8] <- 317
  Annot$MDL$end[ 9] <- 356
  Annot$MDL$end[10] <- 395
  Annot$MDL$end[11] <- 435
  Annot$MDL$end[12] <- 473
  Annot$MDL$end[13] <- 512
  Annot$MDL$end[14] <- 551
  Annot$MDL$end[15] <- 588

  # Confirm the ranges ...
  Annot$MDL$end - Annot$MDL$start

  CX("color sequential /A palette #333339:#9999AA")  # Gradient colour
  myCols <- colorRampPalette(c("#662946", "#B1554D"))(length(Annot$MDL$start))
  myCols[ as.logical((0:9)%%2) ] <- "#DDDDE9"

  CX("~label /A")

  # Show each fragment and label it with its index
  for (i in 1:length(Annot$MDL$start)) {
    CX(sprintf("color /A:%d-%d %s",
               Annot$MDL$start[i],
               Annot$MDL$end[i],
               myCols[i]))
    mid <- round((Annot$MDL$start[i] + Annot$MDL$end[i]) / 2)
    CX(sprintf("label /A:%d text \"%d\"", mid, i))
  }

  # (Note: When I first inspected this result, I noted that many of the repeats
  # have a somewhat open helix at the end that is not annotated to be helical by
  # dssp, presumably because the H-bonds are not formed in their expected way. I
  # have manually adjusted all fragments by extending them to the start of the
  # following helix (without overlap). Note that it is not clear whether the
  # loop between two repeats should be counted with the preceding or the
  # following fragment.

}



# ==   5.02  Split the chain  ==================================================

if (FALSE) {

  # With this, we can split the chain into "submodels, for superposition.

  # re-colour
  CX("color sequential /A palette #333339:#9999AA")  # Gradient colour
  myCols <- colorRampPalette(c("#662946", "#B1554D"))(length(Annot$MDL$start))
  CX("~label /A")

  # Show each fragment and label it with its index
  for (i in 1:length(Annot$MDL$start)) {
    CX(sprintf("color /A:%d-%d %s",
               Annot$MDL$start[i],
               Annot$MDL$end[i],
               myCols[i]))
    mid <- round((Annot$MDL$start[i] + Annot$MDL$end[i]) / 2)
    CX(sprintf("label /A:%d text \"%d\"", mid, i))
  }

  myModels <- sprintf("atoms /A:%d-%d", Annot$MDL$start, Annot$MDL$end)
  myModels <- paste(myModels, collapse = " ")

  CX(sprintf("split %s", myModels))

}


# ==   5.03  Superpositions  ===================================================

if (FALSE) {

  # Finally we can compute superpositions. We suprimpose the submodels 2-15 onto
  # submodel 1. The matchmaker command works from sequnce similarity plus
  # secondary structure patterns - but we already know that sequence similarity
  # is a bit weak, so we raise the relative weight of the sencondary structure
  # component from its default 0.3 to 0.67

  CX("mm #1.2-15 to #1.1 alg sw ssFraction 0.67 showAlignment true")
  CX("cofr #1.1")

  # Unfortunately, this does not compute a superposition-induced sequence
  # alignment. We will need to get that with other means. It is one of the
  # most important tools of structure analysis

  # The superposition is not well visualized in the cartoon representation, we
  # should look at actual backbone atoms instead.

  CX("select #1.1-15")
  CX("cartoon hide sel")
  CX("show (sel-residues & backbone) target ab")
  CX("~select")

  highlightCX("#1.2,4,13")
  # Fragments 2, 4 and 13 do not align well ... we should suspect an indel in
  # the sequence that changes the arrangement of the helices, or perhaps  a
  # proline residue forcing a kink - such sequence signals would be typical for
  # the detailed structuring of interfaces. Let's hide these fragments for the
  # time being:

  CX("select #1.2,4,13")
  CX("hide sel")

  # Speaking of proline: prolines and glycines often play an important role in
  # structural motifs since some backbone turn angles are only accessible to
  # glycine and proline.

  highlightCX("#1.1-15:PRO")

  # Well! Look at that.


}



# ==============================================================================

# To come:
#  - Superposition induced sequence alignment
#  - Multiple sequence alignment
#  - Consensus sequences and sequence logos
#  - Structural sequence signals of this repeat
#  - Functional signals of this protein
#  - ...



# [END]
