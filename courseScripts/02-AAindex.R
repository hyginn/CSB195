# tocID <- "courseScripts/02-AAindex.R"
#
# Purpose: Introduction to the aaindex in R
# Version: 1.1
# Date:    2024-09-28
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#    1.1   Added exploration of values, added TOC
#    1.0   First in-class version, assigned as task
#
# ToDo:
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                           Line
#TOC> -----------------------------------------------
#TOC>   1        BASIC CONCEPTS                    43
#TOC>   2        THE INDEX                         61
#TOC>   3        WORKING WITH THE INDEX           113
#TOC>   4        EXPLORATIONS                     141
#TOC> 
#TOC> ==========================================================================


# Copy this file into your "myScripts" folder. Then open the copy for editing.
# Take notes into your copy - if you save and overwrite a course-script copy,
# it will mess up the downloading updates from GitHub but changes in your own
# copy will not be affected.

# Run the command gAIinit() to place an AI-initialization prompt into your
#  clipboard. Then open a session with your favorite generative
# AI assistant and paste the prompt. Suggested follow-on prompts
# are included below. These prompts cover a bare minimum for guidance, but do
# ask the AI to clarify any terms and concepts that you need more information
# on.


# =    1  BASIC CONCEPTS  ======================================================

# == Prompts:
#
#    What can you tell me about the AAindex collection of
#    amino acid property indices that is distributed as the "aaindex"
#    dataset in the seqinr:: package?.
# --
#    I would need some more information about what an R package is, and
#    how the data() function works.
# --
#    If I run the command in an RStudio session, how do I know that the
#    data set has been successfully made available?
# --
#



# =    2  THE INDEX  ===========================================================

# The aaindex data set was automatically loaded when your class project
# session started up, this is part of the start up procedure which runs
# the file .util.R  - you can open and inspect that file, the command that
# loads the aaindex is in section 3.
#
# If you look at the Environment pane, you will see it among the top of the
# listed elements:
#
#   aaindex   Large list (544 elements, 2.2 Mb).
#
# == Prompts:
#
#    I see that the aaindex is a "Large list". What is that? How is a list
#    different from a scalar, a vector, an array, or a data frame?
# --
#    How can I access a single element of the index to inspect it?
# --
#    Ok: I typed aaindex[[8]] and got:
#      $H
#      [1] "BHAR880101"
#
#      $D
#      [1] "Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)"
#
#      $R
#      [1] "LIT:1414112"
#
#      $A
#      [1] "Bhaskaran, R. and Ponnuswamy, P.K."
#
#      $T
#      [1] "Positional flexibilities of amino acid residues in globular proteins"
#
#      $J
#      [1] "Int. J. Peptide Protein Res. 32, 241-255 (1988)"
#
#      $C
#      [1] "VINM940103    0.869  KARP850102    0.806  WERD780101   -0.803
#           RICJ880111   -0.813"
#
#      $I
#      Ala   Arg   Asn   Asp   Cys   Gln   Glu   Gly   His   Ile
#    0.357 0.529 0.463 0.511 0.346 0.493 0.497 0.544 0.323 0.462
#      Leu   Lys   Met   Phe   Pro   Ser   Thr   Trp   Tyr   Val
#    0.365 0.466 0.295 0.314 0.509 0.507 0.444 0.305 0.420 0.386
#
#    Can you explain this?
# --


# =    3  WORKING WITH THE INDEX  ==============================================


# == Prompts:
#
#    What is the code I need to extract this index from the data set? I
#    would like the result in a "named vector" so I can easily access the
#    results.
# --
#    Nice! Can you show me how to compute whether hydrophilic amino acids
#    are more flexible than hydrophobic amino acids?
# --
#    I see, the values are different - but is the difference significant?
# --
#    I get the following result - but I still can't tell whether the difference
#    is significant (I have never taken a statistics course  :-) .
#
#      data:  flex_hydrophilic and flex_hydrophobic
#      t = 2.4407, df = 17.844, p-value = 0.02532
#      alternative hypothesis: true difference in means is not equal to 0
#      95 percent confidence interval:
#       0.01097766 0.14732537
#      sample estimates:
#      mean of x mean of y
#      0.4703333 0.3911818
# --
#

# =    4  EXPLORATIONS  ========================================================

# Let's use the aaindex to explore some amino acid properties. As part of your
# .util.R start up script, I have defined a function - grepAAindex() - that
# uses regular expressions to scan the Definitions (the $D
# list items) and return the matching list items.
#
# == Prompts:
#
#    What is "grep"?
# --
#    What is a "regular expression"?
# --
#    Please remind me of what a "list" is in R.
# --

# For example, we can serch for the string "hydro". Matches would include
# "hydrophobicity", "hydrophobic", "hydrogen" and more. Execute the next
# expression.

grepAAindex("hydro")

# For example, you should see:
#    $`68`
#    [1] "Consensus normalized hydrophobicity scale (Eisenberg, 1984)"

# ... and these are the corresponding values:

aaindex[[68]]$I

# Let's visualize this. The numbers are a collection of discrete values and
# there is no inherent ordering among them. The most appropriate visualization
# is a "barplot".

# == Prompts:
#
#    What is a "barplot" and when is it used?
# --

# First, we assign the data points to a variable name. Let's call them
# myHyd.

myHyd <- aaindex[[68]]$I

# Then we can plot a barplot with default values. It appears in the plot pane,
# the lower-right of the RStudio window.

barplot(myHyd)

# You see that only some of the names of the vector appear because there is not
#  enough space. To see all of the names you can either make the plot
# window wider (by dragging the vertical separator), or make the font smaller
# (by setting the "cex" parameter to a small value.) I also color the bars by
# a scheme of characteristic colors for the amino acids. (This is contained
# in the vector AACOLS which was also defined in .util.R)
#
# (The expression AACOLS[A2Aaa(names(myHyd))] first converts the three-letter
#  codes of aaindex into one-letter codes, and then uses the one letter codes
#  to retrieve the color values from AACOLS in the order we need.)

barplot(myHyd,
        col = AACOLS[A2Aaa(names(myHyd))],
        main = "Eisenberg Hydrophobicity Scale",
        names.arg = names(myHyd),
        cex.names = 0.4)


# Let's contrast this with an independent property: side chain volume.
# We proceed as above:

grepAAindex("volu")

# For example ...
#    $`150`
#    [1] "Side chain volume (Krigbaum-Komoriya, 1979)"

myVol <- aaindex[[150]]$I
barplot(myVol,
        col = AACOLS[A2Aaa(names(myVol))],
        main = "Side Chain Volume",
        names.arg = names(myVol),
        cex.names = 0.4)

# As you can see, the two scales are different - one way to compare such
# different properties of the same data points is to plot
# one data set against the other in a scatterplot.

# == Prompt:
#
#    What is a "scatterplot" and when is it used?
# --

# Scatterplots are one of the most fundamental paradigms of data visualization.
# In R, they are provided through the function plot().

plot(myHyd, myVol)

# Well - these are literally just scattered points. In order for this to be
# informative, we need to identify which amino acid is being represented
# by each dot. For example, we could color them with the color vector we
# have defined.

plot(myHyd, myVol,
     col = AACOLS[A2Aaa(names(myVol))])

# or, choosing filled points by setting the "pch" parameter:

?pch    # Gives you a help page ... Number 19 is what we need.

plot(myHyd, myVol,
     pch = 19,
     col = AACOLS[A2Aaa(names(myVol))])

# Better ... but even better would be to get the amino acid names in these
# positions. That requires two steps:
# - First we plot() an empty frame with axes, tick markss, labels etc.
# - Then we use the text() function to plot the names.

plot(myHyd, myVol,
     type = "n",  # "n" turns the plotting off
     main = "Amino Acid Hydrophobicity vs. Volume",
     xlab = "Hydrophobicity scale (Eisenberg 1984)",
     ylab = "Side chain volume (Krigbaum 1979)")
text(myHyd, myVol,
     labels = names(myHyd),  # the labels are the names() of either vector
     col = AACOLS[A2Aaa(names(myVol))],
     cex = 1.5)              # make them larger

# You can increase the plot by dragging the separators of the plot window.

# There are many, many more indices in the data set. So the challenge becomes
# how to use these to define "similarity", as a single value with which we
# can evaluate amino acid pairs. Explore some of the indices, extract them,
# plot them, and familiarize yourself with the data and the code to work
# with it.


# ====  DONE  ==================================================================

# If you've arrived here, congratulations. But please remember:
#
#   Have the AI work WITH you, not FOR you.
#
# Don't turn this script into copypasta - but read the answers carefully
# and ask back about anything you don't fully understand, whether biology,
# computing, or statistics.
#
# Have fun!
#





# [END]
