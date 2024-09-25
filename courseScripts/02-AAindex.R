# 02-AAindex.R
#
# Purpose: Introduction to the aaindex in R
# Version: 1.0
# Date:    2024-09-24
# Author:  boris.steipe@utoronto.ca
#
# ToDo:
# Notes:
#
# ==============================================================================

# Copy this file into your "myScripts" folder. Then open the copy for editing.
# Take notes into your copy - if you save and overwrite a course-script copy,
# it will mess up the downloading updates from GitHub.

# Run the command gAIinit(). Then open a session with your favorite generative
# AI assistant and paste the initializing prompt. Suggested follow-on prompts
# are included below.


# ====  BASIC CONCEPTS  ========================================================

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



# ====  THE INDEX  =============================================================

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


# ====  WORKING WITH THE INDEX  ================================================


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
