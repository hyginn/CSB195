# tocID <- "src/R/aminoAcidSimilarity.R"
#
# Purpose: Develop a function to compute a pairwise similarity score for
#          amino acids. The score is based on the distance between two amino
#          acids in a feature spce derived from AAindex indices via their
#          AAontology categories.
#
# Version: 3.0
# Date:    2025-09
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   3.0    2025 standalone version to support Articulate Programming mode. First
#          raw "in class" version.
#   2.1    Some cleanup after class
#   2.0    Work from AAontology considerations, refactor all code
#   0.5    Split file into two portions:
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
# Preconditions:
#   Supplementary Data files 1 and 3 from Breimann et al. (2024)
#   https://www.sciencedirect.com/science/article/pii/S0022283624003267
#
# ToDo:
#
# Notes:  This script contains the rationale for the definition of aaSim().
#         The actual function is created when source()'ing R/aaSim.R
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                          Line
#TOC> --------------------------------------------------------------
#TOC>   01       PROLOGUE                                         94
#TOC>   01.1       AI                                            127
#TOC>   02       BASIC CONCEPTS                                  183
#TOC>   02.1       AAindex                                       209
#TOC>   02.2       Working with AAindex data                     257
#TOC>   02.3       EXPLORATIONS                                  276
#TOC>   03       A COMPUTABLE DEFINITION OF SIMILARITY           448
#TOC>   04       PROBLEMS WITH THE AAINDEX                       455
#TOC>   05       READING PUBLISHED DATA                          478
#TOC>   05.1       Read Source Data                              511
#TOC>   05.2       Scale Raw Values                              552
#TOC>   06       MERGING DATA  FROM TWO SOURCES                  581
#TOC>   06.1       Merge Into Single Data Frame                  639
#TOC>   06.2       Validations                                   676
#TOC>   07       PARALLEL COORDINATES PLOT                       713
#TOC>   08       PREPROCESS THE DATA                             732
#TOC>   08.1       Validate                                      735
#TOC>   08.1.1         Validate Scaling                          736
#TOC>   08.1.2         Check outliers                            740
#TOC>   08.1.3         Validate merging                          746
#TOC>   09       SELECT SCALES                                   755
#TOC>   10       PCA                                             861
#TOC>   11       AMINO ACID FEATURE SPACE                        909
#TOC>   12       AMINO ACID SIMILARITY                           983
#TOC>   13       INITIALIZATION                                 1080
#TOC>   13.1       Packages                                     1083
#TOC>   13.2       Data                                         1096
#TOC>   13.3       Functions                                    1144
#TOC> 
#TOC> ==========================================================================


# =    01  PROLOGUE  ===========================================================

#   This script has two processing layers:
#
#   If you source() the script, packages will be installed that you do not
#   have yet, and some global variables and functions in your
#   "workspace" will be defined. The script only writes files if they do not
#   exist already (we do not overwrite). In that case, the script continues
#   to run. The function will abort however if it doesn't see the standard
#   course directories. In that case you might be doing something unintended
#   and its better not to proceed and to alert you.
#
#   Other parts of the script are interactive. You execute them line by line,
#   and make an attempt to understand the code - even if you did not write it.
#   Or you carefully read the explanations, and invoke your AI alter ego to
#   understand the background. Blocks of this interactive layer are enclosed in
#
#   if (FALSE) {
#
#   }
#   ... statements, that is, they will not be entered in normal execution.
#
#   Therefore:
#
#   In order to get started, source() this script. Sourcing this script skips
#   all the interactive parts, and loads a block of initializations that I have
#   placed at the end of the script.
#
#   Ask your instructor in case you do NOT see the message "Done,
#   successfully initialized." in your console after sourcing.



# ==   01.1  AI  ===============================================================


# We will work through this R script in class, but there are far too many
# concepts to cover in detail during that time. The main goal is to understand
# the computational approach in principle—not to memorize or write everything
# from scratch.

# An AI assistant can be an excellent code tutor for your own study because it
# knows R well and is very patient. However, it’s up to you to ask questions
# whenever you encounter something unfamiliar or that you feel you *should*
# know.

# To make the most of this, it helps to set the context for the AI so that its
# responses align with the principles used in this script. Open your favorite AI
# assistant and paste the initialization prompt below. These are ground rules
# for the AI tutor, not rules for you. Copy the text between the quote marks
# and paste it into the AI input.

# Note: in the CONTEXT section we tell the AI which computer platform we are
# working on. Change that if required.

"
You are my R tutor, and as a novice I need your help with code and concepts. When we discuss, please follow the following guidlines.

CONTEXT
* I am working on a Mac computer.
* I use the RStudio IDE.
* This work is in a first year University course _Computational Biology Foundations_.

STYLE & SCOPE
* End every response with a single # character so I know the response is complete.
* Explain syntax and concepts at a novice level; be concise and stay on topic.
* Prefer base R solutions over additional packages when feasible.
* Prefer package::function() style, over library(package) - it makes the source of functions explicit for me.
* Write function names in your responses with parentheses (e.g., rnorm(), c()) and packages with double colons (e.g., httr::, utils::) - so I am not confused about which is which.
* Use camelCase in examples for variable names. We use snake_case for python.
* Prefer dataFrame[ , col] over dataFrame[[col]] when practical.
* Executing expressions to show their return values is fine (e.g., path.expand("~")); avoid wrapping this in unnecessary print() statements.
* Avoid _magic numbers_: define parameters as variables at the top of a function or before complex calls.

EXPLICITNESS & DEBUGGABILITY
* Favour explicit, stepwise code that is easy to debug; avoid hiding behavior behind implicit conventions. For example, always use return() in functions.
* Minimize deep pipelines; prefer small, named intermediate objects with clear intent, for example use constructs like sel <- is.na(dat); dat <- dat[! sel].

PACKAGE CHOICES
* By default, avoid tidyverse functions and ggplot2.
* If I explicitly ask for tidyverse/ggplot2 (or other packages), you should show me how, but you might briefly note the trade-offs over base-R.

ROBUSTNESS
* My prompts may include leading comment characters (#) when copied from scripts. Please ignore them.
* If my question is ambiguous, ask for clarification; otherwise answer directly and concisely.

Please confirm with one word.
"

# =    02  BASIC CONCEPTS  =====================================================

# In our quest for a computable definition of amino acid similarity, we are
# going to use scales that were determined a few decades ago and collated in
# the AAindex database:
#
#    https://www.genome.jp/aaindex/
#

# == Prompts:
#
#    What can you tell me about the AAindex collection of
#    amino acid property indices that is distributed as the "aaindex"
#    dataset in the seqinr:: package?.
# --
#    I would need some more information about what an R package is, where it
#    is loaded from, and whether it is safe to use.
# --
#    If I run the command:
#      data(aaindex, package = "seqinr")
#    ... in an RStudio session, how do I know that the
#    data set has been successfully made available?
# --
#


# ==   02.1  AAindex  ==========================================================

# The aaindex data set was loaded for this script during initialization
#
# If you look at the Environment pane in RStudio , you will see it at the
# top of the listed elements:
#
#   aaindex   Large list (544 elements, 2.2 Mb).
#
# == Prompts:
#
#    I have loaded the aaindex dataset from the seqinr:: package, now
#    I have the object aaindex in my workspace. How can I access a single
#    element of the aaindex list to inspect it?
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


# ==   02.2  Working with AAindex data  ========================================


# == Prompts:
#
#    What is the code I need to extract this index (number 8) from the
#    data set? I would like the result in a "named vector" so I can easily
#    access the results.
# --
#    Nice! Can you show me how to compute whether hydrophilic amino acids
#    are more flexible than hydrophobic amino acids? Let's do this as
#    concise as reasonably possible without repeating things I have already
#    done, like loading the index.
# --
#    I see a value - but what does it mean? And is this significant?
#   (I have never taken a statistics course  :-) .
# --
#

# ==   02.3  EXPLORATIONS  =====================================================

# Let's use the aaindex to explore some amino acid properties. During
# initialization, we defined a function - grepAAindex() - that uses
# regular expressions to scan the Definitions (the $D list items) and
# return the matching list items.
#
# == Prompts:
#
#    What is "grep"?
# --
#    What is a "regular expression"?
# --
#    Please remind me of what a "list" is in R.
# --

# For example, we can search for the string "hydro". Matches would include
# "hydrophobicity", "hydrophobic", "hydrogen" and more. Execute the next
# expression.
if (FALSE) {

grepAAindex("hydro")


}

# For example, you should see:
#    EISD840101`
#    [1] "Consensus normalized hydrophobicity scale (Eisenberg, 1984)"

# ... and these are the corresponding values:

if (FALSE) {

aaindex[["EISD840101"]]$I


}

# Let's visualize this. The numbers are a collection of discrete values and
# there is no inherent ordering among them. The most appropriate visualization
# is a "barplot".

# == Prompts:
#
#    What is a "barplot" and when is it used?
# --

# First, we assign the data points to a variable name. Let's call them
# myHyd.

if (FALSE) {

myHyd <- aaindex[["EISD840101"]]$I


}

# Then we can plot a barplot with default values. It appears in the plot pane,
# the lower-right of the RStudio window.

if (FALSE) {

barplot(myHyd)


}

# You see that only some of the names of the vector appear because there is not
# enough space. To see all of the names you can either make the plot
# window wider (by dragging the vertical separator), or make the font smaller
# (by setting the "cex" parameter to a small value.)

# It is very helpful to color the bars by
# a scheme of characteristic colors for the amino acids. (This is contained
# in the vector AACOLS which we also defined above).

# (The expression AACOLS[A2Aaa(names(myHyd))] first converts the three-letter
#  codes of aaindex into one-letter codes, and then uses the one letter codes
#  to retrieve the color values from AACOLS in the order we need.)

if (FALSE) {

  barplot(myHyd,
          col = AACOLS[Aaa2A(names(myHyd))],
          main = "Eisenberg Hydrophobicity Scale",
          names.arg = names(myHyd),
          cex.names = 0.4)


  # Let's contrast this with an independent property: side chain volume.
  # We proceed as above:

  grepAAindex("volu")

  # For example ...
  #    KRIW790103
  #    [1] "Side chain volume (Krigbaum-Komoriya, 1979)"

  myVol <- aaindex[["KRIW790103"]]$I
  barplot(myVol,
          col = AACOLS[Aaa2A(names(myVol))],
          main = "Side Chain Volume",
          names.arg = names(myVol),
          cex.names = 0.4)

}

# As you can see, the two scales are different - one way to compare such
# different properties of the same data points is to plot
# one data set against the other in a scatterplot.

# == Prompt:
#
#    What is a "scatterplot" and when is it used?
# --

if (FALSE) {
  # Scatterplots are one of the most fundamental paradigms of data visualization.
  # In R, they are provided through the function plot().

  plot(myHyd, myVol)

  # Well - these are literally just scattered points. In order for this to be
  # informative, we need to identify which amino acid is being represented
  # by each dot. For example, we could color them with the color vector we
  # have defined.

  plot(myHyd, myVol,
       col = AACOLS[Aaa2A(names(myVol))])

  # or, choosing filled points by setting the "pch" parameter:

  ?pch    # Gives you a help page ... Number 19 is what we need.

  plot(myHyd, myVol,
       pch = 19,
       col = AACOLS[Aaa2A(names(myVol))])

  # Better ... but even better would be to get the amino acid names in these
  # positions. That requires two steps:
  # - First we plot() an empty frame with axes, tick marks, labels etc.
  # - Then we use the text() function to plot the names.

  plot(myHyd, myVol,
       type = "n",  # "n" turns the plotting off
       main = "Amino Acid Hydrophobicity vs. Volume",
       xlab = "Hydrophobicity scale (Eisenberg 1984)",
       ylab = "Side chain volume (Krigbaum 1979)")
  text(myHyd, myVol,
       labels = names(myHyd),  # the labels are the names() of either vector
       col = AACOLS[Aaa2A(names(myVol))],
       cex = 1.5)              # make them larger
}

# You can increase the plot by dragging the separators of the plot window.

# There are many, many more indices in the data set. So the challenge becomes
# how to use these to define "similarity", as a single value with which we
# can evaluate amino acid pairs.
#
# TASK:
# Explore some of the indices, extract them,
# plot them, and familiarize yourself with the data and the code to work
# with it.
#
# TASK:
# What if you use three different indices. Can you plot 3D data in R. How?
#
#


# =    03  A COMPUTABLE DEFINITION OF SIMILARITY  ==============================

# Similarity is based on properties that distinguish amino acids. Let us look at two such properties.

# ... barplot, and 2D plot.


# =    04  PROBLEMS WITH THE AAINDEX  ==========================================

# Large collections of observations like the aaindex usually cannot be used
# as-is for data aggregation:
#   - the data may be redundant
#   - features that are similarly important may have been observed with
#     very different frequencies
#   - not all features may be relevant for the question at hand, some may be
#     irrelevant and others even lead to circular reasoning.
#
# What we need is a structured categorization of the data that allows us to
# group, select, and weight the individual observations.

# == Prompts:
#
#    What is an ontology in data science?
# --
#    In that context: what is a category?
# --
#    What could category membership be based on in principle?
# --


# =    05  READING PUBLISHED DATA  =============================================


# Serendipitously, last year a paper was published in the Journal of
# Molecular Biology in which the authors construct an ontology of
# categories and subcategories for the AAindex. You can access the paper here:
#
# https://www.sciencedirect.com/science/article/pii/S0022283624003267
#
# Two of the supplementary data files should have been downloaded when you
# initialized this script. So, first, we confirm that the files exist in the
# expected location:

# Supplementary Table 1:
#   https://docs.google.com/spreadsheets/d/1xtoQ6dggalABsJEe248C74Zw1HTiMe7gLqdp3tGO3_M
#
#   Open this, and go to the "Normalized" sheet.
#
# Supplementary Table 3:
#   https://docs.google.com/spreadsheets/d/1f14pDvEP5cvIM-H6ghY9OyH9v8ark0yT5h_BDmKzPP8
#
#   Open this, and go to the "Scales" sheet.
#
# You notice that in the table that contains the actual values, each dataset
# is contained in a COLUMN, but in the table that contains the category and
# subcategory information, the datasets are listed in ROWS. For simplicity,
# we need all of this data in one data frame.

# == Prompt:
#
#    Please remind me of what a "data frame" is in R.
# --

# ==   05.1  Read Source Data  =================================================


# First, we will read the data into R. It is now in the MicroSoft Excel
# spreadsheet format.

# == Prompt:
#
#    How can I read an Excel file (.xlsx) into an R dataframe?
# --

if (FALSE) {

  list.files("data/", pattern = "xlsx$")



  # ... This MUST give you both files:
  # [1] "Breimann_2024_Supplementary_Table_1.xlsx"
  #     "Breimann_2024_Supplementary_Table_3.xlsx"
  # (If you don't get those two files listed, you can't continue. Ask for help
  #  to fix the problem.)

  # ... then we read them into R. Note that a data frame can only
  # contain a single "sheet", so if the Excel file contains multiple sheets
  # ( as these two do) we need to specify which sheet we want.

  myScales <- readxl::read_excel("data/Breimann_2024_Supplementary_Table_1.xlsx",
                                 sheet = "Normalized")

  myCats   <- readxl::read_excel("data/Breimann_2024_Supplementary_Table_3.xlsx",
                                 sheet = "Scales")

  # The resulting objects are "tibbles" - special data structures that are
  # central to the "tidyverse" programming paradigm. Since I want to avoid this
  # in our course, I convert the objects to normal R data frames.
  myScales <- as.data.frame(myScales)
  myCats   <- as.data.frame(myCats)

}

# ==   05.2  Scale Raw Values  =================================================

# Since the dimensionality reduction method that we will use later (PCA) is
# sensitive to the variance of the data, we will not use the min/max scaling
# between 0 and 1 that the authors used (and called "Normalization"), but a
# proper scaling that transforms the data to a mean of 0 and a
# standard deviation of 1.

if (FALSE) {


  x <- myScales   # make a temporary copy for validation

  for (i in 2:ncol(myScales)) {
    myScales[ , i] <- scale(as.matrix(as.numeric(myScales[ , i])))
  }

  # Validate
  idx <- 63                           # Pick a random dataset
  mean(myScales[ , idx])              # Must be (nearly) zero
  sd(myScales[ , idx])                # Must be (nearly) one
  plot(myScales[ , idx], x[ , idx])   # Points exactly on the diagonal

  # Clean up
  rm(x)
}



# =    06  MERGING DATA  FROM TWO SOURCES  =====================================


# I have produced a summary of the data to share it with the AI. Try this:

if (FALSE) {

  str(myScales, list.len = 5)
  str(myCats)

}

# == Prompt:
#
#    I have read two tables and converted them into data frames. Here is
#    what they look like:
#  > str(myScales, list.len = 5)
#  'data.frame':	20 obs. of  587 variables:
#    $ AA        : chr  "A" "C" "D" "E" ...
#  $ ANDN920101: num  0.494 0.864 1 0.42 0.877 0.025 0.84 0 0.506 0.272 ...
#  $ ARGP820101: num  0.23 0.404 0.174 0.177 0.762 0.026 0.23 0.838 0.434 0.577 ...
#  $ ARGP820102: num  0.355 0.579 0 0.019 0.601 0.138 0.082 0.44 0.003 1 ...
#  $ ARGP820103: num  0.504 0.387 0 0.032 0.67 0.17 0.053 0.543 0.004 0.989 ...
#  [list output truncated]
#
#  > str(myCats)
#  'data.frame':	586 obs. of  5 variables:
#   $ scale_id         : chr  "LINS030110" "LINS030113" "JANJ780101" ...
#   $ category         : chr  "ASA/Volume" "ASA/Volume" "ASA/Volume" ...
#   $ subcategory      : chr  "Accessible surface area (ASA)"  ...
#   $ scale_name       : chr  "ASA (folded coil/turn)"  ...
#   $ scale_description: chr  "Total median accessible [...]"  ...
#
#  Now I would like to combine this information into one single data frame, by
#  appending the amino acid values in the myScales columns to the respective
#  myCats rows. It looks like the column headers of the myScales data frame, and
#  the scale_id of the myCats data frame can be used to combine the data.
#  But they are not in the same order, and I don't know whether all ids are
#  contained in both data frames exactly once. How can I check that before I
#  continue? Can you please suggest code and break this down step by step?
# --


# I get clearly commented and quite reasonable code. You need to do this
# yourself however. What is the result - are the ids identical or not?

# == Prompt:
#    Now I need to take the twenty amino acid values from the myScales
#    data frame and put those into additional comlumns of myCats. Can you
#    please suggest more code and again break this down step by step? Please
#    make sure that the column headers of the new columns have the correct
#    amino acid one letter codes from column one of myScales.
# --

# The code I get here again looks reasonable. Please copy the suggested code
# into your script and execute it. Since it changes the structure of the
# original myCats data frame, make a copy with the name "aaontology".

# ==   06.1  Merge Into Single Data Frame  =====================================

if (FALSE) {


  # Doing this in a "pedestrian" way, without special functions ...

  aaontology <- myCats
  rownames(aaontology) <- aaontology$scale_id

  # 1: Create the extra columns and give them a value of zero
  for (i in 1:nrow(myScales)) {
    aaontology[ , ncol(aaontology) + 1] <- 0
    colnames(aaontology)[ncol(aaontology)] <- myScales$AA[i]
  }

  # Check
  ncol(aaontology)  # Must now be 25
  head(aaontology)

  # 2: for each scale in myScales, put the values in the appropriate row
  for (i in 2:ncol(myScales)) {
    id <- colnames(myScales)[i]
    idx <- which(aaontology$scale_id == id)
    stopifnot(length(idx) == 1)
    stopifnot(sum(aaontology[idx, 6:25]) == 0)
    aaontology[idx, 6:25] <- myScales[ , i]
  }

  # Confirm that this is what we expect
  id <- colnames(myScales)[63]       # Pick a random scale_id
  all(myScales[ , id] == aaontology[id, 6:25])

}



# ==   06.2  Validations  ======================================================


# Of course, we now have to validate the result.
#
# As a first sanity check, do the following:
if (FALSE) {

  aaontology[nrow(aaontology), ]  # prints out the contents of the last row

}

# My result is:
#          scale_id category subcategory
#    586 ZIMJ680105   Others        PC 2
#                               scale_name
#    586 Principal Component 1 (Zimmerman)
#                       scale_description     A C D     E F     G
#    586 RF rank (Zimmerman et al., 1968) 0.444 0 0 0.025 1 0.175
#            H     I     K     L     M     N    P     Q     R     S
#    586 0.338 0.894 0.044 0.925 0.756 0.162 0.75 0.388 0.112 0.256
#            T     V     W     Y
#    586 0.419 0.719 0.894 0.762
#
# ... your result should look the same.

# == Prompt:
#
#    I have called the resulting data frame "aaontology". Can you suggest
#    some tests to validate that the result is complete and correct?
# --

# When I executed the tests suggested by ChatGPT, one of the tests actually
# failed. But I gave the AI the failing output and it was able to correctly
# infer the problem and suggest a fix.


# =    07  PARALLEL COORDINATES PLOT  ==========================================

# It should be straightforward to e.g. reproduce a part of Figure 9 (Shape
# category, subcategory "Side chain length").

# == Prompt:
#
#    Great. I am ready to look at some data. Could you help me with code for a
#    "parallel coordinates plot" of all indices in the subcategory
#    "Side chain length"? I would like to have the lines all in one color -
#    let's use turquoise - and for some reason the amino acids should be ordered
#    as AILMPVFWYCGNQSTDEMKR. The labels should appear in the x-axis. And please
#    add a dashed black line at the value 0.5. Thanks!
# --


# Pretty neat - that works for me.


# =    08  PREPROCESS THE DATA  ================================================


# ==   08.1  Validate  =========================================================
# ===   08.1.1  Validate Scaling                     
#
# Confirm that all means are zero and all SDs are 1.

# ===   08.1.2  Check outliers                       
#
# Identify any value outside the range (-3, 3)

if (FALSE) {

# ===   08.1.3  Validate merging                     
for (i in 2:ncol(myScales)) {
  id <- colnames(myScales)[i]
  stopifnot(all(myScales[ , id] == aaontology[id, 6:25]))
}

}


# =    09  SELECT SCALES  ======================================================

if (FALSE) {

  # Prepare a data frame of subcategory IDs and 0 / 1 selectors

  mySubCats <- data.frame(cat = aaontology$category,
                          subcat = aaontology$subcategory,
                          inc = numeric(nrow(aaontology)))
  mySubCats <- mySubCats[! duplicated(mySubCats$subcat), ]

  for (i in 1:nrow(mySubCats)) {
    cat(sprintf("mySubCats$inc[%d] <- 0    # %s : %s\n",
                i,
                mySubCats$cat[i],
                mySubCats$subcat[i]))
  }

}


if (FALSE) {


  # We select all subcategories that depend on biophysical properties, but
  # not on the propensity for an amino acid to be selected in a specific
  # biological context. That might depend on the genetic code (among other
  # factors).

  mySubCats$inc[ 1] <- 1    # ASA/Volume : Accessible surface area (ASA)
  mySubCats$inc[ 2] <- 1    # ASA/Volume : Buried
  mySubCats$inc[ 3] <- 1    # ASA/Volume : Hydrophobic ASA
  mySubCats$inc[ 4] <- 1    # ASA/Volume : Partial specific volume
  mySubCats$inc[ 5] <- 1    # ASA/Volume : Volume
  mySubCats$inc[ 6] <- 0    # Composition : AA composition
  mySubCats$inc[ 7] <- 0    # Composition : AA composition (surface)
  mySubCats$inc[ 8] <- 0    # Composition : Membrane proteins (MPs)
  mySubCats$inc[ 9] <- 0    # Composition : Mitochondrial proteins
  mySubCats$inc[10] <- 0    # Composition : MPs (anchor)
  mySubCats$inc[11] <- 0    # Composition : Unclassified (Composition)
  mySubCats$inc[12] <- 0    # Conformation : Coil
  mySubCats$inc[13] <- 0    # Conformation : Coil (C-term)
  mySubCats$inc[14] <- 0    # Conformation : Coil (N-term)
  mySubCats$inc[15] <- 0    # Conformation : Linker (>14 AA)
  mySubCats$inc[16] <- 0    # Conformation : Linker (6-14 AA)
  mySubCats$inc[17] <- 0    # Conformation : Unclassified (Conformation)
  mySubCats$inc[18] <- 0    # Conformation : α-helix
  mySubCats$inc[19] <- 0    # Conformation : α-helix (C-cap)
  mySubCats$inc[20] <- 0    # Conformation : α-helix (C-term)
  mySubCats$inc[21] <- 0    # Conformation : α-helix (C-term, out)
  mySubCats$inc[22] <- 0    # Conformation : α-helix (left-handed)
  mySubCats$inc[23] <- 0    # Conformation : α-helix (N-cap)
  mySubCats$inc[24] <- 0    # Conformation : α-helix (N-term)
  mySubCats$inc[25] <- 0    # Conformation : α-helix (N-term, out)
  mySubCats$inc[26] <- 0    # Conformation : α-helix (α-proteins)
  mySubCats$inc[27] <- 0    # Conformation : β/α-bridge
  mySubCats$inc[28] <- 0    # Conformation : β-sheet
  mySubCats$inc[29] <- 0    # Conformation : β-sheet (C-term)
  mySubCats$inc[30] <- 0    # Conformation : β-sheet (N-term)
  mySubCats$inc[31] <- 0    # Conformation : β-strand
  mySubCats$inc[32] <- 0    # Conformation : β-turn
  mySubCats$inc[33] <- 0    # Conformation : β-turn (C-term)
  mySubCats$inc[34] <- 0    # Conformation : β-turn (N-term)
  mySubCats$inc[35] <- 0    # Conformation : β-turn (TM helix)
  mySubCats$inc[36] <- 0    # Conformation : π-helix
  mySubCats$inc[37] <- 1    # Energy : Charge
  mySubCats$inc[38] <- 1    # Energy : Charge (negative)
  mySubCats$inc[39] <- 1    # Energy : Charge (positive)
  mySubCats$inc[40] <- 1    # Energy : Electron-ion interaction pot.
  mySubCats$inc[41] <- 1    # Energy : Entropy
  mySubCats$inc[42] <- 1    # Energy : Free energy (folding)
  mySubCats$inc[43] <- 1    # Energy : Free energy (unfolding)
  mySubCats$inc[44] <- 1    # Energy : Isoelectric point
  mySubCats$inc[45] <- 1    # Energy : Non-bonded energy
  mySubCats$inc[46] <- 1    # Energy : Unclassified (Energy)
  mySubCats$inc[47] <- 0    # Others : Mutability
  mySubCats$inc[48] <- 0    # Others : PC 1
  mySubCats$inc[49] <- 0    # Others : PC 2
  mySubCats$inc[50] <- 0    # Others : PC 3
  mySubCats$inc[51] <- 0    # Others : PC 4
  mySubCats$inc[52] <- 0    # Others : PC 5
  mySubCats$inc[53] <- 0    # Others : Unclassified (Others)
  mySubCats$inc[54] <- 1    # Polarity : Amphiphilicity
  mySubCats$inc[55] <- 1    # Polarity : Amphiphilicity (α-helix)
  mySubCats$inc[56] <- 1    # Polarity : Hydrophilicity
  mySubCats$inc[57] <- 1    # Polarity : Hydrophobicity
  mySubCats$inc[58] <- 1    # Polarity : Hydrophobicity (interface)
  mySubCats$inc[59] <- 1    # Polarity : Hydrophobicity (surrounding)
  mySubCats$inc[60] <- 1    # Polarity : Unclassified (Polarity)
  mySubCats$inc[61] <- 0    # Shape : Graph (1. eigenvalue)
  mySubCats$inc[62] <- 0    # Shape : Graph (2. eigenvalue)
  mySubCats$inc[63] <- 1    # Shape : Reduced distance
  mySubCats$inc[64] <- 1    # Shape : Shape and Surface
  mySubCats$inc[65] <- 1    # Shape : Side chain length
  mySubCats$inc[66] <- 1    # Shape : Steric parameter
  mySubCats$inc[67] <- 0    # Shape : Unclassified (Shape)
  mySubCats$inc[68] <- 1    # Structure-Activity : Backbone-dynamics (-CH)
  mySubCats$inc[69] <- 1    # Structure-Activity : Backbone-dynamics (-NH)
  mySubCats$inc[70] <- 1    # Structure-Activity : Flexibility
  mySubCats$inc[71] <- 1    # Structure-Activity : Flexibility (2 rigid neighbors)
  mySubCats$inc[72] <- 1    # Structure-Activity : Stability
  mySubCats$inc[73] <- 1    # Structure-Activity : Stability (helix-coil)
  mySubCats$inc[74] <- 0    # Structure-Activity : Unclassified (Structure-Activity)

}

# =    10  PCA  ================================================================

# Prompt:
#    Please explain PCA for data analysis.
# --

# Create a selection vector for the data columns we want to include

if (FALSE) {

  sel <- logical(nrow(aaontology))
  for (i in 1:nrow(aaontology)) {
    idx <- which(aaontology$subcategory[i] == mySubCats$subcat)
    if (mySubCats$inc[idx] == 1) {
      sel[i] <- TRUE
    } else {
      sel[i] <- FALSE
    }
  }

}

# Validate: check all selected subcategories
# ...

if (FALSE) {


  # Apply Principal Component Analysis (PCA) to the selected data, to
  # compute components that have reduced dimensions, are orthogonal, and
  # are ordered by importance (contribution to variance).
  aaPCA <- prcomp(aaontology[sel, 6:25], scale. = FALSE)

  # Identify the first principal components that capture 95% of
  # variance as a basis for computing a similarity measure.
  cumulativeProportions <- cumsum(aaPCA$sdev^2) / sum(aaPCA$sdev^2)
  barplot(cumulativeProportions,
          main="Cumulative Proportion of Variance Explained",
          xlab="Principal Component",
          ylab="Cumulative Proportion",
          ylim=c(0, 1),
          border="blue",
          col="lightblue")
  abline(h=0.95, col="red", lty=2)  # A horizontal line to indicate the 95% mark

  numPCsToRetain <- which(cumulativeProportions >= 0.95)[1]  # 12
}

# =    11  AMINO ACID FEATURE SPACE  ===========================================

if (FALSE) {


  # Extract these components into a suitable data structure.
  aaFeatureSpace <- aaPCA$rotation[, 1:numPCsToRetain]
  rownames(aaFeatureSpace) <- rownames(aaPCA$rotation)
  colnames(aaFeatureSpace) <- paste0("PC", 1:numPCsToRetain)
  str(aaFeatureSpace)  # inspect

  # Let's look at the numerical distribution of values.
  for (i in 1:ncol(aaFeatureSpace)) {
    print(mean(aaFeatureSpace[, i]))
  }
  # All the means are "zero" (or at least very, very, very close to zero).

  for (i in 1:ncol(aaFeatureSpace)) {
    print(sd(aaFeatureSpace[, i]))
  }
  # All the standard deviations are 0.229. But that's actually not good. If the
  # values of the higher PCs are distributed the same as those of the first or
  # second PCs, then we give them an inappropriately large contribution
  # to the distance calculation. We should re-scale the columns according
  # to each PCs contribution to the total variance.

  # We can get the variances from aaPCA$sdev - the variances are the squares
  # of the standard deviations:

  aaPCA$sdev[1:ncol(aaFeatureSpace)]    # The standard deviations
  aaPCA$sdev[1:ncol(aaFeatureSpace)]^2  # The variances

  # To scale the PCs that make up the dimensions of our feature space, we can
  # simply multiply them by the variances. (Remember, the standard deviations of
  # the PCs are already all the same. If the standard deviations would have been
  # different, we would have divided by the standard deviation first.)

  # Here we go:
  for (i in 1:ncol(aaFeatureSpace)) {
    aaFeatureSpace[, i] <- aaFeatureSpace[, i] * ( aaPCA$sdev[i]^2 )
  }

  # Save the resulting feature space so we don't have to re-compute it when
  # we need it. saveRDS() saves a single R object in a compressed format. It is
  # safer to use than save() (see below).
  saveRDS(aaFeatureSpace, file = "data/aaFeatureSpace.2025.Rds")

  # saveRDS() / readRDS() save and _recreate_ single compressed R objects.
  # You can assign them  to a new variable (or the same variable name).
  # This is in contrast to save() / load(), which _restores_ an R object to its
  # original name, and overwrites existing objects if you are not careful.
  #
  # Execute this to recreate aaFeatureSpace:
  # aaFeatureSpace <- readRDS("data/aaFeatureSpace.3.0.Rds")
  #

}

# == MILESTONE: Amino Acid Feature Space =======================================

# At this point
#   - we have taken a large set of biophysical and statistical observations
#     of amino acids;
#   - we have scaled them to be numerically comparable among each other;
#   - we have selected catgories of observations that we believe are independent
#     of the genetic code itself, and thus can be used to evaluate the
#     genetic code;
#   - we have used PCA to remove any sampling bias and correlations between
#     the various sets of observations;
#   - we have constructed a "feature space" constructed from the "principal
#     components" of the dataset, which represents all of the information
#     contained in the original observations.


# =    12  AMINO ACID SIMILARITY  ==============================================

# We are now ready to define:
#
#    The similarity between two amino acids is the Euclidian distance between
#    their locations in an amino acid feature space.
#
# This maps the diverse properties of amino acids to a single number.


# Based on this, I have supplied a function - aaSim() -
# which was defined in the initialization block,
# which computes the actual similarities.

#    DEFINE THE FUNCTION HERE ...



if (FALSE) {

  aaSim("Q", "Q") #  0         Identical amino acids: distance is zero.
  aaSim("F", "I") #  0.8016213 Similar amino acids: distance is small.
  aaSim("Q", "F") #  4.763314  Dissimilar amino acids: distance is large.
  aaSim("F", "Q") #  4.763314  Distance between x and y is the same as
  #              the distance between y and x.
  aaSim("Q", "*") #  8.727359  Distance between any amino acid and a stop
  #              codon is large.
  aaSim("*", "*") #  0:        Two stop codons: distance is zero.

}




# Let's do a quick check of which amino acid is the most distinct, and which
# one is the most "plain vanilla" among the twenty. We develop this from
# pseudocode:

# Prompt:
#    Please explain "pseudocode".
# --

#> For each amino acid
#>   For each of the 19 other amino acids
#>     Compute the distance of the amino acid pair
#>   Compute the mean

if (FALSE) {

  # We can use the AADAT dataset to get a vector of one-letter codes.
  myA <- sort(unique(AADAT$A))[-1]
  mySims <- numeric(20)  # make a 20-elelment vector of numbers
  names(mySims) <- myA


  #> For each amino acid
  for (thisA in myA) {
    for (otherA in myA) {
      #>   For each of the 19 other amino acids
      if (thisA == otherA) {
        next
      }
      #>     Sum the distance of the amino acid pair
      mySims[thisA] <- mySims[thisA] + aaSim(thisA, otherA)
    }
    #>   Compute the mean
    mySims[thisA] <- mySims[thisA] / 19
  }

  (mySims <- sort(mySims) ) # sort the vector and show the results

  # The most "distinct" amino acid is "R" (arginine), the most average amino
  # acid is "T" (threonine).

  # Plot this:

  barplot(mySims,
          main = "Amino Acid Similarity",
          ylab = "Mean pairwise distances in feature space (AU)",
          ylim = c(0.0, 3.5),
          col = AACOLS[names(mySims)],
          cex.names = 0.5)

}
# (To follow: interpretation of principal components ...)


# ==============================================================================
#
#    I N I T I A L I Z T I O N   B L O C K   B E L O W
#
# ==============================================================================

# The initialization block starts here. You do not need to execute this
# part of the script again if you are working with it interactively.


# =    13  INITIALIZATION  =====================================================


# ==   13.1  Packages  =========================================================

# My preferred idiom to install packages is:
if (! requireNamespace("readxl", quietly=TRUE)) {
  install.packages("readxl")
}
# ... this installs only if the package has not been installed previously.

# More packages ...
if (! requireNamespace("seqinr", quietly=TRUE)) {
  install.packages("seqinr")
}

# ==   13.2  Data  =============================================================

# Amino acid data
# https://docs.google.com/spreadsheets/d/1R50YeqoplcA-IZ8sFNW9J7JgCzc2WGjTC9-mjPCZ5gY
AADAT <- read.csv("dat/SGC.csv")
rownames(AADAT) <- AADAT$Codon

A <- AADAT$A        # The one-letter amino acid symbols
Aaa <- AADAT$Aaa    # The three-letter amino acid symbols

# Load the AAINDEX data set from the seqinr package
data(aaindex, package = "seqinr")   # creates the "list" aaindex in the current
# workspace

# Download two Excel tables into the dat/ folder (if those don't exist yet) ...

if (! file.exists("dat/Breimann_2024_Supplementary_Table_1.xlsx")) {
  ghURL  <- paste0("https://raw.githubusercontent.com/hyginn/CSB195/main/dat/",
                   "Breimann_2024_Supplementary_Table_1.xlsx")
  myPath <- paste0("dat/",
                   "Breimann_2024_Supplementary_Table_1.xlsx")
  download.file(url = ghURL, destfile = myPath, mode = "wb")
}

if (! file.exists("dat/Breimann_2024_Supplementary_Table_3.xlsx")) {
  ghURL  <- paste0("https://raw.githubusercontent.com/hyginn/CSB195/main/dat/",
                   "Breimann_2024_Supplementary_Table_3.xlsx")
  myPath <- paste0("dat/",
                   "Breimann_2024_Supplementary_Table_3.xlsx")
  download.file(url = ghURL, destfile = myPath, mode = "wb")
}

# Download a precomputed amino acid feature spcae file
if (! file.exists("dat/aaFeatureSpace.3.0.Rds")) {
  ghURL  <- paste0("https://raw.githubusercontent.com/hyginn/CSB195/main/dat/",
                   "aaFeatureSpace.3.0.Rds")
  myPath <- paste0("dat/",
                   "aaFeatureSpace.3.0.Rds")
  download.file(url = ghURL, destfile = myPath, mode = "wb")
}

# Define a consistent set of amino acid colours
AACOLS <- c(G = "#B9C2CD", P = "#D4CF82", C = "#F1DD38", A = "#D2DF40",
            V = "#B4E149", I = "#96E351", L = "#78E65A", M = "#6DCB6E",
            F = "#63B182", W = "#599797", Y = "#4F7CAB", H = "#4562BF",
            R = "#3B48D4", K = "#5F6BD8", Q = "#838FDC", N = "#A8B3E0",
            T = "#AA90BA", S = "#AD6D95", D = "#B04B70", E = "#B3294B")

# ==   13.3  Functions  ========================================================

Aaa2A <- function(txt) {
  # bare-bones - no checking done
  AAmap <- c(Ala = "A", Arg = "R", Asn = "N",
             Asp = "D", Cys = "C", Gln = "Q",
             Glu = "E", Gly = "G", His = "H",
             Ile = "I", Leu = "L", Lys = "K",
             Met = "M", Phe = "F", Pro = "P",
             Ser = "S", Thr = "T", Trp = "W",
             Tyr = "Y", Val = "V")
  return(AAmap[txt])
}

if (FALSE) {
  # Test the function
  Aaa2A(names(aaindex[["EISD840101"]]$I))
  Aaa2A(c("Trp", "His", "Tyr", "?"))
}

grepAAindex <- function(key, idx = aaindex, el = "D", val = "txt", ...) {
  # Search for strings in an aaindex element
  # key:   a regular expression pattern
  # idx:   the list to search in. Default to aaindex (which should exist in the
  #        workspace).
  # el:    which element to search. Default "D" (Definition), also useful might
  #        be "T" (Titles).
  # val:   format of return value: "txt" (default) returns the element contnets,
  #        "idx" return the numeric indices of hits.
  # ...:   you can pass additional parameters into the grep() call.
  # value: a character vector, length 0 if nothing found.

  dat <- character(length(idx))        # set up an empty vector
  for (i in seq_along(idx)) {          # iterate over aaindex ...
    dat[i] <- idx[[i]][[el]]           # ... load element contents
    names(dat)[i] <- idx[[i]][["H"]]   # ... and name for the result vector
  }

  sel <- grep(key, dat, ...)           # a vector of matching index numbers

  if (val == "txt") {                  # "txt" ...
    return(dat[sel]) }                 # ... subset the result vector and return
  else if (val == "idx") {             # "idx" ...
    return(sel) }                      # ... return only the indices
  else {
    stop(sprintf("Unknown option: val = \"%s\"", val))
  }
}

if (FALSE) {
  # Examples of usage:
  grepAAindex("[Hh]ydroph")
  grepAAindex("flex", el = "T", val = "idx")
  aaindex[[8]]
}

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

  SPACEFILE <- "dat/aaFeatureSpace.3.0.Rds"

  # Stop codon distance:
  # The distance of an amino acid to a stop codon is STOPDIST times the
  # maximum distance in the distance matrix.
  STOPDIST <- 1.5

  # Recreate a 12-dimensional feature space of amino acids
  # (see aminoAcidSimilarity.R for details). Rownames should be single
  # letter symbols.
  AASPACE <- readRDS(SPACEFILE)

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
  aaSim("*", "*")   # same

  # Sanity check: Find the smallest and largest pairwise distance
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
  for (i in 1:19) {      # real amino acids only, no stop-codons
    for (j in i:20) {
      if (x[i, j] == dMin) {
        cat(sprintf("min: (%s, %s) = %.3f\n", aa[i], aa[j], x[i, j]))
      } else if (x[i, j] == dMax) {
        cat(sprintf("max: (%s, %s) = %.3f\n", aa[i], aa[j], x[i, j]))
      }
    }
  }
  # This may vary with the distance matrix details, but will probably turn out
  # to make (I, L) the minimum distance pair, (I, K) the maximum distance pair.

  x <- as.vector(x[1:19, 1:19])  # collapse the non-stop elements into a vector
  x <- x[x > 0]                  # remove all zero values
  summary(x)                     # print the summary

}  # end if (FALSE) ...

# In source() mode: announce that you arrived at the end successfully.
#



cat("Done. successfully initialized.\n")

# [END]
