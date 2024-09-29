# tocID <- "./courseScripts/03-AAontology.R"
#
# Purpose: A: Introducing the idea of structuring data through an ontology.
#          B: Reading supplementary data for AAontology (Breimann et. al 2024)
#          C: Combining data for ease of use
#          D: A plot example
#
# Version: 0.9
# Date:    2024-09-28
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   0.1    in tutorial class script
#
# ToDo:
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                   Line
#TOC> -------------------------------------------------------
#TOC>   1        PROBLEMS WITH THE AAINDEX                 40
#TOC>   2        READING PUBLISHED DATA                    63
#TOC>   3        MERGING DATA  FROM TWO SOURCES           135
#TOC>   3.1        Validations                            191
#TOC>   4        PARALLEL COORDINATES PLOT                224
#TOC>
#TOC> ==========================================================================


# Copy this file into your "myScripts" folder and edit your copy with notes,
# comments, and code.

# Run the command gAIinit(), open a session with your favorite generative
# AI assistant, and paste the prompt. More prompts are included below,
# ask the AI to clarify any additional terms and concepts. Some of the prompts
# are quite long, make sure you copy the entire prompt.


# =    1  PROBLEMS WITH THE AAINDEX  ===========================================

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
#    In that context, what is a category?
# --
#    What could category membership be based on?
# --


# =    2  READING PUBLISHED DATA  ==============================================


# Serendipitously, a few days ago a paper was published in the Journal of
# Molecular Biology in which the authors construct an ontology of
# categories and subcategories for the AAindex. You can access the paper here:
#
# https://www.sciencedirect.com/science/article/pii/S0022283624003267
#
# I have placed two of the supplementary data files into the course project
# data folder, and you can also inspect them here:
#
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

# First, we will read the data into R. It is now in the MicroSoft Excel
# spreadsheet format.

# == Prompt:
#
#    How can I read an Excel file (.xlsx) into an R dataframe?
# --

# My preferred idiom to install packages is:
if (! requireNamespace("readxl", quietly=TRUE)) {
  install.packages("readxl")
}
# ... this installs only if the package has not been installed previously.
# First, we check that the files exist in the expected location:

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


# =    3  MERGING DATA  FROM TWO SOURCES  ======================================


# I have produced a summary of the data to share it with the AI. Try this:

str(myScales, list.len = 5)
str(myCats)

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



# ==   3.1  Validations  =======================================================


# Of course, we now have to validate the result.
#
# As a first sanity check, do the following:
aaontology[nrow(aaontology), ]  # prints out the contents of the last row

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


# =    4  PARALLEL COORDINATES PLOT  ===========================================

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



# ====  DONE  ==================================================================

# If you've arrived here, congratulations - you did a quite respectable
# piece of work.
#



# [END]
