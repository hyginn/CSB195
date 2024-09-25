# tocID <- "./courseScripts/tutorial_03.R"
#
# Purpose: A: aaindex properties
#          B: reading supplementary data for AAontology (Breimann et. al 2024)
#
# Version: 0.1
# Date:    2024-09
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#
# ToDo:
# Notes:
#
# ==============================================================================


# =    1  AAindex  =============================================================

# aaindex properties

grepAAindex("hydro")

# For example ...
# $`68`
# [1] "Consensus normalized hydrophobicity scale (Eisenberg, 1984)"

myHyd <- aaindex[[68]]$I
barplot(myHyd,
        col = AACOLS[A2Aaa(names(myHyd))],
        names.arg = names(myHyd),
        cex.names = 0.4)

grepAAindex("volu")

# For example ...
# $`150`
# [1] "Side chain volume (Krigbaum-Komoriya, 1979)"

myVol <- aaindex[[150]]$I
barplot(myVol,
        col = AACOLS[A2Aaa(names(myVol))],
        names.arg = names(myVol),
        cex.names = 0.4)

# These scales are different - one way to compare them is to plot
# one against the other in a scatterplot.

plot(myHyd, myVol)

plot(myHyd, myVol,
     type = "n",
     main = "Amino Acid Hydrophobicity vs. Volume",
     xlab = "Hydrophobicity scale (Eisenberg 1984)",
     ylab = "Side chain volume (Krigbaum 1979)")
text(myHyd, myVol,
     labels = names(myHyd),
     col = AACOLS[A2Aaa(names(myVol))],
     cex = 1.5)

# There are many more properties ...



# =    2  AAontology  ==========================================================

# == Prompt:
#    What is an "ontology"?
# --
#    How can I read an Excel file (.xlsx) in R?
# --

install.packages("readxl")


# My preferred idiom here is:
if (! requireNamespace("readxl", quietly=TRUE)) {
  install.packages("readxl")
}

# Do we have the data files?
list.files("data/", pattern = "xlsx$")

# Read a sheet ...
myScales <- readxl::read_excel("data/Breimann_2024_Supplementary_Table_1.xlsx",
                               sheet = "Normalized")




# =    6  VALIDATION  ==========================================================

#
if (FALSE) {
  # This block is never executed. Place all tests and experiments here that
  # should not be run when this file is source()'ed.





}

# [END]
