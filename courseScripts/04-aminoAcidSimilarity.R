# tocID <- "courseScripts/04-aminoAcidsimilarity.R"
#
# Purpose: Develop a function to compute a pairwise similarity score for
#          amino acids. The score is based on the distance between two amino
#          acids in a feature spce derived from AAindex indices via their
#          AAontology categories.
#
# Version: 2.1
# Date:    2024-10
# Author:  boris.steipe@utoronto.ca; ChatGPT-4 and 4o
#
# Versions:
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
#                                                                              #
#   THIS SCRIPT IS NOT "SAFE TO SOURCE". SOURCING IT WILL HAVE SIDE EFFECTS.   #
#                                                                              #
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                 Line
#TOC> -----------------------------------------------------
#TOC>   1        INITIALIZATIONS                         83
#TOC>   1.1        Parameters                            85
#TOC>   1.2        Packages                              91
#TOC>   2        PREPROCESS THE DATA                     98
#TOC>   2.1        Read Source Data                     100
#TOC>   2.2        Scale Raw Values                     114
#TOC>   2.3        Merge Into Single Data Frame         138
#TOC>   2.4        Validate                             168
#TOC>   2.4.1          Validate Scaling                 169
#TOC>   2.4.2          Check outliers                   173
#TOC>   2.4.3          Validate merging                 177
#TOC>   3        SELECT SCALES                          184
#TOC>   4        PCA                                    282
#TOC>   5        AMINO ACID FEATURE SPACE               322
#TOC>   6        AMINO ACID SIMILARITY                  392
#TOC> 
#TOC> ==========================================================================


# =    1  INITIALIZATIONS  =====================================================

# ==   1.1  Parameters  ========================================================

A <- AADAT$A        # The one-letter amino acid symbols
Aaa <- AADAT$Aaa    # The three-letter amino acid symbols


# ==   1.2  Packages  ==========================================================

if (! requireNamespace("readxl", quietly=TRUE)) {
  install.packages("readxl")
}


# =    2  PREPROCESS THE DATA  =================================================

# ==   2.1  Read Source Data  ==================================================

# Scale values (indices) are in Breimann et al. (2024) Supplementary Table 1
SCALESHEET <- "Raw"  # Not "Normalized" here.
myScales <- readxl::read_excel("data/Breimann_2024_Supplementary_Table_1.xlsx",
                               sheet = SCALESHEET)
myScales <- as.data.frame(myScales)

# Scale categories and subcategories are in Supplementary Table 3
myCats   <- readxl::read_excel("data/Breimann_2024_Supplementary_Table_3.xlsx",
                               sheet = "Scales")
myCats   <- as.data.frame(myCats)


# ==   2.2  Scale Raw Values  ==================================================

# Since the dimensionality reduction method we will use later (PCA) is
# sensitive to the variance of the data, we will not use the min/max scaling
# between 0 and 1 that the authors used (and called "Normalization"), but a
# proper scaling that transforms the data to a mean of 0 and a
# standard deviation of 1.

x <- myScales   # make a temporary copy for validation

for (i in 2:ncol(myScales)) {
  myScales[ , i] <- scale(as.matrix(as.numeric(myScales[ , i])))
}

# Validate
idx <- 57                           # Pick a random dataset
mean(myScales[ , idx])              # Must be (nearly) zero
sd(myScales[ , idx])                # Must be (nearly) one
plot(myScales[ , idx], x[ , idx])   # Points exactly on the diagonal

# Clean up
rm(x)


# ==   2.3  Merge Into Single Data Frame  ======================================

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

# Check
id <- colnames(myScales)[57]       # Pick a random scale_id
all(myScales[ , id] == aaontology[id, 6:25])

# ==   2.4  Validate  ==========================================================
# ===   2.4.1  Validate Scaling            
#
# Confirm that all means are zero and all SDs are 1.

# ===   2.4.2  Check outliers              
#
# Identify any value outside the range (-3, 3)

# ===   2.4.3  Validate merging            
for (i in 2:ncol(myScales)) {
  id <- colnames(myScales)[i]
  stopifnot(all(myScales[ , id] == aaontology[id, 6:25]))
}


# =    3  SELECT SCALES  =======================================================

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


# =    4  PCA  =================================================================

# Prompt:
#    Please explain PCA for data analysis.
# --

# Create a selection vector for the data columns we want to include

sel <- logical(nrow(aaontology))
for (i in 1:nrow(aaontology)) {
  idx <- which(aaontology$subcategory[i] == mySubCats$subcat)
  if (mySubCats$inc[idx] == 1) {
    sel[i] <- TRUE
  } else {
    sel[i] <- FALSE
  }
}

# Validate: check all selected subcategories
# ...

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

# =    5  AMINO ACID FEATURE SPACE  ============================================

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
saveRDS(aaFeatureSpace, file = "data/aaFeatureSpace.3.0.Rds")

# saveRDS() / readRDS() save and _recreate_ single compressed R objects.
# You can assign them  to a new variable (or the same variable name).
# This is in contrast to save() / load(), which _restores_ an R object to its
# original name, and overwrites existing objects if you are not careful.
#
# Execute this to recreate aaFeatureSpace:
# aaFeatureSpace <- readRDS("data/aaFeatureSpace.3.0.Rds")
#


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


# =    6  AMINO ACID SIMILARITY  ===============================================

# We are now ready to define:
#
#    The similarity between two amino acids is the Euclidian distance between
#    their locations in an amino acid feature space.
#
# This maps the diverse properties of amino acids to a single number.


# Based on this, I have supplied a function - aaSim() -
# which computes the actual similarities.

aaSim("Q", "Q") #  0         Identical amino acids: distance is zero.
aaSim("F", "I") #  0.8016213 Similar amino acids: distance is small.
aaSim("Q", "F") #  4.763314  Dissimilar amino acids: distance is large.
aaSim("F", "Q") #  4.763314  Distance between x and y is the same as
                #              the distance between y and x.
aaSim("Q", "*") #  8.727359  Distance between any amino acid and a stop
                #              codon is large.
aaSim("*", "*") #  0:        Two stop codons: distance is zero.


# You can read the construction of the function in the file R/aaSim.R

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


# (To follow: interpretation of principal components ...)

# [END]
