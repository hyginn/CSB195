# tocID <- "courseScripts/04-aminoAcidsimilarity.R"
#
# Purpose: Develop a function to compute a pairwise similarity score for
#          amino acids. The score is based on the distance between two amino
#          acids in a feature spce derived from AAindex indices via their
#          AAontology categories.
#
# Version: 2.0
# Date:    2024-09
# Author:  boris.steipe@utoronto.ca; ChatGPT-4 and 4o
#
# Versions:
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
#TOC>   1        INITIALIZATIONS                         81
#TOC>   1.1        Parameters                            83
#TOC>   1.2        Packages                              89
#TOC>   2        PREPROCESS THE DATA                     96
#TOC>   2.1        Read Source Data                      98
#TOC>   2.2        Scale Raw Values                     112
#TOC>   2.3        Merge Into Single Data Frame         136
#TOC>   2.4        Validate                             166
#TOC>   2.4.1          Validate Scaling                 167
#TOC>   2.4.2          Check outliers                   171
#TOC>   2.4.3          Validate merging                 175
#TOC>   3        SELECT SCALES                          184
#TOC>   4        1                                      186
#TOC>   5        1                                      187
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
# biological context. That might depend on the code (among other factors).

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


# =    4  PCA ==================================================================

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


#       Apply Principal Component Analysis (PCA) to the selected data, to
#       compute components that have reduced dimensions, are orthogonal, and
#       are ordered by importance (contribution to variance).
idxPCA <- prcomp(idxScaled, scale. = FALSE)
#       PYP: Identify the first principal components that capture 95% of
#         variance as a basis for computing a similarity measure.
cumulativeProportions <- cumsum(idxPCA$sdev^2) / sum(idxPCA$sdev^2)
barplot(cumulativeProportions,
        main="Cumulative Proportion of Variance Explained",
        xlab="Principal Component",
        ylab="Cumulative Proportion",
        ylim=c(0, 1),
        border="blue",
        col="lightblue")
abline(h=0.95, col="red", lty=2)  # A horizontal line to indicate the 95% mark

numPCsToRetain <- which(cumulativeProportions >= 0.95)[1]  # 14

#       JOS: Extract these components into a suitable data structure.
aaFeatureSpace <- idxPCA$rotation[, 1:numPCsToRetain]
rownames(aaFeatureSpace) <- rownames(idxPCA$rotation)
colnames(aaFeatureSpace) <- paste0("PC", 1:numPCsToRetain)
str(aaFeatureSpace)  # inspect

# Let's look at the numerical distribution of values.
for (i in 1:14) {
  print(mean(aaFeatureSpace[, i]))
}
# All the means are "zero" (or at least very, very, very close to zero).

for (i in 1:14) {
  print(sd(aaFeatureSpace[, i]))
}
# All the standard deviations are 0.229. But that's actually not good. If the
# values of the higher PCs are distributed the same as those of the first or
# second PCs, then we give them an inappropriately large contribution
# to the distance calculation. We should re-scale the columns according
# to each PCs contribution to the total variance.

# We can get the variances from idxPCA$sdev - the variances are the squares
# of the standard deviations:

idxPCA$sdev[1:14]    # The standard deviations
idxPCA$sdev[1:14]^2  # The variances

# To scale the PCs that make up the dimensions of our feature space, we can
# simply multiply them by the variances. (Remember, the standard deviations of
# the PCs are already all the same. If the standard deviations would have been
# different, we would have divided by the standard deviation first.)

# Here we go:
for (i in 1:14) {
  aaFeatureSpace[, i] <- aaFeatureSpace[, i] * ( idxPCA$sdev[i]^2 )
}

# Save the resulting feature space so we don't have to re-compute it when
# we need it. saveRDS() saves a single R object in a compressed format. It is
# safer to use than save() (see below).
saveRDS(aaFeatureSpace, file = "data/aaFeatureSpace.2.0.Rds")

# saveRDS() / readRDS() save and _recreate_ single compressed R objects.
# You can assign them  to a new variable (or the same variable name).
# This is in contrast to save() / load(), which _restores_ an R object to its
# original name, and overwrites existing objects if you are not careful.
#
# Execute this to recreate aaFeatureSpace:
# aaFeatureSpace <- readRDS("data/aaFeatureSpace.2.0.Rds")
#

#       VAU: Validate PCA Components
#         DUJ: Validation - Compute correlations with existing indices to
#           support biological interpretation of the components.
# Number of original indices and number of retained PCs
numIndices <- nrow(idxScaled)
numPCs <- numPCsToRetain

# Initialize an empty matrix to store the correlations
corMatrix <- matrix(0, nrow=numIndices, ncol=numPCs)
rownames(corMatrix) <- rownames(idxScaled)
colnames(corMatrix) <- colnames(aaFeatureSpace)

# Loop to compute correlations
for (i in 1:numIndices) {
  for (j in 1:numPCs) {
    corMatrix[i, j] <- cor(idxScaled[i, ], aaFeatureSpace[, j])
  }
}  # The matrix 'corMatrix' now contains the desired correlations

# Find a good correlation of a PC with the index ...
# For the first PC, get the index with the maximum absolute correlation
maxIndexForPC1 <- which.max(abs(corMatrix[, 1]))  # column 1

# Retrieve the corresponding row name (index name) and correlation value
bestIndexName <- rownames(corMatrix)[maxIndexForPC1]
bestCorrelation <- corMatrix[maxIndexForPC1, 1]

cat("The index most correlated with PC1 is:", bestIndexName,
    "with a correlation of:", bestCorrelation, "\n")

# Verify with a scatterplot
# Extract PC1 values and NISK860101 index values
x <- idxScaled["NISK860101", ]
y <- aaFeatureSpace[, "PC1"]

# Create a scatterplot
plot(x, y,
     main="Scatterplot of PC1 vs NISK860101",
     xlab="NISK860101 Index Values",
     ylab="PC1 Values",
     pch=19,
     col="skyblue")

# Add a linear regression line for reference
abline(lm(y ~ x), col="#ff33aa")


#         DYS: Validation - Examine a scatterplot of the amino acids
#           PC1 against PC2
plot(aaFeatureSpace[, "PC1"], aaFeatureSpace[, "PC2"],
     main="Scatterplot of PC1 vs PC2",
     xlab="PC1", ylab="PC2", pch=19, col="blue")

# Adding amino acid labels to the points for clarity
text(aaFeatureSpace[, "PC1"], aaFeatureSpace[, "PC2"],
     labels=rownames(aaFeatureSpace), pos=3, cex=0.8)


#   YTI: Define "similarity" as the Euclidean distance between the
#     two input amino acids in the PCA based feature space stored in JOS: .
#     (Such Euclidean distances in feature space are intuitive to interpret.)
#     (Alternatives to consider are:
#       - Average of one-dimensional differences.
#       - Weighted average of one-dimensional differences if some indices
#         carry more importance.
#       - Sum of one-dimensional differences
#       - Clustering values and assigning categories to pairwse differences
#       - Cosine similarity
#       - Any other statistical or machine learning method that can be
#         applied to vector differences.
#       )


# Note: the function below - myAAsim() - which was originally developed in
# conversation with the AI, has been superseded by aaSim() (defined in
# ./R/aaSim.R). I have kept it in this file to keep the development process
# complete.

myAAsim <- function(a1, a2, space = aaFeatureSpace) {
  # a1, a2: amino acid one-letter symbols

  # Extract the feature vectors for the two amino acids
  vector1 <- space[a1, ]
  vector2 <- space[a2, ]

  # Compute the Euclidean distance
  distance <- sqrt(sum((vector1 - vector2)^2))

  return(distance)
}

# Try this:
myAAsim("M", "W")  # should be "similar"
myAAsim("I", "D")  # should be "dissimilar"
myAAsim("Q", "Q")  # should be identical

# Examine heatmap

aminoAcids <- rownames(aaFeatureSpace)
distanceMatrix <- matrix(0, nrow=length(aminoAcids), ncol=length(aminoAcids))
rownames(distanceMatrix) <- aminoAcids
colnames(distanceMatrix) <- aminoAcids

for(i in 1:length(aminoAcids)) {
  for(j in 1:length(aminoAcids)) {
    distanceMatrix[i, j] <- myAAsim(aminoAcids[i], aminoAcids[j])
  }
}

heatmap(distanceMatrix, main="Amino Acid Pairwise Distances",
        Colv=NA, Rowv=NA, scale="none", col=rev(heat.colors(256)),
        xlab="Amino Acid", ylab="Amino Acid")

# Examine ALL the heatmaps
# Function to compute distance matrix for a given feature space subset
computeDistanceMatrix <- function(spaceSubset) {
  aminoAcids <- rownames(spaceSubset)
  distanceMatrix <- matrix(0, nrow=length(aminoAcids), ncol=length(aminoAcids))
  rownames(distanceMatrix) <- aminoAcids
  colnames(distanceMatrix) <- aminoAcids

  for(i in 1:length(aminoAcids)) {
    for(j in 1:length(aminoAcids)) {
      distanceMatrix[i, j] <- myAAsim(aminoAcids[i], aminoAcids[j], spaceSubset)
    }
  }

  return(distanceMatrix)
}

# Create a series of heatmaps

for(i in 1:14) {
  distanceMatrixSubset <- computeDistanceMatrix(aaFeatureSpace[, 1:i, drop=FALSE])
  heatmap(distanceMatrixSubset, main=paste("Distances using PC1 to PC", i),
          Colv=NA, Rowv=NA, scale="none", col=rev(heat.colors(256)),
          xlab="Amino Acid", ylab="Amino Acid")
}

# Do the same thing with dendrogram heatmaps

for(i in 1:14) {
  distanceMatrixSubset <- computeDistanceMatrix(aaFeatureSpace[, 1:i, drop=FALSE])

  # Using heatmap to plot the distance matrix with dendrogram clustering
  heatmap(distanceMatrixSubset, main=paste("Distances using PC1 to PC", i),
          Colv=TRUE, Rowv=TRUE, scale="none", col=rev(heat.colors(256)),
          xlab="Amino Acid", ylab="Amino Acid")
}

# run t-SNE for more exploration:
install.packages("Rtsne")

# Load the Rtsne package
library(Rtsne)

# Apply t-SNE on aaFeatureSpace
set.seed(123)  # Set seed for reproducibility
tsneResults <- Rtsne::Rtsne(aaFeatureSpace,
                     dims=2,
                     perplexity=5,            # was: 10 (default)
                     check_duplicates=FALSE)

# Plot the t-SNE results
plot(tsneResults$Y, t='n', main="t-SNE of aaFeatureSpace", xlab="", ylab="")
text(tsneResults$Y, labels=rownames(aaFeatureSpace), cex=0.8, col="blue")


# === CONTINUE AI CO-DEVELOPED PSEUDOCODE =================

#     DMD: Validate Euclidean Distance computation
#       FKB: Validation - Inspect a distance matrix, ensure metricity:
#         - Non-negativity and identity of indiscernibles:
#           d(a,a) == 0 and d(a,b) >= 0
#         - Symmetry: d(A, B) == d(B, A)
#         - Triangle inequality: d(a,b) + d(b,c) ≥ d(a,c)

# PROCESS:

#   YLC: Apply the similarity function
#     YUC: Accept and validate input
#       Ensure the input is a vector with exactly two amino acid symbols that
#       correspond to the labels used for the feature space.
#       GLI: Represent the two amino acids as a vector in the defined
#         feature space
#       DZF: Compute and RETURN the Euclidean distance between those vectors.
#         The result is the similarity. Lower values are more similar.
#         A value of 0.0 means: they are identical.

# OUTPUT: similarity of an amino acid pair as a single scalar value

# VALIDATION:

#   LIL: Validation - Construct Valid Test Cases for the full workflow
#     UTE: Utilize known similarities and dissimilarities between amino acid
#       pairs for eyeball testing
#     ESX: Construct synthetic data defining artificial amino acids and
#       indices, that give known correct results.
#       XQI: Run tests of the workflow and confirm that the known correct
#         results are obtained

#   WRF: Validation - Analyze Result Distribution
#     WKS: Visualize the distribution of similarity scores using histograms
#     KZL: Obtain statistical summaries of similarity scores to gauge
#       distribution properties

#   CPZ: Validation - Ensure Biological/Chemical Relevance
#     AHQ: Check results for logical consistency with biological/chemical
#       intuitions
#     WVG: Compare results with known sequence alignment matrices like BLOSUM

#   NAJ: Validation - Visualize and Intuitively Assess Results
#     AKK: Visualize PCA-transformed data, exploring clusters or patterns
#       among amino acids
#     OLX: Create a heatmap of the distance matrix to visually explore amino
#       acid pair similarities



# ====  TESTS AND VALIDATION ===================================================

if (FALSE) {



}



# [END]
