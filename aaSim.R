# tocID <- "aaSim.R"
#
# Purpose: define a function to compute a pairwise similarity score for
#          amino acids.
#
# Version: 0.3
# Date:    2023-10
# Author:  boris.steipe@utoronto.ca; CSB195 2023 Class; ChatGPT-4
#
# Versions:
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
# Input:   a vector with two amino acid one-letter codes
# Output:  a scalar distance in a feature space between the amino acids
#
# ToDo:   Code cleanup: make sure the file is safe-to source and defines
#         the aaSim() function upon source()'ing.
#
# Notes:  This script will contain BOTH the code for the function development
#         and the actual function definition. Thus all relevant information
#         is sanely kept in one single file.
#
# ==============================================================================


# ====  PARAMETERS  ============================================================

URL <- paste(c(
  "https://docs.google.com/spreadsheets/d/",       # Prefix
  "11WtjhJYufNqkEjbHkHvtGlqLocXkXUFG-00W1BGPUG4",  # File ID
  "/edit?usp=sharing"),                            # Suffix
  collapse = "")                                   # one string, no blanks

sheet <- "Curation table"                          # Sheet name

A <- AADAT$A        # The one-letter amino acid codes
Aaa <- AADAT$Aaa    # The three-letter amino acid codes

# ====  PACKAGES  ==============================================================


# ====  FUNCTIONS  =============================================================

A2Aaa <- function(a) {
  # Helper function to convert one-letter codes to three letter codes
  if (a %in% A) {   # if the provided amino acid code is a one-letter code
    return(Aaa[which(A == a)])
  } else {          # it is already three-letter. Return it as is.
    return(a)
  }
}
# NOTE: This is how ChatGPT wrote it, I myself would write this more
#   defensively (and more efficiently). But it will do ...



# ====  PROCESS  ===============================================================

# OBJECTIVE: To compute pairwise amino acid similarity from
#   distance in a feature space

# INPUTS:
#   URL, sheet - a Google sheet of amino acid property indices in rows,
#   all 20 amino acids in columns

# PREPARE: Construct a similarity function

#   BEZ: Fetch index data from a curated Google sheet
#     QCR: Use readGsheet(URL, sheet) to fetch the data from Google Sheets.
curatedTable <- readGsheet(URL, sheet)

#     CCH: Validate the fetched data
#       GYE: Ensure data structure is a data frame.
if(!is.data.frame(curatedTable)){
  stop("Fetched data is not a data frame.")
}

#       SXI: Confirm presence of required columns (e.g., amino acid codes).
if(!all(Aaa %in% colnames(curatedTable))){
  stop("Not all required amino acid columns (three-letter codes) are present.")
}

requiredColumns <- c("ID", "Cur.derived")
if(!all(requiredColumns %in% colnames(curatedTable))){
  stop("ID and/or Cur.derived column is missing.")
}
if(!is.numeric(curatedTable$Cur.derived)){
  stop("Cur.derived is not numeric.")
}
if(length(unique(curatedTable$ID)) != nrow(curatedTable)){
  stop("Duplicate ID entries found.")
}
if(any(is.na(curatedTable$ID))){
  stop("The ID column contains NA values.")
}

#       PEM: Check for missing or NA values and remove indices that
#         contain them. (Note: if a lot are missing we might consider
#         imputation of missing data.)
# complete.cases() flags rows in which the requested columns have NA
sel <- complete.cases(curatedTable[ ,c(Aaa, "Cur.derived")])
curatedTable <- curatedTable[sel, ]

#     HCG: Filter the data to include or exclude specific rows
#       Ignore all rows in which the curators have assigned a > 0 probability
#       that the index is influenced by the structure of the genetic code.
sel <- curatedTable$Cur.derived == 0
curatedTable <- curatedTable[sel, ]

#   YLV: Construct a feature space
#     LND: Extract amino acid data only, but keep the ID information
#       as rownames
idxData <- curatedTable[ , Aaa]
rownames(idxData) <- curatedTable$ID
colnames(idxData) <- Aaa

#     EFL: Pre-process the index data
#       JJT: scale() the indices to have mean of 0 and standard deviation of 1
#         in each row. (Normalizing between 0 and 1 could be an alternative.)
idxScaled <- t(scale(t(idxData))) # Transpose, scale, and then transpose back

#       HPR: Identify outlier values to manually review them.
#         We consider outliers to deviate more then 3 sigma from the mean.
#         Since the data is already scaled, we can just look for values
#         that are smaller than -3 or larger than 3.
outliers <- abs(idxScaled) > 3
outlierIndices <- which(abs(idxScaled) > 3, arr.ind = TRUE)
outlierRows <- rownames(idxScaled)[outlierIndices[, 1]]
outlierCols <- colnames(idxScaled)[outlierIndices[, 2]]

# There are indeed typos and anomalies in the data. Let's inspect them all:
for (i in 1:length(outlierRows)) {
  cat("ID:", outlierRows[i], " - Amino Acid:", outlierCols[i], "\n")
}

# ID: RACS820104  - Amino Acid: Cys - Plausible
# ID: FINA910101  - Amino Acid: Asp - Plausible
#    ID: MONM990201  - Amino Acid: Asp - Looks like typo: 15 corrected to 1.5
# ID: AURR980106  - Amino Acid: Glu - Plausible
#    ID: WERD780103  - Amino Acid: Phe - Can't tell, but the Met value is NA in the
#       original paper ... IDK why the aaindex has 0.11. Changed to NA.
# ID: NAKH900113  - Amino Acid: Gly - Plausible
# ID: ISOY800108  - Amino Acid: Gly - Plausible
# ID: FAUJ880105  - Amino Acid: Gly - Plausible
# ID: RACS820109  - Amino Acid: Gly - Plausible
# ID: MAXF760104  - Amino Acid: Gly - Plausible
# ID: KIMC930101  - Amino Acid: Gly - Plausible
#
# ID: CHAM810101  - Amino Acid: Gly - Plausible
# ID: RICJ880115  - Amino Acid: Gly - Plausible
# ID: KOEP990101  - Amino Acid: Gly - Plausible
# ID: NAKH900112  - Amino Acid: Leu - Plausible
# ID: MAXF760103  - Amino Acid: Asn - Plausible
# ID: WERD780102  - Amino Acid: Asn - Plausible
# ID: QIAN880108  - Amino Acid: Pro - Plausible
# ID: QIAN880134  - Amino Acid: Pro - Plausible
# ID: ONEK900101  - Amino Acid: Pro - Plausible
# ID: MUNV940101  - Amino Acid: Pro - Plausible
#
# ID: CHOP780213  - Amino Acid: Pro - Plausible
#    ID: FINA910102  - Amino Acid: Pro - anomalous, but that's in the paper
# ID: RACS820114  - Amino Acid: Pro - Plausible
# ID: TANS770104  - Amino Acid: Pro - Plausible
# ID: BLAM930101  - Amino Acid: Pro - Plausible
# ID: ROBB760103  - Amino Acid: Pro - Plausible
# ID: ONEK900102  - Amino Acid: Pro - Plausible
# ID: GEOR030109  - Amino Acid: Pro - Plausible
# ID: MUNV940105  - Amino Acid: Pro - Plausible
#    ID: FAUJ880113  - Amino Acid: Pro - wrong: value was 0 but should be NA
# ID: AURR980119  - Amino Acid: Pro - Plausible
#
#    ID: BUNA790101  - Amino Acid: Pro - wrong: value was 0 but should be NA
#    ID: YUTK870103  - Amino Acid: Arg - wrong: value was 0 but should be NA
# ID: EISD860102  - Amino Acid: Arg - Plausible
# ID: JOND750102  - Amino Acid: Arg - Plausible
# ID: FAUJ880109  - Amino Acid: Arg - Plausible
#    ID: YUTK870102  - Amino Acid: Arg - wrong: value was 0 but should be NA
# ID: JACR890101  - Amino Acid: Arg - Plausible
# ID: OOBM850102  - Amino Acid: Trp - Plausible
# ID: WEBA780101  - Amino Acid: Trp - Plausible
# ID: GARJ730101  - Amino Acid: Trp - Plausible



#     JKQ: Apply Principal Component Analysis (PCA) to the scaled data, to
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
# we need it:
save(aaFeatureSpace, file = "data/aaFeatureSpace.2.0.RData")

# load("data/aaFeatureSpace.2.0.RData")  # This would reload the save()'d object

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

aaSim <- function(a1, a2, space = aaFeatureSpace) {
  # Convert to three-letter codes if necessary
  a1 <- A2Aaa(a1)
  a2 <- A2Aaa(a2)

  # Extract the feature vectors for the two amino acids
  vector1 <- space[a1, ]
  vector2 <- space[a2, ]

  # Compute the Euclidean distance
  distance <- sqrt(sum((vector1 - vector2)^2))

  return(distance)
}

# Try this:
aaSim("M", "W")  # should be "similar"
aaSim("Met", "Trp")  # should be the same
aaSim("I", "D")  # should be "dissimilar"
aaSim("Q", "Q")  # should be identical

# Examine heatmap

aminoAcids <- rownames(aaFeatureSpace)
distanceMatrix <- matrix(0, nrow=length(aminoAcids), ncol=length(aminoAcids))
rownames(distanceMatrix) <- aminoAcids
colnames(distanceMatrix) <- aminoAcids

for(i in 1:length(aminoAcids)) {
  for(j in 1:length(aminoAcids)) {
    distanceMatrix[i, j] <- aaSim(aminoAcids[i], aminoAcids[j])
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
      distanceMatrix[i, j] <- aaSim(aminoAcids[i], aminoAcids[j], spaceSubset)
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
tsneResults <- Rtsne(aaFeatureSpace,
                     dims=2,
                     perplexity=5,            # was: 10 (default)
                     check_duplicates=FALSE)

# Plot the t-SNE results
plot(tsneResults$Y, t='n', main="t-SNE of aaFeatureSpace", xlab="", ylab="")
text(tsneResults$Y, labels=rownames(aaFeatureSpace), cex=0.8, col="blue")


#     DMD: Validate Euclidean Distance computation
#       FKB: Validation - Inspect a distance matrix, ensure metricity:
#         - Non-negativity and identity of indiscernibles:
#           d(a,a) == 0 and d(a,b) >= 0
#         - Symmetry: d(A, B) == d(B, A)
#         - Triangle inequality: d(a,b) + d(b,c) ≥ d(a,c)

# PROCESS:

#   YLC: Apply the similarity function
#     YUC: Accept and validate input
#       Ensure the input is a vector with exactly two amino acid codes that
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

# ====  APPENDIX: ORIGINAL PSEUDOCODE ==========================================

# OBJECTIVE: To compute pairwise amino acid similarity from
#   distance in a feature space

# INPUTS:
#   URL, sheet - a Google sheet of amino acid property indices in rows,
#   all 20 amino acids in columns

# PREPARE: Construct a similarity function

#   BEZ: Fetch index data from a curated Google sheet
#     QCR: Use readGsheet(URL, sheet) to fetch the data from Google Sheets.
#     CCH: Validate the fetched data
#       GYE: Ensure data structure is a data frame.
#       SXI: Confirm presence of required columns (e.g., amino acid codes).
#       PEM: Check for missing or NA values and remove indices that
#         contain them. (Note: if a lot are missing we might consider
#         imputation of missing data.)
#     HCG: Filter the data to include or exclude rows.
#       ZQL: Determine criteria to exclude specific indices.
#       UOK: Filter out rows of the data frame based on these criteria.

#   YLV: Construct a feature space
#     EFL: Pre-process the index data
#       JJT: scale() the indices to have mean of 0 and standard deviation of 1
#         in each row. (Normalizing between 0 and 1 could be an alternative.)
#       HPR: Flag outlier values to manually review them.
#     JKQ: Apply Principal Component Analysis (PCA) to the scaled data, to
#       compute components that have reduced dimensions, are orthogonal, and
#       are ordered by importance (contribution to variance).
#       PYP: Identify the first principal components that capture 95% of
#         variance as a basis for computing a similarity measure.
#       JOS: Extract these components into a suitable data structure.
#       VAU: Validate PCA Components
#         DUJ: Validation - Compute correlations with existing indices to
#           support biological interpretation of the components.
#         YDB: Validation - Ensure chosen principal components capture
#           significant variance through a scree plot
#         ZQB: Validation - Use biplot to visualize the correlation between
#           original variables and principal components

#   YTI: Define "similarity" as the Euclidean distance between the
#     two input amino acids in the PCA based feature space stored in JOS: .
#     (Such Euclidean distances in feature space are intuitive to interpret.)
#     (Alternatives to consider are:
#       - Average of one-dimensional differences.
#       - Weighted average of one-dimensional differences if some indices
#         carry more importance.
#       - Sum of one-dimensional differences
#       - Clustering values and assigning categories to pairiwse differences
#       - Cosine similarity
#       - Any other statistical or machine learning method that can be
#         applied to vector differences.
#       )
#     DMD: Validate Euclidean Distance computation
#       FKB: Validation - Inspect a distance matrix, ensure metricity:
#         - Non-negativity and identity of indiscernibles:
#           d(a,a) == 0 and d(a,b) >= 0
#         - Symmetry: d(A, B) == d(B, A)
#         - Triangle inequality: d(a,b) + d(b,c) ≥ d(a,c)

# PROCESS:

#   YLC: Apply the similarity function
#     YUC: Accept and validate input
#       Ensure the input is a vector with exactly two amino acid codes that
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


# [END]
