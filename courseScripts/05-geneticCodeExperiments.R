# tocID <- "courseScripts/05-geneticCodeExperiments.R"
#
# Purpose: Summarize rationale and process of performing experiments
#          on the quality of the Universal Genetic Code.
#
# Version: 0.1
# Date:    2024-10
# Author:  boris.steipe@utoronto.ca; ChatGPT-4 and 4o
#
# Versions:
#   0.1    New script with code contributions from
#
# Preconditions:
#    AADAT and GCdf are defined via .util.R
#
# ToDo:
#
# Notes:
#
# ==============================================================================
#                                                                              #
#   THIS SCRIPT IS NOT "SAFE TO SOURCE". SOURCING IT WILL HAVE SIDE EFFECTS.   #
#                                                                              #
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                         Line
#TOC> -------------------------------------------------------------
#TOC>   1        INITIALIZATIONS                                 48
#TOC>   1.1        Parameters                                    50
#TOC>   1.2        Packages                                      53
#TOC>   2        RATIONALE AND PROCESS                           60
#TOC>   3        DATA STRUCTURE FOR A GENETIC CODE               99
#TOC>   4        A RANDOM GENETIC CODE                          128
#TOC>   4.1        Constructing Randomized Codes                170
#TOC>   4.1.1          Redundancy in the Code                   178
#TOC>   4.1.2          The amino acid "alphabet"                225
#TOC>   4.1.3          Codes with the same redundancy           238
#TOC>   4.2        The randomized code generator                362
#TOC>   5        NEIGHBOURING CODONS                            372
#TOC>   6        EVALUATING THE STANDARD GENETIC CODE           385
#TOC>
#TOC> ==========================================================================


# =    1  INITIALIZATIONS  =====================================================

# ==   1.1  Parameters  ========================================================


# ==   1.2  Packages  ==========================================================

# Package loader for packages and data in the Bioconductor project.
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# =    2  RATIONALE AND PROCESS  ===============================================

# In our search to answer the question: "Is the genetic code the best it can
# be?", we observed that identical amino acids are encoded predominantly in
# adjacent positions of the genetic code, and we speculated that the same might
# be true for similar amino acids. If this conjecture holds up, we would
# conclude that the genetic code has been optimized through its evolution to
# mitigate the effects of point mutations, or errors in translating gene
# transcripts into proteins. This gives us (i) an idea of what "the best" could
# mean in this context, and (ii) an approach to make (computational)
# experimental observations.

# Roughly, we could evaluate the quality of any given code by summing over the
# result of all possible point mutations. Each point mutation would change a
# codon XXX into a codon YYY, and if we compute the similarity between the
# encoded amino acid x and amino acid y - via the function aaSim(x, y) - we have
# a measure for the severity of the mutation's consequences. We do this 64
# times, and we take the sum of all mutations to be the overall quality of the
# code.

# Going forward, for clarity we will refer to the genetic code we use as "sGC",
# shorthand for the "Standard Genetic Code".

# The quality measure we just proposed is just a single number for the sGC. But
# we can  compare this with alternative genetic codes that we create through
# some stochastic procedure. Are those alternative codes overwhelmingly more
# robust then sGC? Then we conclude that our conjecture was wrong and sGC is not
# the best for what we think it ought to do. Are those alternative codes usually
# about the same as sGC? Then we conclude that as far as mitigation againt point
# mutations is concerned, sGC performs on only just about a random level. But if
# sGC is generally or overwhelmingly better than our random codes, we can
# conclude that it is really the best (or nearly so), with respect to our
# assumption. And we may even be able to quantify how much better it is than if
# it would have just been constructed through random chance.

# To put this into practice, there are two programming tasks to solve:
# generating a random genetic code, and evaluating its quality.


# =    3  DATA STRUCTURE FOR A GENETIC CODE  ===================================

# There are many ways to store the genetic code as a data structure. An
# authoritative reference is the Genetic Codes repository at the NCBI:
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi In R programming,
# bioinformatics data, tools, and packages are made available for the very
# large, international Bioconductor project (https://www.bioconductor.org/).
# Let's retrieve the Bioconductor version of sGC (or NCBI code 1). Accessing
# Bioconductor packages goes through two steps: (a) loading a package called the
# Bioconductor manager, and (b) retrieving whatever we need. Step (a) was done
# in our initializations.

# The Genetic code table we most commonly use is part of the Biostrings package.
if (! requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")  # NOT install.packages(...) ...
}

?Biostrings::getGeneticCode         # See documentation on the help-page
Biostrings::getGeneticCode("1")     # "1" corresponds to NCBI code 1

# The R/ folder contains a script for a function that prints any genetic code
# in a standard table format: printGCtable(). It gets loaded from the .util.R
# script. We can print the standard code with it:

printGCtable(Biostrings::getGeneticCode("1"))
# Or ...
printGCtable(Biostrings::getGeneticCode("1"), format = "Name")


# =    4  A RANDOM GENETIC CODE  ===============================================

# The basic function to randomize collections in R is sample(). Without sample()
# takes a vector of things as input - that's what is being sampled from. Then it
# takes a number that determines how many elements should be sampled. Another
# parameter determines whether sampling is with replacement, or without
# replacement. And finally non-equal probabilities can be set for each of the
# inoput elements. Here are some examples:
#

# A six-faced die:
sample(1:6, 1)

# Rolling 13 times
sample(1:6, 13, replace = TRUE)

# ... etc. but the default parameters ensure that sampling a vector creates
# a random permutation:

sample(1:6)

# A random genetic code is one that randomly assigns amino acids to the 64
# codons. But there are some additional things to consider:

#  1. All twenty amino acids have to appear, plus at least one stop codon.
#  2. If we simply shuffle the letters, we keep eg. 6 Leu and one Met in the
#     code, and only put them into different positions. That seems overly
#     restrictive.
#  3. If  we allow the level of redundancy to change freely, we will find
#     that the "best" code has one stop codon, one each of every amino acid,
#     and Thr for the remaining 44 codons.

# Coming up with a balanced and realistic scheme for constructing random codes
# is not trivial. One way to achieve randomness but conserve some features is to
# shuffle elements around, and in our case a two-fold shuffling process appears
# appropriate:

# A: Put each amino acid into a random position.
# B: Replace each amino acid with a random other.



# ==   4.1  Constructing Randomized Codes  =====================================

# In the class R-project, we keep the genetic code information in a table that
# we load from the data/ directory during startup. This is the data frame
# "GCdf". Here are the 64 one-letter codes taken from that table:
cat(GCdf$A)


# ===   4.1.1  Redundancy in the Code


# We can sort the amino acids, for clarity
cat(sort(GCdf$A))

# R has the very useful function table() to summarize how many of each we
# had:

table(GCdf$A)  # this means "*" was repeated 3 times, "A" four times, etc.

# Sorted, in descending order of frequency:
sort(table(GCdf$A), decreasing = TRUE)

#    L R S A G P T V * I C D E F H K N Q Y M W
#    6 6 6 4 4 4 4 4 3 3 2 2 2 2 2 2 2 2 2 1 1



# This table means: the genetic code has ...
# L L L L L L   # Three amino acids with six codons each ...
# R R R R R R
# S S S S S S
#
# A A A A       # Five amino acids with four codons each ...
# G G G G
# P P P P
# T T T T
# V V V V
#
# * * *         # Three with three codons ...
# I I I
#
# C C           # Nine with two ...
# D D
# E E
# F F
# H H
# K K
# N N
# Q Q
# Y Y
#
# M             # ... and 2 singletons
# W


# ===   4.1.2  The amino acid "alphabet"

# We can extract the amino acid letters from column "A" of GCdf. If we
# take only the unique ones (no repetitions) we get the "alphabet"
# we are working with.
sort(unique(GCdf$A))


# A different approach to get the same result would have been:
myAA <- sort(names(table(GCdf$A)))
cat(myAA)


# ===   4.1.3  Codes with the same redundancy

# Our goal is to take elements from this alphabet and distribute them so
# that we have three sets of six of the same amino acids, five sets
# of four, ... etc.

# To achieve this, we can make a vector of numbers that we use as indices for
# our myAA vector. The vector should specify: take six of the first amino acid,
# six of the second, six of the third. That makes three groups of six. Then take
# four of the fourth, and the fifth, etc. If we spell this out literally, these
# instructions look like this:
x <- c( 1,  1,  1,  1,  1,  1,
        2,  2,  2,  2,  2,  2,
        3,  3,  3,  3,  3,  3,
        4,  4,  4,  4,
        5,  5,  5,  5,
        6,  6,  6,  6,
        7,  7,  7,  7,
        8,  8,  8,  8,
        9,  9,  9,
        10, 10, 10,
        11, 11,
        12, 12,
        13, 13,
        14, 14,
        15, 15,
        16, 16,
        17, 17,
        18, 18,
        19, 19,
        20,
        21)              # Every number "points" at a position of our alphabet.

length(x)               # Verify that these are 64 numbers
length(unique(x))       # Verify that these are 21 different numbers
all(unique(x) == 1:21)  # Verify that these are the numbers 1:21 in order
cat(x)                  # Visual check

# This vector represents the redundancy of the natural genetic code. We
# can construct it more compactly using the rep() function - but this really
# just makes the same vector ... and writing it out by hand as above is
# also oK: writing 64 numbers is not too many. But here is how we can do this
# using rep():

y <- sort(c(rep( 1:3,  6),
            rep( 4:8,  4),
            rep( 9:10, 3),
            rep(11:19, 2),
            20, 21))

length(y)                # must be 64 numbers
length(unique(myAA[y]))  # must be 21 unique numbers
all(unique(y) == 1:21)   # Verify that these are the numbers 1:21 in order
cat(y)                   # Visually confirm the sets of numbers
all(x == y)              # must be TRUE

# We can write this even more compactly, using both table() and rep():
idxAA <- rep(1:21, times = sort(table(GCdf$A), decreasing = TRUE))

length(idxAA)                # must be 64 numbers
length(unique(idxAA))        # must be 21 unique numbers
all(unique(idxAA) == 1:21)   # Verify that these are the numbers 1:21 in order
cat(idxAA)                   # Visually confirm the sets of numbers
all(x == idxAA)              # must be TRUE

# This last solution is not only more compact, but it is entirely data driven -
# it produces the correct "map" whatever is in the column GCdf$A, no manual
# input is required. But, again, this index vector is just the same that we
# wrote by hand above, just created more efficiently.

# Using such a "map" is a very generic way of producing patterns. The vector
# myAA defines WHAT we want to pattern, the vector idxAA defines HOW we want
# to arrange it.

# So how do we use this map?
# We simply extract elements from our vector according to the map:
cat(myAA[idxAA])  # we use cat() to suppress the quotation marks.

# Now, if we shuffle our "alphabet" with sample(), we get something like this:
randAA <- sample(myAA)
cat(randAA)   # the randomized alphabet

# And if we apply our map to this randomized alphabet:
cat(randAA[idxAA])   # the code-pattern with different amino acids

# The only thing left to do is to shuffle the result again, to break up the
# contiguous sets of letters:
cat(sample(randAA[idxAA]))  # a vector of random letters from our alphabet,
# but with the same redundancy as the original
# code.

# We can write this into a single, compact expression that will give us
# a different random code every time we execute the expression.
# First we shuffle the alphabet (1),
# then we use our map (2),
# then we shuffle the result (3):
#             (3)    (1)       (2)
randAA <- sample(sample(myAA)[idxAA])
cat(randAA)

# If this nested expression is confusing, just execute it step by step.
# Or break it apart with intermediate assignments:
cat(x      <- sample(myAA))  # (1)  shuffle the alphabet
cat(x      <- x[idxAA])      # (2)  use the map
cat(randAA <- sample(x))     # (3)  shuffle the result

# Confirm that this has the same structure as the original genetic code:
sort(table(randAA), decreasing = TRUE)

# With this approach, we can make random codes. Let's call a random code GCrand:

# First: extract the rows we need from GCdf ...
GCrand <- GCdf[ , c("First", "Second", "Third", "Codon", "A")]

# then replace column "A" with a randomized alphabet
GCrand$A <- sample(sample(myAA)[idxAA])

# Verify by printing the random code table in standard order:
printGCtable(nuc = c("T", "C", "A", "G"),
             order = c(3, 2, 1),
             dat = GCrand,
             format = "Aaa")


# ==   4.2  The randomized code generator  =====================================

# I have placed a randomized code generator in R/rGC.R . It
# is sourced from .util.R and thus automatically available in the course project.

# Example:
rGC()
printGCtable(rGC())


# =    5  NEIGHBOURING CODONS  =================================================

# "Neighbouring Codons" can be derived from a given codon by a single point
# mutation.

# Computing all nine neighbouring codons to a given codon just requires a bit of
# bookkeeping. See the function code R/neighCodons.R (sourced from .util.R) if
# you are interested.

# Example
neighCodons("TTT")


# =    6  EVALUATING THE STANDARD GENETIC CODE  ================================

# Define the sGC
# Initialize the sum of distances to zero
# For each codon in sGC
#   Get the encoded amino acid
#   For all nine neighbors
#      get the encoded amino acid
#      compute distance in feature space
#      add to the sum of distances


thisCode <- GCdf$A                      # Define the sGC
names(thisCode) <- rownames(GCdf)
sumDist <- 0                            # Initialize the sum of distances

for (codonX in names(thisCode)) {       # For each codon
  aaX <- thisCode[codonX]               # Get the encoded amino acid

  for (codonY in neighCodons(codonX)) { # For all nine neighbors
    aaY <- thisCode[codonY]             # get the encoded amino acid
    dist <- aaSim(aaX, aaY)             # compute distance in feature space
    sumDist <- sumDist + dist           # add to the sum of distances
  }
}

sumDist
sumDist / (64 * 9)

# There we have it. 1272.9 and a bit. That's the quality of the standard Genetic
# Code.
#
# So what do we do with that number? Is the universal genetic code any good?
#
# Is it?



# [END]
