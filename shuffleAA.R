# shuffleAA.R
#
# Starter code to shuffle 20 Amino acids and a stop codon to get the same
# level of redundancy as the original code but with randomized letters.
#
# CSB195 - 2023-10
# boris.steipe@utoronto.ca
# ==============================================================================

# The universal genetic code maps the 64 possible trinucleotide codons to
# 20 different amino acids (or a "stop" signal). Here is the code (using
# the one-letter notation for amino acids, and "*" for "stop").

#        TTT F	TCT S	 TAT Y	 TGT C
#        TTC F	TCC S	 TAC Y	 TGC C
#        TTA L	TCA S	 TAA *	 TGA *
#        TTG L	TCG S	 TAG *	 TGG W
#
#        CTT L 	CCT P  CAT H	 CGT R
#        CTC L	CCC P	 CAC H	 CGC R
#        CTA L	CCA P	 CAA Q	 CGA R
#        CTG L	CCG P	 CAG Q	 CGG R
#
#        ATT I 	ACT T  AAT N	 AGT S
#        ATC I	ACC T	 AAC N	 AGC S
#        ATA I	ACA T	 AAA K	 AGA R
#        ATG M	ACG T	 AAG K	 AGG R
#
#        GTT V 	GCT A  GAT D	 GGT G
#        GTC V	GCC A	 GAC D	 GGC G
#        GTA V	GCA A	 GAA E	 GGA G
#        GTG V	GCG A	 GAG E	 GGG G

# In the class R-project, we keep the genetic code information in a table that
# we load from the data/ directory during startup. This is the data frame
# "GCdf". Here are the 64 one-letter codes taken from that table:
cat(GCdf$A)

# We can sort them, for clarity
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

# We can extract the amino acid letters from column "A" of GCdf. If we
# take only the unique ones (no repetitions) we get the "alphabet"
# we are working with.
sort(unique(GCdf$A))


# A different approach to get the same result would have been:
myAA <- sort(names(table(GCdf$A)))
cat(myAA)

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
length(unique(myAA[idxAA]))  # must be 21 unique numbers
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
             order = c(1, 2, 3),
             dat = GCrand,
             aaColName = "A")


# That's it.



# [END]
