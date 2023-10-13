# shuffleAA.R
#
# Starter code to shuffle 20 Amino acids and a stop codon to get the same
# level of redundancy as the original code but with randomized letters.
#
# CSB195 - 2023-10
# boris.steipe@utoronto.ca
# ==============================================================================

# Here are the 64 one-letter codes in our table of the genetic code:
cat(GCdf$A)

# We can sort them, for clarity
cat(sort(GCdf$A))

# R has the very useful function table() to summarize how many of each we
# had:

table(GCdf$A)  # this means "*" was repeated 3 times, "A" four times, etc.


# Sorted, in descending order of frequency:
sort(table(GCdf$A), decreasing = TRUE)

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

# The amino acid letters are in column "A" of GCdf. Extract them and
# take only the unique ones (no repetitions). This is the "alphabet"
# we are working with.
myAA <- sort(unique(GCdf$A))

cat(myAA)

# A different approach for the same result would have been:
cat(sort(names(table(GCdf$A))))

# Our goal is to distribute this alphabet so that we have three sets of
# six of the same amino acids, five sets of four, ... etc.

# We can make a vector of numbers that we use as indices for our myAA vector.
# It should define: take six of the first amino acid, six of the second,
# six of the third. That makes three groups of six. Then take four of the
# fourth, and the fifth, etc. if we spell this out literally it would
# lokk like this.
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
       21)

length(x)               # Verify that these are 64 numbers
length(unique(x))       # Verify that these are 21 different numbers
all(unique(x) == 1:21)  # Verify that these are the numbers 1:21 in order
cat(x)                  # Visual check

# This vector represents the redundancy of the natural genetic code. We
# can construct it more compactly using the rep() function - but this really
# just makes the same vector ... and writing it out by hand as above is
# also oK: writing 64 numbers is not too many. But here is how we can do this
# more compactly:

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

# Or, even more compactly, using table()
idxAA <- rep(1:21, times = sort(table(GCdf$A), decreasing = TRUE))

length(idxAA)                # must be 64 numbers
length(unique(myAA[idxAA]))  # must be 21 unique numbers
all(unique(idxAA) == 1:21)   # Verify that these are the numbers 1:21 in order
cat(idxAA)                   # Visually confirm the sets of numbers
all(x == idxAA)              # must be TRUE

# This last solution is not only more compact, but it is entirely generic -
# it produces the correct "map" whatever is in the column GCdf$A, no manual
# input is required. But, again, this index vector is the same that we wrote
# by hand above.

# Using such a "map" is a very generic way of producing patterns. The vector
# myAA defines WHAT we want to pattern, the vector idxAA defines HOW we want
# to arrange it.

# So how do we use this?
# We simply extract elements from our vector according to the map:
cat(myAA[idxAA])

# Now, if we shuffle our "alphabet" with sample(), we get something like this:
randAA <- sample(myAA)
cat(randAA)

# And if we apply our map to this randomized alphabet:
cat(randAA[idxAA])

# The only thing left to do is to shuffle the result again, to break up the
# contiguous sets of letters:
cat(sample(randAA[idxAA]))

# We can write this into a single, compact expression that will give us
# a different random code every time we execute the expression. First we
# shuffle the alphabet, then we use our map, then we shuffle the result:

randAA <- sample(sample(myAA)[idxAA])
cat(randAA)

# Confirm that this has the same structure as the original genetic code:
sort(table(randAA), decreasing = TRUE)

# With this approach, we can make random codes:

# First: extract the rows we need from GCdf ...
GCrand <- GCdf[ , c("First", "Second", "Third", "Codon", "A")]

# then replace column "A" with shuffled set of letters
GCrand$A <- sample(sample(myAA)[idxAA])

# Verify by printing the random code table in standard order:
printGCtable(nuc = c("T", "C", "A", "G"),
             order = c(1, 2, 3),
             dat = GCrand,
             aaColName = "A")


# That's it.



# [END]
