# tocID <- "courseScripts/06-significance.R"
#
# Purpose: Evaluating the significance of an observation.
#
# Version: 1.0
# Date:    2024-10-09
# Author:  boris.steipe@utoronto.ca; ChatGPT-4 and 4o
#
# Versions:
#   1.1    More explanations in comments
#   1.0    In class - Class session 05: new script with code contributions
#          from FND-STA-Significance.R and 05-GeneticCodeExperiments.R
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
#TOC>   Section  Title                                              Line
#TOC> ------------------------------------------------------------------
#TOC>   1        INITIALIZATIONS                                      48
#TOC>   1.1        Parameters                                         50
#TOC>   1.2        Packages                                           53
#TOC>   2        RATIONALE AND PROCESS                                63
#TOC>   3        Significance and p-value                            125
#TOC>   3.1        Significance levels                               145
#TOC>   3.2        probability and p-value                           161
#TOC>   3.2.1          p-value illustrated                           193
#TOC>   4        One- or two-sided                                   250
#TOC>   5        Significance by integration                         297
#TOC>   6        Significance by simulation or permutation           307
#TOC>   7        Finally: The genetic code ...                       430
#TOC>
#TOC> ==========================================================================


# =    1  INITIALIZATIONS  =====================================================

# ==   1.1  Parameters  ========================================================


# ==   1.2  Packages  ==========================================================

# Package loader for packages and data in the Bioconductor project.
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")  # NOT install.packages(...) ...
}

# =    2  RATIONALE AND PROCESS  ===============================================

# We have previously
#  (A) defined an amino acid similarity function
#  (B) defined a way to create "random" genetic codes
#  (C) constructed R code that constructs adjacent codons to any given codon.

#  We can connect these elements together to create and evaluate random
#  genetic codes.

# Consider how we evaluate the standard genetic code:
#

thisCode <- GCdf$A                      # Define the sGC
names(thisCode) <- rownames(GCdf)

sumDist <- 0                            # Initialize a variable to contain
                                        # the sum of distances

for (codonX in names(thisCode)) {       # For each codon in the code
  aaX <- thisCode[codonX]               # ... get the encoded amino acid

  for (codonY in neighCodons(codonX)) { # For all nine neighbors of the codon
    aaY <- thisCode[codonY]             # ... get the encoded amino acid
    dist <- aaSim(aaX, aaY)             # ... compute distance in feature space
    sumDist <- sumDist + dist           # ... add to the sum of distances
  }
}

# This is the result: a measure of how the genetic code preserves amino acid
# similarity if it is subjected to point mutations
sumDist

# But this is just a single number, so we decided to compare this to random
# genetic codes. Essentially this works with the same R code, except we
# start from a random code.

myRandomSeed <- sample(1:.Machine$integer.max, 1)
thisCode <- rGC(myRandomSeed)           # Define a random genetic code
names(thisCode) <- rownames(GCdf)       # Name with the 64 codons

sumDist <- 0                            # Initialize the sum of distances

for (codonX in names(thisCode)) {       # For each codon
  aaX <- thisCode[codonX]               # Get the encoded amino acid

  for (codonY in neighCodons(codonX)) { # For all nine neighbors
    aaY <- thisCode[codonY]             # get the encoded amino acid

    dist <- aaSim(aaX, aaY)             # compute the distance
    sumDist <- sumDist + dist           # add to the sum of distances
  }
}

sumDist
# This number is different. But how different? Are we looking at a difference
# that we should expect simply from random fluctuations, is this a
# high-probability difference? Or is the difference significant? We need to
# talk about statistical probability and significance.



# =    3  Significance and p-value  ============================================

# The idea of the probability of an event has a precise mathematical
# interpretation, but how is it useful to know the probability? Usually we are
# interested in whether we should accept or reject a hypothesis based on the
# observations we have. A rational way to do this is to say: if the probability
# of observing the data is very small under some hypothesis of how it was
# generated, then we will assume the observation is due to something other than
# that hypothesis. Of course there are many ways in which the data could have
# been generated, and we can have many hypotheses about the process. But in
# statistics we are particularly interested in a process that is boring, i.e.
# uninformative - essentially a process that is merely due to random chance.
# Because, if we are able to reject that an observation is due to random chance,
# something else is going on, something that is more interesting and will tell
# us something worth knowing  about the system we are observing. We call this
# "uninformative" hypothesis the "null-hypothesis".

# But what do we mean by "the probability of observing the data is very small"?
# What is the "probability we are speaking of here? And what is "very small"?

# ==   3.1  Significance levels  ===============================================

# A "very small" probability is purely a matter of convention - a cultural
# convention. In the biomedical field we usually call probabilities of less than
# 0.05 (5%) small enough to reject the null-hypothesis. We call observations
# with a probability of less than 0.05 "significant" and if we want to highlight
# this in text or in a graph, we often mark them with an asterisk (*). Also we
# often call observations with a probability of less than 0.01 "highly
# significant" and mark them with two asterisks (**). But there is no special
# significance in these numbers, the cutoff point for significance could also be
# 0.0498631, or 0.03, or 1/(pi^3). 0.05 is just the value that the British
# statistician Ronald Fisher happened to propose for this purpose in 1925.
# Incidentally, Fisher later recommended to use different cutoffs for different
# purposes (cf. https://en.wikipedia.org/wiki/Statistical_significance).


# ==   3.2  probability and p-value  ===========================================

# But what do we even mean by the probability of an observation?
# Assume I am drawing samples from a normal distribution with a mean of 0 and a
# standard deviation of 1. The sample I get is ...

set.seed(sqrt(5))
x <- rnorm(1)
set.seed(NULL)

print(x, digits = 22)
# [1] -0.8969145466249813791748

# So what's the probability of that number? Obviously, the probability of
# getting exactly this number is very, very, very small. But also obviously,
# this does not mean that observing this number is in any way significant - we
# always observe some number. That's not what we mean in this case. There are
# several implicit assumptions when we speak of the probability of an
# observation:

# 1: the observation can be compared to a probability distribution;
# 2: that distribution can be integrated between any specific value
#      and its upper and lower bounds (or +- infinity).

# Then what we really mean by the probability of an observation in the context
# of that distribution is: the probability of observing that value, or a value
# more extreme than the one we have. We call this the p-value. Note that we are
# not talking about an individual number anymore, we are talking about the area
# under the curve between our observation and the upper (or lower) bound of the
# curve, as a fraction of the whole.


# ===   3.2.1  p-value illustrated

# Let's illustrate. First we draw a million random values from our
# standard, normal distribution:

N <- 1e6                             # one million
set.seed(112358)                     # set RNG seed for repeatable randomness
r <- rnorm(N)                        # N values from a normal distribution
set.seed(NULL)                       # reset the RNG

# Let's see what the distribution looks like:

(h <- hist(r, breaks = 40))

# The histogram details are now available in the list h -  e.g. h$counts

# Where is the value we have drawn previously?
abline(v = x, col = "#EE0000")

# How many values are smaller?
sum(r < x)

# Let's color the bars:
#    first, make a vector of red and green colors for the bars with breaks
#    smaller and larger then x, white for the bar that contains x ...
hCol <- rep("#EE000044", sum(h$breaks < x) - 1)
hCol <- c(hCol, "#FFFFFFFF")
hCol <- c(hCol, rep("#00EE0044", sum(h$breaks > x) - 1))
# ... then plot the histogram, with colored bars ...
hist(r, col = hCol, breaks = 40)
# ... add two colored rectangles into the white bar ...
idx <- sum(h$breaks < x)
xMin <- h$breaks[idx]
xMax <- h$breaks[idx + 1]
y <- h$counts[idx]
rect(xMin, 0, x, y, col = "#EE000044", border = TRUE)
rect(x, 0, xMax, y, col = "#00EE0044", border = TRUE)
# ... and a red line for our observation.
abline(v = x, col = "#EE0000", lwd = 2)

# The p-value of our observation is the red area as a fraction of the
# whole histogram (red + green).


# Task:
#    Explain how the expression sum(r < x) works to give us a count of values
#    with the property we are looking for. E.g., examine -4:4 < x

# Task:
#    Write an expression to estimate the probability that a value
#    drawn from the vector r is less-or-equal to x. The result you get
#    will depend on the exact values that went into the vector r but it should
#    be close to 0.185  That expression is the p-value associated with x.

sum(r <= x) / length(r)


# =    4  One- or two-sided  ===================================================

# The shape of our histogram confirms that the rnorm() function has returned
# values that appear distributed according to a normal distribution. In a normal
# distribution, readily available tables tell us that 5% of the values (i.e. our
# significance level) lie 1.96 (or approximately 2) standard deviations away
# from the mean. Is this the case here? How many values in our vector r are
# larger than 1.96?

sum(r > 1.96)
# [1] 24589

# Wait - that's about 2.5% of 1,000,000, not 5% as expected. Why?

# The answer is: we have to be careful with two-sided distributions. 2 standard
# deviations away from the mean means either larger or smaller than 1.96 . This
# can give rise to errors. If we are simply are interested in outliers, no
# matter larger or smaller, then the 1.96 SD cutoff for significance is correct.
# But if we are specifically interested in, say, larger values, because a
# smaller value is not meaningful, then the significance cutoff, expressed as
# standard deviations, is relaxed. We can use the quantile function to see what
# the cutoff values are:

quantile(r)
quantile(r, probs = c(0.025, 0.975)) # for the symmetric 2.5% boundaries
# close to ± 1.96, as expected
quantile(r, probs = 0.95) # for the single 5% boundary
# close to 1.64 . Check counts to confirm:
sum(r > quantile(r, probs = 0.95))
# [1] 50000
# which is 5%, as expected.

# Task:
# Use abline() to add the p = 0.05 boundary for smaller values to the histogram.

abline(v = quantile(r, probs = c(0.05)), col = "#8800BB", lwd = 0.5)

# This means: we consider outliers to the left of this line to be significant observations because they occur less than 5% of the total.

# To summarize: when we evaluate the significance of an event, we divide a
# probability distribution into two parts at the point where the event was
# observed. We then ask whether the integral over the more extreme part is less
# or more than 5% of the whole. If it is less, we deem the event to be
# significant.
#


# =    5  Significance by integration  =========================================

# If the underlying probability distribution can be analytically or numerically
# integrated, the significance of an observation can be directly computed. This
# is the basis of the many, many types of statistical tests that exist, and the
# underlying mathematics can be somewhat challenging. It is generally not
# recommended to publish work in which the conclusions depend on accurate
# statistics if you are not able to involve someone with expertise in the field.


# =    6  Significance by simulation or permutation  ===========================

# Whether the integration is correct, or relies on assumptions that may not
# be warranted for biological data, can be a highly technical question.
# Fortunately, we can often simply run a simulation, a random resampling, or a
# permutation and then count the number of outcomes, just as we did with our
# rnorm() samples. We call this an empirical p-value. (Actually, the "empirical
# p-value" is defined as (Nobs + 1) / (N + 1).  )

# Here is an example. Assume you have a protein sequence and
# you speculate that positively charged residues are close to negatively charged
# residues to balance charge locally. A statistic that would capture this is the
# mean minimum distance between all D,E residues and the closest R,K,H
# residue. Let's compute this for the sequence of yeast Mbp1.

MBP1 <- paste0("MSNQIYSARYSGVDVYEFIHSTGSIMKRKKDDWVNATHILKAANFAKAKRTRILEKEVLK",
               "ETHEKVQGGFGKYQGTWVPLNIAKQLAEKFSVYDQLKPLFDFTQTDGSASPPPAPKHHHA",
               "SKVDRKKAIRSASTSAIMETKRNNKKAEENQFQSSKILGNPTAAPRKRGRPVGSTRGSRR",
               "KLGVNLQRSQSDMGFPRPAIPNSSISTTQLPSIRSTMGPQSPTLGILEEERHDSRQQQPQ",
               "QNNSAQFKEIDLEDGLSSDVEPSQQLQQVFNQNTGFVPQQQSSLIQTQQTESMATSVSSS",
               "PSLPTSPGDFADSNPFEERFPGGGTSPIISMIPRYPVTSRPQTSDINDKVNKYLSKLVDY",
               "FISNEMKSNKSLPQVLLHPPPHSAPYIDAPIDPELHTAFHWACSMGNLPIAEALYEAGTS",
               "IRSTNSQGQTPLMRSSLFHNSYTRRTFPRIFQLLHETVFDIDSQSQTVIHHIVKRKSTTP",
               "SAVYYLDVVLSKIKDFSPQYRIELLLNTQDKNGDTALHIASKNGDVVFFNTLVKMGALTT",
               "ISNKEGLTANEIMNQQYEQMMIQNGTNQHVNSSNTDLNIHVNTNNIETKNDVNSMVIMSP",
               "VSPSDYITYPSQIATNISRNIPNVVNSMKQMASIYNDLHEQHDNEIKSLQKTLKSISKTK",
               "IQVSLKTLEVLKESSKDENGEAQTNDDFEILSRLQEQNTKKLRKRLIRYKRLIKQKLEYR",
               "QTVLLNKLIEDETQATTNNTVEKDNNTLERLELAQELTMLQLQRKNKLSSLVKKFEDNAK",
               "IHKYRRIIREGTEMNIEEVDSSLDVILQTLIANNNKNKGAEQIITISNANSHA")

# First we split this string into individual characters:
v <- unlist(strsplit(MBP1, ""))

# ... and find the positions of our charged residues

ED  <- grep("[ED]", v)   # Vector containing the location of E and D
RKH <- grep("[RKH]", v)  # Vector containing the location of R, K, and H

sep <- numeric(length(ED)) # Vector to hold the distances
for (i in seq_along(ED)) {
  sep[i] <- min(abs(RKH - ED[i]))
}

# Task: read and explain this bit of code

# Now that sep is computed, what does it look like?

table(sep)  # these are the minimum distances
# 24 of D,E residues are adjacent to R,K,H;
# the longest separation is 28 residues.

# What is the mean separation?
mean(sep)

# The value is 4.1 . Is this significant? Honestly, I would be hard pressed
# to solve this analytically. But by permutation it's soooo easy.

# First, we combine what we have done above into a function:

chSep <- function(v) {
  # computes the mean minimum separation of oppositely charged residues
  # Parameter: v (char) a vector of amino acids in the one-letter code
  # Value: msep (numeric) mean minimum separation

  ED  <- grep("[EDed]", v)
  RKH <- grep("[RKHrkh]", v)

  sep <- numeric(length(ED))
  for (i in seq_along(ED)) {
    sep[i] <- min(abs(RKH - ED[i]))
  }
  return(mean(sep))
}

# Execute the function to define it.

# Confirm that the function gives the same result as the number we
# calculated above:
chSep(v)

# Now we can produce a random permutation of v (a "shuffle"), and recalculate

set.seed(112358)    # set RNG seed for repeatable randomness
w <- sample(v)      # This shuffles the vector v. It is equivalent to
                    # sample(v, length(v), replace = FALSE).
                    # Memorize this code paradigm. It is very useful.

chSep(w)            # Recompute the separation with the shuffled vector:
# 3.773 ... that's actually less than what we had before.

# Let's do this 10000 times and record the results (takes a few seconds):

N <- 10000
chs <- numeric(N)
for (i in 1:N) {
  chs[i] <- chSep(sample(v, length(v))) # charge
}

hist(chs, breaks = 50)
abline(v = chSep(v), col = "#EE0000")

# Contrary to our expectations, the actual observed mean minimum charge
# separation seems to be larger than what we observe in randomly permuted
# sequences. Now, what would we mean by: "this is a significant result", or
# "this is not a significant result"? Where does our 5% bound appear? The shape
# of the probability distribution is no longer symmetric, so our assumptions
# about the normal distribution no longer hold...

# So here is where our empirical p-value come to the rescue. Instead of
# integrating over this unknown distribution, we simply count how many values
# are equal or larger than the one we have observed, and divide that by the
# total number of observations.

# Task:
# Calculate the empirical p-value for chsep(v)
x <- (sum(chs >= chSep(v)) + 1) / (length(chs) + 1)

# The result says: the minimum separation between oppositely charged residues in
# this sequence is a bit larger than what we would expect by random chance, but
# the difference is not significant: 12% of random observations are equal or
# even larger.


# =    7  Finally: The genetic code ...  =======================================

# Task: Map out how to calculate an empirical p-value for genetic code
# similarity. You can try simply handing the code snippet that evaluates the
# random code similarity to your AI tutor and asking them how to proceed, to
# explain the steps, and to plot the results. Also ask them to overlay the
# histogram you receive with the curve of a normal distribution.



# [END]
