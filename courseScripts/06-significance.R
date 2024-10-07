# tocID <- "courseScripts/06-significance.R"
#
# Purpose: Evaluating the segnificance of an observation.
#
# Version: 1.0
# Date:    2024-10-07
# Author:  boris.steipe@utoronto.ca; ChatGPT-4 and 4o
#
# Versions:
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
#TOC>   1        INITIALIZATIONS                                      47
#TOC>   1.1        Parameters                                         49
#TOC>   1.2        Packages                                           52
#TOC>   2        RATIONALE AND PROCESS                                62
#TOC>   3        Significance and p-value                            115
#TOC>   3.1        Significance levels                               126
#TOC>   3.2        probability and p-value                           143
#TOC>   3.2.1          p-value illustrated                           175
#TOC>   4        One- or two-sided                                   232
#TOC>   5        Significance by integration                         277
#TOC>   6        Significance by simulation or permutation           283
#TOC>   7        Finally: The genetic code ...                       395
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
#  - defined an amino acid similarity function
#  - defined a way to create "random" genetic codes
#  - constructed code that identifies the adjacent codons to any given codon.

#  We can connect these elements together to create and evaluate random
#  genetic codes.

# The standard genetic code
# =========================
#

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


# A random code
# =============


thisCode <- rGC(112358)                 # Define a random genetic code
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



# =    3  Significance and p-value  ============================================

# The idea of the probability of an event has a precise mathematical
# interpretation, but how is it useful to know the probability? Usually we are
# interested in whether we should accept or reject a hypothesis based on the
# observations we have. A rational way to do this is to say: if the probability
# of observing the data is very small under the null-hypothesis, then we will
# assume the observation is due to something other than the null-hypothesis. But
# what do we mean by the "probability of our observation"? And what is "very
# small"?

# ==   3.1  Significance levels  ===============================================

# A "very small" probability is purely a matter of convention - a cultural
# convention. In the biomedical field we usually call probabilities of less then
# 0.05 (5%) small enough to reject the null-hypothesis. Thus we call
# observations with a probability of less than 0.05 "significant" and if we want
# to highlight this in text or in a graph, we often mark them with an asterisk
# (*). Also we often call observations with a probability of less than 0.01
# "highly significant" and mark them with two asterisks (**). But there is no
# special significance in these numbers, the cutoff point for significance could
# also be 0.0498631, or 0.03, or 1/(pi^3). 0.05 is just the value that the
# British statistician Ronald Fisher happened to propose for this purpose in
# 1925. Incidentally, Fisher later recommended to use different cutoffs for
# different purposes (cf.
# https://en.wikipedia.org/wiki/Statistical_significance).


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

(h <- hist(r))

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
hist(r, col = hCol)
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

abline(v = quantile(r, probs = c(0.05)))

# To summarize: when we evaluate the significance of an event, we divide a
# probability distribution into two parts at the point where the event was
# observed. We then ask whether the integral over the more extreme part is less
# or more than 5% of the whole. If it is less, we deem the event to be
# significant.
#


# =    5  Significance by integration  =========================================

# If the underlying probability distribution can be analytically or numerically
# integrated, the siginificance of an observation can be directly computed.


# =    6  Significance by simulation or permutation  ===========================

# But whether the integration is correct, or relies on assumptions that may not
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

# first we split this string into individual characters:
v <- unlist(strsplit(MBP1, ""))

# and find the positions of our charged residues

ED  <- grep("[ED]", v)
RKH <- grep("[RKH]", v)

sep <- numeric(length(ED)) # this vector will hold the distances
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

# Now we can produce a random permutation of v, and recalculate

set.seed(pi)                       # set RNG seed for repeatable randomness
w <- sample(v, length(v))          # This shuffles the vector v. Memorize this
# code paradigm. It is very useful.
set.seed(NULL)                     # reset the RNG



chSep(w)
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
# sequences. Now, what would we mean by: this is a significant result, or
# this is not a significant result? Where does our 5% bound come in?

# Task:
# Calculate the empirical p-value for chsep(v)
x <- (sum(chs >= chSep(v)) + 1) / (length(chs) + 1)

# =    7  Finally: The genetic code ...  =======================================





# Note: sample(1:.Machine$integer.max, 1)

# [END]
