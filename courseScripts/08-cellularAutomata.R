# tocID <- "courseScripts/08-cellularAutomata.R"
#
# Purpose: Demonstrate the concept of cellular automata (CA); visualize CA
#          evolution; experiment with parameters.
#
# 2022-10 - 2024-10
# Boris Steipe (boris.steipe@utoronto.ca)
#
# Version:  1.3
#
# Versions:
#           1.3   Maintenance 2024
#           1.2.2 Maintenance 2023
#           1.2.1 Move utility scripts to R/
#           1.2   Move most code out ... teaching the code is not a critical
#                   objective here. Convert to "interactive" textbook" paradigm.
#           1.1   Updates after demo in course. More extensive explorations.
#           1.0   First course version
#
#  To Do:
#    Colour patterns (e.g. Gliders)
#    Fade out background ether pattern
#    Preload full window of evolving code
#
#  References:
#   Cook, Matthew (2004) Universality in Elementary Cellular Automata.
#        Complex Systems  15:1-40
#        https://wpmedia.wolfram.com/uploads/sites/13/2018/02/15-1-1.pdf
#
#    Wikipedia (2022) Cellular Automaton.
#        https://en.wikipedia.org/wiki/Cellular_automaton
#
#    Wikipedia (2022) Rule 110.
#        https://en.wikipedia.org/wiki/Rule_110
#
#    Wolfram, Stephen (1984) Universality and Complexity in Cellular Automata.
#        Physica  10D: 1-35
#        https://doi.org/10.1016/0167-2789(84)90245-8
#
#
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                  Line
#TOC> ----------------------------------------------------------------------
#TOC>   1        INTRODUCTION                                             62
#TOC>   2        FUNCTIONS                                                75
#TOC>   3        FIRST STEPS                                             108
#TOC>   4        EXPLORATIONS                                            189
#TOC>   4.1        First exploration. Step through the rules ...         192
#TOC>   4.2        Digression: Fibonacci words                           268
#TOC>   4.3        Four classes                                          297
#TOC>   4.3.1          Digression: initializations                       322
#TOC>   5        EVERY SINGLE CA                                         353
#TOC>   6        LONG EVOLUTIONS                                         364
#TOC>
#TOC> ==========================================================================


# =    1  INTRODUCTION  ========================================================

# Taking a perspective on "life" as a process of choice (and computation), we
# introduce models that are able to **sense** the environment and make
# decisions, they perform - in the simplest, most fundamental way - information
# processing. The simplest such model consists of a single "cell" in a
# one-dimensional world, which can sense the state of its left and right
# neighbour, and then decide how its own state will change according to a
# compact set of rules.

# Such systems are called "Cellular Automata" (CA).


# =    2  FUNCTIONS  ===========================================================
#

# The key functions we need here are in the source()d files below. You are
# welcome to look at them in detail, but there is a lot of bookkeeping going on
# which may obscure that they are really very simple. Just be sure to understand
# what they do in principle. I have added documentation to the functions so can
# always use the help-function to recall the meaning of the parameters: just
# type ?<function name>. And where you find a function call of the form
# plotFigure() below, for code that creates the "figures" for this script, you
# can view the source by calling the function with the parameter "showSource" in
# the arguments. e.g. plotFigure("CA.1", showSource)


source("R/plotFigure.R")    # Figures
# == plotFigure()       === compute and plot a figure
# == catSection()       === optionally cat() the script source

source("R/imPlot.R")        # plot the contents of a matrix
# == imPlot()           === plots a raster image
# == imText()           === adds text to a raster image

source("R/CAtools.R")       # sundry utilities to work with cellular automata
# == printRule()         === prints the rule set of a cellular automaton
# == plotRule()         === plots a graphic of a cellular automaton's rule set
# == applyRule()        === applies a rule to three input elements
# == CA()               === run a cellular automaton

source("R/fibWord.R")       # compute a Fibonacci word
# == fibWord()          === a non-repeating string of 1s and 0s



# =    3  FIRST STEPS  =========================================================
#

# A cell that senses its own state and that of its two neighbors, given that
# the states can be 0 or 1, has a state-space of eight possibilities:

plotFigure("CA.1")

# The "Wolfram code" associates with each of the eight 0s or 1s a zero or one.
# Consider, for example, Rule 18.
printRule(18)
plotFigure("CA.2")

# This means: if a cell is in state 0, and its LEFT neighbour is "1", then it
# switches to state 1. And if a cell is in state 0, and its RIGHT neighbour is
# "1", then it also switches to state 1. But for all other states, it either
# switches to or remains in state 0.

# Note that the input states are always the same, and describe all possible
# states of a "cell" and its immediate left and right neighbours.
#
#    111  110  101  100  011  010  001  000
#
# Incidentally, these descriptive strings of 0's and 1's can be interpreted
# as the eight binary numbers from 7 (111) to zero (000).
#     7    6    5    4    3    2    1    0
# ... but don't be confused: that's just a mapping and has nothing to do with
#  what the state of the cells mean. All of these input states must be accounted
#  for, and each one can result in a 1 or 0 of the middle element - sometimes
#  this changes the state of the element, and sometimes it keeps it the same.
#
# Thus the possible output states are a sequence of 8 binary digits. Just as
# we can associate the three digits of the input state with decimal numbers,
# we can do the same with the output digits.
#
# The rule for converting binary to decimal states that each binary digit
# stands for one of the powers of 2 from left to right these, from zero to
# seven in our case

print( 2^(7:0) )

# So the binary digits \0b 00010010 are converted to decimal digits by
# computing: 0*128 + 0*64 + 0*32 + 1*16 + 0*8 + 0*4 + 1*2 + 0*1 or ...
#    16 + 2, which is 18. This means the Wolfram code is "self-describing"!
# A code number contains the instruction set for its application.

# When we say a CA "evolves", we mean that the CA iterates a 1D input "world"
# through a number of applications of its transformation rules:

plotFigure("CA.3")

# Note that this process is undefined at the edges. We could always pad the
# edges with 0, or we could "wrap the edges around" ... the right-most neighbour
# of the last cell is the first cell, and the left-most neighbour of the first
# cell is the last cell. If you imagine evolving the CA on a piece of graphing
# paper, this is like turning  the paper into a cylinder by gluing the edges
# together. What goes out on the right comes back in on the left. We also call
# this "periodic boundary conditions".

plotFigure("CA.4")

# Let's evolve a few cellular automata. We start from a single point in the
# middle and apply a given rule step by step. After every line we pause;
# hit <return> in the console to advance one more iteration, or hit <esc> to
# stop.

plotFigure("CA.5", rule = 18)

# Try a different rule: Rule 2  - do you see how this works?

printRule(2)
plotFigure("CA.5", rule = 2)

# Can you make it go in the other direction ?

printRule(16)
plotFigure("CA.5", rule = 16)

# To look for interesting rules, we should consider a larger matrix ...


# =    4  EXPLORATIONS  ========================================================


# ==   4.1  First exploration. Step through the rules ...  =====================

# It is tedious to enter the rule number by hand. We can apply a little trick:
# if we assign a value inside parentheses, the result of the assignment is
# available outside the parentheses. We can use this to increment a number, step
# by step. Note what gets printed:

(i <- 0)       # 0
(i <- i + 1)   # 1
(i <- i + 1)   # 2
(i <- i + 1)   # 3

# ... etc.

# We can use the same syntax to plot CAs one after another:

iRule <- 0

# Select the entire funcxtion call below.  When you hit <cmd>-<return>, iRule is
# incremented by one, and the evolution of the CA is displayed. You can step
# through the different rules, simply by hitting <cmd>-<return> again and again.
# (You can also change the parameters ...)

plotFigure("CA.6", iRule = (iRule <- iRule + 1),
           nx = 50, ny = 71, drawGrid = TRUE)

# You will quickly notice that the behaviour of our CAs can be **very**
# different. Rules 0 and 8 do nothing. (Why?) Rules such as 2, 13, 23, 24, and
# 50 show very simple behaviour. Rules 18, 22, and 26 evolve fractals...

plotFigure("CA.6", iRule = 26, nx = 300, ny = 500)

# And rules like 30, and 45? They develop into some disordered state from which
# no order seems to be recoverable.

plotFigure("CA.6", iRule = 30, nx = 300, ny = 500)
plotFigure("CA.6", iRule = 45, nx = 300, ny = 500)

# ... and then there is rule 110 - which appears to sit somewhere between order
# and disorder:

plotFigure("CA.6", iRule = 110, nx = 300, ny = 500)

# Stephen Wolfram surmised in 1984 that CAs would fall into four different
# classes, depending on their long-term behaviour.

# Class I CAs    converge from almost any initial state to a single point
#                in state-space. Rules 4, 8, 36 ... are examples.

plotFigure("CA.6", iRule = 4, nx = 300, ny = 500)

# Class II CAs:  enter into some periodic behavior over long time, similar to
#                the phase-space cycles we had observed in Lotka-Volterra
#                systems.  Rules 1, 23, 25, 26, and 27 are examples.

plotFigure("CA.6", iRule = 1, nx = 300, ny = 500)

# Class III CAs: exhibit aperiodic, or "chaotic" behavior over time. Such
#                systems may appear similar at various times
#                in their trajectory, but they are never exactly the same.
#               Frequently, these have a fractal geometry. Rules 18, 45, 60,
#               75, 86 and 105 are  striking examples.
#
plotFigure("CA.6", iRule = 18, nx = 300, ny = 500)

# Class IV CAs:  are undecidable. The quintessential example is the
#                famous rule 110.

plotFigure("CA.6", iRule = 110, nx = 300, ny = 500)

# Try some of these by setting the parameter iRule to their respective value,
# for example:

plotFigure("CA.6", iRule = 193, nx=300, ny=500)  # Have you seen this before?


# ==   4.2  Digression: Fibonacci words  =======================================

# To illustrate these CAs, we'll not start with a random initial state, but with
# a "Fibonacci Word" - a binary sequence that has locally similar regions but
# that never repeats This sequence is related to the famous Fibonacci sequence
# ...

# 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89 ...                    (OEIS A000045)

# ... i.e. each term is the sum of its two predecessors. The sequence has
# many, many beautiful and fascinating properties, e.g. the ratio of two terms
# approximates the golden ratio, (1 + sqrt(5)) / 2 , as the terms grow larger,
# and Fibonacci numbers come up frequently in the geometry of living systems.
# But that's an exploration for another day. Here we use the binary Fibonacci
# sequence:

fibWord()

# We can see its non-repeating property nicely when we "fold" the sequence
# along one of the Fibonacci numbers:

plotFigure("CA.7", i = 55, j = 34)

# The sequence arranges into columns, more or less, but it is not exactly
# repetitive, there are always locations where it shifts. This is useful for us,
# since it reduces the artefacts of evolving CAs that are introduced
# by wrapping our space and creating a periodic boundary.


# ==   4.3  Four classes  ======================================================

# Here we plot only one of the four classes each with the Fibonacci word initial
# state...

# Class 1
plotFigure("CA.8", iRule = 36,
           nx = 233, ny = 233, vInit = "fib")

# Class 2
plotFigure("CA.8", iRule = 45,
           nx = 233, ny = 233, vInit = "fib")

# Class 3
plotFigure("CA.8", iRule = 86,
           nx = 233, ny = 233, vInit = "fib")
plotFigure("CA.8", iRule = 105,
           nx = 233, ny = 233, vInit = "fib")
# (That's two examples - hard to resist.)

# Class 4
plotFigure("CA.8", iRule = 110,
           nx = 233, ny = 233, vInit = "fib")


# ===   4.3.1  Digression: initializations

# Nb. the code I wrote for this exploration can generate a number of different
# initializations:

# a single cell ...
plotFigure("CA.8", iRule=110, nx=233, ny=233, vInit = 1   )

# several isolated cells ...
plotFigure("CA.8", iRule=110, nx=233, ny=233, vInit = 3   )

# a repeated vector
plotFigure("CA.8", iRule=110, nx=233, ny=233, vInit = c(0,0,0,0,1,1,1,1)   )

# a repeated vector, but created from a string
plotFigure("CA.8", iRule=110, nx=233, ny=233, vInit = c("01001110000011111111"))

# random: a floating point number that defines the fraction of coloured cells.
plotFigure("CA.8", iRule=110, nx=233, ny=233, vInit = 0.5)

# random with p(0) != p(1)
plotFigure("CA.8", iRule=110, nx=233, ny=233, vInit = 0.1)
plotFigure("CA.8", iRule=110, nx=233, ny=233, vInit = 0.9)

# Fibonacci word
plotFigure("CA.8", iRule=110, nx=233, ny=233, vInit = "fib")

# You are invited to explore varying sizes and initializations with the samples
# of evolving CAs we discuss below.


# =    5  EVERY SINGLE CA  =====================================================

# Let's look at all of the CA's to get a sense of what the space of
# possibilities looks like. We can plot all CAs from a loop. After each plot we
# sleep for a second or so (you can adjust that) and you can interrupt the cycle
# by pressing the red STOP-sign.

plotFigure("CA.9", start=0, end=255, sleep=1.0, nx=377, ny=610, vInit = "fib")



# =    6  LONG EVOLUTIONS  =====================================================

# We may be interested to let the evolution run for a longer time to see whether
# the CA ultimately settles into a repeating state. Here is a scrolling
# evolution:

# Note: this loop plots A LOT. The default plot window doesn't like to be
# updated so frequently. It "buffers" output until the process ends. I am sure
# there is a sane reason for this, but for what we are doing here, it breaks the
# process. That's bad. Therefore we will open an external graphics window with
# dev.new(). This is reasonably responsive on the Mac, but I don't know how well
# it works on Windows. Let me know.

dev.new() # open a fresh plot window

# ... then evolve CA 165 over 800 steps
plotFigure("CA.10", iRule=165, nx=233, ny=377, vInit="fib", nSteps = 800)

# Here too, you should be able to interrupt the scrolling by pressing the red
# STOP-sign.
#
# Does this CA ultimately stabilize? Does the chaotic evolution remain chaotic?

# I find this particular evolution intriguing for two reasons:
#
# (1) the rule is symmetric ...

printRule(165)

# ... and although many local axes of symmetry arise, there is one major axis at
# ~ 0.75 nx that remains unperturbed throughout the entire evolution. Indeed, if
# you look closely, you will find that there is a second global axis of symmetry
# that is wrapped around, where you would expect it ... at ~0.25 nx. Which is
# remarkable since nx (=233) is an odd number! If we run this with nx = 234 the
# axis shifts to the center, and the overall feeling is very different. Over
# time, it feels like the character of the evolution is slowly changing towards
# smaller and smaller features... until the large triangles reappear!

plotFigure("CA.10", iRule=165, nx=234, ny=377, vInit="fib", nSteps = 800)

# ... (2) it does not get boring.

# How about this one ...
plotFigure("CA.10", iRule=165, nx=64*8, ny=300, vInit="fib", nSteps = 300)
# Wait for it ...

# There is much to be explored. For this course, we are pursuing this in a
# rather qualitative fashion; I would be happy if you can develop something
# of an intuition what is and what is not possible with these simple automata,
# and get a feel for what "the edge of chaos" looks like.

# Of course, rule 110 is particularly intriguing because of its
# undecidability: will it evolve into a trivial, locally repetitive state? Will
# it remain chaotic? Will it attain some large-scale ordered structure? Here is
# a nice example:

plotFigure("CA.10", iRule=110,
           nx=98, ny=294,
           vInit=0.5,seed=2748,
           nSteps=1200)

# This evolution of CA 110 starts off from a jumbled array of triangles, formed
# and vanishing without order or purpose. Soon a crystalline background ether
# congeals from the tempest, it harbours thicker and thinner slanted fibrils
# that manifest, then interact, then annihilate each other. A much larger
# triangle appears unexpectedly, a singular occurrence, soon fading from view as
# a distinct vertical track of cells runs down the middle for a hundred fifty
# steps maybe. It too soon ceases to be when after some five hundred steps the
# patterns converge on two doubled, repetitive disturbances. The CA has reached
# its final state as those two self-sustaining perturbations spawn and
# annihilate dotted feathers locally, with a periodicity of 30 steps, and
# forever continue on their path through their world - side-by-side, yet
# oblivious of each other. 30-step cycles arise - from the behaviour of cells
# whose horizon is not more than three cells across.

# And sometimes (e.g. seed = 72) they go on, and on - and die.


dev.off()

# [END]
