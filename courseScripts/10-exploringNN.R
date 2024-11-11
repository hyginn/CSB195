# tocID <- "courseScripts/10-exploringNN.R"
#
#
# Purpose: Define simple neural networks. Run them. Experiment with parameters.
#            Find "interesting" variants.
#
#
# Version: 1.10
# Date:    2023-10 - 2024-11
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.10  Add some tools to generate random weight matrices. Add a
#         network plot function with weight-colored edges. Add a phase-space
#         density plot as the default option for the third plotting panel.
#   1.9   Update iGraph network sketch. Changed name from "evolvingNN.R" to
#         "exploringNN.R".
#   1.0   Posted for course assignment
#   0.5   More details about spectrum analysis.
#   0.4   Make sure input is from a column vector!
#   0.3   Refactored for better separation of concerns - Defining the network,
#           constructing input, plotting time-course, analyzing output.
#   0.2   Added theoretical explanations
#   0.1   Starter Code
#
# ToDo:
#   2.0  Turn demo code into .Rmd
#   Use fields::image.plot() for the color legend.
#
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                   Line
#TOC> -----------------------------------------------------------------------
#TOC>   1        Preparation: packages                                     86
#TOC>   2        Introduction                                             107
#TOC>   2.1        Terminology:                                           112
#TOC>   2.2        What goes on in a neuron                               176
#TOC>   2.3        A computational neuron                                 190
#TOC>   2.4        Activation functions                                   234
#TOC>   3        Computational Neurons in a Network                       304
#TOC>   3.1        Defining a Neural Network with a Weight matrix         376
#TOC>   3.1.1          Digression: colors                                 416
#TOC>   3.2        Plot the network with igraph::                         461
#TOC>   4        Running the network                                      529
#TOC>   5        Input                                                    638
#TOC>   5.1        Feeding the network: single pulse                      697
#TOC>   5.2        Feeding the network: harmonic oscillator               710
#TOC>   6        Workbench                                                740
#TOC>   6.1        Running experiments ...                                854
#TOC>   7        "Interesting" Networks                                   889
#TOC>   7.1        Random weight matrices                                 923
#TOC>   7.2        Putting everything together                           1013
#TOC>
#TOC> ==========================================================================


################################################################################
#                                                                              #
#                            S O U R C E   S A F E                             #
#                                                                              #
#     This script is  "source safe".  source()'ing the  code will  define      #
#     (refresh) the global parameters, and define the functions. It won't      #
#     actually run code. Execute the statements in the if (FALSE) { ... }      #
#     blocks to run the interactive parts of this script.                      #
#                                                                              #
#     To make changes in the code and experiment with it, make your own        #
#     copy in your myScripts folder.                                           #
#                                                                              #
################################################################################


# NOTE: You need to read and explore this script, execute the functions but also
# try what happens if you change things.

# OPEN A NEW R SCRIPT into which you write your code experiments. Then you can
# easily reproduce what you did, and the experiments are automatically
# documented for your report.



# =    1  Preparation: packages  ===============================================

# A package to support plotting. We need this  to plot lines with a colour
# gradient which we can't do in base R.
if (! requireNamespace("plotrix", quietly=TRUE)) {
  install.packages("plotrix")
}
# Package information:
#  library(help   = plotrix)     # basic information
#  data(package  = "plotrix")    # available datasets

# A package to compute with graphs. iGraph is the standard package for working
# with all manners of graphs in R, here I only use it to plot the topology
# of some simple neural networks.
if (! requireNamespace("igraph", quietly=TRUE)) {
  install.packages("igraph")
}
# Package information:
#  library(help   = igraph)     # basic information on available commands


# =    2  Introduction  ========================================================

# In this script we explore how to construct simple neural networks and try to
# evolve complex behavior.

# ==   2.1  Terminology:  ======================================================

# Neuron - the elementary unit:
# ======-----------------------
# The term "neuron" comes from an analogy of computational "Artificial Neural
# Networks" (ANN, or simply NN) with biological networks of interconnected nerve
# cells. In biology, a nerve cell is called a neuron. In some contexts one might
# speak of "nodes", especially in the context of graph theory, or when
# discussing non-neural network machine learning models. In our context, we
# always call the elementary unit of computation a "neuron".

# Input - what goes into a neuron:
# =====---------------------------
# Information is passed between neurons, and what a neuron receives is its
# "input". The terms "signal" or "activity" is sometimes used, however, input is
# the most generic and easily understood.

# Output - what a neuron produces:
# ======--------------------------
# The result of a neuron's information processing is its "output". In some
# contexts, "activation" is used, especially when referring to the value a
# neuron produces after applying its activation function. But output is the
# generally accepted term.

# Connection - what allows neurons to exchange information:
# ==========-----------------------------------------------
# In general we speak of two neurons having a "connection" when they exchange
# information. "Edge" is a term borrowed from graph theory (sometimes also
# "arc") and might be encountered when discussing networks from a topological or
# graphical perspective. "Dendrite" and "axon" are biological terms for the
# cellular structures that form connections, and a "Synapse" is the actual point
# of contact and interface for passing information, these can be used to
# emphasize biological parallels. When speaking of the system of connections
# we refer to the"topology", or "architecture" of the network.

# Source - A neuron from which input is received:
# ======-----------------------------------------
# A neuron receives input from a "source". In network theory, "predecessor"
# might be used. In graph theory, we sometimes encounter "tail" to refer to the
# node of origin for an edge in a directed graph. In graph theory a "source" is
# a node that has no incoming edges, so our usage here is more general.

# Target - A neuron to which output is sent:
# ======------------------------------------
# The neuron that receives output from another is its "target". In network
# theory, "successor" might be used. In graph theory, we sometimes encounter
# "head" to refer to the node that an a directed edge connects to, or the
# related term  "sink" - a node that has no outgoing edges.

# So: A NEURON receives INPUT through a CONNECTION to a SOURCE, and transmits
# OUTPUT to a TARGET.

# Iteration - updating all neurons in turn:
# =========--------------------------------
# We call the process of receiving input for all neurons and updating
# the output, one "iteration" of the network. The term "step" is sometimes used
# equivalently, but we will use it only (if at all) to refer to the update of a
# single neuron. Related are "ticks" - elementary units of time, for example
# in the Game of Life; "cycles" is sometimes used to refer to cyclical
# or recurrent processes. "Epoch", and "batch" are machine learning terms: one
# full optimization cycle is an "epoch", a subset of data that is used in
# learning is a "batch".


# ==   2.2  What goes on in a neuron  ==========================================

# A neuron processes its input in the following way:
#   1. Each input is multiplied with a weight that the neuron stores. This
#      is the "strength" of the connection. The weighted inputs are all
#      summed together.
#   2. To this is added a "bias". This determines the baseline level, or
#      sensitivity of the neuron.
#   3. The result is passed through an "activation function". This determines
#      the "output" of the neuron.
#   4. The output is made available as input for every connected neuron.
#   5. The process is iterated.


# ==   2.3  A computational neuron  ============================================

# If we implement such a computational neuron, we could proceed as follows:

if (FALSE) {

  # Assume our neuron receives input from three other neurons. Let's just
  # pick -1, 3, and 5 as those inputs:
  (input <- c(-1, 3, 5))     # NOTE: wrapping a statement with parentheses
  #       performs the assignment AND prints the
  #       result in one step.

  # When the neuron receives that input, it multiplies each input by a weight.
  # So let's define weights, again: arbitrarily
  (weights <- c(-3, 1, -2))

  # Finally, we define a bias ...:
  bias <- 0.1

  # And an activation function:
  fAct <- function(x) { return(tanh(x)) }

  # With that, we can compute one step of information processing.

  # First we multiply the inputs with the weights:
  (x <- input * weights)
  # ... then we add them all together:
  (x <- sum(x))
  # ... then we add the bias:
  (x <- x + bias)
  # ... then we apply the activation function.
  (x <- fAct(x))

}

# That's it! That's how our neuron responds to the input it received
# from three other neurons. That's the value it will send to all neurons
# it connects to.


# Let's talk about the activation function - the function that takes an input
# and computes the actual value of the output.


# ==   2.4  Activation functions  ==============================================

# There are many ways to define activation functions.
# (cf. https://en.wikipedia.org/wiki/Activation_function)
#
# Here we use three common ones but it should be easy for you to define
# more. We will define:
#   - the logistic function (or sigmoid function),
#   - the hyperbolic tangent
#   - the ReLU function

# Let's define them first, then compare them.

# The logistic function returns very small values over most of its negative
# domain then softly switches to positive values. I have transformed its
# location and scale to approximate the behaviour of the ReLU() function.
fLogistic <- function(x, loc = 2.3, scale = 1, A = 5) {
  return(A * plogis(x, location = loc, scale = scale))
}

# The hyperbolic tangent looks similar to the logistic function, but it is
# symmetric around the x-axis. It is the only of the three functions that will
# return negative values (although the other functions too can result in
# negative outputs if the neuron has negative weights on some inputs.)
fTanh <- function(x, A = 5, scale = 0.3) {
  return(A * tanh(x * scale))
}

# The ReLU (Rectifier Linear Unit) activation function is very often used
# in machine learning applications. It is very simple:
# return x, if x is positive, 0 otherwise.
ReLU <- function(x) {
  return(ifelse(x > 0, x, 0))
}

# Plot them to compare them:
if (FALSE) {

  inputs <- seq(-8, 8, by = 0.1)  # Some sequence of x-values for which
                                  # we want our functions plotted.

  # first plot: fTanh()
  plot(inputs, fTanh(inputs, A = 5, scale = 0.3), type = "l", lwd = 2,
       main = "Activation functions", xlab = "input", ylab = "output")
  abline(v = 0 , col = "#CCDDEE")
  abline(h = 0 , col = "#CCDDEE")

  # overlay fLogistic()
  points(inputs, fLogistic(inputs), type = "l", lwd = 2, col = "#55CC99")

  # overlay ReLU()
  points(inputs, ReLU(inputs), type = "l", lwd = 2, col = "#BB0000")

  # Add a legend. Always add legends to plots!
  legend("bottomright",
         legend = c("fTanh()", "fLogistic()", "ReLU()"),
         col = c("#000000", "#55CC99", "#BB0000"),
         lty = 1, lwd = 2, bty = "n")

}


# (For everything that follows below, I will use ReLU() as the default
# activation function. But I will write the code so it is easy to change the
# activation function for experiments.)

# With that, we have everything in place to define individual neurons. Ready for
# the next step: lets connect several neurons together as a network.


# =    3  Computational Neurons in a Network  ==================================

# A neural network is a set of connected neurons. Here is an example:
#
#               +--->[N1]---+
#               |           |
#   >--+--[IN]--+--->[N2]-+-+->[OUT]--+--->     #
#      ^        |         | |         |         # We define a simple network of
#      |        | +<------+ |         |         # computational neurons ...
#      |        | |         |         |         #
#      |        +-+->[N3]---+         |
#      +------------------------------+
#

NN1 <- igraph::make_graph(edges = c("IN",  "N1",
                                    "IN",  "N2",
                                    "IN",  "N3",
                                    "N1",  "OUT",
                                    "N2",  "N3",
                                    "N2",  "OUT",
                                    "N3",  "OUT",
                                    "OUT", "IN"),
                          directed = TRUE)

myLayout <- matrix(c(-1,  0,  # "IN" on the left
                      0,  1,  # "N1" in the center-top
                      0,  0,  # "N2" in the center-middle
                      0, -1,  # "N3" in the center-bottom
                      1,  0), # "OUT" on the right
                    ncol = 2, byrow = TRUE)


if (FALSE) { # plot the network

  plot(NN1,
       layout = myLayout,
       vertex.color = "#DDDDDD",
       vertex.shape = "crectangle",
       vertex.size = 31.5,
       vertex.size2 = 22.5,
       vertex.label.color = "#000000",
       vertex.label.cex = c(1.2, 0.9, 0.9, 0.9, 1.2),
       vertex.label.family = "Helvetica",
       edge.curved = c(rep(0, 7), 0.7))
}


# What happens if we pass a single unit of input into the network?
#
# Assume all weights are the same, and they are 1.0 .
#
# Steps:
# 1. [IN] multiplies it with weight 1, and generates 1 as the output.
# 2. [N1], [N2], and [N3] receive this input, and generate an output of 1.
# 3. [OUT] receives 3 inputs. [N3] receives 1 input.
# 4. [IN] receives the three outputs from [OUT]. [OUT] receives input from [N3].
# 5. The inputs are further circulating through the network. Note that the
#    values get larger and larger.

# To represent this network on a computer, we can list the
# connections that each neuron receives :
# [IN]  <- [OUT]
# [N1]  <- [IN]
# [N2]  <- [IN]
# [N3]  <- [IN], [N2]
# [OUT] <- [N1], [N2], [N3]

# We can then set up nested if-else statements to determine which input
# to multiply with what weight and what output to compute to propagate
# information through the network.


# ==   3.1  Defining a Neural Network with a Weight matrix  ====================

# But there is a simpler way to do this, which is much more flexible and faster
# to work with: we represent the network in a matrix. We keep the weights in
# this matrix. Wherever the weights are zero, this means those neurons have no
# connection. An input that is multiplied with a weight of zero does not
# have an effect.

# Weights don't have to be zero or one. Any non-zero weight will have an effect
# on the output.

# Non-zero weights in a matrix cell mean: the neuron that owns the row of that
# cell is connected to and receives input from the neuron that owns the column.
# Let's construct such a matrix:

NEUNAM <- c("IN", "N1", "N2", "N3", "OUT")   # We define some neuron names
nNeurons <-  length(NEUNAM)                  # and that implies the number
weightMat <- c(
#  IN    N1    N2    N3   OUT
    0,    0,    0,    0,  1.0,  # IN         # Next I just write the weights
  1.0,    0,    0,    0,    0,  # N1         # into a linear vector, with
  1.0,    0,    0,    0,    0,  # N2         # enough cells to form a square
  1.0,    0,  1.0,    0,    0,  # N3         # matrix (with nNeurons^2 cells).
    0,  1.0,  1.0,  1.0,    0)  # OUT

weightMat <- matrix(weightMat,               # We reshape the vector into a 2D
                    nrow = nNeurons,         # square matrix
                    byrow = TRUE)            # row-by-row.
rownames(weightMat) <- NEUNAM                # And we assign the neuron names
colnames(weightMat) <- NEUNAM

vBias <- rep(0.1, nNeurons)                  # We also define a bias vector.
names(vBias) <- NEUNAM                       # You can later experiment with it.


# Lets Visualize this network: we use the igraph package to convert the matrix
# into an adjacency matrix ( a way to represent a mathematical graph, and
# we use the igraph plotting methods to plot the network.)


# ===   3.1.1  Digression: colors

# Never underestimate the importance of color schemes. Colors help to form
# interpretative associations. Using colors well enhances the story your results
# are telling. Poorly chosen colors obscure the story. Inconsistent colors tell
# the wrong story. Taking time and care with your colors is always worth the
# effort!

# Let's compute a consistent set of colors for plotting the network. A green hue
# for input, a red hue for output, and shades of slate grey for intermediate
# neurons. I'll provide a function for this, so we are flexible to change the
# number of neurons if we wish.

neuCols <- function(N,
                    cIN = "#01cb73",
                    cN = c("#274886","#6b89a3","#8caab1"),
                    cOUT = "#d91c68",
                    alpha = "55") {
  # Make a named vector of N colors. The first element is "IN", the last
  # element is OUT, and the intermediate elements are named N1, N2, ...
  # Parameter "alpha" sets the transparency. It is a two-digit hexadecimal
  # value: "ff" is fully opaque, "00 is fully transparent. alpha = 55 is
  # 33% opaque. Using transparent colors in dense plots allows to see
  # overlaps.

  cIN  <- paste0(cIN,  alpha)
  cN   <- paste0(cN,   alpha)
  cOUT <- paste0(cOUT, alpha)
  fColN <- colorRampPalette(cN, alpha = TRUE)
  vCol <- c(cIN, fColN(N - 2), cOUT)
  names(vCol) <- c("IN", paste0("N", 1:(N-2)), "OUT")

  return(vCol)

}

if (FALSE) {
  # check the colors
  barplot(rep(1, 5), col = neuCols(5))
  barplot(rep(1, 9), col = neuCols(9))
}

NEUCOL <- neuCols(nNeurons)


# ==   3.2  Plot the network with igraph::  ====================================

# For future use, we turn plotting our network into a function.

plotNN <- function(wMat){
  # plot a neural network given an adjacency (weight) matrix

  thisNN <- igraph::graph_from_adjacency_matrix(t(wMat),
                                      mode = "directed",
                                      weighted = TRUE)

  # define the layout:
  myLayout <- igraph::layout_in_circle(thisNN, order = c(5,2,1,3,4))
  buffer <- 0.15 * max(abs(myLayout))  # 15% padding

  # define the edge colors
  myWeights <- igraph::E(thisNN)$weight
  wMin <- min(myWeights)
  wMax <- max(myWeights)
  if (wMin == wMax) {
    wNorm <- rep(1, length(myWeights))
    myEdgeLabels <- rep("", length(myWeights))
  } else {
    wNorm <- (myWeights - wMin) / (wMax - wMin)
    myEdgeLabels <- sprintf("%3.2f", myWeights)
  }
  myPal <- colorRampPalette(c("#D9083F", "#878787", "#2CCFB4"))
  myEdgeColors <- myPal(100)[as.numeric(cut(wNorm, breaks = 100))]


  # Reduce margins, make it a single-panel layout, save graphics state
  oPar <- par(mar = rep(1, 4), mfrow = c(1,1))

  plot(thisNN,
       layout = myLayout,
       xlim = range(myLayout[, 1]) + c(-buffer, buffer),
       ylim = range(myLayout[, 2]) + c(-buffer, buffer),
       vertex.size = 40,
       vertex.color = NEUCOL,
       vertex.label.family = "sans",
       vertex.label.cex = 1.6,
       edge.label = myEdgeLabels,
       edge.label.family = "sans",
       edge.label.cex = 0.8,
       edge.curved = 0.1,
       edge.color = myEdgeColors,
       edge.width = 2,
       edge.arrow.size = 0.7,        # Reduce arrowhead size
       edge.arrow.width = 0.8)       # Set arrow mode for a sharper look


  par(oPar)                          # Restore the graphics state


  return(invisible(NULL))            # nothing to return
}

# ... and use this to plot our network:
if (FALSE) {

    plotNN(weightMat)

}

# So far so good, we have a way to represent our neural network, and to plot it.
# Now we can think of running it!


# =    4  Running the network  =================================================

# Running the network means propagating inputs through the network and observing
# what happens. Given the architecture we have defined, that is very easy to do.
# We write a simple function for that - and we make it flexible. For
# example, we pass in the activation function in a parameter, so we can easily
# use different activation functions. Yes, functions in R can be arguments of
# other functions. We simply assign the function name.

if (FALSE) {
  # Check your understanding of the terminology: paste the following into
  # your favorite generative AI.
  t2c("I am new to programming. Please explain the difference between the argument and the parameter of a function. An R code example would be great.")
}


# Take care to not assign the return value,
# but the function! Here is an example of this syntax:

if (FALSE) {

  Sys.time()         # Gives me the current date and time.
  now <- Sys.time()  # "now" is now a time object.
  print(now)         # Printing "now" gives you the time when it was assigned.

  now <- Sys.time    # "now is now a function - Sys.time() by another name.
  print(now)         # Printing it gives me information about the function.
  now()              # But I can use it like a function.

}


iterateNN <- function(v, mat, bias, fAct = ReLU, nrm = TRUE) {
  # Purpose: Iterate the input vector v through the network
  # defined in mat and return the resulting output of all neurons.
  #
  # Parameters:
  #   v:     the input vector, i.e. the output of each neuron after the
  #            previous step. v is converted to a column vector, in case
  #            it isn't one already.
  #   mat:   the matrix of weights of all neuron's inputs. Also implies
  #            the topology of the network.
  #   bias:  biases are added to the weighted sum of inputs - i.e to
  #            the result of the matrix multiplication.
  #   fAct:  the activation function to use: ReLU by default.
  #   nrm:   if TRUE, rescale the result to keep the sum of output values
  #            constant. Only the relative inputs change. This keeps the
  #            internal state of the network commensurate with the constant
  #            levels of the bias vector and the input.

  v <- vc(v)                  # Convert v to a column vector, just in case
  vOut <- mat %*% v           # Multiply weight matrix with input vector
  vOut <- vOut + bias         # Add biases
  vOut <- fAct(vOut)          # Apply the activation function to each neuron
  if(nrm && sum(vOut) != 0) {
    vOut <- vOut / sum(vOut)  # Normalize, to prevent runaway activation
  }
  return(vOut)

}

# Let's give this a try, using the weight matrix we defined above and a
# bias vector of 0.

if (FALSE) {

  vIn <- vc(c(1, 0, 0, 0, 0))  # We pass a column vector with a single
                               # value of 1.0 to the network. This means: [IN]
                               # has an output of 1, all other neurons have
                               # an output of 0.
  # Iterate each neuron for one step
  (vOut <- iterateNN(v = vIn, mat = weightMat, bias = rep(0,5), nrm = FALSE))

  # This returns:

  #      [,1]
  #  IN     0
  #  N1     1
  #  N2     1
  #  N3     1
  #  OUT    0

  # Now we can use this output as the new input, and see how this activation
  # spreads
  vNext <- vOut

  # Select the following line, and execute it a few times. See how the neuron
  # states change. Can you predict the next step after each one? (Refer to the
  # plot.)
  (vNext <- iterateNN(v = vNext, mat = weightMat, bias = rep(0,5), nrm = FALSE))

  # Iteration:    0   1   0   0   0   0
  # ------------------------------------
  #        IN     1   0   0   3   1   0
  #        N1     0   1   0   0   3   1  ... etc.
  #        N2     0   1   0   0   3   1
  #        N3     0   1   1   0   3   4
  #        OUT    0   0   3   1   0   9

  # Don't gloss over this: take your time to understand what's happening here
  # and to see how we are passing input and output through the network, by
  # actually just performing matrix multiplication. Make sure that you see
  # how this happens: it's really one of the key insights for this unit.

}

# At this point we are ready to make our network do some actual work.


# =    5  Input  ===============================================================

# As you have seen above, the values from a single unit of input grow quite
# quickly. This is why we added the option to normalize the output at each step.
# Let's turn normalization back on, iterate for a few cycles and see what
# happens.

if (FALSE) {
  nIter <- 30
  vNext <- vc(c(1, 0, 0, 0, 0))
  out <- vNext
  for (i in 1:nIter) {
    vNext <- iterateNN(v = vNext, mat = weightMat, bias = rep(0, 5))
    out <- cbind(out, vNext)
  }
  plot(0:nIter, out[1,],
       ylim = c(min(out), max(out)),
       ylab="Output", xlab = "iterations",
       type = "l", col = NEUCOL[1])
  for (i in 2:5) {
    points(0:nIter, out[i, ], type = "l", col = NEUCOL[i])
  }
  abline(v = c(10, 14) , col = "#CCDDEE")

  # Let's print the values in the slice between the two vertical lines.
  print(out[ , 10:14])

  # You see the following:
  #  1: An activation pulse cycles through the network along the cyclical
  #       connections we have defined,
  #  2: Over time, the output appears to converge to a steady state.
  #  3: The values for [N1] and [N2] are always the same.

  # Now try running this again for more iterations.
  nIter <- 150
  # ... and select and execute the code above, beginning from
  #   vNext <- vc(...

  # At this steady-state, we have an output that we can feed into the network,
  # and then the network returns the same vector of numbers.

  (vSteady <- vc(out[ , ncol(out)]))
  iterateNN(v = vSteady, mat = weightMat, bias = rep(0,5))

  # Note that this does not mean the network does nothing. Every neuron receives
  # input, and produces output - it is just that the values are perfectly
  # balanced. Which is actually pretty interesting, when you think about it.

  # What happens if you give the network a non-zero bias? Try it!

}


# What we learn is: in this scenario, a single activation excites the network
# and then converges to a non-zero steady state of activations.

# But we're here to explore complexity. Can we get more interesting behaviour?
# What happens when we "feed" an outside signal into the netwoork?

# ==   5.1  Feeding the network: single pulse  =================================

# For convenience, we write our simple single impulse vector into a function:

pulse <- function(N, lev = 1.0) {
  # return a single input at point t = 1, and then zero for the rest of the
  # N time points.
  out <- numeric(N)    # vector pre-filled with zeros
  out[1] <- lev        # first element become "lev"
  return(out)          # return the vector
}


# ==   5.2  Feeding the network: harmonic oscillator  ==========================

# For sustained activity, we can feed a harmonic oscillator into the network
# input. That's just a sine wave over time.

oscillator <- function(i, A = 0.2, f = 1/40, damp = 1.0) {
  # return a sine-function at point t = i with a frequency of f
  val <- A * sin((i - 1) * 2 * pi * f)
  val <- val * damp^i
  return(val)
}

if (FALSE) { # Plot the output of oscillator function

  vFeed <- oscillator(1:200)
  plot(0:(length(vFeed) - 1),
       vFeed,
       type = "b", cex = 0.7,
       xlab = "Iteration", ylab = "Input value",
       main = "Harmonic Oscillator Input")
  abline(h = 0, col = "#AA0000")
  abline(v = c(0, 40, 80, 120, 160, 200), col = "#DDDDFF")
}

# Check and see what happens when you set the damping factor to 0.98

# To use such input to drive our network, we simply add it to the input vector
# in the row of the [IN] neuron.


# =    6  Workbench  ===========================================================

# To experiment more thoroughly, we set up a helper function that takes care of
# the following:
#
# 1: It runs the iterations. For each of N required steps, we add some input
#      into the network and record the output of all nodes.
# 2: It stores the results.
# 3: It creates basic plots on the results.
# 4: It returns the arguments we used and the results in a list, for
#      reference and re-use.


runIterations <- function(vF,                 # input feed vector
                          mat = weightMat,    # the weight matrix
                          bias = vBias,       # the bias vector
                          fAct = ReLU,        # the activation function
                          ...,                #   ... additional parameters
                          doPlot = TRUE,
                          panel3 = "density") {  # plot the results

  # Iterate the values in vFeed over the network specified in mat using
  # "bias" as a bias vector.
  N <- length(vF)                             # number of iterations
  nNodes <- nrow(mat)                         # number of nodes in the network
  vIn <- vc(numeric(nNodes))                  # initialize input: all zeros
  vIn[1,1] <- vF[1]                           # initialize [IN] with vFeed[1]
  results <- matrix(numeric(nNodes * N),      # results matrix: N output vectors
                    ncol = N)                 # ... kept in columns
  results[, 1] <- vIn                         # initialize: first input

  # Now run the actual iteration
  for (i in 2:N) {                            # run for N steps:
    results[,i] <- iterateNN(results[,(i-1)], # Take the previous step's output
                             mat = mat,       # and compute the output for the
                             bias = bias,     # current step.
                             fAct = fAct)

    # Add the value in vF[i] to the input vector for [IN], to be
    # used in the next step.
    results[1, i] <- results[1, i] + vF[i]
  }

  # That's it. The rest is plotting and bookkeeping.

  # Plotting block
  if (doPlot) {

    # Define layout
    lay <- matrix(c(1, 1, 2, 3), nrow = 2)
    layout(lay, widths = c(1, 1), heights = c(1, 1))

    # Plot #1: Neuron states over time
    plot(results[1, ], type = "l",
         ylim = c(min(results), max(results)),
         xlab = "iterations",
         ylab = "neuron output",
         col = NEUCOL[1])

    for (i in 2:(nNodes - 1)) {
      points(results[i, ], type = "l", col = NEUCOL[i])
    }
    points(results[nNodes, ], type = "l", col = NEUCOL[nNodes])

    # Plot #2: Phase space of [IN] and [OUT]
    normIn <- mmScale(results[1, ])
    normOut <- mmScale(results[nNodes, ])
    myCols <- rev(hcl.colors(N, palette = "Red-Blue", alpha = 0.33))
    plot(normIn, normOut,
         main = "Phase Space (trajectory)",
         xlab = "[IN]", ylab = "[OUT]",
         type = "n")

    plotrix::color.scale.lines(normIn, normOut, col = myCols, lwd = 2)
    legend("bottomright",
           legend = c("start", sprintf("%d iterations", N)),
           fill = c(myCols[1], myCols[N]),
           cex = 0.75,
           bty = "n")

    # Panel 3: Conditional plots
    if (panel3 == "density") {
       smoothScatter(normIn, normOut,
                    nbin = 200,
                    nrpoints = Inf,
                    pch = "19", cex = 0.3, col = "#FFFFCC",
                    bandwidth = 0.018,
                    main = "Phase Space (density)",
                    xlab = "[IN]", ylab = "[OUT]",
                    colramp = colorRampPalette(hcl.colors(100, "viridis")))
    } else if (panel3 == "spectrum") {
      stats::spectrum(results[5, ],
                      main = "Spectral density of [OUT]",
                      xlab = "frequency",
                      ylab = "density",
                      col = NEUCOL[nNodes])
    } else {
      stop("Requested mode for \"panel3\" parameter is not recognized.")
    }

  } # end if (doPlot)

  # Create output list
  argList <- list()
  argList$call <- match.call()      # retrieves this run's arguments
  argList$args$NN$weightMat <- mat  # The network's actual weight matrix ...
  argList$args$NN$vBias <- bias     # ... and bias vector
  argList$feed <- vF                # The feed vector
  argList$results <- results        # The results

  return(invisible(argList))
}


# ==   6.1  Running experiments ...  ===========================================


if (FALSE) {   # Experiments

  dev.new()          # set up a nice, large plotting window:

  # Run an iteration of the network. Store the results in a list.

  # 200 steps, default arguments
  outList <- runIterations(vF   = oscillator(1:200),
                           mat  = weightMat,
                           bias = c(0.1, 0.1, 0.1, 0.1, 0.1),
                           fAct = ReLU)

  # You can review the arguments you used:
  outList$call
  outList$args

  # You can also save the list, or print() it in a way that you can
  # store it reproducibly in a text document (like a report):

  dput(list(call = outList$call, args = outList$args)) # Note: we don't really
                                                       # want to include the
                                                       # actual feed or output
                                                       # in such a textual
                                                       # representation.

  # And, we can read the arguments from a file, or a character string after
  # assigning the output, e.g.: newList <- dget(textConnection(txt))


 }


# =    7  "Interesting" Networks  ==============================================

# Now the key question is: what are the properties of neural networks - as
# computational devices - that support "interesting" bevhaviour, i.e.
# behaviour "At the Edge of Chaos"?

# Your task is to use the tools you have here (or construct additional ones)
# to answer this question.

# Here is a brief list of the kind of things we could think of:
#
#  - Can we get continuous cycling output, or does the network always
#    equilibrate?
#  - Can we get truely chaotic behaviour?
#
# - More weights? Fewer weights? Random weights?
# - More cells, fewer cells?
#
# - What is a possible objective function for "interesting"?
# - Different feeds:
#   - pulsed?
#   - white noise?
#   - Fibonacci word (aperiodic)?
#   - ...
# - Can we do computations?
#   - Logic functions?
#   - Period doubling?
#   - Rectifier?
#   - Frequency analyzer?
#   - Noise suppression?
#   - Threshold activation?
#   - Bistable behaviour?
# - ...

# ==   7.1  Random weight matrices  ============================================

# One way to explore the potential behaviour of network topologies is to
# trial random networks. Here is a tool to create random weight matrices.

# rsample() helper function
rsample <- function(n, v = c(0, 1)) {
  # random vector of v with length n. Default to a binary vector.
  return(sample(v, n, replace = TRUE))
}

makeRandomNetwork <- function(nodeNames,
                              topology,
                              func = rnorm,
                              ...,
                              seed = NULL) {
  # Make a random neural network from nodeNames using random network weights
  # which are distributed as per func(). The function accepts a seed to
  # initialize the R RNG (Random Number Generator), so you can get the same
  # matrix (and reproduce your results) by remembering only the seed you
  # used.
  #
  # nodeNames   char      The names of the network nodes, also implies the
  #                         number of nodes.
  # topology    num       A matrix that defines the topology of the network.
  #                         If missing, all weights are kept unchanged. If
  #                         present, all weights are multiplied by the
  #                         corresponding value in the topology matrix
  #                         (Hadamard product).
  #                       Typically the topology matrix would consist of
  #                         1's wherever a node named in the column connects
  #                         to a node named in the row, and 0's elsewhere,
  #                         but non-integral values might also make sense.
  # func        function  Name of the random-value generating  function.
  # ...                   Optional additional arguments passed to func().
  # seed                  Seed for the RNG, for reproducible randomness.
  #                         Default NULL puts the RNG in a fresh state.
  #

  set.seed(seed)                   # Initialize the RNG
  nN <- length(nodeNames)          # number of nodes

  # Call function func() with parameters. The first argument is the square
  # of nN, i.e. the number of values needed to fill the square matrix.
  val <- do.call(func, c(list(n = nN^2), list(...)))

  rMat <- matrix(val, nrow = nN)   # put val into the matrix
  rownames(rMat) <- nodeNames      # name the matrix rows ...
  colnames(rMat) <- nodeNames      # ... and columns

  if (! missing(topology)) {
    rMat <- rMat * topology       # Hadamard product
  }

  return(rMat)
}

# Since makeRandomNetwork() accepts a topology matrix to produce networks
# that are not fully connected, here is the topology for the network we have
# used above:

myTopo <- matrix(c(0, 0, 0, 0, 1,
                   1, 0, 0, 0, 0,
                   1, 0, 0, 0, 0,
                   1, 0, 1, 0, 0,
                   0, 1, 1, 1, 0), nrow = length(NEUNAM), byrow = TRUE)

# When we use this topology mask, we can explore the effects of changing
# weights, while leaving the actual connections unchanged.


# usage examples
if (FALSE) {

  makeRandomNetwork(NEUNAM, func = rsample)  # unit weights
  makeRandomNetwork(NEUNAM, func = runif)    # uniform distribution
  makeRandomNetwork(NEUNAM, func = rnorm)    # normal distribution
  makeRandomNetwork(NEUNAM, func = rsample, v = c(-1, 0, 1))
  makeRandomNetwork(NEUNAM, func = runif, min = -1, max = 1)
  makeRandomNetwork(NEUNAM, func = rnorm, mean = 0.5, sd = 1.5)

  # Reproducible random, uniform distributed weights between -1 and 1, with
  # the topology we defined above:
  makeRandomNetwork(NEUNAM,
                    topology = myTopo,
                    func = runif, min = -1, max = 1,
                    seed = 271828182)

}

# ==   7.2  Putting everything together  =======================================

# Here is one mode of using these tools for explorations:

if (FALSE) {

# select the code between here ...
(mySeed <- sample(1:.Machine$integer.max, 1))  # a large number
outList <- runIterations(vF   = pulse(200),
                         mat  = makeRandomNetwork(NEUNAM,
                                                  topology = myTopo,
                                                  func = rnorm, mean=0.1, sd=2,
                                                  seed = mySeed),
                         bias = c(0.01, 0.01, 0.01, 0.01, 0.01),
                         fAct = ReLU)
# ... and here.

# Then execute it with <command><return> again and again. Observe the plots. If
# you find an interesting combination of values, copy the last seed from
# the console and note it down. This allows you to reproduce your network at
# any time later.
#
# And inspect your network:

plotNN(outList$args$NN$weightMat)

# For an example of something "interesting", you could try:
# mySeed <- 423686718
#
# Note: not all random networks will work - some parameter values will cause
# errors. Don't worry about that. Just try again.

}


# Have fun!



# [END]
