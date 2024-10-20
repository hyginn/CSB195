# tocID <- "courseScripts/07-Optimization.R"
#
# Purpose: Optimizing objects with an objective function
#
# Version: 1.0
# Date:    2024-10-18
# Author:  boris.steipe@utoronto.ca; ChatGPT-4 and 4o
#
# Versions:
#   1.0    New script discussing optimization and a worked example for
#          optimizing a linear ordering of the standard amino acids.
#
# Preconditions:
#
# ToDo:
#
# Notes:
#
# ==============================================================================
#                                                                              #
#   THIS SCRIPT IS "SAFE TO SOURCE" WITHOUT SIDE EFFECTS.                      #
#                                                                              #
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                  Line
#TOC> ------------------------------------------------------
#TOC>   1        RATIONALE                                43
#TOC>   2        ALGORITHM                                72
#TOC>   3        IMPLEMENTATION                          113
#TOC>   3.1        State Representation                  115
#TOC>   3.2        State Modification                    121
#TOC>   3.3        Objective function                    138
#TOC>   3.4        Simulated Annealing Function          175
#TOC>   4        RUN THE OPTIMIZATION                    264
#TOC>   5        SOME TECHNICAL CONSIDERATIONS           356
#TOC> 
#TOC> ==========================================================================


# =    1  RATIONALE  ===========================================================

# There are many approaches to optimization ... Let's discuss.
gAIinit()

# Task: Navigate to your AI tutor of choice and paste the initialization prompt.
# Then continue:

# Prompts:
#
# ---
#    Could you summarize for me different concepts for optimization and mention
#    what they are most useful for and some pros and cons?
# ---

#    I would like to run a simulated annealing of a linear ordering of the
#    twenty amino acid one-letter codes. I have a function
#    aaSim(aa1, aa2) that returns the distance in feature space of two
#    amino acids aa1 and aa2; amino acid pairs with a smaller distance are
#    considered more similar. The linear ordering I would like, minimizes the
#    sum of similarities of all twenty letters, where adjacent neighbours
#    each count for 0.31, neighbours separated by one each account for 0.125,
#    and neighbours separated by two amino acids each account for  0.065
#    of the total. Also, I would like to consider the ordering to have
#    periodic boundaries.
#
#   Can you design how this  can be done?


# =    2  ALGORITHM  ===========================================================

# ChatGPT's response is a good example of the essential components of a
# simulated annealing optimization algorithm. (Edited for clarity):

# Here are the necessary components:

#   State Representation:
#     The state of our system will be a permutation (ordering) of the 20 amino
#     acid one-letter codes. To optimize our system, we need to represent it in
#     the computer. We'll store the state as a vector of characters and call
#     it our State Vector.

#   Objective Function:
#     This defines the objective of our optimization. Here, we want to
#     minimize the sum of the similarities, while accounting for neighbors at
#     different distances. Thus the objective function will compute the
#     weighted sum of distances for adjacent pairs, pairs separated by one
#     position, and pairs separated by two positions. We'll also treat the
#     ordering as having periodic boundaries, meaning the last and first
#     position of our State Vector will also be considered neighbors.
#
#   State Modification:
#     We will modify the current State Vector in small steps to evaluate the
#     result with our objective function. The simplest approach is to swap the
#     position of letters.
#
#   Simulated Annealing Algorithm (SA):
#     SA itself involves initializing the State Vector with a random ordering
#     and then iterating for many steps of modification, evaluation, and
#     acceptance or rejection of the modification. If the modified state is
#     better than the current state, we always accept it. If it is worse
#     than the current state, we accept it with a small probability that
#     depends on how much worse it is, and a "temperature". We slowly lower
#     the temperature over the course of the iterations. This allows the system
#     to explore its state space broadly at the beginning, and find a minimum
#     as the optimization progresses.
#

# With this, we can start implementing.

# =    3  IMPLEMENTATION  ======================================================

# ==   3.1  State Representation  ==============================================

aaAlphabet <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")


# ==   3.2  State Modification  ================================================

modifyState <- function(state) {
  # Modify the input state by swapping the letters in two random positions

  # Select two distinct positions to swap
  idxSwap <- sample(1:length(state), 2)

  # Swap the contents of the two positions
  x <- state[idxSwap[1]]
  state[idxSwap[1]] <- state[idxSwap[2]]
  state[idxSwap[2]] <- x

  return(state)
}


# ==   3.3  Objective function  ================================================

evalObjective <- function(state) {
  # Calculate the weighted sum of distances for a given state of amino
  # acid letters

  n <- length(state)

  # define some weights decreasing in the golden ratio
  phi <- (1 + sqrt(5)) / 2
  weight    <- 1 / phi                    # For adjacent neighbors
  weight[2] <- (1 - weight[1]) / phi      # For neighbors separated by 1
  weight[3] <- 1 - weight[1] - weight[2]  # For neighbors separated by 2

  sumDist <- 0

  for (i in 1:n) {   # Loop through each position

    # Neighbors at distance 1 (wrap around using modulo arithmetic)
    sumDist <- sumDist +
      weight[1] * (aaSim(state[i], state[ (i      %% n) + 1]) +
                   aaSim(state[i], state[((i - 2) %% n) + 1]))

    # Neighbors at distance 2
    sumDist <- sumDist +
      weight[2] * (aaSim(state[i], state[((i + 1) %% n) + 1]) +
                   aaSim(state[i], state[((i - 3) %% n) + 1]))

    # Neighbors at distance 3
    sumDist <- sumDist +
      weight[3] * (aaSim(state[i], state[((i + 2) %% n) + 1]) +
                   aaSim(state[i], state[((i - 4) %% n) + 1]))
  }

  return(sumDist)
}

# ==   3.4  Simulated Annealing Function  ======================================

# A Simulated Annealing function with Progress Bar and Score Tracking
simulatedAnnealing <- function(state = sample(aaAlphabet),
                               maxIterations = 1000,
                               initialTemp = 100,
                               coolingRate = 0.99,
                               showPlot = TRUE,
                               showProgress = TRUE) {
  # Initial setup
  currentState <- state
  currentScore <- evalObjective(currentState)

  bestState <- currentState
  bestScore <- currentScore
  temp <- initialTemp

  # Vector to track scores over time
  scoreHistory <- numeric(maxIterations)

  if (showProgress) {
    # Set up a progress bar
    pb <- utils::txtProgressBar(min = 0, max = maxIterations, style = 3)
  }

  # Iteration loop
  for (i in 1:maxIterations) {
    # Modify the current state
    newState <- modifyState(currentState)
    newScore <- evalObjective(newState)

    # Check if the new solution is better or accept with a temperature-
    # dependent probability
    if (newScore < currentScore)  { # Better, therefore accept
      bestState <- newState
      bestScore <- newScore
      acceptThisState <- TRUE
    } else if (runif(1) < exp((currentScore - newScore) / temp)) {
      # Worse, but acceptable at this temperature
      acceptThisState <- TRUE
    } else {
      acceptThisState <- FALSE
    }

    if (acceptThisState) {             # Update new best state
      currentState <- newState
      currentScore <- newScore
    }
    scoreHistory[i] <- currentScore

    # Cool down the temperature
    temp <- temp * coolingRate

    if (showProgress) {
      # Update progress bar
      utils::setTxtProgressBar(pb, i)
    }
  }

  # Done iterations.
  if (showProgress) {
    # Close progress bar
    close(pb)
  }

  # Result: Plot scores of accepted states over time
  if (showPlot) {
    plot(scoreHistory,
         type = "l",
         main = "Simulated Annealing Score Improvement",
         xlab = "Iteration",
         ylab = "Objective Score",
         col  = "#8c3861")
  }

  # Rotate best solution to always begin with "A"
  idxA <- which(bestState == "A")
  if (idxA != 1) {
    bestState <- c(bestState[idxA:length(bestState)],
                   bestState[1:(idxA - 1)])
  }

  # Return the results
  return(list(bestState = bestState,
              bestScore = bestScore,
              scoreHistory = scoreHistory))
}


# =    4  RUN THE OPTIMIZATION  ================================================

# That's it for the definitions. Now we can run the actual optimization. I put
# this code in a "guard block" if (FALSE) { ... }, that does not get executed if
# you source() the file. source()'ing the file thus has no side effects, other
# than defining the functions. But of course you can select and execute the
# code from the script itself.

if (FALSE) {

  # Try this a few times to see how the optimization goes:

  x <- simulatedAnnealing(maxIterations = 10000, coolingRate = 0.9991)
  x$bestState
  x$bestScore


  # Then run it 100 times to get a sampling of possible results
  myRes <- list()
  myMin <- 999999
  for (i in 1:100) {
    myRes[[i]] <- simulatedAnnealing(maxIterations = 10000,
                                     coolingRate = 0.9991,
                                     showPlot = FALSE,
                                     showProgress = FALSE)
    if (myRes[[i]]$bestScore < myMin) {
      myMin <- myRes[[i]]$bestScore
    }
    cat(sprintf("i: %3i, score: %3.8f solution: %s\n", i,
                myRes[[i]]$bestScore,
                paste(myRes[[i]]$bestState, collapse = "")))
  }

  # When you run the optimization many times, you should see a few varying non-
  # optimal solutions, where the optimization gets stuck in a local minimum,
  # and several instances of the same optimal solution with a score of
  # 51.94000448. Since this solution is reached from a large number of starting
  # states taken from 20! = 2.4e+18 distinct permutations, we can have some
  # confidence that this is a global optimal solution (in the context of
  # the given weights and similarity function). However if you look closely,
  # you will notice that these are actually a family of  similar solutions
  # with almost indistinguishable scores, that are less than 1e-14 different.
  # Thus for all practical purposes, there is not a single global optimum
  # for this problem. The lowest score I find is :
  #
  # $bestScore
  # 51.94000447780499740702
  #
  # $bestState
  # "A" "C" "V" "L" "I" "F" "M" "W" "Y" "H"
  # "R" "K" "E" "D" "Q" "N" "S" "P" "T" "G"
  #
  # Let me know in case you find a better one  :-)


  # Here are just a few ways to visualize the results.

  # We can extract some of the near optimal solutions like so:
  myEps <- 1e-14

  for (i in 1:length(myRes)) {
    if (abs(myRes[[i]]$bestScore - myMin) < myEps) {
      cat(sprintf("i: %3i, score: %3.8f solution: %s\n", i,
                  myRes[[i]]$bestScore,
                  paste(myRes[[i]]$bestState, collapse = "")))
    }
  }

  # ... and we can plot the ordering like so:
  for (i in 1:length(myRes)) {                  # retrieve the best solution
    if (myRes[[i]]$bestScore == myMin) {
      bestIdx <- i
    }
  }
  bestState <- myRes[[bestIdx]]$bestState

  oPar <- par(mar=rep(0, 4))                    # turn off plotting margins
  plot(c(0, 20), c(0, 20),                      # empty, square plot window
       type = "n",
       axes = FALSE,
       xlab = "", ylab = "",
       asp = 1)
  for (i in 1:20) {                             # plot coloured rectangles
    rect(i - 0.4, 9.6, i + 0.4, 10.4,
         col = AACOLS[bestState[i]])
  }
  text(1:20, rep(9.0, 20), labels = bestState)  # write AA codes

  par(oPar)                                     # restore plotting margins

}

# =    5  SOME TECHNICAL CONSIDERATIONS  =======================================

# If you ever try your own optimizations, the following notes may be helpful:


# Choosing the right parameters for simulated annealing, especially the starting
# and ending temperatures (T_start and T_end), can be problem-specific and may
# require some experimentation. However, there are some general guidelines and
# considerations that might help you find a good starting point:

# Starting Temperature (T_start)
# ==============================
# High enough to accept worse solutions: Initially, you want T_start to be high
# enough that the algorithm is likely to accept worse solutions with reasonable
# probability, which helps to explore the solution space. You might set it so
# that exp(-ΔE/T_start) is reasonably high (e.g., > 0.5) for typical values of
# ΔE (the difference in objective function values between successive solutions).

# Adjustment:
# It may be worth experimenting with a few different starting temperatures
# to see how it affects the search and adjust based on your observations.

# Ending Temperature (T_end)
# ==========================
# You could define a maximal number of steps, or use an ending temperature at
# which you would consider your system "frozen". An ending temperature should be
# low enough to refine solutions, i.e. low enough that towards the end of the
# search, the algorithm is primarily exploring the areas of the solution space
# it has found to be good, and is unlikely to accept worse solutions. In
# practice, you might set T_end so that exp(-ΔE/T_end) is quite low (e.g.
# < 0.05) for the smallest ΔE values you expect to encounter.

# Be wary of setting T_end too high, as it might result in premature convergence.

# A Practical Approach
# ====================

# Rule of thumb: Sometimes, researchers choose T_start and T_end based on their
# knowledge of the problem domain or via preliminary experiments.

# Automatic tuning: In some cases, automatic parameter tuning methods or
# adaptive simulated annealing approaches can be used, where the temperature
# schedule is adjusted dynamically based on the search's progress.


# [END]
