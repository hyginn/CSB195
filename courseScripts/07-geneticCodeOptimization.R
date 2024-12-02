# tocID <- "courseScripts/07-geneticCodeOptimization.R"
#
# Purpose: Optimizing the Genetic Code
#
# Version: 1.0
# Date:    2024-10-09
# Author:  boris.steipe@utoronto.ca; ChatGPT-4 and 4o
#
# Versions:
#   1.0    In tutorial - Tutorial session 05: new script taking the
#          insights from 05-geneticCodeExperiments.R further to arrive at
#          better-than-natural codes.
#
# Preconditions:
#    AADAT and GCdf are defined via .util.R, functions aaSim(), neighCodons()
#    and rGC() are available.
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
#TOC>   Section  Title                    Line
#TOC> ----------------------------------------
#TOC>   1        INITIALIZATIONS            42
#TOC>   1.1        Parameters               44
#TOC>   1.2        Packages                 47
#TOC>   2        RATIONALE                  50
#TOC>   3        PROCESS                    71
#TOC>
#TOC> ==========================================================================


# =    1  INITIALIZATIONS  =====================================================

# ==   1.1  Parameters  ========================================================


# ==   1.2  Packages  ==========================================================


# =    2  RATIONALE  ===========================================================

# We have previously constructed code that evaluates the "quality" of a non-
# natural genetic code. We have also seen that the standard genetic code is
# much better than the randomized codes. But how could nature have found such a
# highly performant code, if not by random chance? And can we do even better?
# can we do even better? Let's try our hand at optimizing genetic codes.

# There are many approaches to optimization ... Let's discuss.
gAIinit()

# Task: Navigate to your AI tutor of choice and paste the initialization prompt.
# Then continue:

# Prompt:
#    I would like to optimize the genetic code under some objective function.
#    What options for such a multidimensional optimization problem do I have
#    in principle?
# --

# We usually get a list of options and you can have the AI explain the conceptual background, and the strengths and weaknesses of each one separately:
# -  Gradient descent
# -  Genetic algorithms
# -  Simulated annealing
# -  Monte Carlo sampling
# -  Bayesian Optimization
# -  Swarm intelligence

# Prompt:
#    (1) My genetic code is a vector named with the 64 codons, each element
#    represents one of the 20 amino acids in one-letter code or a stop
#    codon "*".
#    (2) I have an objective function qGC() that takes as input
#    two amino acids (or stop codons) and outputs a quality measure. For
#    example  qGC("A", "Q"). This should be minimized.
#    (3) I also have a function neighCodons() that generates the nine
#    neighbouring codons of the input.
#    (4) I want to allow two types of changes to the code:
#      (a) the identity of two amino acids (but not stop codons) are swapped.
#      (b) the position of two amino acids in the code are swapped.
#    Together this guarantees that the redundancy (information) of the code
#    does not change, and that that the number of stop codons stays the same.
#
#    Given this context, which optimization algorithm would you choose?
# --

modifyCode <- function(code, identitySwapProb = 0.5) {
  # Randomly choose one of the two modification rules based on identitySwapProb
  rule <- sample(c("identitySwap", "positionSwap"), 1, prob = c(identitySwapProb, 1 - identitySwapProb))

  if (rule == "identitySwap") {
    # Get all indices of amino acids (excluding stop codons)
    aminoAcidIndices <- which(code != "*")

    # Rule (a): Swap the identities of two amino acids
    swapIndices <- sample(aminoAcidIndices, 2)  # Choose two distinct amino acids
    temp <- code[swapIndices[1]]
    code[swapIndices[1]] <- code[swapIndices[2]]
    code[swapIndices[2]] <- temp
  } else {
    # Rule (b): Swap positions of two amino acids in the code
    swapPositions <- sample(1:length(code), 2)  # Choose any two positions in the code
    temp <- code[swapPositions[1]]
    code[swapPositions[1]] <- code[swapPositions[2]]
    code[swapPositions[2]] <- temp
  }

  return(code)
}




validateModifyCode <- function(originalCode, modifiedCode) {
  # 1. Check that the number of stop codons remains the same
  stopCodonCheck <- sum(originalCode == "*") == sum(modifiedCode == "*")
  if (!stopCodonCheck) {
    return("Fail: Stop codon count changed.")
  }

  # 2. Identify positions where the codes differ
  changedPositions <- which(originalCode != modifiedCode)

  # 3. If exactly two positions are different, it's a position swap
  if (length(changedPositions) == 2) {
    if (originalCode[changedPositions[1]] == "*" || originalCode[changedPositions[2]] == "*") {
      return("Fail: Stop codons were swapped.")
    } else {
      return("Pass: Valid position swap.")
    }
  }

  # 4. For more than two changes, check for a valid identity swap
  if (length(changedPositions) > 2) {
    # Get unique elements in the changed positions
    originalElements <- unique(originalCode[changedPositions])
    modifiedElements <- unique(modifiedCode[changedPositions])

    # There should be exactly two unique amino acids swapped
    if (length(originalElements) == 2 && length(modifiedElements) == 2 && all(originalElements %in% modifiedElements)) {
      return("Pass: Valid identity swap.")
    } else {
      return("Fail: Invalid identity swap pattern.")
    }
  }

  # Catch any other unexpected outcomes
  return("Fail: Unexpected modification pattern.")
}

# [...]

testModifyCode <- function(originalCode, identitySwapProb, numTrials = 10) {
  stopCodonCount <- sum(originalCode == "*")  # Initial number of stop codons

  # Run multiple trials
  for (trial in 1:numTrials) {
    modifiedCode <- modifyCode(originalCode, identitySwapProb)

    # Check if the number of stop codons remains the same
    if (sum(modifiedCode == "*") != stopCodonCount) {
      cat("Test failed at trial:", trial, " - Stop codon count mismatch!\n")
      return(FALSE)
    }

    # Print some trials to show changes
    if (trial <= 3) {
      cat("Original Code:\n", originalCode, "\nModified Code:\n", modifiedCode, "\n---\n")
    }
  }

  cat("All tests passed for identitySwapProb =", identitySwapProb, "\n")
  return(TRUE)
}

# Example usage with a starting code
myGC <- rGC(seed = 112358)

# Run separate tests for each rule
cat("Testing Identity Swap Rule:\n")
testModifyCode(myGC, identitySwapProb = 1.0)

cat("Testing Position Swap Rule:\n")
testModifyCode(myGC, identitySwapProb = 0.0)

# Note: a better test would have been to directly validate the resulting code
# using table()

validateModifyCode <- function(originalCode, identitySwapProb, numTrials = 10) {
  originalTable <- table(originalCode)  # Get frequency counts of the original code

  for (trial in 1:numTrials) {
    modifiedCode <- modifyCode(originalCode, identitySwapProb)
    modifiedTable <- table(modifiedCode)  # Get frequency counts of the modified code

    # Compare frequency counts for amino acids (excluding stop codons)
    if (!all(modifiedTable[names(modifiedTable) != "*"] == originalTable[names(originalTable) != "*"])) {
      cat("Test failed at trial:", trial, " - Frequency mismatch!\n")
      cat("Original frequencies:\n", originalTable, "\nModified frequencies:\n", modifiedTable, "\n")
      return(FALSE)
    }

    # Print sample results for inspection
    if (trial <= 3) {
      cat("Original Code:\n", originalCode, "\nModified Code:\n", modifiedCode, "\n---\n")
    }
  }

  cat("All tests passed for identitySwapProb =", identitySwapProb, "\n")
  return(TRUE)
}

# Example usage with a starting code
myGC <- rGC(seed = 112358)

# Run separate tests for each rule
cat("Validating Identity Swap Rule:\n")
validateModifyCode(myGC, identitySwapProb = 1.0)

cat("Validating Position Swap Rule:\n")
validateModifyCode(myGC, identitySwapProb = 0.0)


# [...]

simulatedAnnealing <- function(initialCode, initialTemp = 100, maxIter = 1000, identitySwapProb = 0.5) {
  # Initialize variables
  currentCode <- initialCode
  currentEnergy <- energyFunction(currentCode)
  bestCode <- currentCode
  bestEnergy <- currentEnergy
  energyHistory <- numeric(maxIter)  # Store energy at each iteration

  # Print progress bar setup
  cat("Simulated Annealing Progress:\n")
  pb <- txtProgressBar(min = 0, max = maxIter, style = 3)  # Create a progress bar

  # Simulated Annealing main loop
  for (iteration in 1:maxIter) {
    # Update the temperature using the cooling schedule
    temperature <- coolingSchedule(iteration, initialTemp)

    # Generate a new code by making a random modification
    newCode <- modifyCode(currentCode, identitySwapProb)
    newEnergy <- energyFunction(newCode)

    # Determine if the new code should be accepted
    if (newEnergy < currentEnergy || runif(1) < acceptanceProbability(currentEnergy, newEnergy, temperature)) {
      currentCode <- newCode
      currentEnergy <- newEnergy
    }

    # Update best solution found so far
    if (currentEnergy < bestEnergy) {
      bestCode <- currentCode
      bestEnergy <- currentEnergy
    }

    # Store energy for plotting
    energyHistory[iteration] <- currentEnergy

    # Update progress bar
    setTxtProgressBar(pb, iteration)
  }
  close(pb)  # Close the progress bar

  # Plot the energy history
  plot(1:maxIter, energyHistory, type = "l", col = "blue", lwd = 2,
       xlab = "Iteration", ylab = "Energy (Quality Measure)",
       main = "Simulated Annealing Optimization")

  cat("Optimization Complete!\n")
  cat("Best energy found:", bestEnergy, "\n")

  return(list(bestCode = bestCode, bestEnergy = bestEnergy))
}

# Supporting functions

# 1. Energy function: Calculate energy of the genetic code
energyFunction <- function(code) {
  energy <- 0
  codonList <- names(code)

  # Loop through all codons and calculate the quality for neighboring codons
  for (i in 1:length(codonList)) {
    neighbors <- neighCodons(codonList[i])  # Get neighboring codons for the current codon
    for (neighbor in neighbors) {
      if (neighbor %in% codonList) {
        energy <- energy + qGC(code[codonList[i]], code[neighbor])  # Calculate pairwise energy
      }
    }
  }

  return(energy)
}

# 2. Acceptance probability function
acceptanceProbability <- function(oldEnergy, newEnergy, temperature) {
  exp(-(newEnergy - oldEnergy) / temperature)  # Metropolis acceptance criterion
}

# 3. Cooling schedule function
coolingSchedule <- function(iteration, initialTemp) {
  initialTemp / (1 + log(1 + iteration))  # Example logarithmic cooling schedule
}

# Example usage with a starting code:
qGC <- aaSim
myGC <- rGC(seed = 112358)

result <- simulatedAnnealing(myGC, initialTemp = 100, maxIter = 500, identitySwapProb = 0.7)



# ```

### Key Features:
1. **Progress Bar**: `txtProgressBar()` is used to show the progress of iterations.
2. **Energy Tracking**: `energyHistory` stores the energy value at each iteration, which is plotted after the optimization completes.
3. **Cooling Schedule**: `coolingSchedule()` function controls the temperature to ensure gradual exploration.
4. **Acceptance Probability**: Uses the Metropolis criterion to decide whether to accept a new code.

### Example Usage:
1. **Define `myGC`** using your custom `rGC()` function.
2. Run `simulatedAnnealing(myGC, initialTemp = 100, maxIter = 500, identitySwapProb = 0.7)` to test it.

This structure should provide a good foundation for further tuning. Let me know if you would like to discuss any specific parts or add features!#


# (running 500000 iterations: best energy 1304.83)

# Great - that runs well. The optimization is quite noisy. Can we modify it to use the more "invasive" identity swaps early on, and then later favour only position swaps?

  simulatedAnnealing <- function(initialCode, initialTemp = 100, maxIter = 1000) {
    # Initialize variables
    currentCode <- initialCode
    currentEnergy <- energyFunction(currentCode)
    bestCode <- currentCode
    bestEnergy <- currentEnergy
    energyHistory <- numeric(maxIter)  # Store energy at each iteration

    # Print progress bar setup
    cat("Simulated Annealing Progress:\n")
    pb <- txtProgressBar(min = 0, max = maxIter, style = 3)  # Create a progress bar

    # Simulated Annealing main loop
    for (iteration in 1:maxIter) {
      # Update the temperature using the cooling schedule
      temperature <- coolingSchedule(iteration, initialTemp)

      # Adjust identitySwap probability based on the current iteration
      identitySwapProb <- swapProbabilitySchedule(iteration, maxIter)

      # Generate a new code by making a random modification
      newCode <- modifyCode(currentCode, identitySwapProb)
      newEnergy <- energyFunction(newCode)

      # Determine if the new code should be accepted
      if (newEnergy < currentEnergy || runif(1) < acceptanceProbability(currentEnergy, newEnergy, temperature)) {
        currentCode <- newCode
        currentEnergy <- newEnergy
      }

      # Update best solution found so far
      if (currentEnergy < bestEnergy) {
        bestCode <- currentCode
        bestEnergy <- currentEnergy
      }

      # Store energy for plotting
      energyHistory[iteration] <- currentEnergy

      # Update progress bar
      setTxtProgressBar(pb, iteration)
    }
    close(pb)  # Close the progress bar

    # Plot the energy history
    plot(1:maxIter, energyHistory, type = "l", col = "blue", lwd = 2,
         xlab = "Iteration", ylab = "Energy (Quality Measure)",
         main = "Simulated Annealing Optimization with Dynamic Swap Probability")

    cat("Optimization Complete!\n")
    cat("Best energy found:", bestEnergy, "\n")

    return(list(bestCode = bestCode, bestEnergy = bestEnergy))
  }

# Function to adjust the identitySwap probability dynamically
swapProbabilitySchedule <- function(iteration, maxIter) {
  # Start with high identity swap probability and gradually reduce it
  startProb <- 0.9  # Initial probability for identity swaps
  endProb <- 0.0    # Final probability for identity swaps

  # Linear interpolation: Decrease from startProb to endProb over iterations
  prob <- startProb - (startProb - endProb) * (iteration / maxIter)
  return(prob)
}

# Supporting functions (unchanged)
energyFunction <- function(code) {
  energy <- 0
  codonList <- names(code)
  for (i in 1:length(codonList)) {
    neighbors <- neighCodons(codonList[i])
    for (neighbor in neighbors) {
      if (neighbor %in% codonList) {
        energy <- energy + qGC(code[codonList[i]], code[neighbor])
      }
    }
  }
  return(energy)
}

acceptanceProbability <- function(oldEnergy, newEnergy, temperature) {
  exp(-(newEnergy - oldEnergy) / temperature)
}

coolingSchedule <- function(iteration, initialTemp) {
  initialTemp / (1 + log(1 + iteration))
}

# run
myGC <- rGC(seed = 112358)
result <- simulatedAnnealing(myGC, initialTemp = 100, maxIter = 50000)






# =    3  PROCESS  =============================================================


# Prompt:
#    This is code I have for my objective function:

qGC <- function(GC) {

  sumDist <- 0                      # Initialize the objective variable

  for (codonX in names(GC)) {       # For each codon in the code
    aaX <- GC[codonX]               # ... get the encoded amino acid

    for (codonY in neighCodons(codonX)) { # For all nine neighbors of the codon
      aaY <- GC[codonY]             # ... get the encoded amino acid
      dist <- aaSim(aaX, aaY)       # ... compute distance in feature space
      sumDist <- sumDist + dist     # ... add to the sum of distances
    }
  }
  return(sumDist)
}




#    Also, I want to optimize according to the following two constraints:
#      (a) the number of stop codons "*" should not change;
#      (b) the level of redundancy in the code should not change
#    Given that, which of the optimization options is the most straightforward
#    to implement, and how would you go about writing the code?
# --

generate_neighbor <- function(GC) {
  # Randomly select two codons to swap their encoded amino acids
  codon_pair <- sample(names (GC), 2)
  # Swap their associated amino acids
  new GC < - GC
  new_GC[codon_pair[1]] < - GC[codon_pair[2]] new_GC[codon_pair[2]] <- GC [codon_pair[1]]
  return(new GC)
}





# Note: The most important part of the task is to ensure that the produced
# code is logically sound. We would also like to visualize the progress of
# the optimization, and, of course, to be able to see the resulting code.


# [END]
