# tocID <- "./R/gameOfLife.R"
#
# Demo code
#
# 2022-10  -  2023-10
#
# Boris  Steipe (boris.steipe@utoronto.ca)
#
#
# Version:  1.2
#
# Versions:
#           1.2   Changed logic of the next-step computation to support non-
#                 integral rules. Added programming exercises.
#           1.1.1 Moved src and pattern files to ./R/ resp. ./data/ directories
#           1.1   Updated, expanded and added examples section
#           1.0   In class demo 2022
#
# ==============================================================================
#


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                             Line
#TOC> -------------------------------------------------
#TOC>   1        Structures and Functions            44
#TOC>   1.1        wrap()                            50
#TOC>   1.2        getFate()                         63
#TOC>   1.3        runGOL()                          88
#TOC>   2        Explorations                       144
#TOC>   2.1        The Glider                       167
#TOC>   2.2        A Glider-gun                     178
#TOC>   2.3        Max Fill                         190
#TOC>   2.4        Pentomino                        205
#TOC>   3        Next ?                             216
#TOC>
#TOC> ==========================================================================


# This is a demo implementation of John Horton Conway's "Game of Life".
# cf. https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life

# =    1  Structures and Functions  ============================================

source("R/imPlot.R")
source("R/seedGOL.R")


# ==   1.1  wrap()  ============================================================

wrap <- function(w) {
# wrap a 2D matrix of size x*y into a torus (x-2) * (y-2)
#
  w[1, ]       <- w[nrow(w)-1, ]    # wrap last actual row to row 1
  w[nrow(w), ] <- w[2,         ]    # wrap row 2 to last row
  w[ ,1]       <- w[ ,ncol(w)-1]    # wrap last actual column to column 1
  w[ ,ncol(w)] <- w[ ,2        ]    # wrap column 2 to last column

  return(w)
}

# ==   1.2  getFate()  =========================================================

getFate <- function(cellState, cellNeigh, rL, rD) {
  # determine the next state of a cell in cellState given the number of
  # cell neighbors, and the rule set for live and dead cells.

  if (cellState == 1) {
    ruleSet <- rL
  } else {
    ruleSet <- rD
  }

  pAlive <- ruleSet[cellNeigh + 1]  # get the probability of a cell being
                                    # alive in the next step, given the number
                                    # of neighbors it has and the rule set.
  if (runif(1) < pAlive) {
    nextState <- 1
  } else {
    nextState <- 0
  }
  return(nextState)

}


# ==   1.3  runGOL()  ==========================================================
runGOL <- function(w0, rL, rD, nTick = 100, drawGrid = TRUE) {

  # "0" means: "dead in the next step.
  # "1" means: "alive" in the next step.
  if (missing(rL)) { # Rules what to do with "live" cells
    # Neighbours:    0  1  2  3  4  5  6  7  8
    rL <-          c(0, 0, 1, 1, 0, 0, 0, 0, 0)
  }
  if (missing(rD)) { # Rules what to do with "dead" cells
    # Neighbours:    0  1  2  3  4  5  6  7  8
    rD <-          c(0, 0, 0, 1, 0, 0, 0, 0, 0)
  }

  w0 <- wrap(w0)  # w0: the current "world"
  w1 <- w0        # w1: the "world" of the next step

  for (tick in seq_len(nTick)) {    # for each tick of the universe clock
    pBar(tick, nTick)               # progress bar

    # compute the world of the next step
    for (i in 2:(nrow(w0)-1)) {     # for each actual row
      for (j in 2:(ncol(w0)-1)) {   # for each actual column

        # Compute the number of neighboring cells: the sum of live cells in
        # the 3 x 3 neighbourhood minus the value of the center cell.
        nNeigh <- sum(w0[(i-1):(i+1), (j-1):(j+1)]) - w0[i, j]

        # Apply the rules for the next step
        w1[i, j] <- getFate(cellState = w0[i, j],
                            cellNeigh = nNeigh,
                            rL = rL,
                            rD = rD)
      } # end column
    } # end row

    # wrap contents of world edge into virtual cells
    w1 <- wrap(w1)

    w0 <- w1  # update the world ...
    # ... and plot the new world state.
    imPlot(w0[2:(nrow(w1)-1), 2:(ncol(w1)-1)], drawGrid = drawGrid)

    safeSleep(0.1)
    dev.flush()

  } # end tick

  return(invisible(w1))
}




if (FALSE) {

# =    2  Explorations  ========================================================

  graphics.off()    # close all open graphics devices
  dev.new()         # re-open the standard device
  dev.new()         # open a second, floating graphics window

  # Initialize dimensions and probabilities
  X <- 100
  Y <- 100
  p <- 0.5

  # Fill and display a starting state ...
  myWorld <- sample(0:1, X * Y, prob=c(1-p, p), replace = TRUE)
  dim(myWorld) <- c(Y,X)
  imPlot(myWorld)

  # Run this evolution for maximally 1000 steps. Interrupt it with the red
  # STOP sign if needed ... will it freeze earlier? Undecided.
  runGOL(myWorld, nTick = 1000)

  # Try changing the starting probability to 0.1, or 0.9. Surprising?


# ==   2.1  The Glider  ========================================================
  # People have been searching for interesting patterns for a long time. This
  # one is a "glider", it can be used to send information around ...

  X <- 50; Y <- 50; myWorld <- matrix(integer(X * Y), ncol = X)
  myWorld[25:27,27]<-1;myWorld[26,25]<-1;myWorld[27,26]<-1   # glider
  imPlot(myWorld)

  runGOL(myWorld, nTick = 1000)


# ==   2.2  A Glider-gun  ======================================================
  # This "glider gun", spawns gliders. Regularly, like clockwork. Or like the
  # clock cycles of a computer. But this world is small. What goes around comes
  # around ... and around  ... and around

  X <- 57; Y <- 50; myWorld <- matrix(integer(X * Y), ncol = X)
  myWorld <- placeCell("data/GosperGun.cell", myWorld, 5, 20)
  imPlot(myWorld, drawGrid = FALSE)

  runGOL(myWorld, nTick = 1000, drawGrid = FALSE)


# ==   2.3  Max Fill  ==========================================================
  # Meet Max. Max is hungry. He devours the void, filling it with non-void
  # at a breathtaking rate. But did I mention, this world is small? So nothing
  # good can come of that ... and that results in a rather spectacular
  # conflagration as non-void turns back into void, and only glowing embers
  # remain. My, my. Nb. it all looks rather random ... but do you notice the
  # two-fold symmetry axis?

  X <- 150; Y <- 150; myWorld <- matrix(integer(X * Y), ncol = X)
  myWorld <- placeLif("data/maxFill.lif", myWorld, 60, 60)
  imPlot(myWorld, drawGrid = FALSE)

  runGOL(myWorld, nTick = 1000, drawGrid = FALSE)


# ==   2.4  Pentomino  =========================================================
  # What is this? A fly speck? a pimple? Five cells of barely nothing? Ah,
  # but perhaps a sesame seed? Or a riot?

  X <- 150; Y <- 150; myWorld <- matrix(integer(X * Y), ncol = X)
  myWorld <- placeCell("data/pentomino.cell", myWorld, 75, 75)
  imPlot(myWorld, drawGrid = FALSE)

  runGOL(myWorld, nTick = 1000, drawGrid = FALSE)


# =    3  Next ?  ==============================================================

  # What's next?

  # 3.1 Look for patterns?
    # There is a lot of literature on that ... common cell formats are .lif and
    # .cell and I have included code that reads them both. Perhaps have a look
    # at the "catagolue" https://catagolue.appspot.com/home with its dictionary
    # and embedded simulations...

  # 3.2 Change the rules?
    # How would probabilistic outcomes change the behaviour? Try it out, the
    # sample code should make that easy to implement. Just choose
    # non-integral rule sets.

  # 3.3 Speed up the code?
    # Can the code be significantly faster if we don't consider every cell, but
    # only every live cell? But how would we do gthat and get the same results?


  # How far can one take this? Check out this video!
  # https://www.youtube.com/watch?v=4lO0iZDzzXk
  #
  # Is this a universal computer?
  # (You bet!)

}

# [END]
