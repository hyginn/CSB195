# tocID <- "courseScripts/08.1-gameOfLife.R"
#
# Demo code
#
# 2022 - 10  -  2024 - 10
#
# Boris  Steipe (boris.steipe@utoronto.ca)
#
#
# Version:  1.2
#
# Versions:
#           1.2   2024 Updates. Extensively commented, added URL's and examples
#           1.1.1 Moved src and pattern files to R/ resp. data/ directories
#           1.1   Updated, expanded and added examples section
#           1.0   In class demo
#
# Notes:
#
# To Do:
#
# ==============================================================================
#                                                                              #
#   THIS SCRIPT IS "SAFE TO SOURCE" WITHOUT SIDE EFFECTS.                      #
#                                                                              #
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                     Line
#TOC> -----------------------------------------
#TOC>   1        THE GAME OF LIFE            47
#TOC>   2        IMPLEMENT                   82
#TOC>   2.1        wrap()                    91
#TOC>   2.2        runGOL()                 111
#TOC>   3        EXPLORATIONS               168
#TOC>   3.1        The Glider               193
#TOC>   3.2        A Glider-gun             204
#TOC>   3.3        Max Fill                 217
#TOC>   3.4        Pentomino                232
#TOC>   4        NEXT ?                     243
#TOC>
#TOC> ==========================================================================


# =    1  THE GAME OF LIFE  ====================================================

# Cellular automata can be generalized from one dimension to to two dimensions,
# and the best known implementation of this is John Horton Conway's "Game of
# Life" (GOL). cf. https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life This is
# a demo implementation of GOL in R for you to explore it. R and RStudio are
# however not the ideal environments to actually evolve GOL scenarios, due to
# the computational overhead of plotting. There are dedicated sites with
# Javascript implementations that are much faster, but they generally don't
# allow you to poke around in the actual code. See for example:

#   https://playgameoflife.com/
#   https://conwaylife.com/

# For an intriguing vision of what the game is capable of, have a look
# at this youtube video:

#   https://www.youtube.com/watch?v=4lO0iZDzzXk  (you don't have to watch all
#                                                 of it  :-)


# But what are all these dots, and why do they behave the way they do?

# Simple: These dots are Cellular Automata, but of course not just the white
# dots, the black dots too are automata. They are white if they are in one of
# two states and black if they are in the other. In GOL terminology we call the
# white cells "alive" and the black cells "dead". These automata sense their
# environment and evolve step after step (also "generation", or "tick")
# according to very simple rules:

#   0: All cells remain in their state except:
#   1: If a live cell has less than two or more than three neighbours, it dies.
#   2: If a dead cell has three neighbours, it becomes alive.


# =    2  IMPLEMENT  ===========================================================

# Source some utility functions:
source("R/imPlot.R")    # Code that plots the contents of a matrix in a way that
                        # is useful for us here.
source("R/seedGOL.R")   # Code that can load some sample scenarios from the
                        # data folder into a GOL matrix.


# ==   2.1  wrap()  ============================================================

# GOL can evolve in an infinite, open world, or in a world with periodic
# boundary conditions. This means a particle that moves off the top comes back
# in at the bottom, and a particle that moves of to the left comes back in at
# the right. Topologically speaking, this world is a "torus". To implement this,
# we use a function we call wrap().

wrap <- function(w) {
# wrap the 2D matrix "world" w of size x*y into a torus (x-2) * (y-2)

  w[1, ]       <- w[nrow(w)-1, ]    # wrap last actual row to row 1
  w[nrow(w), ] <- w[2,         ]    # wrap row 2 to last row
  w[ ,1]       <- w[ ,ncol(w)-1]    # wrap last actual column to column 1
  w[ ,ncol(w)] <- w[ ,2        ]    # wrap column 2 to last column

  return(w)
}


# ==   2.2  runGOL()  ==========================================================

# This function takes a matrix and evolves it step by step. It plays the basic
# game by default, but you can define your own rules for living and dead cells.

runGOL <- function(w0, rL, rD, nTick = 100, drawGrid = TRUE) {

  # Evolve world w0 according to rule rL and rD for nTick generations

  # Rules: the rL rules mean: if a living cell (1) has 4 living neighbours, it
  #        turns into a dead cell (0) in the next generation. And if it
  #        has 3 living neighbors, it stays alive (1). Etc.
  if (missing(rL)) {
    # Neighbours:    0  1  2  3  4  5  6  7  8
    rL <-          c(0, 0, 1, 1, 0, 0, 0, 0, 0)
  }
  if (missing(rD)) {
    # Neighbours:    0  1  2  3  4  5  6  7  8
    rD <-          c(0, 0, 0, 1, 0, 0, 0, 0, 0)
  }

  w0 <- wrap(w0)   # "wrap" the input world
  w1 <- w0         # initialize the next-generation world

  for (tick in seq_len(nTick)) {    # for each tick
    pBar(tick, nTick)               # progress bar

    for (i in 2:(nrow(w0)-1)) {     # for each actual row
      for (j in 2:(ncol(w0)-1)) {   # for each actual column

        # sum of neighbours in 3 x 3 neighbourhood - self
        nNeigh <- sum(w0[(i-1):(i+1), (j-1):(j+1)]) - w0[i, j]
        # Rules
        if (w0[i, j] == 1) { w1[i, j] <- rL[nNeigh + 1]
        } else             { w1[i, j] <- rD[nNeigh + 1] }
      } # end column
    } # end row
    w1 <- wrap(w1)  # "wrap" the result

    # visualize the next stage
    imPlot(w1[2:(nrow(w1)-1), 2:(ncol(w1)-1)], drawGrid = drawGrid)

    w0 <- w1        # set the new world to be the current world

    safeSleep(0.1)  # smooth
    dev.flush()     # display

  } # end tick

  return(invisible(w1))
}




if (FALSE) {

# =    3  EXPLORATIONS  ========================================================

  # Execute this code from the script ...

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


# ==   3.1  The Glider  ========================================================
  # People have been searching for interesting patterns for a long time. This
  # one is a "glider", it can be used to send information around ...

  X <- 50; Y <- 50; myWorld <- matrix(integer(X * Y), ncol = X)
  myWorld[25:27,27]<-1;myWorld[26,25]<-1;myWorld[27,26]<-1   # glider
  imPlot(myWorld)

  runGOL(myWorld, nTick = 1000)


# ==   3.2  A Glider-gun  ======================================================
  # This "glider gun", spawns gliders. Regularly, like clockwork. Or like the
  # clock cycles of a computer. But this world is small. What goes around comes
  # around ... and around  ... and around - with explosive and surprising
  # results. Will it all stop? Or go on forever?

  X <- 57; Y <- 50; myWorld <- matrix(integer(X * Y), ncol = X)
  myWorld <- placeCell("data/GosperGun.cell", myWorld, 5, 20)
  imPlot(myWorld, drawGrid = FALSE)

  runGOL(myWorld, nTick = 1000, drawGrid = FALSE)


# ==   3.3  Max Fill  ==========================================================
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


# ==   3.4  Pentomino  =========================================================
  # What is this? A fly speck? a pimple? Five cells of barely nothing? Ah,
  # but perhaps a sesame seed? Or a riot?

  X <- 150; Y <- 150; myWorld <- matrix(integer(X * Y), ncol = X)
  myWorld <- placeCell("data/pentomino.cell", myWorld, 75, 75)
  imPlot(myWorld, drawGrid = FALSE)

  runGOL(myWorld, nTick = 1000, drawGrid = FALSE)


# =    4  NEXT ?  ==============================================================

  # What's next? Look for patterns?

  # There is a lot of literature on that ... common cell formats are .lif and
  # .cell and I have included code that reads them both and when you find
  # an interesting pattern online you can download it and load it like so:

  X <- 600; Y <- 150; myWorld <- matrix(integer(X * Y), ncol = X)
  myWorld <- placeLif("data/linePuff.lif", myWorld, 30, 25)
  imPlot(myWorld, drawGrid = FALSE)
  runGOL(myWorld, nTick = 1000, drawGrid = FALSE)


  # Perhaps have a look
  # at the "catagolue" https://catagolue.appspot.com/home with its dictionary
  # and embedded simulations...

  # Or change the rules? Probabilistic? Non-integral? The sample code
  # should make that easy.

  # Is this a universal computer?
  # (You bet!)


} # end if (FALSE) { ...}

# [END]
