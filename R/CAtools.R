# tocID <- "CAtools.R"
#
# Utilities for Cellular Automata
#
# 2022 - 10  -  2022 - 11
#
# Boris Steipe (boris.steipe@utoronto.ca)
#
# Version:  2.0
#
# Versions:
#           2.0   Update logic to be use applyRule in standard, probabilistic
#                   or time-reversible mode.
#           1.0   First version
#
#  Notes:
#           ...
#
#  To Do:
#           ...
#
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                Line
#TOC> ------------------------------------
#TOC>   1        showRule()             39
#TOC>   2        plotRule()             83
#TOC>   3        applyRule()           156
#TOC>   4        CA()                  213
#TOC>
#TOC> ==========================================================================


# =    1  showRule()  ==========================================================

showRule <- function(R, quietly = FALSE) {
  #' @title showRule()
  #' @description Print out the rule-set for Wolfram code R
  #' @param R  int.  The Wolfram code for a 1D cellular automaton rule-set;
  #'                 an integer in [0, 255].
  #' @param quietly  If FALSE, print. If TRUE, do not print.
  #'
  #' @return   int.  Invisibly, the rule set as a named integer vector.
  #'                 If "qietly" is false, the rule set is printed as a
  #'                 side-effect.
  #'
  #' @details showRule() prints out the construction of a cellular automaton
  #'          rule-set computed from its Wolfram code. This code is derived from
  #'          the state-configurations of a cell's neighborhood, ordered in
  #'          descending order - in our 1D case: [111][110][101]...[000].
  #'          Each of the eight configurations can be taken as a digit in
  #'          an eight-digit binary number, and thus each of the possible
  #'          255 integers that can be constructed from 8-bit specifies
  #'          a unique rule set.
  #'
  #' @examples   showRule(18)
  #'             myR <- showRule(110)
  #'             myR["010"]
  if (R < 0 || R > 255) {
    stop("R is not in [0,255].")
  }

  v <- as.integer(intToBits(R)[8:1] == "01")  # convert integer to binary array

  names(v) <- c("111", "110", "101", "100",   # name the elements
                "011", "010", "001", "000")

  if (!quietly) {
    cat(sprintf("\n  Rule: %d\n  %s\n   %s  ",
                R,
                paste(names(v), collapse = " "),
                paste0(v, collapse = "   ")))
  }

  return(invisible(v))
}

# =    2  plotRule()  ==========================================================

plotRule <- function(R, main, cap) {
  #' @title plotRule()
  #' @description Plot the rule-set for Wolfram code R as a graphic.
  #' @param R  int  The Wolfram code for a 1D cellular automaton rule-set;
  #'                 an integer in [0, 255].
  #' @param nain char Optional: a figure title passed to imPlot() .
  #' @param cap  char Optional: a figure caption passed to imPlot() .
  #' @return   int  The rule set is plotted as a side-effect.
  #'
  #' @details plotRule() prints out the construction of a cellular automaton
  #'          rule-set computed from its Wolfram code. This code is derived from
  #'          the state-configurations of a cell's neighborhood, ordered in
  #'          descending order - in our 1D case: [111][110][101]...[000].
  #'          Each of the eight configurations can be seen as a digit in
  #'          an eight-digit binary number, and thus each of the 255 numbers
  #'          specifies a unique rule set.
  #'
  #' @examples   plotRule(18)
  #'   plotRule(110,
  #'            main = "Rule 110",
  #'            cap="The Rule 110 CA is a \"universal computer\".")

  if (R < 0 || R > 255) {
    stop("R is not in [0,255].")
  }

  if (missing(main)) { main <- "" }
  if (missing(cap )) { cap  <- "" }

  s <- "
.................................
.III.II0.I0I.I00.0II.0I0.00I.000.
..v...v...v...v...v...v...v...v..
..a...b...c...d...e...f...g...h..
................................."

  v <- showRule(R, quietly = TRUE)  # convert integer to binary array

  s <- gsub("a", v["111"], s)
  s <- gsub("b", v["110"], s)
  s <- gsub("c", v["101"], s)
  s <- gsub("d", v["100"], s)
  s <- gsub("e", v["011"], s)
  s <- gsub("f", v["010"], s)
  s <- gsub("g", v["001"], s)
  s <- gsub("h", v["000"], s)

  imPlot(s,
         colMap = c("." = "#F4FAFE",
                    "v" = "#F4FAFE",
                    "0" = "#dbeaf0",
                    "1" = "#5994a6",
                    "O" = "#b5d8e3",
                    "I" = "#19738f"),
         drawGrid = TRUE,
         main = main,
         cap = cap)

  imText(s,
         charMap = c("."="", "O"="0", "I"="1", "v"= "↓"),
         colMap = c("1"="#FFFFFF", "I"="#FFFFFF"),
         cex = 0.7)

  return(invisible(NULL))
}

if (FALSE) {
  plotRule(94)
}


# =    3  applyRule()  =========================================================

applyRule <- function(R, s, mode = "std", p, prev) {
  #' @title applyRule()
  #' @description The new state of a cellular automaton is computed by applying
  #'   the Wolfram code rule-set R to the neighbourhood state s.
  #'
  #' @param  R    int The Wolfram code for a 1D cellular automaton rule set; an
  #'   integer in [0, 255].
  #' @param  s    int A binary vector of length three.
  #' @param  mode chr The mode of applying the rule. One of ["std", "tRev", "pI"
  #'   and "pO"].
  #' @param p     num Probability of making an error while applying the code
  #'   rule.
  #' @param prev  int The state of the cell at t-1
  #' @details   Computing the output state requires only the rule number, since
  #'   the binary representation of the Wolfram code number IS the rule set. It
  #'   does not need to be separately stored. This is done in the "std" mode
  #'   (default). There are three modifications to this "standard" way of
  #'   updating the cell's state: in mode "pI", each digit of the input state s
  #'   is inverted (1->0 and 0->1) with probability p. In mode "pO" the output
  #'   state is inverted with probability p. In mode tRev the output state is
  #'   combined xor() with the previous state (prev), to produce a
  #'   time-reversible state change
  #' @return      0 or 1
  #' @examples
  #'   applyRule(110, c(1,0,1))                      # standard
  #'   applyRule(239, c(1,0,0), mode="pI", p=0.5)    # probabilistic
  #'   applyRule(  0, c(1,0,0), mode="pO", p=0.5)    # probabilistic
  #'   applyRule( 37, c(1,0,0), mode="tRev", prev=1) # time-reversible

  if (mode == "pI") {    # probabilistic change of the input state
    x <- runif(3) < p    # three random values compared to the threshold
    s[x] <- 1 - s[x]     # invert the randomly selected cells
  }


  v <- as.integer(intToBits(R)[8:1] == "01")  # compute the output integers
                                              # from the binary representation
                                              # of R
  idx <-  8 - (s[1]*4) - (s[2]*2) - (s[3]*1)  # compute the index in v that
                                              # represents s
  S <- v[idx]                                 # fetch the state S from v


  if (mode == "pO") {                # probabilistic change of the output state
    if (p > runif(1)) {              # random value compared to the threshold
      S <- 1-S                       # invert output state
    }
  } else if (mode == "tRev") {       # compute time-reversible output
    S <- as.numeric(xor(S, prev))    # xor() with state at t-1
  }

  return(S)                          # return output state
}


# =    4  CA()  ================================================================

CA <- function(R,
               ruleMode = "std",
               pR,
               nx = 100, ny = 100,
               vInit, seed = NULL,
               wrap = TRUE, N = NULL) {
  #' @title CA()
  #' @description Evolve a cellular automaton with Wolfram code R in a nx * ny
  #'   matrix.
  #'
  #' @param  R   int. The Wolfram code for a 1D cellular automaton rule set; an
  #'   integer in [0, 255].
  #' @param ruleMode  Run this as a standard, probabilistic, or time-reversible
  #' application of rule R. (See documentaation for applyRule.)
  #' @param pR  probability for one of the probabilistic rule applications.
  #' @param  nx  the width of the space in which the automaton evolves
  #' @param  ny  the number of steps to compute
  #' @param  vInit   the initialization vector. See details.
  #' @param  seed    Initializaton of the RNG. Default NULL.
  #' @param  wrap    if TRUE, wrap the vertical edges of the matrix (periodic
  #'   boundary).
  #' @param  N    If this is not NULL, the function plots every updated row
  #'   until the matrix is filled, then scrolls the matrix for a total of N
  #'   steps (or until the process is aborted.).
  #' @return      the evolved matrix as 0s and 1s.
  #' @details If vInit is missing or NULL, a single 1 is placed at the centre of
  #'   the first row. If vInit is in [0, 1[ the vector is filled with random 1s
  #'   of p == vInit. If vInit is >= 1, round(vInit) 1s are equally spaced into
  #'   vInit. If vInit is the keyword "fib", the first row is initialized with a
  #'   Fibonacci word. If length(vInit) > 1 the vector is recycled into the
  #'   length of vInit. In particular, this allows to pass in a custom
  #'   initialization vector of length nx - this vector must be a numeric vector
  #'   of 0s and 1s. If vInit is a character-string of 0s and 1s, it is
  #'   converted into such an initialization vector.
  #' @examples
  #'

  i0 <- 1 # this is the index of the starting row. It is always 1, except for
          # time reversible CAs. Then it is 2.

  m <- matrix(integer(nx * ny), ncol = nx)

  # All the initialization stuff.
  if (missing(vInit) || is.null(vInit)) { # place a single 1 in the middle

    m[1, round((nx) / 2)] <- 1

  } else if (ruleMode == "tRev"){
    stopifnot(ncol(vInit) == nx && nrow(vInit) == 2)
    m[1:2, ] <- vInit
    i0 <- 2

  } else if (length(vInit) > 1){  # recycle this vector
    m[1, ] <- rep(vInit, ceiling(nx/length(vInit)))[1:nx]

  } else if (is.numeric(vInit) && vInit >= 0 && vInit < 1) { # probability
    oldSeed <- .Random.seed
    set.seed(seed)
    m[1, ] <- sample(0:1, nx, replace= TRUE, prob = c(1-vInit, vInit))
    .Random.seed <- oldSeed  # Be considerate, reset the RNG to how you found it

  } else if (is.numeric(vInit) && vInit >= 1) { # place isolated 1s
    idx <- cumsum(c(0.5, rep(1, vInit - 1), 0.5))
    idx <- round((idx / max(idx)) * nx)[-length(idx)]
    m[1, idx] <- 1

  } else if (vInit == "fib") { # Fibonacci word
    m[1, ] <- as.integer(unlist(strsplit(fibWord(nx), "")))

  } else if (grepl("^[01]+$",vInit)) { # string of 0s and 1s
    x <- as.integer(unlist(strsplit(vInit, "")))
    m[1, ] <- rep(x, ceiling(nx/length(vInit)))[1:nx]

  } else {
    stop("PANIC: I don't know what to do with this vInit.")
  }

  if (is.null(N)) {  # No plotting of progress
    N <- 0
  } else {
    Nmax <- N
    cat("\n")
  }

  if (wrap) {
    m <- cbind(m[ ,nx], m, m[ ,1]) # wrap first and last column
    nx <- nx + 2                  # increase nx accordingly
  }


  for (i in i0:(ny-1)) {
    for (j in 2:(nx-1)) {
      m[(i+1), j] <- applyRule(R,
                               m[i,(j-1):(j+1)],
                               ruleMode,
                               pR,
                               prev = m[i-1, j])
    }
    if (wrap) {
      m[i+1, 1] <- m[i+1, j]     # wrap last cell to beginning
      m[i+1, nx] <- m[i+1, 2]    # wrap first cell to end
    }
    if (N > 0) {
      imPlot(m)                  # plot it out row-by-row
      dev.flush()
      N <- N-1
      cat(sprintf("\rStep %d", Nmax - N))
    }
  }

  for (k in seq_len(N)) { # scroll if N is > 0
    m[1:(ny-1), ] <- m[2:ny, ] # shift the matrix one row up
    for (j in 2:(nx-1)) {
      m[(i+1), j] <- applyRule(R,
                               m[i,(j-1):(j+1)],
                               ruleMode,
                               pR,
                               prev = m[i-1, j])
    }
    m[i+1, 1] <- m[i+1, nx-1]  # wrap last cell to beginning
    m[i+1, nx] <- m[i+1, 2]    # wrap first cell to end
    imPlot(m)
    dev.flush()
    cat(sprintf("\rStep %d", Nmax - (N - k)))

  }

  if (wrap) {
    m <- m[ , -c(1, nx)] # remove wrapped columns
  }

  return(m)
}

if (FALSE) {

  # ==  VALIDATE ===============================================================

  # Standard mode, Rule 18, all defaults
  m <- CA(18); imPlot(m)

  # probabilistic CA
  m <- CA(18, ruleMode = "pI", pR = 0.0002);   imPlot(m)
  m <- CA(18, ruleMode = "pO", pR = 0.0002);   imPlot(m)

  # time-reversible CA
  X <- 200
  Y <- 323

  mI <- matrix(integer(2*X), nrow = 2)   # two empty rows
  mI[2, sample(1:X, 13)] <- 1             # thirteen random 1s

  m <- CA(37, nx=X, ny=Y, ruleMode="tRev", vInit=mI);   imPlot(m)

  mI <- matrix(integer(2*X), nrow = 2)   # two empty rows
  mI[2, sample(1:X, 17)] <- 1             # thirteen random 1s
  m <- CA(22, nx=X, ny=Y, ruleMode="tRev", vInit=mI);   imPlot(m)

}


# [END]
