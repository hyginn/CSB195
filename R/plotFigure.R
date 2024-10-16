# tocID <- "plotFigure.R"
#
#   Plotting figures for the cellularAutomata.R script; optionally
#   showing the source code for that section. First attempt at a more
#   interactive-textbook like coding style.
#
# Boris Steipe (boris.steipe@utoronto.ca)
#
# Date:       2022-10 - 2024-10
#
# Version:    1.1
#
# Versions:   1.1  Maintenance
# Versions:   1.0  Published to class
#
# To Do:
#             ...
#
# ==============================================================================
#

# NOTE:
#         Do not change the formatting of if() / else if() statements
#         since catSection() makes assumptions about this.

plotFigure <- function(fig, showSource, ...) {
  #' @title plotFigure()
  myArgs <- list(...)
  thisFile <- normalizePath("./R/plotFigure.R")
  if (FALSE) {


# == BEGIN FIGURE CODE =========================================================

# == 1  CA (Cellular Automata)

  } else if (fig == "CA.1" ) { # ===============================================
    s <- "
.................................
.111.110.101.100.011.010.001.000.
................................."
    imPlot(s, cap = "Figure: A cell and its two neigbours can be in eight distinct binary states.")
    imText(s, charMap= c("."=""), colMap = c("1"="#ffffff"), cex = 0.8)


  } else if (fig == "CA.2" ) { # ===============================================
    plotRule(18,
             main = "Rule 18",
             cap =  "Figure: Input and output states for \"Rule 18\"")

  } else if (fig == "CA.3" ) { # ===============================================
    s <- "
................
....01001010....
................
....OIO.........
.....0..........
.....IOO........
......1.........
......OOI.......
.......1........
.......OIO......
........0.......
........IOI.....
.........0......
.........OIO....
..........0.....
................
....-011000-....
................"
    imPlot(s,
           colMap = c(
             "." = "#F4FAFE",
             "-" = "#F4FAFE",
             "v" = "#F4FAFE",
             "O" = "#efefef",
             "I" = "#bbbbbb",
             "0" = "#b5d8e3",
             "1" = "#19738f",
             NULL),
           drawGrid = TRUE,
           main = "One step of evolving \"01001010\" with rule 18" ,
           cap = "Figure: To \"evolve\" a CA, consider its overlapping segments.")

    imText(s,
           charMap = c("."="", "O"="0", "I"="1", "v"= "↓"),
           colMap = c("1"="#FFFFFF", "I"="#FFFFFF"),
           cex = 0.7)


  } else if (fig == "CA.4" ) { # ===============================================
    s <- "
....01001010....
................
...Q01001010D...
................
...DOO....IOQ...
....1......1....
....|......|....
....10110001....
................"
    imPlot(s,
           colMap = c(
             "." = "#F4FAFE",
             "|" = "#F4FAFE",
             "-" = "#F4FAFE",
             "v" = "#F4FAFE",
             "I" = "#bbbbbb",
             "O" = "#efefef",
             "Q" = "#efefef",
             "D" = "#efefef",
             "1" = "#19738f",
             "0" = "#b5d8e3",
             "q" = "#b5d8e3",
             "d" = "#b5d8e3",
             NULL),
           drawGrid = TRUE,
           main = "Periodic boundary conditions for \"0100101\"" ,
           cap = "Figure: First and last cells are wrapped around.")

    imText(s,
           charMap = c("."="",
                       "v"= "↓",
                       "O"="0","q"="0","d"="0","Q"="0","D"="0",
                       "I"="1"),
           colMap = c( "q"="#FF0000", "Q"="#FF0000",
                       "d"="#007700", "D"="#007700",
                       "1"="#FFFFFF", "I"="#FFFFFF"),
           cex = 0.7)



  } else if (fig == "CA.5" ) { # ===============================================

    nx <- 11  # width
    ny <- 9   # height
    m <- matrix(integer(nx * ny), nrow = ny)  # integer() fills the matrix with 0s

    m[1,6] <- 1  # Initialize the first row
    imPlot(m, drawGrid = TRUE,
           main = sprintf("Rule %d", myArgs$rule))

    R <- myArgs$rule
    printRule(R)
    xl <- c(nx, 1:(nx-1))
    xr <- c(2:nx, 1)
    for (i in 1:(ny-1)) {
      for (j in 1:nx) {
        m[(i+1), j] <- applyRule(R, c(m[i,xl[j]], m[i,j], m[i,xr[j]]))
      }
      imPlot(m, drawGrid = TRUE)
      readline("Hit <Return> for next row >")
    }



  } else if (fig == "CA.6" ) { # ===============================================

    m <- CA(myArgs$iRule, nx = myArgs$nx, ny=myArgs$ny)
    imPlot(m,
           drawGrid = myArgs$drawGrid,
           cap = sprintf("Rule %d", myArgs$iRule))
    printRule(myArgs$iRule)



  } else if (fig == "CA.7" ) { # ===============================================

    i <- myArgs$i
    j <- myArgs$j

    w <- fibWord(i * j, mode = "num")
    imPlot(matrix(w, ncol = i, byrow = TRUE),
           drawGrid = TRUE,
           main = "Fibonacci Word",
           cap = sprintf("length %d * %d (%d per row)", i, j, i))



  } else if (fig == "CA.8" ) { # ===============================================

    m <- CA(myArgs$iRule, nx = myArgs$nx, ny=myArgs$ny, vInit = myArgs$vInit,
            seed = myArgs$seed)
    if (is.null(myArgs$vInit)) { myArgs$vInit <- "NULL" }
    if (length(myArgs$vInit) > 1 ) { myArgs$vInit <- "vector" }
    imPlot(m,
           cap = sprintf("CA %d (initialization mode: %s)",
                         myArgs$iRule, as.character(myArgs$vInit)))


  } else if (fig == "CA.9" ) { # ===============================================


    if (is.null(myArgs$vInit)) { myArgs$vInit <- "NULL" }
    if (length(myArgs$vInit) > 1 ) { myArgs$vInit <- "vector" }

    for (i in myArgs$start:myArgs$end) {
      myCol <- sample(hcl.colors(256, "Temps"),1)  # random color
      m <- CA(i, nx = myArgs$nx, ny=myArgs$ny, vInit = myArgs$vInit,
              seed = myArgs$seed)
      printRule(i)
      cat("\n--------------------------------------\n")
      imPlot(m, colMap = c("0" = "#f3f3f3", "1"=myCol))
      safeSleep(myArgs$sleep)
    }


  } else if (fig == "CA.10" ) { # ===============================================

    m <- CA(myArgs$iRule,
            nx = myArgs$nx,
            ny=myArgs$ny,
            vInit = myArgs$vInit,
            seed = myArgs$seed,
            N = myArgs$nSteps)



    # == END OF FIGURE CODE ========================================================
  } else { stop(sprintf("Figure # \"%s\" not defined in %s", fig, thisFile)) }
  if (! missing(showSource)) { catSection(thisFile, fig) }
  return(invisible(NULL))
}


# == catSection() ==============================================================
catSection <- function(fN, sec) {
  txt <- readLines(fN)
  pat <- sprintf("^\\s*\\}\\s*else\\s+if\\s*\\(fig ==") # else-if statements
  idx <- grep(pat, txt)
  # index the end of the last statement
  idx <- c(idx, (grep("^\\s*\\}\\s*else\\s*\\{\\s*stop\\s*\\(", txt) - 1))
  # find the if-satement with the section
  iI <- grep(sprintf("fig\\s*==\\s*\\\"%s\\\"", sec), txt[idx])

  cat(sprintf(" == Script source for Figure \"%s\" ===================\n", sec))
  cat(txt[(idx[iI]+1):(idx[iI + 1]-1)], sep = "\n")
  cat(" ======================================================\n")

  return(invisible(NULL))
}



# [END]
