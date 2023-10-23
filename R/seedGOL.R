# tocID <- "seedGOL.R"
#
# Seeding a GOL world with patterns
#
# 2022-10-31
# Boris Steipe (boris.steipe@utoronto.ca)
#
# Version:  1.0
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                Line
#TOC> ------------------------------------
#TOC>   1        placeCell()            21
#TOC>   2        placeLif()             55
#TOC> 
#TOC> ==========================================================================


# =    1  placeCell()  =========================================================
#

placeCell <- function(inFile, W, X, Y) {
# Place a cell formatted pattern from inFile into W at X, Y

  s <- readLines(inFile)
  s <- s[! grepl("^! ", s)]   # remove comments
  s <- gsub("O", "1", s)
  s <- gsub("\\.", "0", s)
  dx <- nchar(s[1])
  dy <- length(s)
  s <- as.integer(unlist(strsplit(s, "")))
  s <- matrix(s, ncol = dx, byrow = TRUE)


  # write the matrix s into W at position [Y, X]
  W[Y:(Y + dy - 1), X:(X + dx - 1)] <- s

  return(W)
}



if (FALSE) {

  myWorld <- matrix(numeric(60 * 40), ncol = 60)
  myWorld <- placeCell("GosperGun.cell", myWorld, 2, 3)

  imPlot(myWorld, drawGrid = TRUE)

}


# =    2  placeLif()  ==========================================================

placeLif <- function(inFile, W, X, Y) {
# Place a lif formatted pattern from inFile into W at X, Y
#
# For file format see: https://conwaylife.com/wiki/Life_1.05
# For a catalogue see http://www.radicaleye.com/lifepage/patterns/contents.html

  s <- readLines(inFile)
  iBlocks <- grep("^#P ",s)
  iBlocks <- c(iBlocks, length(s)+1)  # terminate
  blocks <- vector(mode = "list", length = length(iBlocks) - 1) # init list

  # read each block
  for (i in seq_len(length(iBlocks) - 1) ) {  # all but last index
    blocks[[i]]$xy    <- as.integer(strsplit(s[iBlocks[i]], "\\s")[[1]][2:3])
    blocks[[i]]$block <- s[(iBlocks[i] + 1):(iBlocks[i+1] - 1)]
  }

  # determine bounding box coordinates
  bbox <- integer(4)
  names(bbox) <- c("xi", "xa", "yi", "ya")  # min / max, counting from left-top

  for (i in seq_len(length(blocks)) ) {
    bbox["xi"] <- min(bbox["xi"], blocks[[i]]$xy[1])
    bbox["xa"] <- max(bbox["xa"],
                      max(blocks[[i]]$xy[1] + nchar(blocks[[i]]$block)))
    bbox["yi"] <- min(bbox["yi"], blocks[[i]]$xy[2])
    bbox["ya"] <- max(bbox["ya"], blocks[[i]]$xy[2] + length(blocks[[i]]$block))
  }

  # make a matrix to hold the patterns
  dx <- bbox["xa"] - bbox["xi"]
  dy <- bbox["ya"] - bbox["yi"]
  m <- matrix(integer(dx * dy), ncol = dx)

  # add each block to the matrix
  for (i in seq_len(length(blocks)) ) {
    dLeft <- -bbox["xi"]
    dTop  <- -bbox["yi"]
    for (j in seq_len(length(blocks[[i]]$block))) { # for each row
      s <- blocks[[i]]$block[j]                     # convert to vector of (0,1)
      s <- gsub("\\.", "0", s)
      s <- gsub("\\*", "1", s)
      s <- as.integer(strsplit(s, "")[[1]])
      x <- blocks[[i]]$xy[1] + dLeft + 1            # first column
      y <- blocks[[i]]$xy[2] + dTop + j             # row
      m[y, (x:(x + length(s) - 1))] <- s

    }
  }

  # write the matrix into W at position [Y, X]
  W[Y:(Y - 1 + dy), X:(X - 1 + dx)] <- m

  return(W)
}


if (FALSE) {
  inFile <- "Rabbits.lif"
  myWorld <- matrix(numeric(20 * 30), ncol = 30)
  myWorld <- placeLif("Rabbits.lif", myWorld, 5, 4)

  imPlot(myWorld, drawGrid = TRUE)


  inFile <- "pentomino.lif"
  myWorld <- matrix(numeric(15 * 20), ncol = 20)
  myWorld <- placeLif(inFile, myWorld, 8, 10)

  imPlot(myWorld, drawGrid = TRUE)


  inFile <- "maxFill.lif"

  myWorld <- matrix(numeric(100 * 100), ncol = 100)
  myWorld <- placeLif("maxFill.lif", myWorld, 30,30)

  imPlot(myWorld, drawGrid = TRUE)


  inFile <- "linePuff.lif"

  myWorld <- matrix(numeric(100 * 100), ncol = 100)
  myWorld <- placeLif("linePuff.lif", myWorld, 30,30)

  imPlot(myWorld, drawGrid = TRUE)


}


# [END]
