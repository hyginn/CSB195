# imPlot.R
#
# Plot a matrix of 0 and 1 as a rasterImage()
# Plot a 2D text image as a rater image, and overlay text
#
# Boris Steipe (boris.steipe@utoronto.ca)
# 2022-10-23
#
# ToDo: curves
#
# ==============================================================================

imPlot <- function(m, colMap, main, cap, drawGrid) {
  #' @title imPlot()
  #' @description Plot the contents of a matrix, color coded according to
  #' colMap, such that the matrix cell [1.1] is at the top left of the
  #' image and the cells are drawn with an aspect ratio of 1.
  #'
  #' @param  m   A numeric or character matrix, or a 2D text that can be
  #' converted to a character matrix by txt2m().
  #' @param  colMap   A named vector of valid colour specifications for each
  #' character in m. If missing, a default map is produced.
  #' @param main Optional. A title above the plot.
  #' @param cap  Optional. A caption below the plot.
  #' @param drawGrid  logical. Whether to draw gridlines on the output.
  #' @return      NULL, invisible.
  #' @details Input matrix m is converted to a character matrix and the values
  #' are replaced by the colors found in colMap. No checks are made that all
  #' values actually appear in colmap
  #' @examples
  #'   # pseudo-order in the non-repeating Fibonacci word of length 8*13
  #'   s <- c("0100101001001010010100100101001001010010100100101001",
  #'          "0100100101001001010010100100101001001010010100100101")
  #'   s <- unlist(strsplit(s, ""))
  #'   s <- matrix(s, nrow=8, byrow=TRUE)
  #'   print(s)
  #'   imPlot(s, colMap = c("0"="#88B8CE", "1"="#1978A5"), drawGrid = TRUE)

  if (is.numeric(m)) { m[] <- as.character(m) } # square brackets preserve dim()
  if (is.null(dim(m)) && is.character(m)) { m <- txt2m(m) } # convert 2D text

  nx <- dim(m)[2]
  ny <- dim(m)[1]

  if (missing(colMap) || is.null(colMap)) {
    val <- as.character(sort(unique(as.vector(m))))
    colMap <- colorRampPalette(c(
      "#f8f8f8",
      "#1978a5",
      NULL))(length(val))
    names(colMap) <- val
  }

  im <- matrix(colMap[m], ncol = nx)

  if (missing(main) || is.null(main) || length(main) == 0 || main == "") {
    mTop <- 0.5
    main <- ""
  } else {
    mTop <- 5
  }

  if (missing(cap) || is.null(cap) || length(cap) == 0 || cap == "") {
    mBot <- 0.5
    cap <- ""
  } else {
    mBot <- 4
  }

  if (missing(drawGrid) || is.null(drawGrid)) {
    drawGrid <- FALSE
  }


  opar <- par("mar" = c(mBot, 0.5, mTop, 0.5))

  plot(c(0, nx),c(0,ny),           # empty plot to define the frame
       ylim = c(-1, nrow(im) + 1),
       main = main,
       type = "n",
       axes = FALSE,
       xaxs = "i", yaxs = "i",
       xlab = "", ylab = "",
       asp = 1)

  if (cap != "") {
    mtext(cap, side = 1, cex = 0.7)
  }

  rasterImage(im,
              xleft = 0,
              ybottom = 0,
              xright = nx,
              ytop = ny,
              main = main,
              interpolate = FALSE)

  if (drawGrid) {
    segments(0:nx, rep(0, nx+1),               # draw vertical gridlines
             0:nx, rep(ny, nx+1),
             col = "#ddddff")
    segments(rep(0, ny+1), 0:ny,               # draw horizontal gridlines
             rep(nx, ny+1), 0:ny,
             col = "#ddddff")
  }

  par(opar)                                    # reset the graphics state

  return(invisible(NULL))                      # return nothing
}


# == imText() ===========
#

imText <- function(txt, charMap, colMap, cexMap, cex = 1.0) {
  #' @title imText()
  #' @description make a text-plot into an exisiting plot frame of a 2D
  #'              text.
  #' @param txt
  #' @param charMap Convert the input characters to output characters
  #'                according to the map. If missing, leave all characters
  #'                as they are. Characters in input that are not included
  #'                here are left as they are.
  #' @param colMap  Plot the characters in the color specified. Default #000000.
  #' @param cex     Global chharacter expansion value (cex). Default 1.0.
  #' @param cexMap  Plot the characters in the size specified. Default 1.0.
  #'                Values override cex.
  #'
  #' @return NULL invvisibly. Text is added to an existing plot as a
  #'         side-effect.

  m <- txt2m(txt)
  alf <- unique(as.vector(m))

  # prepare maps with default, overwrite with passed parametrs
  x <- alf; names(x) <- alf
  if (missing(charMap)) {
    charMap <- x
  } else {
    x[names(charMap)] <- charMap
    charMap <- x
  }

  x <- rep("#000000", length(alf)); names(x) <- alf
  if (missing(colMap)) {
    colMap <- x
  } else {
    x[names(colMap)] <- colMap
    colMap <- x
  }

  x <- rep(cex, length(alf)); names(x) <- alf
  if (missing(cexMap)) {
    cexMap <- x
  } else {
    x[names(cexMap)] <- cexMap
    cexMap <- x
  }


  for (i in seq_len(nrow(m))) {
    for (j in seq_len(ncol(m))) {
      ch <- m[i, j]
      text(j-0.5, (nrow(m) - i) + 0.5,
           adj = 0.5,
           labels = charMap[ch],
           col = colMap[ch],
           cex = cexMap[ch])
    }
  }

  return(invisible(NULL))

}


# == txt2m() ============

txt2m <- function(txt) {
#' @title txt2m
#' @description A utility function to convert a 2D text to a character matrix.
#' @param txt A 2D text.
#' @details   The txt must have exactly as many line breaks as the matrix will
#'            have rows. All whitespace will be removed. Use "." instead and
#'            imgPlot() with a light color.
#' @return    txt converted to character matrix
#' @examples   txt2m("
#' ..........
#' .01001010.
#' ..........")

  dy <- lengths(gregexpr("\\n", txt))
  txt <- gsub("\\s", "", txt)
  txt <- txt[txt != ""]

  m <- matrix(unlist(strsplit(txt, "")), nrow = dy, byrow = TRUE)

  return(m)
}


# == safeSleep()
safeSleep <- function(tSleep) {
#' @title safeSleep
#' @param tSleep Sleep for tSleep seconds
#' @return NULL invisibly
#' @details Alternative to Sys.sleep() that handles interrupts more safely.
#' cf. https://stackoverflow.com/questions/1174799

  then <-Sys.time()
  while((as.numeric(Sys.time()) - as.numeric(then)) < tSleep){
    ;   # twiddle your thumbs
  }
  return(invisible(NULL))
}




if (FALSE) {
  # pseudo-order in the non-repeating Fibonacci word of length 8*13
  s <- c("0100101001001010010100100101001001010010100100101001",
         "0100100101001001010010100100101001001010010100100101")
  s <- unlist(strsplit(s, ""))
  s <- matrix(s, nrow=8, byrow=TRUE)
  print(s)
  imPlot(s, colMap = c("0"="#88B8CE", "1"="#1978A5"), drawGrid = TRUE)


  myTxt <-"
           ..........
           .01001010.
           .........."
  (myM <- txt2m(myTxt))
  imPlot(myM)

  imPlot("
          010010100100101001010")

}

# [END]

