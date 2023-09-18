# tocID <- ".util.R"
#
# CSB195 Class project: utility scripts
#
# 2022-09  - 2022-10
# boris.steipe@utoronto.ca
#
# This file is source()d upon startup by .Rprofile
#
# ==============================================================================
#


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                           Line
#TOC> ---------------------------------------------------------------
#TOC>   1        Remote control of ChimeraX                        29
#TOC>   2        A progress bar for long-running code             104
#TOC>   3        Find Keywords in aaindex                         127
#TOC>   4        A colour palette for amino acids                 165
#TOC>   5        Extracting R code from Google docs               198
#TOC>   6        Reading Google sheets                            268
#TOC>   7        Plotting amino acids as 2D scatterplot           315
#TOC>
#TOC> ==========================================================================


# =    1  Remote control of ChimeraX  ==========================================
CXPORT <- 61803
cat(sprintf("  Defining ChimeraX port (CXPORT) as %d.\n", CXPORT))

cat("  Defining CX() ...\n")
CX <- function(cmd, port = CXPORT, quietly = FALSE) {
  # send a command to ChimeraX listening on port CXPORT via its REST
  # interface.
  # Parameters:
  #   cmd      char     a ChimeraX commandline command
  #   port     int      the portnumber on which ChimeraX is listening
  #   quietly  logical  if FALSE, cat() the contents of the response
  #
  # Value:  the reply by ChimeraX, invisibly.

  # (A) construct the base address, port, and command:
  CXREST <- sprintf("http://127.0.0.1:%s/run?", CXPORT)

  # (B) sanitize the user-entered variable cmd to be properly encoded in
  # a http message. (No blanks, some specially handled characters, ...)
  # Here we use gsub(), which seeks for a pattern defined in its first
  # argument, substitutes the contents of the second argument, in the
  # string that is identified with the third argument.
  #
  # Patterns are defined as "regular expressions".
  #
  cmd <- gsub("(^\\s+)|(\\s+$)", "", cmd)  # trim whitespace
  # percent encode reserved characters
  cmd <- gsub("%",   "%25", cmd)          #   %
  cmd <- gsub("#",   "%23", cmd)          #   #
  cmd <- gsub("/",   "%2F", cmd)          #   /
  cmd <- gsub(":",   "%3A", cmd)          #   :
  cmd <- gsub("@",   "%40", cmd)          #   @
  cmd <- gsub(",",   "%2C", cmd)          #   ,
  cmd <- gsub("\\*", "%2A", cmd)          #   *
  cmd <- gsub("\\?", "%3F", cmd)          #   ?
  cmd <- gsub("!",   "%21", cmd)          #   !
  cmd <- gsub("=",   "%3D", cmd)          #   =
  cmd <- gsub("\\(", "%28", cmd)          #   (
  cmd <- gsub("\\)", "%29", cmd)          #   )
  cmd <- gsub("\\[", "%5B", cmd)          #   [
  cmd <- gsub("\\]", "%5D", cmd)          #   ]
  cmd <- gsub("&",   "%26", cmd)          #   &
  cmd <- gsub("\\+", "%2B", cmd)          #   +

  cmd <- gsub("\\s+", "+", cmd)            # whitespace to "+"
  cmd <- URLencode(cmd)                    # encode other special characters

  # Combine the base-address and the current contents of cmd ...
  cmd <- paste0(CXREST, "command=", cmd, collapse = "")

  # send the command to ChimeraX, and capture the response ...
  r <- httr::GET(cmd)
  # ... the response is a list-object and we can analyze it.

  if (! r$status_code == 200) {  # If we did NOT receive a "success" status code
    stop("ChimeraX returned status code %d", r$status_code)
  }

  if (length(r$content) == 0) {
    reply <- ""
  } else {
    reply <- rawToChar(r$content)  # reformat
  }

  if (quietly == FALSE) {
    cat(reply)                     # print the response
  }

  return(invisible(reply))         # return the reply  but do not also
                                   # print it.

}


# =    2  A progress bar for long-running code  ================================

cat("  Defining pBar() ...\n")
pBar <- function(i, l, nCh = 50) {
  # Draw a progress bar in the console
  # i: the current iteration
  # l: the total number of iterations
  # nCh: width of the progress bar
  ticks <- round(seq(1, l-1, length.out = nCh))
  if (i < l) {
    if (any(i == ticks)) {
      p <- which(i == ticks)[1]  # use only first, in case there are ties
      p1 <- paste(rep("#", p), collapse = "")
      p2 <- paste(rep("-", nCh - p), collapse = "")
      cat(sprintf("\r|%s%s|", p1, p2))
      flush.console()
    }
  }
  else { # done
    cat("\n")
  }
}

# =    3  Find Keywords in aaindex  ============================================

cat("  Defining grepAAindex() ...\n")

grepAAindex <- function(key, el = "D") {
  # Search for strings in an aaindex element
  # key:   a regular expression pattern
  # el:    which element to search. Default "D" (Definition), also useful might
  #        be "T" (Titles).
  # value: a character vector, length 0 if nothing found.

  library(seqinr)
  data("aaindex")

  dat <- character(length(aaindex))
  for (i in 1:length(aaindex)) {
    dat[i] <- aaindex[[i]][el]
  }
  names(dat) <- 1:length(aaindex)

  sel <- grep(key, dat)
  return(dat[sel])
}

# grepAAindex("[Hh]ydroph")
# idx <- 544
# aaindex[[idx]]$D
# aaindex[[idx]]$T
# cat(sprintf("\n%s\t%s", names(aaindex[[idx]]$I), aaindex[[idx]]$I))
#
# grepAAindex("[Vv]olum")
# idx <- 150
# aaindex[[idx]]$D
# aaindex[[idx]]$T
# cat(sprintf("\n%s\t%s", names(aaindex[[idx]]$I), aaindex[[idx]]$I))



# =    4  A colour palette for amino acids  ====================================

cat("  Defining AACOLS ...\n")

# A colour palette for amino acid properties
AACOLS <- character()
AACOLS["R"] <- "#5770ff" # Positive
AACOLS["K"] <- "#4785EE" #
AACOLS["H"] <- "#37a1de" #
AACOLS["E"] <- "#ff6f59" # Negative
AACOLS["D"] <- "#ff7391" #
AACOLS["N"] <- "#C9D4FF" # Hydrophilic
AACOLS["Q"] <- "#CADFFC" #
AACOLS["S"] <- "#CBEAF9" #
AACOLS["T"] <- "#CDF5F7" #
AACOLS["Y"] <- "#FBFFC9" # Hydrophobic
AACOLS["W"] <- "#EDFDC8" #
AACOLS["F"] <- "#DFFCC8" #
AACOLS["I"] <- "#D2FBC8" #
AACOLS["L"] <- "#C4FAC7" #
AACOLS["M"] <- "#B7F9C7" #
AACOLS["V"] <- "#A9F8C7" #
AACOLS["A"] <- "#9CF7C7" #
AACOLS["G"] <- "#d2d2d2" # Glycine
AACOLS["C"] <- "#fff963" # Cysteine
AACOLS["P"] <- "#edc06d" # Proline
AACOLS <- gsub("$", "AA", AACOLS)  # Make the colors 33% transparent
# AACOLS <- gsub("$", "80", AACOLS)  # Make the colors 50% transparent
# AACOLS <- gsub("$", "55", AACOLS)  # Make the colors 67% transparent
# AACOLS <- gsub("AA$", "", AACOLS) # Remove transparency
# barplot(rep(1, 20), col = AACOLS, names.arg = names(AACOLS), cex.names=0.5)


# =    5  Extracting R code from Google docs  ==================================

cat("  Defining fetchGoogleDocRCode ...\n")

fetchGoogleDocRCode <- function (URL,
                                 delimB = "^# begin code",
                                 delimE = "^# end code",
                                 myExt = ".R") {

  # Retrieve text from a Google doc, subset to a delimited range, write to
  # a tempfile() with extension ".R", and open it in the RStudio editor.
  # Parameters:
  #    URL     chr   URL of a Google doc that is open to share or contained
  #                  in a shared folder
  #    delimB  chr   regex pattern for the begin-delimiter
  #    delimE  chr   regex pattern for the end-delimiter
  #    myExt  chr   extension of tempfile. Default ".R"
  # Value:           None. Executed for its side-effect of writing
  #                  text to tempfile() and opening it in the editor
  #

  # Parse out the ID
  ID <- regmatches(URL, regexec("/d/([^/]+)/", URL))[[1]][2]

  # make a retrieval URL
  URL <- sprintf("https://docs.google.com/document/d/%s%s",
                 ID,
                 "/export?format=txt")

  # GET() the data.
  response <- httr::GET(URL)
  if (! httr::status_code(response) == 200) {
    stop(sprintf("Server status code was \"%s\".",
                 as.character(httr::status_code(response))))
  }

  s <- as.character(response)
  s <- strsplit(s, "\r\n")[[1]]   # split into lines, delimited with \r\n
  iBegin <- grep(delimB, s)       # find the two delimiter indices
  iEnd   <- grep(delimE, s)

  # Sanity checks
  if (length(iBegin) == 0) {
    stop("Begin-delimiter was not found in document.")
  } else if (length(iEnd) == 0) {
    stop("End-delimiter was not found in document.")
  } else if (length(iBegin) > 1) {
    stop("More than one Begin-delimiter in document.")
  } else if (length(iEnd) > 1) {
    stop("More than one End-delimiter in document.")
  } else if ((iEnd - iBegin) < 2) {
    stop("Nothing delimited or delimiter tags not correctly ordered.")
  }

  s <- s[(iBegin+1):(iEnd-1)]          # extract delimited text

  myFile <- tempfile(fileext = ".R")   # get name for temporary file
  write(s, myFile)                     # write s into temporary file
  file.edit(myFile)                    # open in editor

  return(invisible(NULL))              # return nothing
}

if (FALSE) {
  fetchGoogleDocRCode("https://docs.google.com/document/d/15qUO3WwKZSqK84gNj8XZIrCe6Ih791oFfGTJ82nuM_w/edit?usp=sharing")

}



# =    6  Reading Google sheets  ===============================================

cat("  Defining read.gsheet() ...\n")

read.gsheet <- function(URL, sheet, ...) {
  # Read a sheet from a Google shets URL.
  # URL: the URL of the Google sheets spreadsheet. Note: the document must
  #      have been opened to read for everyone with the URL.
  # sheet: the name of the sheet
  # ... : other arguments that are passed to read.csv()
  # value: a dataframe

  # We assume the URL was received from the sheet's sharing button, thus
  # we parse out only the ID

  ID <- regmatches(URL, regexec("/d/([^/]+)/", URL))[[1]][2]

  # make a retrieval URL
  URL <- sprintf("https://docs.google.com/spreadsheets/d/%s%s%s",
               ID,
               "/gviz/tq?tqx=out:csv&sheet=",
               gsub(" ", "+", sheet))

  # the GET() function from httr will get the data.
  response <- httr::GET(URL)
  if (!httr::status_code(response) == 200) {
    stop(sprintf("Server status code was \"%s\".",
                 as.character(httr::status_code(response))))
  }

  x <-as.character(response)
  x <- strsplit(x, "\n")[[1]]
  x[1] <-      gsub("\\\"", "", x[1])
  cNames <-  strsplit(x[1], ",")[[1]]  # Restore the original names from
                                       # row 1 so we don't have to use the
                                       # names that R assigns by default.
  tbl <- utils::read.csv(text = x)
  colnames(tbl) <- cNames

  return(tbl)

}

 # X <- "https://docs.google.com/spreadsheets/d/1tRCPhaua5cjcH_0DuZOiv8BVbdr_V6miC2JeKiOYj-o/edit?usp=sharing"
 # Y <- "AA styles"
 # Z <- read.gsheet(X, Y)


# =    7  Plotting amino acids as 2D scatterplot  ==============================

cat("  Reading Google sheet AADAT ...\n")

AADAT <- read.gsheet("https://docs.google.com/spreadsheets/d/1tRCPhaua5cjcH_0DuZOiv8BVbdr_V6miC2JeKiOYj-o/edit?usp=sharing"
                     , "AA styles")

cat("  Defining plotAA() ...\n")

plotAA <- function(x, y, aaDat = AADAT, ...) {
  # Plot amino acids in 2D
  # x, y: numeric vectors named vectors as Aaa three letter code. They can be
  #       passed in any order, but the order must be the same for x and y
  # ... : other arguments that are passed to plot()
  # aaDat: a dataframe with plotting properties.
  #        (Default: AADAT - defined above)
  # value: none. A plot is produced.

  o <- names(x)
  if (! all(o %in% aaDat$Aaa)) {
    stop("Names of first vector do not match three-letter codes.")
  }
  if (! all(o %in% names(y))) {
    stop("Names of first and second vector do not match.")
  }
  if (! all(o == names(y))) {
    stop("Order of first and second vector is not the same.")
  }
  ord <- match(o, aaDat$Aaa) # index vector that maps the order of input
                             # to the order in aaDat

  # rescale aaDat$vol between 0.7 and 5
  CEXMIN <- 0.7
  CEXMAX <- 4.3
  ci <- min(aaDat$vol)
  ca <- max(aaDat$vol)
  cexVol <- ((aaDat$vol - ci) / (ca - ci) * CEXMAX) + CEXMIN
  names(cexVol) <- aaDat$Aaa

  plot(x, y,
       pch = aaDat$pch[ord],
       cex = cexVol[ord],
       col = gsub("..$", "FF", aaDat$col[ord]),
       bg = aaDat$col[ord])

  CEXMIN <- 0.7
  CEXMAX <- 1.3
  cexVol <- ((aaDat$vol - ci) / (ca - ci) * CEXMAX) + CEXMIN
  names(cexVol) <- aaDat$Aaa
  text(x, y, labels = aaDat$A[ord], cex = cexVol[ord], col="#000000AA")
}

# library(seqinr)
# data("aaindex")
# x <- aaindex[[150]]$I
# y <- aaindex[[544]]$I
# plotAA(x, y)


# [END]
