# tocID <- ".util.R"
#
# CSB195 Class project: utility scripts
#
# 2022-09  - 2023-09
# boris.steipe@utoronto.ca
#
# This file is source()d upon startup by .Rprofile
#
# ==============================================================================
#


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                            Line
#TOC> ----------------------------------------------------------------
#TOC>   01       Install missing packages                           34
#TOC>   02       Load required libraries                            50
#TOC>   03       Generative AI                                      57
#TOC>   03.1       t2c - write text to clipboard                    59
#TOC>   03.2       Initialize generative AI initial prompt          74
#TOC>   04       Remote control of ChimeraX                        107
#TOC>   05       A progress bar for long-running code              183
#TOC>   06       Find Keywords in aaindex                          219
#TOC>   07       A colour palette for amino acids                  259
#TOC>   08       Extracting R code from Google docs                300
#TOC>   09       Reading Google sheets                             372
#TOC>   10       Plotting amino acids as 2D scatterplot            422
#TOC> 
#TOC> ==========================================================================


# =    01  Install missing packages  ===========================================

if (!requireNamespace("httr", quietly=TRUE)) {
  utils::install.packages("httr")
}

if (!requireNamespace("clipr", quietly=TRUE)) {
  utils::install.packages("clipr")
}

if (!requireNamespace("seqinr", quietly=TRUE)) {
  utils::install.packages("seqinr")
}



# =    02  Load required libraries  ============================================
# ... only if we must.
#
  library(seqinr)  # need to use library() since we load data from the package



# =    03  Generative AI  ======================================================

# ==   03.1  t2c - write text to clipboard  ====================================
cat("  Defining t2c() ...\n")
t2c <- function(txt) {
  #' Convenience function to send R objects to the clipboard.
  #' @examples
    #' t2c("This string was written to the clipboard by t2c().")
    #'
  if (! clipr::clipr_available()) {
    stop("The clipr:: package is not available. This must be fixed.")
  }
  clipr::write_clip(txt)
  return(invisible(NULL))
}


# ==   03.2  Initialize generative AI initial prompt  ==========================

cat("  Defining AIinit() ...\n")

AIinit <- function() {
  #' Loads a prompt to initialize a Generative AI session into the clipboard.
  #' @usage AIinit()    # then just paste the clipboard contents
  #'                    # into the AI assistant session.

  txt <- "
I would like you to act as an R language tutor and answer my prompts to help me learn R. As a tutor, you do the following:
* You write # as the very last character of your response, so I know that the response is complete.
* You keep in mind that I am a programming beginner and you explain syntax and concepts at a novice level.
* You are concise.
* You avoid using packages when a base R function is trivial to write for the purpose.
* When you must use non-standard packages, you write code as package::function() and do not use library(package) if possible.
* You do not use tidyverse functions.
* You use dataFrame[ , col] notation, not dataFrame[[col]] if possible.
* When a direct code solution can address the question, simply provide the code - do not assign the result and print it. For example:
Good:
path.expand(\"~\")

Poor:
home_directory <- path.expand(\"~\")
print(home_directory)

Please confirm with one word. Then we begin.
"
  t2c(txt)
}



# =    04  Remote control of ChimeraX  =========================================
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



# =    05  A progress bar for long-running code  ===============================

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

# Usage example:
if (FALSE) {

N <- 123
for (i in 1:N) {
  pBar(i, N)
  Sys.sleep(0.05)
}

}



# =    06  Find Keywords in aaindex  ===========================================

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

# Usage:
#
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



# =    07  A colour palette for amino acids  ===================================

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

# Usage:
#   Just use the vector like any named dataset.
# barplot(rep(1, 20), col = AACOLS, names.arg = names(AACOLS), cex.names=0.5)

# You can use gsub() to add a value for the colours' transparency (alpha
# channel):
# AACOLS <- gsub("$", "AA", AACOLS)    # Make the colors 33% transparent
# AACOLS <- gsub("$", "80", AACOLS)    # Make the colors 50% transparent
# AACOLS <- gsub("$", "55", AACOLS)    # Make the colors 67% transparent
# AACOLS <- gsub("..$", "FF", AACOLS)  # Reset the two last digits to "FF"
                                       #   to remove transparency



# =    08  Extracting R code from Google docs  =================================

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

  # Usage example:
  fetchGoogleDocRCode("https://docs.google.com/document/d/15qUO3WwKZSqK84gNj8XZIrCe6Ih791oFfGTJ82nuM_w/edit?usp=sharing")

}



# =    09  Reading Google sheets  ==============================================

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




# =    10  Plotting amino acids as 2D scatterplot  =============================

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
       bg = aaDat$col[ord],
       ...)

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
# plotAA(x, y, xlab = aaindex[[150]]$D, ylab = aaindex[[544]]$D)



# [END]
