# tocID <- ".util.R"
#
# CSB195 Class project: utility scripts
#
# 2022-09  - 2023-10
# boris.steipe@utoronto.ca
#
# This file is source()'d upon startup by .Rprofile
#
# ==============================================================================
#


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                               Line
#TOC> -------------------------------------------------------------------
#TOC>   1        Install missing packages                              50
#TOC>   2        Load required libraries                               70
#TOC>   3        Load datasets                                         74
#TOC>   4        Generic Utilities                                     81
#TOC>   4.1        vr - make a row-vector                              84
#TOC>   4.2        vc - make a column-vector                          101
#TOC>   4.3        mmScale - min-max Scaling                          118
#TOC>   4.4        A progress bar for long-running code               130
#TOC>   4.5        randSeed - large, random seeds                     166
#TOC>   4.6        Random IDs                                         197
#TOC>   5        Generative AI                                        259
#TOC>   5.1        t2c - write text to clipboard                      263
#TOC>   5.2        Initialize generative AI initial prompt            281
#TOC>   6        Working with Google assets                           311
#TOC>   6.1        Extracting R code from Google docs                 314
#TOC>   6.2        Reading Google sheets                              387
#TOC>   7        Remote control of ChimeraX                           460
#TOC>   8        Bioinformatics Utilities                             541
#TOC>   8.1        Find Keywords in aaindex                           544
#TOC>   8.2        A colour palette for amino acids                   583
#TOC>   8.3        Load the genetic code into a data frame            623
#TOC>   8.4        Load an amino acid dataset                         631
#TOC>   8.5        Convert one-letter symbols to three-letter         647
#TOC>   8.6        Amino acid similarity                              729
#TOC>   8.7        Dotplot                                            738
#TOC>   8.8        Plotting amino acids as 2D scatterplot             857
#TOC>   9        Plot Utilities                                       916
#TOC>   9.1        Draw a triangle on an existing plot                918
#TOC> 
#TOC> ==========================================================================


# =    1  Install missing packages  ============================================
#

if (!requireNamespace("httr", quietly=TRUE)) {
  utils::install.packages("httr")
}

if (!requireNamespace("clipr", quietly=TRUE)) {
  utils::install.packages("clipr")
}

if (!requireNamespace("stringr", quietly=TRUE)) {
  utils::install.packages("stringr")
}

if (!requireNamespace("seqinr", quietly=TRUE)) {
  utils::install.packages("seqinr")
}


# =    2  Load required libraries  =============================================
# None needed currently.
#

# =    3  Load datasets  =======================================================
#

cat("  Loading aaindex dataset from sequinr:: ...\n")
utils::data(aaindex, package = "seqinr")


# =    4  Generic Utilities  ===================================================
#

# ==   4.1  vr - make a row-vector  ============================================
#

cat("  Defining vr() ...\n")
vr <- function(v, rowName, colNames) {
  #' Drop dimension of input vector and make a row-vector. Keep names().
  #' @examples
  #' vr(c(1, 1, 2, 3, 5, 8, 13, 21))
  if (missing(rowName))  { rowName  <- NULL }
  if (missing(colNames)) { colNames <- names(v) }
  rr <- matrix(v, nrow = 1)
  rownames(rr) <- rowName
  colnames(rr) <- colNames
  return(rr)
}


# ==   4.2  vc - make a column-vector  =========================================
#

cat("  Defining vc() ...\n")

vc <- function(v, rowNames, colName) {
  #' Drop dimension of input vector and make a column-vector. Keep names().
  #' @examples
  #' vc(c(1, 1, 2, 3, 5, 8, 13, 21))
  if (missing(rowNames)) { rowNames <- names(v) }
  if (missing(colName))  { colName  <- NULL }
  cc <- matrix(v, ncol = 1)
  rownames(cc) <- rowNames
  colnames(cc) <- colName
  return(cc)
}

# ==   4.3  mmScale - min-max Scaling  =========================================
#

cat("  Defining mmScale() ...\n")

mmScale <- function(x) {
  #' min-max scaling of x. The returned vector has min(x) == 0 and max(x) == 1
  #' @examples
  #' summary(mmScale(rnorm(10000)))
  return((x - min(x)) / (max(x) - min(x)))
}

# ==   4.4  A progress bar for long-running code  ==============================
#

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

# ==   4.5  randSeed - large, random seeds  ====================================
#

cat("  Defining randSeed() ...\n")

randSeed <- function(verbose = TRUE) {
  #' Return a random, positive integer between 1 and the largest integer that
  #' can be represented on the system, usually 2,147,483,647. Also print that
  #' number, if desired.
  #' @examples
  #' randSeed()
  oldSeed <- .Random.seed                      # Save the current state
  set.seed(NULL)                               # Randomize the PRNG
  mySeed <- sample(0:.Machine$integer.max, 1)  # Fetch a value
  .Random.seed <- oldSeed                      # Restore the PRNG
  if (verbose) { print(sprintf("Seed: %i",     # print the number
                               mySeed)) }
  return(mySeed)
}

if (FALSE) {
  randSeed()
  set.seed((x <- randSeed()))  # x gets defined, but the value also gets used
  sample(0:9)    # This is the result we want to reproduce
  tmp <- runif(100)   # run the PRNGn for a bit
  sample(0:9)    # Do we now get a different permutation? Yes.
  set.seed(x)    # Reset the PRNG seed
  sample(0:9)    # should be the same as before

}

# ==   4.6  Random IDs  ========================================================
#

cat("  Defining rID() ...\n")

rID <- function(n = 1, l = 5, mode = "alf") {
  # Create n random IDs of length l, from a choice of alphabets
  # Modes:                                         Keyspace for length 5:
  #   dec: decimal                                 (1.00e+05)
  #   hex: hexadecimal                             (1.05e+06)
  #   let: letters                                 (1.19e+07)
  #   alf: alphabetic - LETTERS, letters, 0:9      (9.16e+08)
  #   asc: printable ASCII except for:             (5.90e+09)
  #          34: "
  #          35: #
  #          39: '
  #          96: `
  #         127: DEL

  if (mode == "dec") {          # 0123456789
    a <- as.character(0:9)
  } else if (mode == "hex") {   # dec + ABCDEF
    a <- c(0:9, LETTERS[1:6])
  } else if (mode == "let") {   # abcdefghijklmnopqrstuvwxyz
    a <- letters
  } else if (mode == "LET") {   # ABCDEFGHIJKLMNOPQRSTUVWXYZ
    a <- LETTERS
  } else if (mode == "alf") {   # LET + let + dec
    a <- c(LETTERS, letters, 0:9)
  } else if (mode == "asc") {   # ASCII characters
    a <- 32:127                                # printable ASCII range
    a <- a[! (a %in% c(34, 35, 39, 96, 127))]  # remove problematic characters
    a <- unlist(strsplit(intToUtf8(a), ""))    # convert to character array
  } else {
    stop("Unsupported mode")
  }

  ids <- character(n)
  for (i in 1:n) {
    ids[i] <- paste(sample(a, l, replace = TRUE), collapse = "")
  }
  return(ids)

}

# Usage example:
if (FALSE) {

  rID()
  rID(l=9, mode = "dec")
  rID(l=9, mode = "hex")
  rID(l=9, mode = "let")
  rID(l=9, mode = "LET")
  rID(l=9, mode = "alf")
  rID(l=9, mode = "asc")
  rID(mode = "qrk")
  rID(15)
  barplot(rep(1, 20), col = sprintf("#%s", rID(20, l = 6, mode = "hex")))

}


# =    5  Generative AI  =======================================================
#


# ==   5.1  t2c - write text to clipboard  =====================================
#

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


# ==   5.2  Initialize generative AI initial prompt  ===========================
#

cat("  Defining gAIinit() ...\n")

gAIinit <- function() {
  #' Loads a prompt to initialize a Generative AI session into the clipboard.
  #' @usage gAIinit()    # then just paste the clipboard contents
  #'                     # into the AI assistant session.

  txt <- "
I would like you to act as an R language tutor and answer my prompts to help me learn R. As a tutor, you do the following:

* You write # as the very last character of your response, so I know that the response is complete.
* You keep in mind that I am a programming beginner and you explain syntax and concepts at a novice level.
* You are concise.
* You avoid using packages when a base R function is trivial to write for the purpose.
* when you mention functions, you identify them by typing parentheses after the name: e.g. rnorm(), c()
* when you mention packages, you identify them by typing two colons after the name: e.g. httr::, utils::
* When you must use non-standard packages, you write code with package::function() and do not use library(package) if possible.
* You do not use tidyverse functions.
* You use dataFrame[ , col] notation, not dataFrame[[col]] if possible.
* When the return value of an expression is the required output, you simply provide expression. You do not wrap it in a print() statement, or assign it to a variable. For example: path.expand(\"~\") NOT print(path.expand(\"~\")) and NOT home_directory <- path.expand(\"~\"); print(home_directory).

Please confirm with one word.
"
  t2c(txt)
}


# =    6  Working with Google assets  ==========================================
#

# ==   6.1  Extracting R code from Google docs  ================================
#

cat("  Defining fetchGoogleDocRCode ...\n")

fetchGoogleDocRCode <- function (URL,
                                 delimB = "^\\s*# begin code",
                                 delimE = "^\\s*# end code",
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


# Usage example:
if (FALSE) {

  fetchGoogleDocRCode("https://docs.google.com/document/d/15qUO3WwKZSqK84gNj8XZIrCe6Ih791oFfGTJ82nuM_w/edit?usp=sharing")

}


# ==   6.2  Reading Google sheets  =============================================
#

cat("  Defining readGsheet() ...\n")

readGsheet <- function(URL, sheet, ...) {
  # Read a sheet from a Google sheets URL to a spreadsheet file.
  # URL: the URL of file location. Note: the document must have permissions
  #      set for everyone with the link to be allowed to read.
  # sheet: the name of the sheet
  # ... : other arguments that are passed to utils::read.csv()
  # value: a data frame

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
  if (httr::status_code(response) != 200) {
    stop(sprintf("Server status code was \"%s\".",
                 as.character(httr::status_code(response))))
  }

  # Check whether reading the sheet returned "text/csv". Otherwise it may
  # not be a spreadsheet, or we are getting a sign-in page (i.e. a permissions
  # error).
  if (! grepl("text/csv", response$headers$`content-type`) ) {
    stop("
  Request did not return \"text/csv\" content. This could mean:
    -  Sheet access was denied and we got a sign-in page. Check permissions.
       The sheet needs to readable by \"everyone with the link\".
    -  The URL does not lead to a spreadsheet. Check it.
    -  Something else happened. Check what the URL actually retrieves.")
  }

  x <-as.character(response)
  x <- strsplit(x, "\n")[[1]]
  x[1] <-      gsub("\\\"", "", x[1])
  cNames <-  strsplit(x[1], ",")[[1]]  # Restore the original names from
  # row 1 so we don't have to use the
  # names that R assigns by default.
  tbl <- utils::read.csv(text = x, ...)
  colnames(tbl) <- cNames

  return(tbl)

}

# Usage example
if (FALSE) {
  # This should work ...
  x <- "1tRCPhaua5cjcH_0DuZOiv8BVbdr_V6miC2JeKiOYj-o"
  x <- sprintf("https://docs.google.com/spreadsheets/d/%s/edit?usp=sharing", x)
  y <- "AA styles"
  z <- readGsheet(x, y)

  # This should fail (permissions) ...
  x <- "1E4GHlgRlTU5Z1FdglNF5T_eR4-Guzkukq99stl35_I0"
  x <- sprintf("https://docs.google.com/spreadsheets/d/%s/edit?usp=sharing", x)
  y <- "GC.csv"
  z <- readGsheet(x, y)

}


# =    7  Remote control of ChimeraX  ==========================================
#

CXPORT <- 61803
cat(sprintf("  Defined ChimeraX port (CXPORT) as %d.\n", CXPORT))


cat("  Defining CX() ...\n")

CX <- function(cmd, port = CXPORT, quietly = FALSE) {
  # send a command to ChimeraX listening on port CXPORT via its REST
  # interface.
  # Parameters:
  #   cmd      char     a ChimeraX command line command
  #   port     int      the port number on which ChimeraX is listening
  #   quietly  logical  if FALSE, cat() the contents of the response
  #
  # Value:  the reply by ChimeraX, invisibly.

  # (1) construct the base address, port, and command. This is the prefix
  #     that httr::GET() needs to send the command to the right port on
  #     the right machine:
  CXREST <- sprintf("http://127.0.0.1:%s/run?", CXPORT)

  # (2) sanitize the user-entered variable cmd so it can be properly encoded
  #     in a http message. (No blanks, some specially handled characters, ...)
  #     Here we use gsub(), which seeks for a pattern defined in its first
  #     argument, substitutes the contents of the second argument, in the
  #     string that is identified with the third argument.
  #
  #     Patterns are defined as "regular expressions".
  #
  cmd <- gsub("(^\\s+)|(\\s+$)", "", cmd)  # trim whitespace

  # "percent encode" reserved characters
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

  # send the command to ChimeraX, and capture the response in the variable r:
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


# =    8  Bioinformatics Utilities  ============================================
#

# ==   8.1  Find Keywords in aaindex  ==========================================
#

cat("  Defining grepAAindex() ...\n")

grepAAindex <- function(key, el = "D") {
  # Search for strings in an aaindex element
  # key:   a regular expression pattern
  # el:    which element to search. Default "D" (Definition), also useful might
  #        be "T" (Titles).
  # value: a character vector, length 0 if nothing found.

  data(aaindex, package = "seqinr")

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


# ==   8.2  A colour palette for amino acids  ==================================
#

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

# ==   8.3  Load the genetic code into a data frame  ===========================
#

cat("  Loading dataset GCdf from ./data/GeneticCode.csv ...\n")

GCdf <- utils::read.csv("data/GeneticCode.csv")


# ==   8.4  Load an amino acid dataset  ========================================
#

cat("  Loading reference dataset AADAT from a Google sheet (Course Data) ...\n")

URL <- paste(c("https://docs.google.com/spreadsheets/d/",
               "1tRCPhaua5cjcH_0DuZOiv8BVbdr_V6miC2JeKiOYj-o",
               "/edit?usp=sharing"),
             collapse = "")
sheet <- "Data"

AADAT <- readGsheet(URL, sheet)
rownames(AADAT) <- AADAT$A



# ==   8.5  Convert one-letter symbols to three-letter  ========================
#

cat("  Defining A2Aaa ...\n")

A2Aaa <- function(aa, out = ifelse(nchar(aa[1]) == 3, "A", "Aaa") ) {
  # aa:   a vector of amino acid one letter or three letter symbols,
  #       or IUPAC names, to be converted.
  # out:  default: convert three-letter symbols to one-letter,
  #          everything else to three-letter.
  #       "A"    - convert all input to IUPAC one letter symbols
  #       "Aaa"  - convert all input to IUPAC three letter symbols
  #       "Name" - convert all input to IUPAC  name
  # Note: Whereas IUPAC uses Aspartic acid and Glutamic acid, we always use
  #       the names of the ionized forms aspartate and glutamate.
  #
  #       IUPAC uses capitalized names in their table, but lowercase
  #       in text. We use lower-case throughout for names, upper case
  #       for one-letter symbols, title case for three-letter symbols.
  #
  #       If your input is well controlled and your code needs to be
  #       fast, use the following idioms directly:
  #       (Assuming your variable name is x)
  #       Initialize: DAT <- GCdf[! duplicated(GCdf$A), c("A", "Aaa")]
  #       x <- DAT$A[match(x, DAT$Aaa)]   # for three-to-one conversion
  #       x <- DAT[x, "Aaa"]                # for one-to-three conversion

  stopifnot(out %in% c("A", "Aaa", "Name"))

  DAT <- GCdf[! duplicated(GCdf$A), c("A", "Aaa", "Name")]
  results <- character(length(aa))

  if (all(nchar(aa) == 1)) {
    source <- "A"
  } else if (all(nchar(aa) == 3)) {
    source <- "Aaa"
  } else {
    source <- "Name"
  }
  target <- out

  results <- DAT[match(aa, DAT[ , "A"]), target]

  sel <- is.na(results)
  if (any(sel)) {    # repeat for Aaa
    results[sel] <- DAT[match(aa[sel], DAT[ , "Aaa"]), target]

    sel <- is.na(results)
    if (any(sel)) {  # repeat for Name
      results[sel] <- DAT[match(aa[sel], DAT[ , "Name"]), target]
    }
  }

  if (any(is.na(results))) {
    iNA <- which(is.na(results))[1]
    stop(sprintf("Symbol %i (\"%s\") is not in the symbol table.",
                 iNA,
                 aa[iNA]))
  }

  return(results)

}

if (FALSE) {

  A2Aaa("")                              # Error
  A2Aaa("Q")                             # "Gln"
  A2Aaa("Leu")                           # "L"
  A2Aaa("*")                             # "***"
  A2Aaa("K", out = "Name")               # "lysine"
  A2Aaa("Phe", out = "Name")             # "phenylalanine"
  A2Aaa("tyrosine")                      # "Tyr"
  A2Aaa("tyrosine", out = "A")           # "Y"
  A2Aaa(c("K", "L", "H"))                # "Lys" "Leu" "His"
  A2Aaa(c("K", "L", "H"), out = "A")     # "K" "L" "H"
  A2Aaa(c("Trp", "Y", "alanine", "Q"))   # "W" "Y" "A" "Q"
  A2Aaa(c("Trp", "Tyr", "Qrk", "Ala"))   # Error
  A2Aaa(c("Quackophane"))                # Error

}

# ==   8.6  Amino acid similarity  =============================================
#

cat("  Defining aaSim() ...\n")

source("./R/aaSim.R")



# ==   8.7  Dotplot  ===========================================================
#

cat("  Defining dotPlot2() ...\n")


dotPlot2 <- function(A, B,        # sequence vectors
                     f,           # filter
                     MDM,         # A Mutation Data Matrix
                     palette,     # a function that returns color values
                     xlab = "",
                     ylab = "") {
  # Purpose:
  #     Create a dotplot to measure sequence similarity between
  #     two amino acid sequences
  # Version:  1.0
  # Date:     2016-09
  # Author:   Boris Steipe
  #
  # Parameters:
  #     A, B: vectors that contain no letters that are not found in MDM
  #     f: filter matrix to weight an average around the neighborhood of
  #        an amino acid pair. Default to the identity matrix if missing.
  #        Average over a window of length f if length(f) is 1.
  #     MDM: A mutation Data matrix. If missing, create one from aaSim()
  #     palette: rainbow(), cm.colors(), or another function that returns
  #              a palette of color hexcodes. If missing, make our own
  #              palette.
  # Value:
  #     none. creates a dotplot.


  if (missing(f)) {
    f <- matrix(1) # default
  } else if (length(f) == 1) {
    if (! f %% 2) {stop("Sorry: f must be odd.")}
    w <- f
    f <- matrix(numeric(w * w), nrow = w)
    for (i in 1:w) { f[i, i] <- 1 }  # identity matrix
  }

  if (missing(MDM)) {
    aa <- c(AADAT$A, "*")
    MDM <- matrix(numeric(length(aa)^2), nrow = length(aa))
    for (i in 1:length(aa)) {
      for (j in 1:length(aa)) {
        MDM[i,j] <- aaSim(aa[i], aa[j])
      }
    }
    rownames(MDM) <- aa
    colnames(MDM) <- aa
  }


  if (missing(palette)) {
    palette <- colorRampPalette(c("#FE333D", # red
                                  "#F09371",  # orange
                                  "#CCCCCC",  # grey
                                  "#D7D7D7",  # grey
                                  "#E1E1E1",  # grey
                                  "#ECECEC",  # grey
                                  "#F6F6F6"), # white
                                bias = 0.8)
  }
  lA <- length(A)
  lB <- length(B)

  m <- matrix(numeric(lA * lB), nrow = lA, ncol = lB)
  for (i in 1:lA) {
    for (j in 1:lB) {
      m[i, j] <- MDM[A[i], B[j]]
    }
  }
  m2 <- m
  wr <- floor((dim(f)[1] - 1) / 2)  # half-window size for rows
  wc <- floor((dim(f)[2] - 1) / 2)  # half-window size for columns

  for (i in (wr + 1):(lA - wr)) {
    for (j in (wc + 1):(lB - wc)) {
      # apply the filter to each value in m by weighting and summing
      # over its wr x wc neighborhood. Put the new value in m2
      m2[i, j] <- sum(f * m[(i-wr):(i+wr), (j-wc):(j+wc)])
    }
  }
  image(1:lA, 1:lB, m2,
        col = palette(24),
        ylim=c(lB,1), xlim=c(1,lA),
        xlab = xlab,
        ylab = ylab,
        axes = FALSE)
  box()

  # find good values for axis ticks and gridlines
  steps <- c(1, 2, 5, 10, 20, 50, 100, 200, 500,
             1000, 2000, 5000, 10000, 20000, 50000)
  gridStep <- sum(steps < max(lA, lB))

  # draw axes
  axis(1, at = c(1, seq(steps[gridStep - 3], lA, by=steps[gridStep - 3])))
  axis(2, at = c(1, seq(steps[gridStep - 3], lB, by=steps[gridStep - 3])))
  axis(3, at = c(1, seq(steps[gridStep - 3], lA, by=steps[gridStep - 3])))
  axis(4, at = c(1, seq(steps[gridStep - 2], lB, by=steps[gridStep - 3])))

  # draw grid with thin, transparent lines
  for (pos in seq(steps[gridStep - 2], lA, by = steps[gridStep - 2])) {
    abline(v=pos, col = "#FFFFFF44", lwd = 0.5)
  }
  for (pos in seq(steps[gridStep - 2], lB, by = steps[gridStep - 2])) {
    abline(h=pos, col = "#FFFFFF44", lwd = 0.5)
  }

  return(invisible(m2))

}





# ==   8.8  Plotting amino acids as 2D scatterplot  ============================
#

cat("  Defining plotAA() ...\n")

plotAA <- function(x, y, aaDat = AADAT, ...) {
  # Plot amino acids in 2D
  # x, y: numeric vectors named vectors as Aaa three letter symbol. They can be
  #       passed in any order, but the order must be the same for x and y
  # ... : other arguments that are passed to plot()
  # aaDat: a dataframe with plotting properties.
  #        (Default: AADAT - defined above)
  # value: none. A plot is produced.

  o <- names(x)
  if (! all(o %in% aaDat$Aaa)) {
    stop("Names of first vector do not match three-letter symbols.")
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

# Usage example
if (FALSE) {

  data(aaindex, package = "seqinr")
  x <- aaindex[[150]]$I
  y <- aaindex[[544]]$I
  plotAA(x, y, xlab = aaindex[[150]]$D, ylab = aaindex[[544]]$D)

}

# =    9  Plot Utilities  ======================================================

# ==   9.1  Draw a triangle on an existing plot  ===============================


cat("  Defining triangle() ...\n")

triangle <- function(xT, yT,          # tip coordinates
                     frac = 0.05,     # size as fraction of x-axis width
                     asp = 1.5,       # aspect ratio: > 1 for taller triangles
                     tip = "below",   # tip above or below base?
                     ...){            # additional parameters for polygon()

  # Adds an equilateral triangle to the current plot. Tip is at (xT, yT).
  # Side length is frac fraction of x-axis width in user coordinates.
  # Triangle can point up or down.


  usr <- par("usr")                           # Get user coordinates
  dx <- usr[2] - usr[1]                       # Calculate x-axis range
  dy <- usr[4] - usr[3]                       # Calculate y-axis range
  aT <- frac * dx                             # Calculate side length
  hT <- (sqrt(3) / 2) * aT * (dy / dx) * asp  # Calculate the height
  ht <- ifelse(tip == "above", -hT, hT)       # Point triangle up or down

  vT <- matrix(nrow = 3, ncol = 2)            # Store three vertex coordinates
  vT[1,] <- c(xT, yT)                         # Tip (middle)
  vT[2,] <- c(xT - aT / 2, yT + hT)           # Left
  vT[3,] <- c(xT + aT / 2, yT + hT)           # Right
  rownames(vT) <- c("tip", "left", "right")   # Add row names, and
  colnames(vT) <- c("x", "y")                 # ... column names

  polygon(vT, ...)                            # Draw the triangle
  return(invisible(vT))                       # Return the vertices

}




# [END]
