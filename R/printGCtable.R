# printGCtable.R
#
# Purpose: Define the function printGCtable() to print alternate tables for
#          presenting the genetic code.
# Version: 1.1
# Date:    2023-09-25
# Author:  CSB195 with ChatGPT-4
#
#  THIS SCRIPT CAN BE SOURCE()'D WITHOUT SIDE EFFECT TO LOAD THE FUNCTION.
#
# Versions:
#      1.1  Change to work with both one-letter and three-letter codes
#      1.0  Basic running code
#
# ToDo:
#    - Comment the code
#    - Refactor for coding style (variable names!)
#    - Refactor to use alternative column names
#    - Return the reordered data in a format that is suitable for further
#      processing / printing / analysis
#    - Write test
#
# Notes:
#   The v 1.0 code was developed in class, in a ChatGPT dialogue. The dialogue
#   is here:
#     https://chat.openai.com/share/53e3b709-7067-4842-b61e-5a7010181f28
#
#   The function does not change any external objects (no side effects).
#   That is good practice.
#
#   The function validates user input. We did not ask ChatGPT to do that.
#   Doing so is good practice.
#
#   ChatGPT likes to use pothole_style for names. We prefer camelCase. But
#   when we implicitly suggested printGCtable() for our function in
#   the example function signature we gave it, it used that name
#   going forward. That is good dialogue.
#
#   ChatGPT added an empty line after each fourth row of the table.
#   That structures the output and makes it more readable. We did not ask
#   for that, it came up with that on its own. Doing so is a good idea.
#
#   The variable names are not very good. For example: don't use "c",
#   it is easy to misread for the inbuilt function "c()".
#
# ==============================================================================


# ====  PARAMETERS  ============================================================

GCDAT <- "data/GeneticCode.csv"


# ====  FUNCTIONS  =============================================================

# ChatGPT code, unaltered:
printGCtable <- function(nuc, order, dat, aaColName = "Aaa") {
  # Parameters:
  #   nuc: a permutation of {A, C, G, T} in the order they should appear
  #        in the table
  #   order: a permutation of {1, 2, 3} corresponding which nucleotide
  #        should vary in the tables, rows, columns, and blocks, resspectively
  #   dat: a data frame with the genetic code: required columns are
  #        "First", "Second", "Third", "Codon" and the column name for the
  #        encoded amino acid - either "Aaa" (default) or "A".

  # Check input consistency
  if (length(nuc) != 4 ) {
    stop("Ensure \"nuc\" has length 4.")
  }
  if (length(order) != 3 || !all(order %in% 1:3)) {
    stop("Ensure \"order\" has three unique values between 1 and 3.")
  }
  if (! any(aaColName == colnames(dat))) {
    stop(sprintf("Columname \"%s\" does not appear in the dataframe \"dat\"",
                 aaColName))
  }

  row_pos <- order[1]
  col_pos <- order[2]
  cell_pos <- order[3]

  # Define function to extract value based on position
  extract_val <- function(df_row, pos) {
    if (pos == 1) return(df_row$First)
    if (pos == 2) return(df_row$Second)
    return(df_row$Third)
  }

  # Create nested loops to go through nucleotide combinations based
  # on the given order
  for (r in nuc) {
    for (cell in nuc) {
      line <- sapply(nuc, function(c) {
        # Find corresponding Codon and AAA based on r, c, and cell
        matching_row <- dat[extract_val(dat,  row_pos) == r &
                            extract_val(dat,  col_pos) == c &
                            extract_val(dat, cell_pos) == cell, ]

        paste0(matching_row[ , "Codon"], " ", matching_row[ , aaColName])
      })
      cat(paste(line, collapse="\t"), "\n")
    }
    cat("\n")  # Add an empty line between groups
  }
}



# ====  PROCESS  ===============================================================
#

if (FALSE) {
  GCdf <- read.csv(GCDAT)
  printGCtable(nuc = c("A", "C", "G", "T"), order = c(2, 3, 1), dat = GCdf)


}

# ====  TESTS  =================================================================
#

if (FALSE) {

  # Validate:

  # This should produce the standard way to print the code:
  # ... three-letter code
  printGCtable(nuc = c("T", "C", "A", "G"), order = c(1, 2, 3), dat = GCdf)

  # ... one-letter code
  printGCtable(nuc = c("T", "C", "A", "G"),
               order = c(1, 2, 3),
               dat = GCdf,
               aaColName = "A")


}

# [END]
