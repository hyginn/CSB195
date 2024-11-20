# tocID <- "./courseScripts/14-parsingAnnotations.R"
#
# Using reglar expressions to parse phenotype annotations from ensembl
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                             Line
#TOC> -------------------------------------------------
#TOC>   1        Preparations                        19
#TOC>   2        Processing ensembl Genes            28
#TOC>   3        Parsing ...                         98
#TOC>
#TOC> ==========================================================================


# =    1  Preparations  ========================================================

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

# =    2  Processing ensembl Genes  ============================================

if (FALSE) {
  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")


  # ====
  STRINGedges <- readRDS("./data/STRINGedges.rds")

  # Get the unique Ensembl identifiers from this list of interactions:
  IDs <- unique(c(STRINGedges$a, STRINGedges$b))

  N <- length(IDs)

  myPhenotypes <- data.frame(sym = character(N),
                             phe = character(N))
  tStart <- Sys.time()
  for (i in 1:N) {
    pBar(i, N) # ... progress bar
    ID <- IDs[i]
    tmp <- biomaRt::getBM(filters = "ensembl_peptide_id",
                          attributes = c("hgnc_symbol",
                                         "phenotype_description"),
                          values = ID,
                          mart = myMart)
    if (length(tmp$hgnc_symbol) == 0) {
      myPhenotypes$sym[i] <- NA
      myPhenotypes$phe[i] <- NA
    } else {
      myPhenotypes$sym[i] <- unique(tmp$hgnc_symbol)
      myPhenotypes$phe[i] <- paste(unique(tmp$phenotype_description),
                                   collapse = ", ")
    }
  }
  tEnd <- Sys.time()
  cat(sprintf("%d genes processed. Time taken: %s\n", N, tEnd - tStart))
  # ... this took about 2 and a half hours.

  # saveRDS(myPhenotypes, "./data/ensemblPhenotypes.rds")

  # To parse the phenotypes, we need to split them apart, and
  # concatenate them together

  # Let's first remove all rows from the dataframe that have NA values in the
  # phenotypes. Those are genes whose function has not been well characterized.

  tmp <- myPhenotypes
  sel <- is.na(tmp$phe) | (tmp$phe == "NA")
  sum(sel)
  tmp <- tmp[! sel, ]

  head(tmp)

  # Next we strsplit()
  tmp <- strsplit(tmp$phe, ", ")
  tmp[[1]]

  # We can simply use unlist() to turn this into a vector
  tmp <- unlist(tmp) # 68,823 annotations

  # We could study these further, for example by using table() to find the most
  # frequently used annotations.

  # Next we use unique() to get our phenotype list:
  tmp <- unique(tmp)  # 10,252 unique phenotypes.

  # We save those ... and we can start parsing the terms
  # saveRDS(tmp, "./data/uniquePhenotypes.rds")
}

# =    3  Parsing ...  =========================================================

if (FALSE) {
  myPhenotypes <- readRDS("./data/uniquePhenotypes.rds")

  # Now what ... ?

  # Example
  hits <- grep("[Cc]ancer", myPhenotypes)
  myPhenotypes[hits]


  # Task:
  # Find a way to identify many (all?) cancer- associated phenotype
  # annotations in our list, using regular expressions.

}

# [END]
