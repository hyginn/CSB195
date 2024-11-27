# tocID <- "./courseScripts/14-parsingAnnotations.R"
#
#
# Purpose: Using reglar expressions to parse phenotype annotations from
#          ensembl phenotype data.
#
# Version: 1.1
# Date:    2023-10 - 2024-11
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.1   Add original ensembl IDs as rownames to datasets to be able to
#         correlate results back with original STRING graph data. Some code
#         cleanup.
#   1.0   New script, first version used in course.
#
# ToDo:
#
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                                Line
#TOC> --------------------------------------------------------------------
#TOC>   1        Preparations                                           41
#TOC>   2        Processing ensembl Genes                               50
#TOC>   2.1        Define list of IDs in dataset                        59
#TOC>   2.2        Retrieve phenotype information from biomart          66
#TOC>   2.3        Process phenotype information                        98
#TOC>   2.3.1          Remove not-available data                       104
#TOC>   2.3.2          Process annotations into individual terms       128
#TOC>   2.3.3          Keep only unique terms                          145
#TOC>   3        Parsing ...                                           154
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
  # define which biomart dataset to use
  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  # read STRING dataset from course folder
  STRINGedges <- readRDS("./data/STRINGedges.rds")

# ==   2.1  Define list of IDs in dataset  =====================================
  #
  # Get the unique Ensembl identifiers from this list of interactions:
  IDs <- unique(c(STRINGedges$a, STRINGedges$b))

  N <- length(IDs)

# ==   2.2  Retrieve phenotype information from biomart  =======================


  myPhenotypes <- data.frame(ID = character(N),
                             sym = character(N),
                             phe = character(N))
  tStart <- Sys.time()
  for (i in 1:N) {
    pBar(i, N) # ... progress bar
    ID <- IDs[i]
    myPhenotypes$ID[i] <- ID
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


# ==   2.3  Process phenotype information  =====================================

  # To get a vector with all phenotype terms, we first need to split our saved
  # strings apart again, and concatenate the individual terms together.
  pheno <- readRDS("./data/ensemblPhenotypes.rds")

# ===   2.3.1  Remove not-available data                  

  # Let's first remove all rows from the dataframe that have NA values in the
  # phenotypes. Those are genes whose function has not been well characterized.
  #
  # Note: biomart wrote the string "NA" into its output. We used the inbuilt
  #       value of NA (i.e. Not Available) instead wherever the STRING ID
  #       was not found in biomart. Thus: we have NA wherever the gene
  #       symbol was not available - the gene was no longer in ensembl,
  #       probably because of a database update, and we have "NA" wherever
  #       no phenotypes had been annotated to that gene.

  sum(is.na(pheno$phe))                    # 142 missing symbols
  sum(pheno$phe == "NA", na.rm = TRUE)     # 5,632 unannotated genes


  head(pheno)                                   # inspect

  sel <- is.na(pheno$phe) | (pheno$phe == "NA")   # select NA and "NA"
  sum(sel)                                        # count and validate
  pheno <- pheno[! sel, ]                         # remove the selected rows

  head(pheno)                                   # confirm

# ===   2.3.2  Process annotations into individual terms  

  # Next we use strsplit() to break all the phenotype strings that we had
  # paste()'ed together with ", " apart again. strsplit() produces a list
  # with one vector of results per list element.
  pheno <- strsplit(pheno$phe, ", ")
  pheno[[1]]

  # We can simply use unlist() to turn this list into a vector that contains
  # all the terms at once:
  pheno <- unlist(pheno) # 68,823 annotations

  # We could study these further, for example by using table() to find the most
  # frequently used annotations
  tPhe <- sort(table(pheno), decreasing = TRUE)
  tPhe[1:10]  # ... many cancer annotations!

# ===   2.3.3  Keep only unique terms                     

  # Finally we use unique() to get our vector of unique phenotype terms:
  pheno <- unique(pheno)  # 10,252 unique phenotypes.

  # We save those ... and we can start parsing the terms
  # saveRDS(pheno, "./data/uniquePhenotypes.rds")
}

# =    3  Parsing ...  =========================================================

if (FALSE) {
  myPhenotypes <- readRDS("./data/uniquePhenotypes.rds")

  # Now what ... ?
  #
  # We need to find regular expressions that identify all cancer associated
  # phenotypes. Then we can look back into our STRING genes and find out which
  # of those are cancer associated, and compute enrichment scores.

  # Example
  hits <- grep("[Cc]ancer", myPhenotypes)
  myPhenotypes[hits]


  # Task:
  # Find a way to identify many (all?) cancer- associated phenotype
  # annotations in our list, using regular expressions.

}

# [END]
