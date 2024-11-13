# tocID <- "courseScripts/12-proteinInteractionAnalysis.R"
#
#
# Purpose:  R code for analysis of the human functional protein
#           interaction network.
#
# Version:   2.0
#
# Date:     2017-08  -  2024-11
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           2.0    Focus more on graph-topological measures and clearly work
#                  out the enriched properties. Make it source() - safe. Update
#                  STRING data from latest (v12). Add gAI hints.
#           1.4    Update vector ID's for betweenness centrality.
#           1.3    Bugfix: called the wrong function on ENSPsel in l. 220
#           1.2    2020 Updates; Rewrite for new STRINg V11;
#                  Deprecate save()/load() for saveRDS()/readRDS()
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
#           1.0    First live version
#           0.1    First code copied from 2016 material.
#
# TODO:
#
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                           Line
#TOC> ---------------------------------------------------------------
#TOC>   1        Setup and data                                    63
#TOC>   1.1        Processing String database data                 76
#TOC>   1.2        Read processed STRING edges                    168
#TOC>   2        Functional Edges in the Human Proteome           183
#TOC>   2.1        Cliques                                        251
#TOC>   2.2        Communities                                    296
#TOC>   2.3        Betweenness Centrality                         314
#TOC>   3        biomaRt                                          364
#TOC>   3.1        The most central proteins ...                  451
#TOC> 
#TOC> ==========================================================================


################################################################################
#                                                                              #
#                            S O U R C E   S A F E                             #
#                                                                              #
#     This script is  "source safe".  source()'ing the  code will  define      #
#     (refresh) the global parameters, and define the functions. It won't      #
#     actually run code. Execute the statements in the if (FALSE) { ... }      #
#     blocks to run the interactive parts of this script.                      #
#                                                                              #
#     To make changes in the code and experiment with it, make your own        #
#     copy in your ./myScripts folder.                                         #
#                                                                              #
################################################################################


# =    1  Setup and data  ======================================================


# Not surprisingly, the analysis of PPI networks needs iGraph:

if (! requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
# Package information:
#  library(help = igraph)       # basic information
#  browseVignettes("igraph")    # available vignettes
#  data(package = "igraph")     # available datasets

# ==   1.1  Processing String database data  ===================================
#

# In order for you to explore some real, biological networks, I give you a data
# frame of functional relationships of human proteins that I have downloaded
# from the STRING database. The full table has 8.5 million records, but I
# processed it to contain only a selection. You do not have to repeat these
# steps. If you want to skip ahead to the next section, that's fine. Ther code
# below is just for illustration how to do this and what considerations might
# apply.
if (FALSE) {
  # 1: Navigated to https://string-db.org/cgi/download
  # 2: Download "9606.protein.links.v12.0.txt.gz" This is a compressed file of
  #    edges in the human interaction network with confidence scores.
  # 3: Read this into a dataframe. (This won't work unless you download the
  #    (large) file.)
  STRINGedges <- read.table("9606.protein.links.v12.0.txt.gz",
                            header = TRUE,
                            sep = " ")
  # STRINGedges is now a data frame with 13,715,404 edges and 3 variables
  #
  # 4: Adjust columnames
  colnames(STRINGedges) <- c("a", "b", "score")
  head(STRINGedges)
  #  a                    b score
  #  1397 9606.ENSP00000000412 9606.ENSP00000438085   993
  #  1602 9606.ENSP00000000412 9606.ENSP00000349437   991
  #  2319 9606.ENSP00000001008 9606.ENSP00000482075   990
  #  2396 9606.ENSP00000001008 9606.ENSP00000360609   999
  #  2440 9606.ENSP00000001008 9606.ENSP00000302961   998
  #  2491 9606.ENSP00000001008 9606.ENSP00000231509   996

  # 5: Remove taxonomy ID substrings from the IDs
  #
  #   STRING has appended the taxonomy-ID for Homo sapiens - 9606 - to the
  #   Ensemble transcript identifiers that start with ENSP. We'll remove them:

  STRINGedges$a <- gsub("^9606\\.", "", STRINGedges$a)
  STRINGedges$b <- gsub("^9606\\.", "", STRINGedges$b)

  head(STRINGedges)

  # 6: Remove duplicate edges
  #
  #   Are the edges unique? Or are the edges duplicated, once for the
  #   forward direction, once for the backward direction? To find this out, we
  #   do the following:
  #   a: we lexical-sort each row so that the "smaller" ID is in column a, the
  #      "larger" of the pair is in column b.
  sel <- STRINGedges$a > STRINGedges$b      # a > b: these need to be swapped
  sum(sel)                                  # 6857702 pairs
  tmp <- STRINGedges$a[sel]                 # save the IDs from a
  STRINGedges$a[sel] <- STRINGedges$b[sel]  # overwrite a with b
  STRINGedges$b[sel] <- tmp                 # put the saved values into b

  #       confirm:
  sel <- STRINGedges$a > STRINGedges$b      # no wrongly ordered pairs remain
  sum(sel)
  sel <- STRINGedges$a == STRINGedges$b     # ... and there are no self-edges
  sum(sel)
  #    b: we produce a vector of strings that concatenate the IDs from
  #       columns a and b
  tmp <- paste(STRINGedges$a, STRINGedges$b)
  head(tmp, 3)

  # [1] "ENSP00000000233 ENSP00000356607"
  #     "ENSP00000000233 ENSP00000427567"
  #     "ENSP00000000233 ENSP00000253413"
  #
  #    c: We find duplicate strings. These mark the positions of edges that
  #       appear twice, regardless of the order the proteins originally had.
  sel <- duplicated(tmp)
  sum(sel)   # 6,857,702
  #    d: We keep only the rows that do NOT contain duplicate edges.
  STRINGedges <- STRINGedges[! sel, ]

  # 6: Select a high-confidence subset
  #    Select a subset with score >= 990
  sel <- STRINGedges$score >= 960
  sum(sel) # 51,826 edges selected
  #   Keep only those records:
  STRINGedges <- STRINGedges[sel, ]

  # 7: How many unique proteins do we have?
  length(unique(c(STRINGedges$a, STRINGedges$b)))  # 9,858

  # 8: Save for future use.
  # saveRDS(STRINGedges, "./data/STRINGedges.rds")
} # End - preparing STRING data



# ==   1.2  Read processed STRING edges  =======================================

# The selected set of high-confidence edges from the functional human protein
# interaction network is a dataframe with about 50,000 edges and about 10,000
# unique proteins. Incidentaly, that's about the size of a fungal proteome. Load
# the saved dataframe here to continue with the script. (To read more about what
# the scores mean, see http://www.ncbi.nlm.nih.gov/pubmed/15608232 ).

STRINGedges <- readRDS("./data/STRINGedges.rds")

if (FALSE) {
  head(STRINGedges)
  str(STRINGedges)
}

# =    2  Functional Edges in the Human Proteome  ==============================


# There are many possibilities to explore interesting aspects of biological
# networks, we will keep with some very simple procedures here but you have
# to be aware that this is barely scratching the surface of possibilities.
# However, once the network exists in your computer, it is comparatively
# easy to find information online about the many, many options you could explore to analyze this network for interesting aspects. Just ask ChatGPT.


# Make a graph from this dataframe

if (FALSE) {
  ?igraph::graph_from_data_frame
}


gSTR <- igraph::graph_from_data_frame(STRINGedges, directed = FALSE)

# CAUTION you DON'T want to plot a graph with 10,000 nodes and 50,000 edges -
# layout of such large graphs is possible, but requires specialized code. Google
# for <layout large graphs> if you are curious. Also, consider what one can
# really learn from plotting such a graph ...

# Of course simple computations on this graph are reasonably fast:


if (FALSE) {
  compSTR <- igraph::components(gSTR)
  summary(compSTR) # our graph is fully connected!

  hist(log(igraph::degree(gSTR)), col="#FEE0AF", breaks = 50)
  # this actually does look rather scale-free

  (freqRank <- table(igraph::degree(gSTR)))

  # let's put that into a dataframe
  freqDF <- data.frame(logRank = log10(as.numeric(names(freqRank)) + 1),
                       logFreq = log10(as.numeric(freqRank)))

  plot(freqDF$logRank,
       freqDF$logFreq,
       type = "b",
       pch = 21, bg = "#FEE0AF",
       xlab = "log(Rank)", ylab = "log(frequency)",
       main = "9,858 nodes from the human functional interaction network")

  # This looks quite scale-free indeed,

  (regressionLine <- lm(freqDF$logFreq ~ freqDF$logRank))
  abline(regressionLine, col = "#AA0000")


  # ...  except perhaps for a break around the
  # eleventh rank - could these be two distinct populations? This would make our
  # question a "mixture deconvolution problem" - but I can't immediately
  # see a biological reason for this.

  (regA <- lm(freqDF$logFreq[1:12] ~ freqDF$logRank[1:12]))
  abline(regA, col = "#00AAFF")

  (regB <- lm(freqDF$logFreq[13:120] ~ freqDF$logRank[13:120]))
  abline(regB, col = "#00AAFF")


}


# ==   2.1  Cliques  ===========================================================

# Let's find the largest cliques. Remember: a clique is a fully connected
# subgraph, i.e. a subgraph in which every node is connected to every other.
# Biological complexes often appear as cliques in interaction graphs.
if (FALSE) {

  igraph::clique_num(gSTR)
  # The largest clique has 83 members.

  (myClique <- igraph::largest_cliques(gSTR)[[1]])

  # Pick one of the proteins and find out what this fully connected cluster of
  # 83 proteins is (you can simply Google for any of the IDs). Is this expected?

  # Make a subgraph from this set of vertices
  mySubG <- igraph::induced_subgraph(gSTR, myClique)

  # color the vertices along a color spectrum
  vCol <- rainbow(igraph::gorder(mySubG)) # the "order" of a graph is the
                                          # number of nodes it contains

  # color the edges to have the same color as the originating node
  eCol <- character()
  for (i in seq_along(vCol)) {
    eCol <- c(eCol, rep(vCol[i], igraph::gorder(R)))
  }

  oPar <- par(mar= rep(0,4)) # Turn margins off
  plot(mySubG,
       layout = igraph::layout_in_circle(mySubG),
       vertex.size = 3,
       vertex.color = vCol,
       edge.color = eCol,
       edge.width = 0.1,
       vertex.label = NA)
  par(oPar)

  # ... well: remember: a clique means every node is connected to every other
  # node. We have 83 * 83 = 6,889 edges. This is what a matrix model of PPI
  # networks looks like for large complexes.

}


# ==   2.2  Communities  =======================================================

if (FALSE) {

  set.seed(112358)                 # set RNG seed for repeatable randomness
  gSTRclusters <- igraph::cluster_infomap(gSTR)
  set.seed(NULL)                   # reset the RNG

  igraph::modularity(gSTRclusters) # ... measures how separated the different
  # membership types are from each other
  tMem <- table(igraph::membership(gSTRclusters))
  length(tMem)  # About 900 communities identified
  hist(tMem, breaks = 50, col = "skyblue")  # most clusters are small ...
  range(tMem) # ... but one has > 280 members

}


# ==   2.3  Betweenness Centrality  ============================================

# Betweenness Centrality is a very meaningful value to compute for connected
# networks. It is perhaps not the most intutive measure, but once you can wrap
# your head around what it means, you realize that it identifies the crucial
# bottlenecks of information flow in a network.

# Let's find the nodes with the 10 - highest betweenness centralities.
#
BC <- igraph::centr_betw(gSTR)   # Compute the Betweenness Centralities for all
                                 # nodes
if (FALSE) {
  # BC$res contains the results
  head(BC$res)
}


# to get the ten-highest nodes, we simply label the elements of BC with their
# index ...
names(BC$res) <- as.character(1:length(BC$res))

# ... and then we sort:
sBC <- sort(BC$res, decreasing = TRUE)

if (FALSE) {
  head(sBC)
}


# This ordered vector means: node 1928 has the highest betweenness centrality,
# node 7475 has the second highest, etc.

BCsel <- as.numeric(names(sBC)[1:10])  # select the top 10

# We can use the labels of the top 10 to fetch the corresponding IDs from gSTR:
ENSPsel <- names(igraph::V(gSTR)[BCsel])

if (FALSE) {
  ENSPsel
}


#  Next, to find what these proteins are...

# We could now Google for all of these IDs to learn more about them. But really,
# googling for IDs one after the other, that would be lame. Let's instead use
# the very, very useful biomaRt package to translate these Ensemble IDs into
# gene symbols.


# =    3  biomaRt  =============================================================


# IDs are just labels, but for _bio_informatics we need to learn more about the
# biological function of the genes or proteins that graph data mining finds for
# us. biomaRt is the tool of choice. It's a package distributed by the
# bioconductor project. This here is not a biomaRt tutorial (that's for another
# day), simply a few lines of sample code to get you started on the specific use
# case of retrieving descriptions for ensembl protein IDs.

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
# Package information:
#  library(help = biomaRt)       # basic information
#  browseVignettes("biomaRt")    # available vignettes
#  data(package = "biomaRt")     # available datasets

# define which dataset to use ... this takes a while for download
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

if (FALSE) {
  # what filters are defined?
  ( filters <- biomaRt::listFilters(myMart) )


  # and what attributes can we filter for?
  ( attributes <- biomaRt::listAttributes(myMart) )


  # Soooo many options - let's look for the name of filters that are
  # useful for ENSP IDs ...
  filters[grep("ENSP", filters$description), ]

  # ... and the  attribute names for gene symbols and descriptions ...
  attributes[grep("symbol", attributes$description, ignore.case = TRUE), ]
  attributes[grep("description", attributes$description, ignore.case = TRUE), ]


  # ... we can put this together: here is an example:
  (myAnn <- biomaRt::getBM(filters = "ensembl_peptide_id",
                           attributes = c("hgnc_symbol",
                                          "wikigene_description",
                                          "interpro_description",
                                          "phenotype_description"),
                           values = "ENSP00000405330",
                           mart = myMart))

  # This is a very comprehensive set of annotations, but since many of the
  # filters retrieve multiple values, and each row of the returned data frame is
  # supposed to have a unique combination, there is also a lot (!) of
  # repetition. We'll "unique()" those values, but we can't do that in a data
  # frame since every column of a data frame must have the same number of
  # elements. But we can turn this into a list instead.

  (myAlist <- list(hgnc_symbol           = unique(myAnn$hgnc_symbol),
                   wikigene_description  = unique(myAnn$wikigene_description),
                   interpro_description  = unique(myAnn$interpro_description),
                   phenotype_description = unique(myAnn$phenotype_description)))

}

# A simple loop will now get us the information for our 10 most central genes
# from the human subset of STRING.

cenP <- list()  # Since we don't know how many matches one of our queries
                # will return, we'll write the result as a list of lists.

for (i in 1:length(ENSPsel)) {
  pBar(i, length(ENSPsel), nCh = length(ENSPsel)) # ... progress bar
  ID <- ENSPsel[i]
  tmp <- biomaRt::getBM(filters = "ensembl_peptide_id",
                        attributes = c("hgnc_symbol",
                                        "wikigene_description",
                                        "interpro_description",
                                        "phenotype_description"),
                        values = ID,
                        mart = myMart)
  cenP[[i]] <- list(Symbol    = unique(tmp$hgnc_symbol),
                    Wgene     = unique(tmp$wikigene_description),
                    Interpro  = unique(tmp$interpro_description),
                    Phenotype = unique(tmp$phenotype_description))
}

# ==   3.1  The most central proteins ...  =====================================

# So what do the proteins with the ten highest betweenness centralities do?
#  ... are you surprised? (I am! Really.)
if (FALSE) {
  for (i in 1:length(cenP)) {
    cat("Symbol:  ", cenP[[i]]$Symbol, "\n")
    cat("WGene:   ", paste(cenP[[i]]$Wgene,     collapse = ", "), "\n")
    cat("Interpro:", paste(cenP[[i]]$Interpro,  collapse = ", "), "\n")
    cat("Phenotype:",paste(cenP[[i]]$Phenotype, collapse = ", "), "\n")
    cat("\n")
  }
}


# Ponder over this output for a bit. Then think about whast would be next in
# your analysis.


# [END]
