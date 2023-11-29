# tocID <- "./sequenceBLAST.R"
#
#
# Purpose: Explore homologues of a sequence (PPP2R1A).
#
#
# Version: 1.0
# Date:    2023-11
# Author:  boris.steipe@utoronto.ca
#
# Versions:
#   1.0  First lecture version.
#
# ToDo:
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                            Line
#TOC> ----------------------------------------------------------------
#TOC>   1        Preparation: packages                              57
#TOC>   1.1        Preparation: Helper functions                    96
#TOC>   2        Introduction: BLAST                               101
#TOC>   3        Coverage of the Tree of Life                      141
#TOC>   4        Data                                              156
#TOC>   5        Conservation patterns: ChimeraX                   265
#TOC>   5.1        Coloring by Sequence Conservation: 1B3U         328
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
################################################################################



# NOTE: You need to read and explore this script, execute the functions but also
# try what happens if you change things.

# OPEN A NEW R SCRIPT into which you write your code experiments. Then you can
# easily reproduce what you did, and the experiments are automatically
# documented for your report.



# =    1  Preparation: packages  ===============================================

# Many things we do in computational biology require tools of the BioConductor
# project. These are loaded from their own installation framework: the
# "Bioconductor manager" ...
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Package information:
#  library(help = Biostrings)       # basic information
#  browseVignettes("Biostrings")    # available vignettes
#  data(package = "Biostrings")     # available datasets


# Once the BiocManager is installed, we can use it to install packages like
# we would install from CRAN. Biostrings is one of the foundational packages
# that has functions to handle the basic objects of biological sequences.
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

if (! requireNamespace("msa", quietly=TRUE)) {
  BiocManager::install("msa")
}
# Package information:
#  library(help=msa)       # basic information
#  browseVignettes("msa")  # available vignettes


# We need the httr:: package for remote control of ChimeraX
if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
# Package information:
#  library(help = httr)       # basic information
#  browseVignettes("httr")    # available vignettes



# ==   1.1  Preparation: Helper functions  =====================================

# ...


# =    2  Introduction: BLAST  =================================================

# In this script we explore the evolutionary descent  of the PP2A
# (Protein Phosphatase 2) subunit A.

# Ask ChatGPT to explain the principles evolutionary relationships:

if (FALSE) {

  t2c("I am studying PPP2R1A. Please explain the concept of \"homology\", and in particular the difference betwwen \"orthologues\" and \"paralogues\". If sequences are othologous (or paralogous), what difference does it make for their sequence conservation and conservation of function? And what are \"isoforms\"?")


  t2c("Following up on that: what do I need to know about \"neo-functionalization\" and \"sub-functionalization\"? Also, do you know whether either or both apply to human PPP2R1A?")


  t2c("What is the best tool that I can use to find protein sequences that are related to PPP2R1A to explore these issues? How should I use it? Can you walk me through an example?")


  t2c("Once I find sequences related to PPP2R1A, how would I estimate what portion of the \"Tree of Life\" shares the PPP2R1A protein? What would be the most likely sequence(s) present in the the \"Last Common Ancestor\"?")

  # Now perform an actual BLAST search. Make sure you are running the
  # BLASTp program (protein sequence) against the "landmark" database. Use
  # default parameters otherwise. You don't need to enter the actual sequence
  # into the search field, the RefSeq ID (NP_055040) suffices. Once you
  # have a result, define the RID in your script. Then you can ask ChatGPT-4
  # (or Bing) to help you interpret the complex but very comprehensive
  # and informative results.

  myRID <- "P8911YEM013"

  t2c(sprintf("I have run a BLAST search but find the results very complex. Here is the link to the result page: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=%s Can you be so kind and have a look and walk me through the highlights?",
              myRID))


  # Conclusion?

  }



# =    3  Coverage of the Tree of Life  ========================================


if (FALSE) {

  # https://en.wikipedia.org/wiki/Tree_of_life_(biology)
  # http://tolweb.org/tree/
  # https://timetree.org/
  #


  # Conclusion?
}


# =    4  Data  ================================================================


if (FALSE) {

  # A summary of data sources:
  #
  # Symbol:
  #   PPP2R1A
  #     (Authority: HGNC - HGNC9302
  #     https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:9302 )
  # Name:
  #   protein phosphatase 2 scaffold subunit Aalpha
  # Protein:
  #   UniProt:  2AAA_HUMAN
  #             P30153      https://www.uniprot.org/uniprotkb/P30153/entry
  #   NCBI:     NP_055040   https://www.ncbi.nlm.nih.gov/protein/NP_055040
  # Gene:
  #   NCBI:     PPP2R1A     https://www.ncbi.nlm.nih.gov/gene/5518
  #
  # Other:
  #   Genecards: https://www.genecards.org/cgi-bin/carddisp.pl?gene=PPP2R1A
  #
  # Phylogenetic:
  #   OrthoDB:   https://www.orthodb.org/?level=&species=&query=759105at2759
  #   OMA:       https://omabrowser.org/oma/vps/HUMAN45628/
  #

  # Sadly, the BLAST interface is nearly useless, returning conflicting
  # and incomplete results, a plethora of poorly documented choices, and
  # a deluge of completely uninformative variants, based on a misunderstood
  # concept of similarity.

  # To get a reasonable selection of sequences, we need to download
  # sequences basically by hand, using pairwaise BLAST from organism that
  # we select from the OrthoDB species tree:

  # Even so we run into baffling issues: the taxonomy database
  # constructs a "Common Tree", but it does not give us the TaxIDs with the
  # tree.

  # Eukaryota
  # ++Apicomplexa
  # |   \-Theileria annulata (taxid:5874)           THEAN  XP_953138.1
  # ++Apusomonadida
  # |   \-Thecamonas trahens (taxid:529818)         THETR  XP_013756917.1
  # ++Euglenozoa
  # |   \-Leishmania donovani (taxid:5661)          LEIDO  XP_003860412
  # ++Evosea
  # |   \-Dictyostelium discoideum (taxid:44689)    DYCDI  XP_638954.1
  # ++Fornicata
  # |   \- Giardia lamblia (taxid:5741)             GIALA  XP_001704098.1
  # ++Fungi
  # |   \-Saccharomyces cerevisiae (taxid:4932)     SACCE  NP_009386.1
  # ++Heterolobosea
  # |   \-Naegleria gruberi (taxid:5762)            NAEGR  XP_002682119.1
  # +-Metazoa
  # |   +-Chordata
  # |   |   ++Mammalia
  # |   |   |   \-Homo sapiens (taxid:9606)         HOMSA  NP_055040
  # |   |   \+Actinopteri
  # |   |       \-Danio rerio (taxid:7955) 	        DANRE  NP_001299846
  # |   ++Arthropoda
  # |   |   \-Drosophila melanogaster (taxid:7227)  DROME  NP_001027225.1
  # |   ++Mollusca
  # |   |   \-Aplysia californica (taxid:6500)      APLCA  XP_012936226
  # |   ++Nematoda
  # |   |   \-Caenorhabditis elegans (taxid:6239)   CAEEL  NP_498162.2
  # |   \+Cnidaria
  # |       \-Actinia tenebrosa (taxid:6105)        ACTTE  XP_031562500.1
  # ++Oomycota
  # |   \-Saprolegnia diclina (taxid:112098)        SAPDI  XP_008614287
  # ++Rhodophyta
  # |   \-Cyanidioschyzon merolae (taxid:45157)     CYAME  XP_005535422.1
  # \+Viridiplantae
  #     \-Arabidopsis thaliana (taxid:3702)         ARATH  NP_001326959.1


  # N.b. Apicomplexans (Eimeria tenella and Cystoisospora suis) stubbornly did
  # not yield a hit until I did a separate BLAST search against all
  # apicomplexans, which yielded Theileria annulata as a suitable source. Such a
  # finding raises the question: was the gene lost in some apicomplexans, was it
  # gained in T. annulata, or did it just diverge too far to be found with
  # standard BLAST search parameters?

  # I have computed a COBALT alignment of the 65 hits, manually selected
  # only the best hit for each of the 16b reference species, recomputed
  # the alignment, and saved the result ...
  #  - as a multi-FASTA file (with gaps):
  #    ./data/PPP2R1A-16_species_COBALT.mfa
  #  - in CLUSTAL format (aligned blocks):
  #    ./data/PPP2R1A-16_species_COBALT.aln

  # I  then edited the results to replace the refSeq IDs with the five character
  # species codes that I define based on the binomial name of the species. This
  # is essential because you must be able to immediately see the species in
  # which the sequence evolved. Any kind of alignment, tree, superposition etc.
  # that obscures this crucial information about the source of the data is
  # _useless_.

  mfaFile <- "./data/PPP2R1A-16_species_COBALT.mfa"
  PPP2R1Aset <- Biostrings::readAAStringSet(mfaFile)



  # Conclusion?
}


# =    5  Conservation patterns: ChimeraX  =====================================


if (FALSE) {
  t2c("I have a COBALT computed multiple sequence alignment for human PPP2R1A and fifteen other eukaryotic species for comparison. I have it both as a multi-FASTA format, as well as in CLUSTAL format. Now I want to map the degree of conservation onto the human structure in 1B3U using ChimeraX. Can you find me a tutorial that covers the process?")

  #  ChatGPT-4 2023-11-29)
  #
  #  I found a comprehensive tutorial on the UCSF ChimeraX website that guides
  #  you through the process of coloring a structure by sequence conservation
  #  using a multiple sequence alignment. This tutorial is specifically tailored
  #  for ChimeraX and covers various aspects of the process:

  #   1. **Overview of the Tutorial**: The tutorial is structured in three
  #   parts, each focusing on a different aspect of coloring by sequence
  #   conservation
  #   #   ChimeraX Tutorial: Coloring by Sequence Conservation]
  #   https://www.rbvi.ucsf.edu/chimerax/data/conservation-coloring/conservation-coloring.html


  #  2. **Starting with ChimeraX**: It begins with instructions on how to start
  #  ChimeraX and access the tutorial within the ChimeraX browser【12†source】.

  #  3. **Fetching Structures and Sequence Alignments**: The tutorial provides
  #  detailed steps on how to fetch structures from the RCSB Protein DataBank
  #  and sequence alignments from various sources, including directly from the
  #  tutorial's webpage.

  #  4. **Sequence Viewer and Conservation Histogram**: It explains how to use
  #  the Sequence Viewer, automatic association of sequences with structures,
  #  and how conservation values are assigned and visualized in a
  #  histogram.

  #  5. **Coloring by Conservation**: Detailed commands are provided for
  #  coloring the structure based on conservation values, adjusting the coloring
  #  palette, and emphasizing certain conservation ranges.

  #  6. **Handling Gaps and Non-associated Residues**: The tutorial also covers
  #  how to deal with gaps in the sequence alignment and residues that are not
  #  associated with a sequence alignment.

  #  7. **Adjusting Conservation Parameters**: Instructions are included for
  #  changing conservation settings, such as the allowed gap fraction, to cover
  #  more of the structure with conservation values.

  #  8. **Adding a Color Key**: The tutorial also explains how to add a color
  #  key to the structure visualization for better interpretation of the
  #  conservation data.

  #  9. **Other Visualization Options**: Besides coloring by sequence
  #  conservation, the tutorial discusses other visualization options like
  #  molecular surfaces and presets for publication-quality images.

  # This tutorial seems to cover all the necessary steps and provides specific
  # examples and commands to map the degree of conservation onto a structure,
  # making it highly suitable for your needs with human PPP2R1A and other
  # eukaryotic species for comparison. You can access the tutorial directly at
  # this link:
  #   ChimeraX Tutorial: Coloring by Sequence Conservation]
  #   https://www.rbvi.ucsf.edu/chimerax/data/conservation-coloring/conservation-coloring.html

}

# ==   5.1  Coloring by Sequence Conservation: 1B3U  ===========================

if (FALSE) {

  # Ok ... let's follow what we learn from the tutorial.

  # Obviously, the first step is to open the 1B3U structure in ChimeraX:
  #
  # (1) Open a fresh, new session in a recently updated version of ChimeraX
  # (2) Execute this command to place the string into the clipboard:
  t2c("remotecontrol rest start port 61803")
  # (3) PASTE this into the ChimeraX console and hit return. Then you can
  #     send commands to ChimeraX directly from your RStudio session.

  CX("open 1B3U")
  CX("select /B")          # First selecting, then deleting, allows you a sanity
  CX("delete atoms sel")   # check: is this really what you wanted to delete?
  CX("cofr /A")            # cofr: center-of-rotation
  CX("camera sbs")
  CX("lighting soft")
  CX("color sequential #1 & protein target abc palette powderblue:orchid:white")
  CX("label defaultheight 2")
  CX("label /A:1, 588")

  # Now we load the alignment file (CLUSTAL format)
  alnFile <- file.path(getwd(), "./data/PPP2R1A-16_species_COBALT.aln")
  CX(sprintf("open \"%s\"", alnFile))

  # The sequence viewer opens with the multiple alignment, and the first
  # sequence in the alignment automatically associates with the structure.
  # It gets a colored box around it.

  CX("color byattr seq_conservation")  # Default: spread conservation values
                                       # between blue - white - red. Not a good
                                       # choice of colours. Why?

  # Let's define our own colors, based on black-body radiation. "Hot" residues
  # should be glowing orange and yellow, "cool residues should be dark.

  # I start from an image of an erupting volcano at night to give me a basic
  # smapling of how black-body colors appear on the computer screen. Then I use
  # a site that is currently my go-to for constructing color spectra:
  # https://coolors.co to set up a spectrum that matches my ideas: slowly
  # increasing across the low end and a rapid jump for emphasis at the high end.
  # I use the eye-dropper tool from the picker to get the basic spectrum
  # covered, then interpolate between colors to adjust the curve.

  # Finally, there is a tool to export the spectrum to code.
  #
  vBB <- c("#160405",
           "#2f0604",
           "#470803",
           "#5f0a02",
           "#770c00",
           "#982300",
           "#b93a00",
           "#fa6700",
           "#fddd45",
           "#fdfeff")

  # The function colorRampPalette() creates a function that returns a
  # required number of colors along a defined palette. (Bias defines
  # whether the spectrum should be expanded on the left (bias < 1), kept as
  # is (default), or expanded on the right (bias > 1).)
  BBcolors <- colorRampPalette(vBB, bias = 0.67)

  N <- 30
  barplot(rep(1, N), col = BBcolors(N))

  myPal <- paste(BBcolors(20), collapse = ":")
  CX(sprintf("color byattr seq_conservation palette %s novalue #666666", myPal))

  # But what does this mean?
  # Where are the regions of sequence variation relative to the catalytic
  # subunit and the regulatory subunit? Let's superimpose the structure of the
  # entire complex - the "holoenzyme"...

  CX("open 2IAE")
  CX("select #2/D,E,F,M,N")
  CX("delete atoms sel")
  CX("select #2")
  CX("hide sel target ab")
  CX("show sel target c")
  CX("color sequential #2/C target abc palette #202C3B:#2E4057:#707C8C")
  CX("color sequential #2/B target abc palette #406450:#66A182:#A0C5B1")
  CX("color  #2/A #9999BB33")

  CX("mm #2/A to #1/A")
  CX("cofr #1")







  # Conclusion?
}





# [END]
