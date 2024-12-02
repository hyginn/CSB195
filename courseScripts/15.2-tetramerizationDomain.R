# tocID <- "./courseScripts/15.2-tetramerizationDomain.R"
#
#
# Purpose:  R code to display a model of the p53 tetramer / DNA complex
#           produced by AlphaFold 3, and superimpose it on the 3KMD
#           experimental structure.
#
# Notes:    AlphaFold3 server: https://alphafoldserver.com/
#           Model: https://alphafoldserver.com/fold/4255d107f0d38953
#           Paper: https://www.nature.com/articles/s41586-024-07487-w
#
#
# Version:   1.0
#
# Date:     2024-12
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    New script
#
# TODO:
#
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                 Line
#TOC> -----------------------------------------------------
#TOC>   1        PREPARE                                 35
#TOC>   2        LOAD THE REFERENCE STRUCTURE            46
#TOC>   3        LOAD THE MODEL                          65
#TOC> 
#TOC> ==========================================================================


# =    1  PREPARE  =============================================================

# (1) Open a fresh, new session in a recently updated version of ChimeraX
# (2) Execute this command to place the string into the clipboard:
t2c("remotecontrol rest start port 61803")
# (3) PASTE this into the ChimeraX console and hit return. Then you can
#     send commands to ChimeraX directly from your RStudio session using
#     the translations defined through the CX() function in your .util.R
#     file.


# =    2  LOAD THE REFERENCE STRUCTURE  ========================================

# The reference structure

CX("open 3KMD")                   # download a file from PDB and open it
CX("lighting soft")
CX("hide #1 atoms")
CX("show #1 cartoons")

CX("color sequential #1 target abc palette #225522:#CCFFDD")

CX("select #1::name=D*")
CX("show sel atoms")
CX("nucleotides sel tube/slab")  # schematic view of the DNA ...
CX("color sel #55AA99")

CX("hide #1 atoms,bonds,cartoons")


# =    3  LOAD THE MODEL  ======================================================

a3Model <- file.path(getwd(), "./data/alphafold3_p53_tertramer_dna_complex.cif")
CX(sprintf("open \"%s\"", a3Model))
CX("cofr #2")

alnFile <- file.path(getwd(), "./data/p53alignments.aln")  # Define the path
CX(sprintf("open \"%s\"", alnFile))                        # Open in ChimeraX

CX("show #2 atoms")
CX("style #2 ball")
CX("size #2 ballScale 0.4 stickRadius 0.5")
CX("select #2/I,J")
CX("show sel atoms,cartoons")
CX("nucleotides sel tube/slab")  # schematic view of the DNA ...
CX("~select #2")


vBB <- c("#160405", "#2f0604", "#470803", "#5f0a02", "#770c00",
         "#982300", "#b93a00", "#fa6700", "#fddd45", "#fdfeff")
BBcolors <- colorRampPalette(vBB, bias = 1.2)
myPal <- paste(BBcolors(20), collapse = ":")
CX(sprintf("color byattr seq_conservation #2 palette %s novalue %s",
           myPal,
           "#9999AA66"))  # default: items without conservation score

CX("select #2/A,B,C,D")
CX("hide sel atoms")
CX("show sel cartoons")
CX("select #2/A,B,C,D:1-84,358-393")
CX("hide sel cartoons")

CX("show #1 cartoons")
CX("select #1::name=D*")
CX("show sel atoms")
CX("nucleotides sel tube/slab")  # schematic view of the DNA ...

# 1 = SUPERIMPOSE THE MODEL ONTO THE REFERENCE

# Note: chains are differently ordered in the two structures. Ordering is
# arbitrary, but the matching domains need to be indicated to the
# alignment tool.
CX("match #2/A,B,C,D,I,J to #1/B,C,D,A,E,F pair ss")
CX("cofr #2/I,J")






# [END]
