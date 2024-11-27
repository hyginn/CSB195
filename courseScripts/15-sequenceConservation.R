# tocID <- "./courseScripts/15-sequenceConservation.R"
#
#
# Purpose:  R code to analyze sequence conservation in the p53 protein via
#           BLAST searches, multiple sequence alignment, mapping conservation
#           to structure, cross-referencing with genome coordinates, and
#           comparing with intogen data.
#
#
# Version:   1.0
#
# Date:     2024-11
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
#TOC>   Section  Title                                      Line
#TOC> ----------------------------------------------------------
#TOC>   1        MOTIVATION                                   57
#TOC>   2        SEQUENCE SOURCE                             104
#TOC>   3        BLAST                                       143
#TOC>   4        MULTIPLE SEQUENCE ALIGNMENT                 222
#TOC>   5        3D SEQUENCE CONSERVATION: 3KMD              259
#TOC>   5.1        Visualize the protein/DNA complex         274
#TOC>   5.2        Read the alignment                        296
#TOC>   5.3        Color by conservation score               318
#TOC>   5.4        Cancer mutation hotspots                  396
#TOC>   6        NEXT ?                                      497
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


# =    1  MOTIVATION  ==========================================================

# We have established that the proteins with the highest betweenness centrality
# in the human functional interaction network are significantly more likely to
# be associated with a cancer phenotype than we would expect by chance: they are
# "enriched" for an association with cancer. Among those, the protein with the
# highest betweenness centrality of all is p53, the protein encoded by the TP53
# gene. Indeed, p53 responds to DNA damage by blocking cell-division to give the
# cell time for DNA repair, or by inducing cellular suicide - "apoptosis" - if
# the damage is too extensive and repair is not possible.

# p53 performs its function through several distinct domains of its structure:
# there is a tetramerization domain, a DNA binding domain, a domain that
# interacts with other transcription factors... If we thus ask "what parts are
# important?", this question takes us on an exploration of sequence -function
# relationships.

# But how to determine "importance" of a sequence-feature in a protein in the
# first place? This somewhat diffuse question can be phrased precisely by
# recognizing that "important" features ought to be conserved throughout
# evolution, since those parts of the protein sequence are under selective
# pressure not to be lost.

# The recipe for assessing importance thus could look like this:
#
#   1.  Define a well-balanced set of species on the evolutionary tree
#       that covers relatives throughout the evolutionary history of
#       p53.
#
#   2.  For each of those species, find their own p53 sequence.
#
#   3.  Perform a multiple sequence alignment, a computational procedure
#       that arranges groups of sequence in a table so that each column
#       contains those residues that are equivalent to each other, i.e.
#       they have the same position as their shared ancestral sequence.
#
#   4.  For each of these columns, compute a conservation score: a measure
#       how much variation has been observed  at this position over the
#       course of evolution.
#
#   5.  Map these conservation scores back to the 3D-structure of our
#       protein of interest - the human p53 protein - to identify where
#       the conserved important residues are located.
#
#   6.  Interpret this in terms of known (or unknown) functions


# =    2  SEQUENCE SOURCE  =====================================================

# To obtain a well-balanced set of genes for assessing residue-level
# conservation scores, we first define a set of species that are distributed
# across the evolutionary tree. We find sample species at ensembl ...

# https://useast.ensembl.org/info/about/speciestree.html

# ... and we can pick about 20 - 25 species to work with. I have collected basic
# species information - names, taxonomy IDs etc. into a CSV file: read it in
# from here:
if (FALSE) {
  (p53species <- read.csv("./data/p53species.csv"))
}

# TimeTree (https://timetree.org/) is a well-designed site that evaluates the
# evolutionary distance of species over time. For example you can find out when
# the last common ancestor between you yourself (Homo sapiens) and the
# indomitable honey badger (Mellivora capensis) lived. 94 million years ago! Who
# knew. But here we use the tool to evaluate a species tree. That requires to
# upload a file with species names.

if (FALSE) {

  writeLines(p53species$species, "./myScripts/p53species.txt")

  # Upload the file to TimeTree and compute a tree of evolutionary distances -
  # the time that has passed since two species shared a common ancestor. The
  # resulting tree shows that the most closely related species in our list are
  # Felis catus and Vulpes vulpes - and their last common ancestor lived about
  # 55 million years ago.

  # In case you are wondering what these Latin names are ... here you go:

  cat(sprintf("%s - %s\n", p53species$species, p53species$name))

}


# =    3  BLAST  ===============================================================

# BLAST - the Basic Local Alignment Search Tool - is the single most frequently
# used tool for bioinformatics and computational biology, by a margin. BLAST
# searches for similar sequences in very large sequence databases. Both the NCBI
# as well as the EBI offer BLAST services that are integrated with their
# database holdings: BLAST searches can be accessed directly from the page for a
# protein entry.

#   https://www.uniprot.org/uniprotkb/P04637/entry   UniProt page
#   https://www.ncbi.nlm.nih.gov/protein/NP_000537   NCBI RefSeq page

# Running a BLAST search is trivially easy, but the challenge is to perform a
# BLAST search for p53, restricted to our reference species. If we don't do
# that, there is simply no way to pick out the information we are looking for
# from the tens of thousands of hits that are returned. Unfortunately the NCBI
# interface for reference species is brittle and does not allow upload of lists
# of restrictions. By contrast, this task works out of the box at the EBI: just
# upload the following comma-separated list of taxonomy IDs:

if (FALSE) {

  cat(paste(p53species$taxid, collapse = ","))

}

# The search sequence is the 393 amino acid long sequence of isoform A of
# the human p53 protein. You can copy it from here:

if (FALSE) {
P53_HOMSA <- "
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI
EQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ
KTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDST
PPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGN
LRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRP
ILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALEL
KDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
"

# Turn it into a vector of single letters...
P53_HOMSA <- unlist(strsplit(gsub("[^A-Z]", "", P53_HOMSA), ""))

}

# The BLAST search returns OTO 600 hits. Most of these are various isoforms and
# fragments of the same thing over and over again. Finding the information you
# need in such large datasets is often not easy. You can download the data and
# process it with a script - that's the preferred way if you need to repeat this
# kind of analysis. Or - if it's a one-off thing - you can just look at the
# results by hand and choose the first (best) hit that appears for every
# species. That's what I did.

# Retrieval of the actual sequences is straightforward via the EBI "DB fetch"
# interface (https://www.ebi.ac.uk/Tools/dbfetch/). Just choose the database
# (UniProtKB) upload the database identifiers, white-space separated ...

if (FALSE) {

  cat(paste(p53species$uniprotID, collapse = " "))

}

# ... choose "fasta" as the output format and "raw" as the style, and retrieve.

# The results still needs a little bit of TLC, I always edit the descriptions
# lines in multi FASTA files to begin with a protein ID and a five-letter
# species code. This is what will be displayed by tree-creation programs or in
# multiple sequence alignments, and it is crucially important to know what you
# are looking at during your analysis. And we put the human sequence as the
# first entry - ChimeraX will need it there, later.

# Anyway - the resulting file is in "./data/p53sequences.mfa" (.mfa is the
# customary extension for multi-fasta formatted files). Have a look.

# We're ready to compute an alignment.


# =    4  MULTIPLE SEQUENCE ALIGNMENT  =========================================

# Several algorithms for MSAs are available on the EBI website
# (https://www.ebi.ac.uk/jdispatcher/msa), that immediately tells you that this
# is not a "solved problem". It is not, because the solutions are ambiguous, and
# which variant (and algorithm) to choose depends on the precise question we are
# trying to answer with our MSA. Are we trying to reproduce the evolutionary
# steps? Are we trying to maximise similarity overall? Perhaps we need to align
# to conserved patterns. Or we need to ensure that structurally-equivalent
# positions are aligned. Each of these options leads to a slightly different
# preferred outcome.

# Here we choose ClustalOmega (https://www.ebi.ac.uk/jdispatcher/msa/clustalo).

#  1.  Select the button for "Protein" input.
#  2.  Copy the contents of ./data/p53sequences.mfa and paste it into
#      the text box.
if (FALSE) {
  # copy file contents to clipboard
  t2c(paste(readLines("./data/p53sequences.mfa"), collapse = "\n"))
}
#  3.  Choose an OUTPUT FORMAT: ClustalW
#  4.  Open the "More Options" menu and choose ORDER: input.
#      This keeps the order of sequences as they appear in the input file,
#      and does not rearrange them to put more closely related sequences
#      next to each other. We need this option because the human p53
#      sequence has to appear as the first sequence, when we use the alignment
#      later to be read by ChimeraX.
#  5.  Submit the alignment. The job will be queued on the EBI's compute
#      cluster, but usually takes less than a minute to complete.
#  6.  When it is done, click on "View Results". Check that the alignment
#      appears in blocks, each line is prepended by the first word of the
#      input sequences - i.e. our species codes, and that P53_HOMSA is the
#      first line. Then click on Download. (Or don't - I have placed a copy
#      in the ./data directory: "p53alignments.aln". )


# =    5  3D SEQUENCE CONSERVATION: 3KMD  ======================================

if (FALSE) {

  # The first step is to open the 3KMD structure in ChimeraX:

  # (1) Open a fresh, new session in a recently updated version of ChimeraX
  # (2) Execute this command to place the string into the clipboard:
  t2c("remotecontrol rest start port 61803")
  # (3) PASTE this into the ChimeraX console and hit return. Then you can
  #     send commands to ChimeraX directly from your RStudio session using
  #     the translations defined through the CX() function in your .util.R
  #     file.

}
# ==   5.1  Visualize the protein/DNA complex  =================================


if (FALSE) {

  CX("open 3KMD")                   # download a file from PDB and open it
  CX("turn Y 180")                  # turn it around to see the DNA
  CX("lighting soft")
  CX("hide atoms")
  CX("show cartoons")               # view the general fold
  CX("select /E,F")                 # chains E and F contain the DNA strands
  CX("show sel atoms")
  CX("hide /E,F:HOH atoms")
  CX("nucleotides sel tube/slab")  # schematic view of the DNA ...
  CX("surface sel")
  CX("transparency 60")            # ... with a transparent surface
  CX("select clear")
  CX("color sequential #1 & protein target abc palette powderblue:orchid:white")

}


# ==   5.2  Read the alignment  ================================================

if (FALSE) {

  # ChimeraX can directly read and interpret .aln files (ClustalW format). The
  # sequence corresponding to the structure model that the alignment is applied
  # to is expected at the top.
  alnFile <- file.path(getwd(), "./data/p53alignments.aln")  # Define the path
  CX(sprintf("open \"%s\"", alnFile))                        # Open in ChimeraX

  # The log should look like this:
  # Alignment identifier is p53alignments.aln
  # Associated 3kmd chain A to P53_HOMSA with 0 mismatches
  # Associated 3kmd chain B to P53_HOMSA with 0 mismatches
  # Associated 3kmd chain D to P53_HOMSA with 0 mismatches
  # Associated 3kmd chain C to P53_HOMSA with 0 mismatches

  # The sequence viewer opens with the multiple alignment, and the first
  # sequence in the alignment automatically associates with the structure.
  # It gets a colored box around it.
}

# ==   5.3  Color by conservation score  =======================================

if (FALSE) {

  CX("color byattr seq_conservation")
  # ChimeraX's default is to spread conservation values between blue - white -
  # red. This is simple to write, but not a good choice of colours. Why? The
  # spectrum is "divergent", and we are looking at a "sequential" property. (Cf.
  # https://blog.datawrapper.de/diverging-vs-sequential-color-scales/ for more
  # discussion.)

  # Let's define our own colors, based on black-body radiation. "Hot" residues
  # should be glowing orange and yellow, "cool residues should be dark.

  # I start from an image of an erupting volcano at night to give me a basic
  # sampling of how black-body colors appear on the computer screen. Then I use
  # a site that is currently my go-to for constructing color spectra:
  # https://coolors.co to set up a spectrum that matches my ideas: slowly
  # increasing across the low end and a rapid jump for emphasis at the high end.
  # I use the eye-dropper tool from the picker to get the basic spectrum
  # covered, then interpolate between colors to adjust the curve.

  # Once I am satisfied with the colors, there is a tool to export the
  # spectrum to code.
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
  BBcolors <- colorRampPalette(vBB, bias = 1.2)

  N <- 30
  barplot(rep(1, N), col = BBcolors(N))

  # How is this not a divergent spectrum? Well - in a sense it is. But
  # psychologically we are more likely to associate such a spectrum with values
  # of a common property than with qualitative differences.

  myPal <- paste(BBcolors(20), collapse = ":")
  CX(sprintf("color byattr seq_conservation palette %s novalue %s",
             myPal,
             "#9999AA66"))  # default: items without conservation score

  # Cartoons are not well suited for interpreting conservation scores, since
  # conservation is side-chain property. A better view is some variant of
  # atom-depiction since it adds information about residue size and contacts,
  # which is more relevant than just looking at the approximate location of the
  # variable spots.

  CX("hide cartoons")
  CX("select /A,B,C,D")
  CX("select subtract :HOH")
  CX("show sel atoms")
  CX("style sel ball")
  CX("size sel ballScale 0.4 stickRadius 0.5")

  # But what does this mean?

  # Variation is predominantly found at the DNA-binding interface (functional
  # variation), whereas the core of the protein is highly conserved (structural
  # conservation). However, where interfaces are conserved, there is some
  # characteristic variation in the second layer, behind the interface. This
  # is a common feature, it causes subtle repositioning of interface residues
  # without outright disruption of important interactions.

}

# ==   5.4  Cancer mutation hotspots  ==========================================

# We had visited the intogen database (https://www.intogen.org/search?gene=TP53)
# and established that there are three hotspots that together account for 16% of
# all mutations in TP53, in cancer genomes (1,726 of 10,550). Where are those
# residues? The "lolliplot" - um .. needle plot - of cancer mutations identifies
# them as positions 175, 248 and 273, but we can NEVER be sure that sequence
# numbers in PDB files are identical to sequence numbers for any transcript,
# protein, or isoform. In fact, more often than not there is significant
# mismatch, and this is due to a large number of experimental reasons and
# community conventions, which are not trivial to identify. There ought to be a
# standard for sequence numbering but AFAIK there is no such thing.

if (FALSE) {
  # Let's check with the N- and C- terminal residues of 3KMD:

  CX("select /D")
  CX("hide sel atoms")
  CX("show sel cartoons")
  CX("color sequential sel target c palette powderblue:orchid:white")
  CX("select /D:92-96,287-291 ")
  CX("hide sel cartoons")
  CX("show sel atoms")
  CX("label sel")
  CX("cofr /D")

  paste(P53_HOMSA[ 92:96 ], collapse = "")
  paste(P53_HOMSA[287:291], collapse = "")

  # This identifies the N-terminal five residues 92-96 to be PLSSS and the
  # C-terminal five residues 287-291 to be ENLRK. Does this match the sequence?
  paste(P53_HOMSA[ 92:96 ], collapse = "")
  paste(P53_HOMSA[287:291], collapse = "")
  # In this case we are lucky, and the mapping is correct.


  # Unfortunately, this information is hard to find at intogen, because the
  # database only publishes genome coordinates, and does not translate the old
  # and new codons:
  # 17:7675088:C
  # 17:7674220:C
  # 17:7673803:G

  # So we navigate to the NCBI gene information page
  # (https://www.ncbi.nlm.nih.gov/gene/7157), click on the link for the current
  # human genome assembly (GRCh38.p14), hoping that intogen used the same
  # assembly for its coordinates, we click on the link to chromosome 17 and
  # follow the link to the graphics display
  # (https://www.ncbi.nlm.nih.gov/nuccore/NC_000017.11?report=graph). We enter
  # the coordinate 17:7675088 into the search box, reconfigure the tracks to
  # show us a six-frame translation (the only option, there is no option to show
  # the annotated protein sequence, and even that is OFF by default), and we
  # make a guess as to which of the six frames we might be looking at.

  # (https://www.ncbi.nlm.nih.gov/nuccore/568815581?report=graph&tkey=ZU8HYmptYWhvZmtkQFF0aFIaCAcvFBQsOAcYMhGGMYMmojf_IdrNsv_pCyMRF2InEUM&assm_context=GCF_000001405.40&mk=7675088|7675088|blue|9&v=7675048:7675127&c=FFCC00&select=null&slim=0)

  # Perhaps the correct track is the one shown in bold? Maybe. (Spoiler:
  # sometimes, but not always). This would identify 17:7675088 to be the first
  # position in a CGC arginine codon. What sequence do we expect there?

  idx <- 175
  paste(paste(P53_HOMSA[(idx-5):(idx-1)], collapse = ""),
        P53_HOMSA[idx],
        paste(P53_HOMSA[(idx+1):(idx+5)], collapse = ""))

  # The first nucleotide of codon 175 is mutated most frequently to T, i.e. TGC,
  # i.e. cysteine.

  # Incidentally, the one-letter annotation is printed at the middle nucleotide,
  # and the sequence runs from right to left, on the bottom strand. That takes
  # some effort to figure out! But once we understand this, we can check all
  # three hot-spots:

  # 17:7675088:C -> T  : CGC -> TGC :  EVVR-R175C-CPHH
  # 17:7674220:C -> T  : CGG -> TGG :  GGMN-R248W-RPIL
  # 17:7673803:G -> A  : CGT -> CAT :  SFEV-R273H-VCAC

  # In this case we are lucky and the annotated sequence numbers appear to match
  # our PDB structure. But, seriously, that is NOT always the case.

  # Here we go, identifying those three residues in the structure:

  CX("select /D:92-96,287-291 ")
  CX("~label sel")
  CX("hide sel atoms")
  CX("show sel cartoons")
  CX("select :175, 248,273")
  CX("show sel atoms")
  CX("style sel sphere")
  CX("label sel")
  CX("cofr sel")

}

# Bottom line? It looks like cancer phenotypes for p53 arise from mutations in
# functionally important residues, not necessarily structurally important
# residues. Rather than outright destroying the protein, they instead corrupt
# detailed aspects of p53 function in subtle ways. But it takes a combined view
# of sequence, structure, and evolution to be able to tell.


# =    6  NEXT ?  ==============================================================

# Finally - what about the residues that are not visible? Can we model them?
# According to intogen, a large number of residues 1:91 and 292:392 are also
# involved in cancer phenotypes. They are found in the TAD domains
# (Transactivating domains - domains that recruit other transcription factors to
# a specific location in the DNA), and in the tetramerization domain, but except
# for one position (337), those are virtually all truncating mutations that
# either destroy the protein (upstream of the structural domain), or eliminate
# the tetramerization domains that is crucial for the proper assembly of the p53
# complex.

# (But that's for another day :-)




# [END]
