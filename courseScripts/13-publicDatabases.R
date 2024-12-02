# tocID <- "/courseScripts/13-publicDatabases.R"
#
#
# Purpose:  Public databases for the annotation of gene features and functions.
#
# Version:   0.1
#
# Date:     2024-11
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           0.1    First code 2024
#
# TODO: Review, annotate, retrieve and analyze results
#
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                       Line
#TOC> ---------------------------------------------------------------------------
#TOC>   1        Our Example: TP53                                             35
#TOC>   2        PubMed - The Literature database                              41
#TOC>   3        GO - The Gene Ontology                                        48
#TOC>   4        KEGG - The Kyoto Encyclopedia of Genes and Genomes            56
#TOC>   5        Reactome - The Pathway database                               62
#TOC>   6        The UCSC Genome Browser                                       68
#TOC>   7        Ensembl - Annotated Genomes at the EBI                        74
#TOC>
#TOC> ==========================================================================


# =    1  Our Example: TP53  ===================================================

# = 1  HGNC
#
# HGNC
# https://www.genenames.org/
#
# https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:11998


# =    2  PubMed - The Literature database  ====================================

# https://pubmed.ncbi.nlm.nih.gov/
#
# Example:
# https://pubmed.ncbi.nlm.nih.gov/?term=(TP53[title]+OR+p53[title])+AND+review[PT]
#

# A bookmarklet for faster access to licensed sources via the UofT library.
# javascript:(function(){var url=window.location.href;var re=/\/([\w.]+)\/(.*$)/;var match=url.match(re);var newURL="http://"+match[1]+".myaccess.library.utoronto.ca/"+match[2];window.location.href=newURL;})();void 0

# =    3  GO - The Gene Ontology  ==============================================

# https://geneontology.org/
#
# Note: GO and GOA
#
# GO contains three independent Ontologies (memorize their names):
#      F: molecular function
#      C: cellular component
#      P: biological process
#
#
# Search in GOA:
#
# Example:
# https://amigo.geneontology.org/amigo/search/annotation?q=%22TP53%22
# - Put TP53 in quotation marks for exact matches.
# - Filter for Organism: Homo sapiens
# - Filkter for Ontology (aspect): P for biological process


# =    4  KEGG - The Kyoto Encyclopedia of Genes and Genomes  ==================

# https://www.genome.jp/kegg/pathway.html
#
# Example:
# Enter prefix: "hsa", then serch for "TP53" ...
# https://www.kegg.jp/pathway/map=hsa04110&keyword=tp53




# =    5  Reactome - The Pathway database  =====================================

# https://reactome.org/
#
# Example:
# https://reactome.org/content/detail/R-HSA-69488
# https://reactome.org/PathwayBrowser/#/R-HSA-69620&PATH=R-HSA-1640170&FLG=UniProt:P04637




# =    6  The UCSC Genome Browser  =============================================

# https://genome.ucsc.edu/
#
# Example:
# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A7668421%2D7687489&hgsid=2382624373_uWNIHeKw7msVY5ZuSH1e3vAJfJA7

# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A7668421%2D7687489



# =    7  Ensembl - Annotated Genomes at the EBI  ==============================

# https://useast.ensembl.org/index.html
#
# Note: BioMart
#
# Example:
# https://useast.ensembl.org/Homo_sapiens/Gene/Phenotype?db=core;g=ENSG00000141510;r=17:7661779-7687546



# [END]
