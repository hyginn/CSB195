# tocID <- "courseScripts/Structure_analysis-Scripting_ChimeraX.R"
#
# Purpose:  Computational Biology Foundations:
#              R code demonstrating amino acid properties with remote
#              scripted ChimeraX visualization
#
# Version:  1.2
#
# Date:     2022-09 - 2023-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    2023 Updates and enhancements. Demo H-bonds of a residue
#                    and zone selections.
#           1.1    Clean up post lecture.
#           1.0    First version for lecture.
#
#
# PREREQUISITES:
#   Installation of a recent version of ChimeraX
#   Some familiarity with ChimeraX command syntax
#
# TODO:
#    ...
#
#
# ==============================================================================
#


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                                  Line
#TOC> ----------------------------------------------------------------------
#TOC>   1        ChimeraX REMOTE SCRIPTING                                51
#TOC>   1.1        Defining a Port                                        69
#TOC>   1.2        Prepare ChimeraX                                       92
#TOC>   1.3        A function to communicate with ChimeraX               123
#TOC>   2        Visualize the 7ACN fold.                                136
#TOC>   2.1        Composition of the protein's structure                158
#TOC>   2.1.1          Visualize a "slice" through the protein           189
#TOC>   2.1.2          Continuous translation  via a for-loop ...        207
#TOC>   2.2        Study a single amino acid in context:                 247
#TOC>   2.2.1          Side chain packing and steric complementarity     256
#TOC>   2.2.2          Hydrogen bonds                                    314
#TOC> 
#TOC> ==========================================================================


#
# =    1  ChimeraX REMOTE SCRIPTING  ===========================================


# One of the cool features of ChimeraX is that it can be driven by Python code,
# both within a running session and through Python scripts. What I find even
# cooler though is that ChimeraX can be driven from any programming language via
# its remote control function that can listen to commands sent from any other
# application. The interface that is used here is the standard REST (method) -
# the GET and POST verbs that underly the communication of clients and servers
# everywhere on the Web.

# In order to establish the communication between this script and ChimeraX, all
# we need to do is:
#  - open ChimeraX;
#  - tell it to listen on a specific "port";
#  - send commands to that port via httr:: package functions.


# ==   1.1  Defining a Port  ===================================================

# The httr:: package needs to be installed. The code is here, just for
# completeness but if the start-up script have executed correctly, the
# installation has already been done (from the .utils.R script).

# if (! requireNamespace("httr", quietly = TRUE)) {
#   install.packages("httr")
# }


# We need to define a port - that's something like a mailbox number which
# processes agree on to share information. Any available port number between
# 49152-65535 is fine. We'll choose 61803. Actually, this too has already been
# taken care if in the installation script. You do not need to execute this
# again.

# CXPORT <- 61803

# Check that our current version of R supports sockets (default since V 3.3):
capabilities("sockets")   # MUST be TRUE. If not, don't continue.


# ==   1.2  Prepare ChimeraX  ==================================================

#  - Open a fresh, new session of a recently updated version of ChimeraX.
#  - In the ChimeraX commandline type:
#
#
#       remotecontrol rest start port 61803
#
#
#    (... where 61803 is the port number that you have assigned to CXPORT.)

# Now watch what happens in ChimeraX when you execute the following line in
# RStudio:

( x <- httr::GET("http://127.0.0.1:61803/run?command=open+7ACN") )

# This is a pretty big deal: you can remote-control the ChimeraX viewer
# from R scripts!

# Notes:

#   (127.0.0.1 is the "localhost" IP address. Only this one IP address is
#   supported for security reasons, but ChimeraX can be opened: cf.
#   https://www.rbvi.ucsf.edu/pipermail/chimerax-users/2017-July/000113.html .)

#   (7ACN is a structure of Aconitase - an enzyme of the citric acid cycle. Cf.
#     - https://pdb101.rcsb.org/motm/154
#     - https://pdb101.rcsb.org/motm/89
#     - https://en.wikipedia.org/wiki/Citric_acid_cycle )


# ==   1.3  A function to communicate with ChimeraX  ===========================
#
# To make it convenient to send commands to Chimera, I have provided a
# function that abstracts away the technical details. The function CX()
# is loaded from the .util.R script on start up:
#   - it accepts ChimeraX commands as a parameter;
#   - it sends the command to Chimera, correctly formatted;
#   - it captures output that may be sent from ChimeraX.
#
# The function is loaded in the environment, you can click its name and
# examine the code. Well make extensive use of it below.


# =    2  Visualize the 7ACN fold.  ============================================
#
# We can use the CX() function to set up a script for the visualization of
# the 7ACN fold. Execute these commands one by one to walk through the
# demonstration.

CX("close")        # close all models, start in a defined state.
CX("open 7ACN")
CX("camera sbs")   # practice your stereo vision
# CX("camera mono") # This would turn stereo off again
CX("lighting soft")
CX("lighting shadows true intensity 0.8")

# Define a color gradient to be emphasize the fold of the protein
# cf. https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html
CX("color sequential #1 & protein target abc palette #B0C4DE:rgb(44,50,56):orchid:maroon")

# Study the fold. Note where it begins, how different elements of "secondary
# structure" (helices, strands) form compact domains, and how the domains
# assemble to a structural whole.


# ==   2.1  Composition of the protein's structure  ============================

# Now change the view to show the actual atoms (as little spheres.)
#
CX("hide cartoons")
CX("show atoms")
CX("style sphere")

# Notice the little red spheres? These are water molecules (H2O); what you see
# is the position of each oxygen atom. We hide them so we don't clutter
# the screen too much.
CX("hide ::name='HOH'")

# First we color all atoms gray ...
CX("color #1 & protein #434c56")

# Then we select two categories of amino acids:
#   First, we select side-chains of hydrophilic amino acids ...
#   (these have atoms in their side-chains that make favourable
#   interactions with water)
CX("select :SER, ASN, ASP, GLN, GLU, ARG, LYS, HIS & sidechain")
CX("color sel #173599")

#   Then we select hydrophobic amino acids ...
#   (these have mostly only carbon atoms in their side chains, and they cannot
#    form hydrogen-bonds, which gives them an "oily" nature).
CX("select :VAL, LEU, ILE, MET, PHE, TRP & sidechain")
CX("color sel #dced21")
CX("select clear")  # clear the selection so it doe not clutter our view


# ===   2.1.1  Visualize a "slice" through the protein      

# set clipping planes to show a slice of the protein
CX("clip near -2.0 far 2.0 position #1")

# set the lighting model to close the slice
CX("lighting flat")
CX("graphics silhouettes false")

# move the structure away from us and towards us...
CX("move z 1.0")
CX("move z 1.0")
CX("move z 1.0")
CX("move z -1.0")
CX("move z -1.0")
CX("move z -1.0")


# ===   2.1.2  Continuous translation  via a for-loop ...   
# ... instead of manually moving along a a view step by step:

CX("clip near -32.0 far -28.0 position #1")

delta <- 0.06     # step-size
N     <- 1000     # number of steps
for (i in c(1:N, -(1:N))) { # Do this N times forward and N times back
#  print(sprintf("move z %f", sign(i) * delta))
  CX(sprintf("move z %f", sign(i) * delta))
  CX("wait 1")
}

# Watch the protein slowly move through the slice-zone. Pay particular attention
# to where the blue and yellow side chains are located. Do you see them in the
# interior of the protein, or on the surface? Are they clustered, or evenly
# distributed? Is their arrangement random, or non- random?

# (If this takes too long and you want to abort it, just click on the red
# "Stop sign" in the menu bar of the RStudio console.)

# As an alternative to a sliding motion, we can rotate the protein slowly
# around its active site:
CX("cofr #1:755")
CX("clip near -2.5 far 2.5 position #1:755")
CX("roll y 0.1")   # (Documentation for all commands is available via the
                   #  ChimeraX help functions. Easiest access is
                   #  by clicking the hyperlinked command in the log window
                   #  of Chimerax. It is either open, or can be opened via
                   #  the menu: Tools > Log.
                   #
CX("stop") # stop the rotation

# Task:
# =====
#   What have you observed?
#   Why is this so?
#   What does this imply?


# ==   2.2  Study a single amino acid in context:  =============================

# Reset our view ...
CX("clip off")
CX("hide atoms")
CX("show cartoons")
CX("color sequential #1 & protein target abc palette #B0C4DE:rgb(44,50,56):orchid:maroon")


# ===   2.2.1  Side chain packing and steric complementarity

# Let's focus on a single amino acid and examine its context: how does it fit
# into the protein that surrounds it.

# Create a second copy of the protein and delete everything but this single
# amino acid, to work with a single residue only. (Isn't there a way to
# duplicate just one residue or other molecule? IDK.)
CX("open 7acn")  # Opens the same structure as model #2

# Open the Sequence viewer window:
CX("sequence chain #2/A")

# To choose an amino acid of a certain type hover over it.
# You will get information like "7acn #1 \A ARG 114 CG" - which is:
#
#   7acn: the PDB ID
#   #1:   model number 1
#   \A:   chain A
#   ARG:  amino acid type (ARG -> arginine)
#   114:  position in sequence  ("sequence number")
#   CG:   atom topology, "Cgamma", third atom of the side-chain

# Let's display a residue together with its solvent accessible surface and see
# how it packs into the rest of the protein. We'll choose a histidine because
# histidine has a pretty shape. For example HIS 298. You could repeat this with
# your "adopted" amino acid. We define HIS 298 as a named selection "myAA":
CX("name myAA #2/A:298")              # define
CX("select #2/A")                     # select the entire chain
CX("select subtract myAA")            # subtract the named selection
CX("delete sel")                      # delete what is selected, myAA remains
CX("view myAA")                       # center the view
CX("show myAA atoms")                 # show atoms, and ...
CX("style myAA ball")                 # style as ball-and-stick
CX("color myAA byelement target ab")  # color by element
CX("surface myAA")                    # compute surface
CX("color myAA #3674b3 target s")     # color wispy pale blue
CX("transparency myAA 50 target s")   # set 50% transparency

# Now lets calculate a surface of the original protein _without_ myAA, to see
# how the amino acid is packed against the rest of the structure.
CX("name myAA")                       # list the selection again to remind us
CX("select #1/A:298")                 # select it in model #1 ...
CX("delete sel")                      # ... and delete it
CX("surface #1")                      # compute a surface of model #1 (without
                                      # HIS 298 of course, that was deleted.)
CX("color #1 #3d5c58 target s")       # color it
CX("transparency #1 70 target s")     # set transparency
CX("cofr myAA")                       # set center of rotation (cofr)
CX("clip near -10.0 far 20.0 position myAA")  # open up the clipping slice a bit

# Now rotate, scale, and inspect ... note how snugly the shape of the histidine
# fits into the matrix of surrounding amino acids. We call this "steric
# complementarity" and it is one of the most important factors that make
# proteins fold into a specific structure. This promotes atom-atom contacts,
# which enable the cohesive "Van der Waals forces" to hold everything together.


# ===   2.2.2  Hydrogen bonds                               

# Another important source of fold-specificity are "hydrogen bonds". These
# are  "non-covalent bonds". They are weaker than a chemical bond
# but functionally important because they can be formed and released relatively easily.

# Hydrogen bonds (H-bonds) can be formed between hydroxyl- (-O-H), amino-
# (-N-H), or thiol (-S-H) groups, and acceptor oxygen (O-) or nitrogen (N-)
# atoms. So a typical hydrogen bond could be indicated like (-O-H⋯O-). In a
# typical protein, MOST H-bond donors and acceptors are actually involved in a
# bond and participate in extended, mutually stabilizing networks of
# cross-linking interactions.

# Let's close everything and reload 7ACN to bring our model back
# into a sane state:

CX("close")
CX("open 7ACN")
CX("camera sbs")
CX("cofr #1:298@CE1")
CX("name myAA #1:298")
CX("color sequential #1 target abc palette #213754:#482154:#542153:#542121")
CX("color myAA #09ad97")
CX("cartoon hide")
CX("show target ab")
CX("style ball")
CX("color #1@N,N* #0858a3")  # Nitrogen atoms blue
CX("color #1@O,O* #a3081d")  # Oxygen atoms red
CX("lighting soft")
CX("lighting shadows true")
CX("clip near -10.0 far 20.0 position myAA")

# Now we calculate and display the H-bonds
CX("hbonds")


# Have a look, rotate, scale, and explore, but this scene can get confusing
# because of the number of atoms and bonds that are shown. Let's focus only on
# the environment around His 298 to show how this particular residue is held in
# a network of stabilizing bonds.

CX("~hbonds")  # Turn hbonds off

# First we select a zone of residues around His 298, which have at least
# one atom within 10 Å (1 nm) of His 298.
CX("select zone myAA 10.0 #1 extend true residues true")
CX("name frozen myZone sel") # name the selection and don't change it later
CX("select")                 # select all
CX("~select myZone")         # de-select the atoms in myZone
CX("hide sel target ab")     # undisplay all that remained  selected
CX("~select")                # de-select
CX("clip near -8.0 far 8.0 position myAA")  # tighten the clipping slice

# Now we display only specific H-bonds. We use the distance command, it allows
# us to draw pseudobonds between specific atoms.

# Set the pseudobond display style:
CX("distance style color #c1daee decimalPlaces 2")

# Here are the two H-bonds that the side chain of His 289 forms:
CX("distance #1:298@ND1 #1:1069") # ... with water molecule HOH 1069
CX("distance #1:298@NE2 #1:22@OE2") # ... with a carboxyl-oxygen of Glu22

# Here are the two H-bonds that hold the backbone atoms of His 289 in place:
CX("distance #1:298@N #1:295@O") # ... His298 N with Phe295 O
CX("distance #1:299@N #1:296@O") # ... Leu299 N with Lys296 O

# These bonding partners are themselves held in place.
# The water molecule HOH 1069 forms H-bonds ...
CX("distance #1:1069  #1:41@NZ")  # ... with the amino group of Lys 41
CX("distance #1:1069  #1:29@NH2") # ... with an amino group of Arg 29
CX("distance #1:1069  #1:297@O")  # ... with the backbone carbonyl of Asp 297

# and the carboxylic acid side chain of Glu 22 is fixed by:
CX("distance #1:22@OE1  #1:2@NH1") # ... both amino functions of Arg 2.
CX("distance #1:22@OE1  #1:2@NH2")

# Obviously, any other residue than histidine would not fit so well.
# Non-disruptive replacements would need to keep the shape and bonding
# properties intact.

# More things we could do is to:
#   - label the residue (command "label");
#   - export a high-resolution image for our journal (command "save");
#   - ... and more.


# When we are done ...
CX("remotecontrol rest stop")  # ... release the socket.



# [END]
