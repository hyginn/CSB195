# tocID <- "R-intro-01-Using_the_IDE.R"
#
# Purpose: Get familiar with the IDE. Play with
#          code. Use a generative AI (gAI) assistant.
#
# Version: 2.0
#
# Date:    2016-12  -  2023-10
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 2.0    Renamed from R_Exercise-BasicSetup.R in the eponymous project;
#          added Generative AI topics.
# V 1.2.2  2022 - Maintenance
# V 1.2.1  2021 - Maintenance
# V 1.2    Maintenance, and removing files that are no longer needed
# V 1.1    Add fast phi calculation and complex expression selection
# V 1.0    Golden Spiral example
# V 0.2    Additions to integrate with tutorial
# V 0.1    First code
#
# TODO:
#


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                         Line
#TOC> -------------------------------------------------------------
#TOC>   1        The RStudio layout                              49
#TOC>   2        Executing code                                  81
#TOC>   2.1        Up-arrow recalls previous commands           114
#TOC>   2.2        History tab                                  133
#TOC>   2.3        Command-specific history                     139
#TOC>   2.4        Cmd-Enter to execute selected code           148
#TOC>   2.5        Analyzing expressions inside-out             227
#TOC>   3        Files, and the working directory               323
#TOC>   4        Change some code and save the change           339
#TOC>   5        Quit, and restart where you left off           370
#TOC>   6        Shortcuts and autocomplete                     388
#TOC>   7        Summary                                        455
#TOC>
#TOC> ==========================================================================


# This is our starting point for actually working with R. Please go through
# the code below, read it carefully, and follow the tasks.


# =    1  The RStudio IDE layout  ==============================================

# Script Pane:
# ============
# In the default layout of RStudio's panes, this script was opened in
# the top-left pane; I refer to this as the "Script Pane". This is where you
# read and edit files.

# Console:
# ========
# The pane below (bottom-left) is the "Console". This is where you type
# commands, or execute them from a script. This is also the pane to where
# RStudio sends its its output.

# Environment Pane:
# =================
# The pane at the top-right is the "Environment Pane". This is where you can see
# which variables or functions have currently been defined ("assigned") and what
# their values are. It also contains the Git tab, which we need from time to
# time to update the project from GitHub.

# Files Pane:
# ===========
# The pane at the bottom-right displays different kinds of information in its
# tabs: In the "Files Tab" you see a directory listing of the contents of the
# current working directory. The "Plots Tab" shows plots and graphs that R has
# created. And the "Help Tab" displays the information in R's help system.

# You can drag the dividers between the panes to resize them.



# =    2  Executing code  ======================================================

#  You can execute code by typing it into the console.
#  Type: 3 + 4  and then hit <enter>
#
#  Now type: 1 + sqrt(5) / 2
#  (This is NOT the Golden Ratio, phi, isn't it?)

# Let's get our AI tutor involved right away.

# Task:
# =====
#
# Open an AI assistant session:
#
#  (1) Open a new session with your favorite generative AI assistant (gAI) in
#      a Web browser window (I use ChatGPT-4, but any gAI should work).
#  (2) In the console type: gAIinit() and hit <return>.
#      This places an initialization command into your clipboard.
#  (3) In your gAI dialogue, <paste> the command and send it off. A typical
#      response should be: "Confirmed. #"
#  (4) Now copy the following request and paste it into the dialogue
#
#         What is an R expression for the Golden Ratio (phi)?
#
#  (5) Look at the response. Is it adequate? Ask any follow-up question
#      you might be interested in. Then ask the following question:
#
#         I was told this is NOT phi: 1 + sqrt(5) / 2 ... but why?
#
#  (6) You should get an explanation.


# ==   2.1  Up-arrow recalls previous commands  ================================
#
#  To correct your previous expression and make it compute the Golden Ratio,
#  use the up-arrow key to get the expression back on the command line,
#  then use the left and right arrow keys to edit the expression and
#  correct it, like you were instructed by the gAI. Then press <return>.
#
#
#  Task:
#  =====
#
#  Correct the previous (wrong) expression for phi
#
#  If your correction was correct,  R should respond with
#
#       [1] 1.618034
#
#  If you don't get exactly the correct value, figure out how to fix this.

# ==   2.2  History tab  =======================================================
#
# Another way to access previously typed code is through the History tab in the
# Environment Pane. Scroll to find a line that you want to execute or edit,
# double-click it, and it will appear in the Console. Try this.

# ==   2.3  Command-specific history  ==========================================
#
# Even more convenient is to type a few characters of a previous command, then
# hit <cmd><up-arrow>. A little window with previous commands that begin with
# the characters you just typed pops up and if you select one and hit <enter>,
# it gets loaded into the command line and you can either hit <enter> to execute
# it again, or edit it. This won't make a lot of sense now because you haven't
# entered many commands yet - but it will be quite useful later.

# ==   2.4  Cmd-Enter to execute selected code  ================================
#
# TAKE NOTE! WHAT FOLLOWS IS SUPER USEFUL!
#
# You can ALSO execute code by selecting code in the Script Pane (where you are
# reading right now) and pressing <command><enter> (or <ctrl><enter> on
# Windows). This passes code from the Script Pane to the console and executes it
# automatically.

# This is very convenient - in fact, this is the preferred way to
# work with code in RStudio. Try this immediately: place the cursor anywhere
# into the expression below and type <command><enter>:

cat("\n  /\\_/\\\n ( o.o )\n  > ^ <\n\n")

# (No, the function cat() is not specifically meant to draw images of cats!
# Its name comes from conCATenate.)

# There are several variants to selecting and executing code. Read this
# carefully and try the alternatives - these three crazy tricks will change your
# life:

#  - When nothing is selected in the Script Pane, pressing <command><enter> will
#    execute the current line and move the cursor to the next line. You can
#    "walk" through code line by line in this way.
#
#  - This is the same as if you would have selected an entire line.
#
#  - If you select more than one line, you can execute an entire block of code
#    at once. Try this: select the block of code below, (from the first
#    assignment of FibPrev to the print() statement) and hit <command><enter>
#    to execute it. The code calculates successive approximations of the Golden
#    Ratio from Fibonacci numbers, and then prints the "true" value.

# Select from below here ...
FibPrev <- 1L    # assign integer 1
FibCurr <- 1L
for (i in 1:40) {
    print(FibCurr / FibPrev, digits = 22) # approximate phi as the ratio of
                                          # two successive Fibonacci numbers
    tmp <- FibCurr
    FibCurr <- FibPrev + FibCurr          # calculate the next Fibonacci number
    FibPrev <- tmp                        # swap
}
print((1 + sqrt(5)) / 2, digits = 22)     # the real Golden Ratio for comparison
print(FibCurr)                            # the Fib_40 is pretty large
# ... to above here.

# Nb. This is a slow algorithm - it adds about a digit of accuracy every 3.5
# iterations. While we're here, let's play a bit: Here is a MUCH faster
# algorithm, based on the following identities:
#    Given F[k] and F[k + 1]:
#        F_2k   == F_k * ((2 * F_k+1) - F_k)
#        F_2k+1 == F_k^2 + F_k+1^2

FibPrev <- 1    # assign 1
FibCurr <- 1
for (i in 1:7) { # Only 7 iterations needed
    print(FibCurr / FibPrev, digits = 22) # approximate phi as the ratio of
                                          # two successive Fibonacci numbers
    Fk  <- FibPrev
    Fk1 <- FibCurr
    FibPrev <- Fk * ((2 * Fk1) - Fk)  # calculate F_2k
    FibCurr <- Fk^2 + Fk1^2           # calculate F_2k+1
}
print((1 + sqrt(5)) / 2, digits = 22)   # The real Golden Ratio (approximately)
print(FibCurr, digits = 22)             # Huge number, efficiently found


# This is just playing with numbers (i) for you to practice reading code, (ii)
# experimenting with alternatives (for example you could install the Rmpfr::
# package for arbitrarily accurate computations, and calculate phi to 100 digits
# accuracy), and (iii) understand that efficient algorithms can be _very_ much
# faster than naive approaches to calculate. Often, whether something can be
# practically calculated at all depends on whether you can approach the problem
# in a way that it maps to a problem for which an efficient algorithm is
# available.


# ==   2.5  Analyzing expressions inside-out  ==================================
#
# Selecting and executing code is also really useful to analyze complex, nested
# R expressions from the inside out. To do this, you select _less_ than
# one line of code.

# Here's an example: assume you would like to try an intermittent fasting
# regime, fasting every fifth day, starting tomorrow. What will the next five
# weekdays of fasting be?

format(Sys.Date() + 1 + seq(0, length.out = 5, by = 5), "%a")

# That seems to be the correct result, but how on earth does this work?

# Task:
# =====
#
# Let's ask the gAI: select the following code block, execute it, and paste it
# into the gAI dialogue. (If you need to open a new dialogue, you need to
# initialize it again with the basic R-tutor persona instructions from the
# gAIinit() command.) Execute the following command which places a multi-line
# dialogue on the clipboard.

t2c("
Can you explain what the expression below does and how it works?
format(Sys.Date() + 1 + seq(0, length.out = 5, by = 5), \"%a\")
")

# Now <paste> this into your gAI dialogue. You will get useful information.
# But I have also seen Bing make mistakes. Let's take this step by step.
#
# To understand this, we analyze the this expression from its components.
# Here it is again, for convenience:

format(Sys.Date() + 1 + seq(0, length.out = 5, by = 5), "%a")

# First selectSys.Date() and execute it. You should get today's date.

# Next select Sys.Date() + 1 and execute that. This should give you
# tomorrow's date.
#
# Next select seq(0, length.out = 5, by = 5) and execute it. You should get
# [1]  0  5 10 15 20 ... i.e a vector of numbers, starting from 0, 5 numbers
#                            in all, each incremented by 5.
#
# What do we get when we add one to this vector? Select
# 1 + seq(0, length.out = 5, by = 5)
# and execute it to find out.
#
# And what do we get when we add a date to this resulting vector? Select
# Sys.Date() + 1 + seq(0, length.out = 5, by = 5)
# and execute it. These are the next five fasting days.

# And finally, what weekdays are these? format() is a versatile function that
# can handle many different types of R objects. We call this a "generic"
# function. If the object is of class "Date", the function sends the data to the
# format.Date() function, which knows how to convert this into weekdays (see
# ?strptime for a complete list of format options).

# Get into the habit to decompose R expressions and analyze what the individual
# components do, how they are nested, and how they work together. This is what
# reading code is about. If you just gloss over such complicated expressions
# you won't learn.
#
#    Use your gAI to help you with understanding syntax. However, don't
#    ------------------------------------------------------------------
#
#    assume it understands logic as well. The gAI is a LANGUAGE MODEL, not
#    ---------------------------------------------------------------------
#
#    a logic engine.
#    ---------------
#

# I don't always write code in such nested ways for two reasons: (1) its
# not always obvious what happens, and code MUST be obvious or it is bad code;
# and (2) this is hard to debug if something goes wrong. Thus for course code
# I would probably write the same thing as:

x <- Sys.Date() + 1
x <- x + seq(0, length.out = 5, by = 5)
format(x, "%a")

# In this case we can watch how x changes line by line in the Environment
# Pane (top right) and check that what we do is correct.
#
#    N.b. Some people advocate using the pipe operator "%>%" from the magrittr::
#    package (also part of the "tidyverse"). And they often write very
#    enthusiastically about its benefits. After many years of considering it
#    I still think this is a poor idea. If you like, we can discuss this at some
#    point. But the idiom of using intermediate assignments to some variable
#    x or similar is what you will encounter most frequently in this course.
#    It simply makes the process more explicit. And in software development
#    EVERYTHING is about making things explicit.


# =    3  Files, and the working directory  ====================================
#
# We need to make sure the working directory is set correctly. Execute
list.files()
# ... to get a list of files in the current working directory.
# If this is not the list of files you see in RStudio's  Files Pane, execute
getwd()
# ... to see what the working directory is set to. Then use the menu:
# Session -> Set Working Directory -> To Project Directory
# This executes the correct setwd() command. Confirm by executing
list.files()
# ... again. But this should really not be necessary. If it is, there may
# be a configuration problem. Talk to me to fix this.
#


# =    4  Change some code and save the change  ================================
#
# Here is a bit of code that plots 1.5 turns of the Golden Spiral. Don't worry
# if you don't understand the details - that's not the point right now
# (and you can always ask the gAI to explain it.). Execute the code.

turns <- 1.5                                 # number of turns to plot
p <- seq(0, (turns * 2 * pi),
         length.out = round(turns * 100))    # intervals, in radians
plot(exp(0.3061868562 * p) * cos(p),         # x-coordinates
     exp(0.3061868562 * p) * sin(p),         # y-coordinates
     type = "l", col = "#CC0000",            # draw as a red line
     xaxt = "n", yaxt = "n",                 # turn tick-marks off
     xlab = "", ylab = "",                   # turn axis labels off
     asp = 1.0)                              # plot with 1:1 aspect ratio
abline(v = 0, col = "#D2DAFF")               # add pale-blue vertical ...
abline(h = 0, col = "#D2DAFF")               # ... and horizontal line

# Now change the value of "turns" from 1.5 to 2.75 in the code above. Note that
# the filename in the Script Pane changes colour: it turns red and has a star
# behind it. This means: this script has unsaved changes.

# Execute the modified block of code. The number of the turns increases.

# Now save the change. The filename changes to black again. The Modified date
# and time of the script in the Files Pane is _now_.

# Note that this is a file from the course project. Your changes will get
# overwritten when you necxt "pull" code from the repository.


# =    5  Quit, and restart where you left off  ================================
#
# Finally, we'll quit RStudio and restart the project. It should be obvious how
# to quit. But here's how to reload a project.

# Choose: File -> Recent Projects -> <the project you need>

# Note: Be sure to choose the "Recent Projects" menu option,
#       not "Recent Files".

# Try this out, use the File -> Quit Session... option to quit - and if you are
# being asked if you want to save anything, for the purpose of this course or
# workshop _always_ answer "No". Then reload the project, note that the value of
# turns should now be what you had changed it to.

# Plot another beautiful spiral for good measure.


# =    6  Shortcuts and autocomplete  ==========================================

# RStudio's editor is actually one of the best text-editor available. It has
# has a great number of useful shortcuts and aids for efficient work - not just
# for coding. Try them.
#
#   - Typing an opening parenthesis automatically types the closing parenthesis.
#
#   - Selecting text and typing a quotation mark quotes the selected text.
#
#   - Both bracketing and quoting of selected code work with single- and
#     double quotation marks, parentheses, square brackets and curly braces.
#
#   - Typing a newline character automatically indents the following line.
#
#   - Typing a newline inside a comment automatically prefixes the line
#     with a comment character. (To get a plain line back, type <cmd><return>.)
#
#   - Typing the first three letters of a variable or function name lists
#     possible choices to autocomplete. Use the arrow keys to select, hit
#     <enter> or <tab> to execute. Note that this is not case sensitive.
#     Autocomplete is therefore especially useful if you are not sure about
#     the correct case. Try it: type  Tur  and you should get the variable
#     "turns" selected.
#
#   - typing <alt>-<up-arrow> shifts a line (or selected block) up or down.
#     (I use this a lot :-)
#
#   - holding <alt> while selecting things allows you to select columns.
#     Try this: in the block below select and delete the extra two
#     columns of sevens:
#
#     012345677789
#     012345677789
#     012345677789
#
#   - note how the insert cursor spans all three lines. You can also use this
#     to type into several lines at once! Try this: place the cursor behind the
#     three nines, then type abcdef. Hex.
#
#   - You can customize the keyboard layout and this makes it easy to switch
#     between programs that you are used to without having to remember new
#     every shortcuts every time. One thing I always change in my work is
#     the code to move between open tabs (I usually have quite a few open).
#     For Chrome, on the Mac the code is Alt+Cmd+Left resp. Alt+Cmd+Right to
#     step forward resp. backwards through tabs. On Chrome for Windows it
#     is Ctrl+PgDown resp. Ctrl+PgUp. The RStudio default is Ctrl+Tab resp.
#     Ctrl+Shift+Tab.
#
#     Task:
#     =====
#     Open Tools > Modify Keyboard Shortcuts...
#     Type  tab  into the Filter box.
#     Find "Open Next Tab" and "Open Previous Tab".
#     Click on the current key combination to edit it.
#     Press the same keycode that works for Chrome on your computer.
#     Click apply to save your choices.
#
#     Now open a few files from the  Files  pane and confirm that the
#     new key combo works for you, forward and back.
#
#  -  Sometimes we want to compare files with each other, and switching tabs
#     is not so convenient. You should know that you can UNDOCK tabs too, by
#     dragging them out of the main RStudio window. That's very useful. Try
#     it!
#

# =    7  Summary  =============================================================
#
# Summary:  - The layout of the Panes
#           - Typing code
#           - First steps with a gAI assistant
#           - Editing code in the Console, using the arrow keys
#           - Retrieving code via the History tab
#           - Selecting code in the Script Pane and executing it
#             with <command><enter>
#           - Analyzing R expressions by selecting and executing parts
#           - Setting the Working Directory
#           - Editing code in the script and saving the change
#           - Quitting RStudio and restarting a recent project.
#           - Some amazingly useful typing and selecting shortcuts


# [END]
