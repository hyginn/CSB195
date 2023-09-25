# GCtables.R
#
# Purpose: Make alternate tables for presenting the genetic code
# Version: 0.1
# Date:    2023-09-25
# Author:  CSB195
#
# Input:
# Output:
# Dependencies:
#
# ToDo:
# Notes:
#
# ==============================================================================


# ====  PARAMETERS  ============================================================
# Define and explain all parameters. No "magic numbers" in your code below.

AADAT <- "data/GeneticCode.csv"


# ====  PACKAGES  ==============================================================
# Check that required packages have been installed. Install if needed.



# ====  FUNCTIONS  =============================================================

# Define functions or source external files
source("<myUtilityFunctionsScript.R>")

myFunction <- function(a, b=1) {
	# Purpose:
	#     Describe ...
	# Parameters:
	#     a: ...
	#     b: ...
	# Value:
	#     result: ...

	# code ...

	return(result)
}



# ====  PROCESS  ===============================================================
# Enter the step-by-step process of your project here. Strive to write your
# code so that you can simply run this entire file and re-create all
# intermediate results.


AAdf <- read.csv(AADAT)





# ====  TESTS  =================================================================
# Enter your function tests here...
if (FALSE) {
  # This block is never executed. Place all tests and experiments here that
  # should not be run when this file is source()'ed.





}

# [END]
