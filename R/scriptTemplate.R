# scriptTemplate.R
#
# Purpose:
# Version:
# Date:
# Author:
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



# ====  PACKAGES  ==============================================================
# Check that required packages have been installed. Install if needed.

if (! requireNamespace("here", quietly=TRUE)) {
  install.packages("here")
}

# Note: use package functions with the :: operator - eg.
# here::here("K")



# ====  FUNCTIONS  =============================================================

# Define functions or source external files
source(here::here("<myUtilityFunctionsScript.R>"))

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








# ====  VALIDATION AND TESTS  ==================================================
# Enter validation of your code logic and your function tests here...
if (FALSE) {
  # This block is never executed. Place all tests and experiments here that
  # should not be run when this file is source()'ed.





}

# [END]
