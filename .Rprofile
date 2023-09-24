# .Rprofile
#
# This file is automatically run when your R project starts up.

# Set a CDN backed CRAN mirror
options(repos = c(CRAN="https://cloud.r-project.org/"))

# Run initialization scripts
cat("\n")
cat("Hello - welcome to the CSB195 class project.\n\n")

cat("Sourcing utility scripts ...\n")
source(".util.R")
if (file.exists("../dev/initDev.R")) { source("../dev/initDev.R") }

cat("\n")
cat("Initializations complete.\n")

cat("\n")
cat("Please remember to \"pull\" updated code from the GitHub repository.\n")
cat("Also, from time to time run:\n")
cat("  update.packages(ask = FALSE, type = \"binary\")\n")
cat("... to bring your installed packages up to date.\n")
cat("\n")


# [END]
