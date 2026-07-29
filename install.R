# install.R -- one-shot installer for the R packages used across the IMCS lectures.
#
# Each lecture also loads (and auto-installs) its own packages at knit time via
# pacman::p_load(), so running this first is optional -- it just avoids surprise
# installations in the middle of a knit.

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

cran_packages <- c(
  # document formats
  "rmarkdown", "knitr", "bookdown", "tint", "formatR", "fontawesome",
  # data handling and plotting
  "tidyverse", "data.table", "here", "skimr", "magrittr",
  "gridExtra", "ggExtra", "RColorBrewer", "VIM",
  # distributions and simulation
  "mvtnorm", "EnvStats", "flexsurv", "rsimsum",
  # estimation and inference
  "lmtest", "sandwich", "coin", "boot", "MASS", "Publish",
  # tooling
  "remotes", "devtools"
)

new_packages <- setdiff(cran_packages, rownames(installed.packages()))
if (length(new_packages) > 0) install.packages(new_packages)

# GitHub-only packages
# AIPW is not loaded by any lecture at knit time; it is the worked example in
# the intro's GitHub-installation walkthrough, so pre-installing it here spares
# students the API-limit errors that walkthrough describes.
if (!requireNamespace("AIPW", quietly = TRUE)) {
  remotes::install_github("yqzhong7/AIPW")
}
if (!requireNamespace("looplot", quietly = TRUE)) {
  remotes::install_github("matherealize/looplot")
}

# LaTeX: the lectures compile with xelatex. If you do not have a TeX
# distribution, run:
#   install.packages("tinytex"); tinytex::install_tinytex()
# tinytex installs any missing LaTeX packages automatically on first compile.
