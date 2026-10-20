# Precompute the vignettes that need data downloads, long fits or packages that are not on CRAN builds.
# Run from the package root:  Rscript vignettes/knit-vignettes.R
# Each vignettes/NAME.Rmd.orig is knitted to vignettes/NAME.Rmd with its output in place; the figures are moved from the
# knitr default figure/ to NAME-figs/ (R CMD check flags a directory called figure/ as a knitr leftover), so building the
# package does not run the analysis. Downloaded data are cached in vignettes/cache/ (git-ignored, not redistributed), as are
# the knitr caches of the slow chunks; delete the cache to recompute.
owd <- setwd("vignettes")
for (f in list.files(pattern = "\\.Rmd\\.orig$")) {
  name <- sub("\\.Rmd\\.orig$", "", f)
  out <- paste0(name, ".Rmd")
  knitr::knit(f, output = out)
  figs <- paste0(name, "-figs")
  unlink(figs, recursive = TRUE)
  if (dir.exists("figure")) file.rename("figure", figs)
  writeLines(gsub("figure/", paste0(figs, "/"), readLines(out), fixed = TRUE), out)
}
setwd(owd)
