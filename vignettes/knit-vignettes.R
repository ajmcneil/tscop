# Precompute the vignettes that need data downloads, long fits or packages that are not on CRAN builds.
# Run from the package root:  Rscript vignettes/knit-vignettes.R
# Each vignettes/NAME.Rmd.orig is knitted to vignettes/NAME.Rmd with its output and figures (vignettes/figure/) in place,
# so building the package does not run the analysis. Downloaded data are cached in vignettes/cache/ (git-ignored, not
# redistributed); delete the cache to download afresh.
old <- setwd("vignettes"); on.exit(setwd(old))
for (f in list.files(pattern = "\\.Rmd\\.orig$")) {
  knitr::knit(f, output = sub("\\.orig$", "", f))
}
