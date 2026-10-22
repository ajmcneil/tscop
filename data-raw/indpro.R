# US industrial production: total index, seasonally adjusted, monthly (Federal Reserve Board, release G.17).
# Run from the package root:  Rscript data-raw/indpro.R
# The Board's file lists the monthly values by year (twelve columns per row), for several series; the total index is B50001
# (2017 = 100). The values are revised regularly, so the dataset records when it was downloaded.
url <- "https://www.federalreserve.gov/releases/g17/ipdisk/ip_sa.txt"
lines <- readLines(url, warn = FALSE)
start <- which(lines == "\"B50001: Total index\"")[1]
rows <- lines[(start + 1):length(lines)]
rows <- rows[startsWith(rows, "\"B50001\"")]            # the rows of the total index (stop at the next series)
vals <- lapply(rows, function(l) {
  z <- strsplit(trimws(sub("^\"B50001\"", "", l)), "[[:space:]]+")[[1]]
  list(year = as.integer(z[1]), x = suppressWarnings(as.numeric(z[-1])))
})
year <- rep(vapply(vals, `[[`, 0L, "year"), each = 12)
month <- rep(1:12, length(vals))
index <- unlist(lapply(vals, function(v) c(v$x, rep(NA_real_, 12 - length(v$x)))))
d <- data.frame(date = as.Date(sprintf("%d-%02d-01", year, month)), index = index)
d <- d[!is.na(d$index) & d$date >= as.Date("1960-01-01"), ]
rownames(d) <- NULL
indpro <- d
attr(indpro, "downloaded") <- Sys.Date()
save(indpro, file = "data/indpro.rda", compress = "xz")
cat("indpro:", nrow(indpro), "months,", format(min(indpro$date)), "to", format(max(indpro$date)), "\n")
