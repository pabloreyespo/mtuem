# Copy the enut-i and enut-ii pipeline outputs into data-raw/*.rds, then run
# data-raw/sampleData.R. Run from the package root:
#   Rscript data-raw/update_from_pipelines.R
enut_repos <- Sys.getenv("ENUT_REPOS", "C:/Users/pablo/Documents/GitHub/enut")
sources <- c(
  enut_i = "enut-i/data/enut-i.dta",
  enut_i_raw = "enut-i/data/enut-i-raw.dta",
  enut_ii = "enut-ii/data/enut-ii.dta",
  enut_ii_raw = "enut-ii/data/enut-ii-raw.dta"
)

# Stata stores integers as doubles; keep whole-number columns as integers like
# the previous package data.
as_package_data <- function(dta) {
  data <- as.data.frame(haven::zap_formats(haven::zap_labels(dta)))
  whole <- vapply(data, function(x) {
    is.double(x) && all(is.na(x) | (x == round(x) & abs(x) < .Machine$integer.max))
  }, logical(1))
  data[whole] <- lapply(data[whole], as.integer)
  data
}

paths <- file.path(enut_repos, sources)
names(paths) <- names(sources)
missing <- paths[!file.exists(paths)]
if (length(missing) > 0) stop("Missing pipeline outputs: ", paste(missing, collapse = ", "))
# enut-i/.gitignore excludes *raw.*, so a pull from the cluster can leave an old
# raw file next to a new processed one.
for (survey in c("enut_i", "enut_ii")) {
  raw <- paths[[paste0(survey, "_raw")]]
  if (file.mtime(raw) < file.mtime(paths[[survey]]) - 3600) {
    warning(raw, " is older than ", paths[[survey]],
      "; keeping the current data-raw copy. Copy the new one from the cluster.", call. = FALSE)
    paths <- paths[names(paths) != paste0(survey, "_raw")]
  }
}

for (name in names(paths)) {
  data <- as_package_data(haven::read_dta(paths[[name]]))
  saveRDS(data, file.path("data-raw", paste0(name, ".rds")))
  cat(name, nrow(data), "rows,", ncol(data), "columns\n")
}
source("data-raw/sampleData.R")
