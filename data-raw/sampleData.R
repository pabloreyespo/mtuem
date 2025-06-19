enut.i <- readRDS('data-raw/enut_i.rds')
enut.ii <- readRDS('data-raw/enut_ii.rds')
maed <- readRDS('data-raw/maed.rds')

usethis::use_data(enut.i)
usethis::use_data(enut.ii)
usethis::use_data(maed)
