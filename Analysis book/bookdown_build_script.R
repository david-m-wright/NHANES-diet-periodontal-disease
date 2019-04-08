# Bookdown build script
library(here)
library(bookdown)

setwd(here("Analysis book"))

preview_chapter("01-nutrients.Rmd")

bookdown::render_book("index.Rmd", "bookdown::gitbook")

