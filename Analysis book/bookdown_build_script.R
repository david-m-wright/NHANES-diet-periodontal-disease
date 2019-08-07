# Bookdown build script
library(here)
library(bookdown)
library(servr)
library(checkpoint)
checkpoint("2019-07-31")

setwd(here("Analysis book"))

# Preview live copy of the book for editing
#serve_book()

# Render a single chapter
preview_chapter("01-nutrients.Rmd")
#preview_chapter("02-food_groups.Rmd")
preview_chapter("03-references.rmd")

# Render the entire book
bookdown::render_book("index.Rmd", "bookdown::gitbook")
