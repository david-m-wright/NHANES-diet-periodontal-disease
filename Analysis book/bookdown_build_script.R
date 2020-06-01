# Bookdown build script
library(knitr)
library(here)
library(bookdown)
library(servr)

setwd(here("Analysis book"))

# Preview live copy of the book for editing
#serve_book()

# Render a single chapter
#preview_chapter("01-nutrients.Rmd")
#preview_chapter("02-food_groups.Rmd")
#preview_chapter("03-figures_for_paper.Rmd")
#preview_chapter("04-food_groups_PCA.Rmd")
#preview_chapter("05-figures_for_revised_paper.rmd")
#preview_chapter("06-references.rmd")


# Render the entire book
bookdown::render_book("index.Rmd", "bookdown::gitbook")
