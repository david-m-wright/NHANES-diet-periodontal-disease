# Bookdown build script
library(here)
library(knitr)
library(checkpoint)
checkpoint_date <- "2020-04-02"
checkpoint(checkpoint_date, use.knitr = TRUE, project = here())
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
#preview_chapter("06-socio_economic_status.rmd")
# preview_chapter("07-figures_for_fifth_revision.Rmd")
#preview_chapter("08-references.rmd")


# Render the entire book
bookdown::render_book("index.Rmd", "bookdown::gitbook", new_session = T)

bookdown::render_book("index.Rmd", "bookdown::pdf_book", new_session = T)
