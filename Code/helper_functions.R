# Functions to use when manipulating datasets

library(tidyverse)
library(broom)
library(scales)

# Percentage formatter
perc <- scales::percent_format()

# Returns the inverse of the logit transformation (from the arm package)
invlogit <-   function (x) 
{
  1/(1 + exp(-x))
}

# Format P values so that those < 0.001 are marked as such and all other returned to three digit precision
FormatPValue <- function(pv){ifelse(pv < 0.001, "<0.001", formatC(pv, digits = 3, format = "f"))}

# Format Odds Ratios or Risk Ratios given estimate and standard error from logistic regression
FormatOddsRatio <- function(estimate, std_error){
  paste0(formatC(exp(estimate), digits = 2, format = "f"), 
        " (", 
        formatC(exp(estimate - 1.96*std_error), digits = 2, format = "f"),
        ", ",
        formatC(exp(estimate + 1.96*std_error), digits = 2, format = "f"),
        ")")
}

# Transform a proportion measures so that all values lie within the interval (0,1) for beta regression
# Arguments: x = vector to transform
# Value: vector transformed to (0,1) interval
PropTransform <- function(x){
  (x * (length(x)-1) + 0.5)/length(x)
}
#PropTransform(seq(0,1, 0.01))

# ntile(1:10, n= 4)
# ntile(10:1, n= 4)


# Function to take quantiles of columns of a matrix and then stack
# so that they are returned in a sparse long format
# Arguments:
# mat = input matrix
# fill_val = value to fill empty columns with (default, 0)
# Value: tibble with quantiles in long format
StackQuantiles <- function(mat, fill_val = 0, ...){
  
  # Take quantiles
  qtx <- apply(mat, 2, function(x) quantile(x, ...)) %>% 
    as_tibble(rownames = "qtile")
  
  # Put into long format  
   bind_rows(replicate(ncol(qtx)-1, qtx, simplify=F)) %>% 
     mutate(component = rep(names(qtx)[-1], nrow(qtx))) %>% 
     gather(key, value, -qtile, -component) %>% 
     mutate(value = if_else(key == component, value, fill_val),
            qtile = as.numeric(str_replace(qtile, "%", ""))) %>% 
     spread(key, value) %>% 
     arrange(component, qtile)
  }



# # Function to present a formatted table of risk ratios from a fitted log-binomial model
# TidyLogBinomal <- function(lb_model){
#   tidy(lb_model) %>% 
#     mutate(Term = term, 
#            `Relative Risk` = estimate, 
#            lcl = estimate - std.error * 1.96,
#            ucl = estimate + std.error * 1.96) %>%
#     mutate_at(vars(`Relative Risk`, lcl, ucl), function(x){formatC(exp(x), digits = 2, format = "f")}) %>%
#     transmute(Term = term, `Relative Risk`, `95% CI` = paste0("(", lcl, ", ", ucl, ")"), 
#               P = FormatPValue(p.value)) %>%
#     filter(Term != "(Intercept)")
# }
# 
