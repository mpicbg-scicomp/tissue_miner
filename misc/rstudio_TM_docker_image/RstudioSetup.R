#!/usr/bin/env Rscript

## preselect a cran mirror
r = getOption("repos") # hard code the UK repo for CRAN

r["CRAN"] = "http://ftp5.gwdg.de/pub/misc/cran/"
options(repos = r)
rm(r)

# Install packages required to run the TM Manual as a markdown document
install.packages("knitr")
install.packages("htmltools")
install.packages("caTools")
install.packages("bitops")
install.packages("rmarkdown")

# Install packages required to save graphs as svg

#install.packages("svglite")
