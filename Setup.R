#!/usr/bin/env Rscript

#devtools::source_url("http://dl.dropbox.com/u/113630701/rlibs/base-commons.R")

#devtools::source_url("https://dl.dropboxusercontent.com/u/113630701/datautils/R/core_commons.R")
#devtools::source_url("https://dl.dropboxusercontent.com/u/113630701/datautils/R/ggplot_commons.R")

## Missing dependecies
#install.packages('codetools')

## preselect a cran mirror
r = getOption("repos") # hard code the UK repo for CRAN 
r["CRAN"] = "http://ftp5.gwdg.de/pub/misc/cran/"
options(repos = r)
rm(r)


if (!require("codetools")) install.packages("codetools")
if (!require("devtools")) install.packages("devtools")

scriptsDir=Sys.getenv("TM_HOME")
source(file.path(scriptsDir, "commons/TMCommons.R"))


#devtools::install_github("hadley/lazyeval")
#install.packages("Rcpp")
#install.packages("dplyr")

## spefic libraries needed in som of the workflows

require.auto(zoo) # for na.locf
require.auto(sp)
require.auto(docopt)

require.auto(graph)
require.auto(igraph) ## lineage analysis

## lineage color optimization
require.auto(RBGL) # --> no longer on cran but somehow still in their index
#source("http://bioconductor.org/biocLite.R")
#biocLite("RBGL")

#require.auto(grDevices) ## still needed needed?

## don't install on headnode but worker node because of external binaries
require.auto(rgeos)

require.auto(hash)

## force dplyr 0.4.1 since 0.4.2 is buggy
if(packageVersion("dplyr") != "0.4-1"){
    packageurl = "http://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.4.1.tar.gz"
    install.packages(packageurl, repos=NULL, type="source", dependencies = TRUE)
}


########################################################################################################################
### quit R restart and run



sessionInfo()
## Required minimal resuilt to run the workflow
#> sessionInfo()
#R version 3.1.2 (2014-10-31)
#Platform: x86_64-unknown-linux-gnu (64-bit)
#
#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=de_DE.utf8        LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
#attached base packages:
#[1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base
#
#other attached packages:
# [1] doMC_1.3.3         iterators_1.0.7    foreach_1.4.2      RColorBrewer_1.1-2 png_0.1-7          sqldf_0.4-10       RSQLite_1.0.0
# [8] DBI_0.3.1          gsubfn_0.6-6       proto_0.3-10       data.table_1.9.4   scales_0.2.4       ggplot2_1.0.1      digest_0.6.8
#[15] tidyr_0.2.0        magrittr_1.5       dplyr_0.4.1        reshape2_1.4.1     stringr_1.0.0      plyr_1.8.1
#
#loaded via a namespace (and not attached):
# [1] assertthat_0.1   bitops_1.0-6     chron_2.3-45     codetools_0.2-11 colorspace_1.2-6 devtools_1.7.0   gtable_0.1.2     httr_0.6.1
# [9] MASS_7.3-39      munsell_0.4.2    Rcpp_0.11.5      RCurl_1.95-4.5   stringi_0.4-1    tools_3.1.2
