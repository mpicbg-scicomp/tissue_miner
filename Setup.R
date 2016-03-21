#!/usr/bin/env Rscript

## preselect a cran mirror
r = getOption("repos") # hard code the UK repo for CRAN 
r["CRAN"] = "http://ftp5.gwdg.de/pub/misc/cran/"
options(repos = r)
rm(r)

## needed since igraph will ask for it otherwise
userLibDirectory=lib=Sys.getenv("R_LIBS_USER")
dir.create(normalizePath(userLibDirectory), recursive=T, showWarnings=F)
update.packages(ask=F, lib=userLibDirectory)

if (!require("codetools")) install.packages("codetools")
if (!require("devtools")) install.packages("devtools")

scriptsDir=Sys.getenv("TM_HOME")
if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))


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

## needed for svg export with ggsave2
require.auto(svglite)


########################################################################################################################
### quit R restart and run

sessionInfo()
writeLines(capture.output(sessionInfo()), ".tm_setup_sessionInfo.txt")

#R version 3.2.2 (2015-08-14)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS release 6.3 (Final)
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
# [1] hash_2.2.6         rgeos_0.3-12       RBGL_1.44.0        igraph_1.0.1       graph_1.46.0       docopt_0.4.3.3     sp_1.2-0
# [8] zoo_1.7-12         doMC_1.3.3         iterators_1.0.7    foreach_1.4.2      RColorBrewer_1.1-2 png_0.1-7          sqldf_0.4-10
#[15] RSQLite_1.0.0      DBI_0.3.1          gsubfn_0.6-6       proto_0.3-10       data.table_1.9.6   scales_0.3.0       ggplot2_1.0.1
#[22] digest_0.6.8       tidyr_0.3.1        magrittr_1.5       dplyr_0.4.3        reshape2_1.4.1     stringr_1.0.0      plyr_1.8.3
#[29] devtools_1.9.1     codetools_0.2-14
#
#loaded via a namespace (and not attached):
# [1] Rcpp_0.12.1         tools_3.2.2         memoise_0.2.1       gtable_0.1.2        lattice_0.20-33     curl_0.9.3          httr_1.0.0
# [8] stats4_3.2.2        R6_2.1.1            BiocGenerics_0.14.0 MASS_7.3-43         assertthat_0.1      colorspace_1.2-6    stringi_0.5-5
#[15] munsell_0.4.2       chron_2.3-47


########################################################################################################################
## installation test

# not all are listed here but the tricky ones
req_pckgs <- c(
"dplyr", "ggplot2", "doMC", "tidyr", "zoo", "sp", "docopt", "graph", "igraph", "RBGL", "rgeos", "hash", "svglite", "sp"
)

if(all(req_pckgs %in% .packages(all.available=TRUE))){
}else{
    warning(paste("The following packages failed to install:", paste(req_pckgs[!(req_pckgs %in% .packages(all.available=TRUE))], collapse=", ")))
    warning("Please check the installation logs or submit a ticket at https://github.com/mpicbg-scicomp/tissue_miner\n")
    exit(1)
}