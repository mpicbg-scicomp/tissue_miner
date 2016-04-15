
## make sure that TM_HOME is defined and points to a directory
scriptsDir=Sys.getenv("TM_HOME")
if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME not correctly defined (",scriptsDir ,")"))
}


## silence all the package loading
# (see http://stackoverflow.com/questions/2723034/suppress-one-commands-output-in-r)
sink(file=file("/dev/null", "w"), type="message")

## Source common functions

if (!require("devtools")) install.packages("devtools")

devtools::source_url("https://raw.githubusercontent.com/holgerbrandl/datautils/v1.13/R/core_commons.R")
devtools::source_url("https://raw.githubusercontent.com/holgerbrandl/datautils/v1.13/R/ggplot_commons.R")
devtools::source_url("https://raw.githubusercontent.com/holgerbrandl/datautils/v1.13/R/datatable_commons.R")

require.auto(sqldf)


source(file.path(scriptsDir, "commons/MovieFunctions.R"))
source(file.path(scriptsDir, "commons/RoiCommons.R"))

## restore R output
sink(file=NULL, type="message")

# disabled for now because not needed for main workflow
#source(file.path(scriptsDir, "commons/MultipleQueriesFunctions.R"))
#source(file.path(scriptsDir, "commons/MultiplePlotsFunctions.R"))

## Set up main parameters

# enable plyr parallelization
sink(file=file("/dev/null", "w"), type="message")
  require.auto(doMC);
sink(file=NULL, type="message")


## todo should go into config
isCluster=Sys.getenv("LSF_SERVERDIR")!=""
isEatonPC=Sys.info()[["nodename"]]=="eaton-pc-2"
registerDoMC(cores=ifelse(isCluster, 5, ifelse(isEatonPC,6,3)));
# require.auto(foreach); require.auto(doMC); registerDoMC(cores=20)

require.auto(dplyr)

## allows processing also on headless systems (like when running as cluster job)
#options(device="png")

## force long integers (cell ids etc) to render without using scientific notation
## http://stackoverflow.com/questions/9397664/force-r-not-to-use-exponential-notation-e-g-e10
options("scipen"=100)

####################################################################################################
## ensure sufficient package versions

## see http://stackoverflow.com/questions/9314783/require-minimum-version-of-r-package
check_version = function(pkg_name, min_version) {
    cur_version = packageVersion(pkg_name)

    if(cur_version < min_version){
        stop(sprintf("Package %s needs a newer version, found %s, need at least %s", pkg_name, cur_version, min_version))
    }
}

check_version("dplyr", "0.4-3")
check_version("ggplot2", "2.0.0")
#check_version("magrittr", "1.5")

####################################################################################################
## Helper functions

## Establish DB connection ####
openMovieDb <- function(movieDir){
  
  # Description: open a SQLite database connection
  # Usage: openMovieDb(movieDir)
  # Arguments: movieDir = path to a given movie folder
  
  db_name=basename(movieDir)
  dbFile=file.path(movieDir, paste0(db_name, ".sqlite"))
  
  if(str_detect(dbFile, "project-raphael")) {
    dbSizeBytes=file.info(dbFile)$size
    
    if(is.na(dbSizeBytes) || dbSizeBytes==0) stop(paste0("db '",dbFile,"'in is empty or does not exist"))
    
    tmpDbFile <- paste0("/tmp/",db_name, "__", dbSizeBytes, ".sqlite")
    ## copy db to tmp on madmax because sqlite driver doesn't seem to like lustre    
    if(!file.exists(tmpDbFile)){
      echo("creating database copy under '",tmpDbFile,"' for db: ", db_name)
      system(paste("cp ",dbFile,tmpDbFile))
    }else{
      echo("using cached db:", tmpDbFile)
    }
    
    dbFile=tmpDbFile;
  }
  
  db <- dbConnect(SQLite(), dbname=paste0(dbFile))
  
  return(db)
}


ma <- function(x,n=11, ...){as.numeric(stats::filter(x,rep(1/n,n), sides=2, ...))}


rearrange_cell_bonds <- function(df) arrange(df,  frame, cell_id, bond_order)


addCellShapes <- function(dfWithCellIdAndFrame){
    stopifnot(all(c("frame", "cell_id") %in% colnames(dfWithCellIdAndFrame))) # data must have frame and cell_id

    local(get(load(file.path(movieDir, "cellshapes.RData")))) %>%
        inner_join(dfWithCellIdAndFrame, by=c("frame", "cell_id")) %>%
        rearrange_cell_bonds()
}


########################################################################################################################
### DB Query Utils

## refactor into commons
cdByDaughters <- function(db){
  cdEvents <- dbGetQuery(db, "select cell_id as mother_cell_id, left_daughter_cell_id, right_daughter_cell_id from cell_histories")
  cdEventsLong <- subset(melt(cdEvents, id.vars="mother_cell_id", value.name="cell_id"), !is.na(cell_id), select=-variable)

  ## first occurence to daughter
  firstOcc <- dbGetQuery(db, "select cell_id, first_occ from cell_histories")
  merge(cdEventsLong, firstOcc)
}


########################################################################################################################
### Wing comparision --> MOVED TO MultipleQueriesFunctions.R


addTimeFunc <- function(movieDb, df){
  time <-  dbGetQuery(movieDb, "select * from frames")
  timeInt <- cbind(time[-nrow(time),], timeInt_sec=diff(time$time_sec))
  result <- dt.merge(df, timeInt, by="frame")
  return(result)
}


add_dev_time <- function(df){
    ## Given a data.frame with a frame column, we add the developmental time according to the shift definition from the configuration.

    algnData <- dt.merge(df, get_movie_time_shift(df$movie), by="movie") %>% mutate(dev_time=(time_sec+time_shift)/3600) ## 54000 sec offset, namely 15hAPF

    if(nrow(algnData)==0) stop("could not add development time to data. Make sure that get_movie_time_shift is correclyt implemented!")

    return(algnData)
}


## DEPRECATED: use dplyr::distinct
unique_rows <- function(df, columns){
    warning("unique_rows() is DEPRECATED: use dplyr::distinct instead")
    unique(setkeyv(data.table(df), columns)) %>% as.df()
}


## angle conversion  helpers
mod2pi <- function(x) x - floor(x/(2*pi))* 2*pi
modPi <- function(x) x - floor(x/(pi))*pi
# ggplot(data.frame(x=c(-10,10)), aes(x=x))+ stat_function(fun=modPiOv2, n=1000)
modPiOv2 <- function(x) x - floor(x/(pi/2))*pi/2


########################################################################################################################
#### Integrate User Configuration

configFile <- Sys.getenv("TM_CONFIG")

## or fall back to defaults if not defined
if(configFile==""){
    configFile=file.path(scriptsDir, "config/default_config.R")
}

if(!file.exists(configFile)){
    stop(paste("could not find configFile", configFile))
}

print(paste("using config file", configFile))
source(configFile)

#########################################################################################################################
#### Function to open any file from a system command
open_file <- function(path){
  if (Sys.info()["sysname"]=="Darwin") system(paste0("open ",path))
  if (Sys.info()["sysname"]=="Linux") system(paste0("xdg-open ",path, " > /dev/null 2>&1"))
  if (Sys.info()["sysname"]=="Windows") system(paste0("start ",path))
}

## todo post-process to create meaningful defaults (like image size dependent grid size)