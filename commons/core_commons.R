########################################################################################################################
## set a default cran r mirror  and customize environment

#cat(".Rprofile: Setting Goettingen repository\n")
r = getOption("repos") # hard code the UK repo for CRAN(
r["CRAN"] = "http://ftp5.gwdg.de/pub/misc/cran/"
options(repos = r)
rm(r)


## user browser for help
options(help_type="html")

## plot more characters per line
options(width=150)

# for sqldf to avoid the use of tckl
options(gsubfn.engine = "R")


########################################################################################################################
## automatic package installation

require.auto <-  function(x){
  x <- as.character(substitute(x))
  
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text=paste("require(", x, ",  quietly=T)", sep="")))
  } else {
    #        update.packages(ask=F) # update dependencies, if any.
    eval(parse(text=paste("install.packages('", x, "')", sep="")))
  }
  
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text=paste("require(", x, ",  quietly=T)", sep="")))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #        biocLite(character(), ask=FALSE) # update dependencies, if any.
    eval(parse(text=paste("biocLite('", x, "')", sep="")))
    eval(parse(text=paste("require(", x, ",  quietly=T)", sep="")))
  }
}


########################################################################################################################
## load core packages

require.auto(plyr)
require.auto(stringr)
require.auto(reshape2)
#require.auto(reshape2, quietly=T, warn.conflicts=F)

## load on purpose after plyr
require.auto(dplyr)
require.auto(magrittr)
require.auto(tidyr)

## needed for caching
require.auto(digest)

## moved into datatable_commons because replaced almost everywhere with dplyr
#require.auto(data.table)




########################################################################################################################
#### Convenience aliases


echo <- function(...) cat(paste(...), fill=T)

ac <- function(...) as.character(...)

# string concatenation without space gaps (could/should use paste0 instead)
## Deprecated: use paste0 instead
#concat <- function(...) paste(..., sep="")

unlen <- function(x) length(unique(x))

pp <- function(dat) page(dat, method = "print")

as.df <- function(dt) as.data.frame(dt)


########################################################################################################################
#### data.frame manipulation


shuffle <- function(df) df[sample(nrow(df)),]

first <- function(x, n=1) head(x,n)



#vec2df <- function(namedVec) namedVec %>% {data.frame(name=names(.), value=., row.names=NULL)}
vec2df <- function(namedVec) data.frame(name=names(namedVec), value=namedVec, row.names=NULL)


rownames2column <- function(df, colname){
  df <- as.df(df)
  df$tttt <- rownames(df);
  rownames(df) <- NULL;
  plyr::rename(df, c(tttt=colname))
}


column2rownames<- function(df, colname){
  #browser()
  ## force into df to avoid dplyr problems
  df <- as.df(df)
  
  rownames(df) <- ac(df[,colname])
  df[colname] <- NULL
  return(df)
}

## pushing some columns to the right end of a data.frame
push_right <- function(df, pushColNames){
  df[,c(setdiff(names(df),pushColNames), pushColNames)]
}


push_left <- function(df, pushColNames){
  df[,c(pushColNames, setdiff(names(df),pushColNames))]
}


#http://astrostatistics.psu.edu/datasets/R/html/base/html/formals.html
set_names <- function(df, ...){
  #browser()
  newnames <- as.character(unlist(list(...)))
  
  ## automatically convert matrices to data.frames (otherwise the names set would fail
  if(is.matrix(df)) df %<>% as.data.frame()
  
  names(df) <- newnames;
  return(df)
}
#iris %>% set_names(c("setosa", "hallo")) %>% head
#iris %>% set_names("setosa", "hallo") %>% head


head_html <- function(df, n=5) head(df, n) %>% knitr::kable(format="html") %>% print()

print_head <- function(df, desc=NULL){
  print(head(df))
  print(nrow(df))
  return(df)
}


fac2char <- function(mydata, convert=names(mydata)[sapply(mydata, is.factor)]){
  if(length(convert)==0)
    return(mydata)
  
  inputColOrder <- names(mydata)
  
  convertData <- subset(mydata, select= names(mydata)%in%convert)
  convertData <- as.data.frame(lapply(convertData, as.character), stringsAsFactors = FALSE)
  
  keepData <-  subset(mydata, select=!(names(mydata)%in%convert))
  newdata <- cbind(convertData, keepData)
  newdata <- newdata[,inputColOrder]
  
  return(newdata)
}

## replace R within pipe change just use ... %>% do(replaceNA(0)) %>% ...
replaceNA <- function(x, withValue) { x[is.na(x)] <- withValue; x }


## workaround for biomart
## Deprecated: load dplyr after biomart to avoid this problem
#dselect <- function(...) dplyr::select(...)


########################################################################################################################
#### Result Caching for long running tasks

## related: http://cran.r-project.org/web/packages/R.cache/R.cache.pdf

cache_it <- function(expr, filePrefix="cache"){
  cacheFile <- paste0(filePrefix, "_", substr(digest::digest(deparse(expr)), 1,6)) %>% paste0(".", ., ".RData")
  
  if(file.exists(cacheFile)){
    local(get(load(cacheFile)))
  } else {
    result <- eval(expr)
    save(result, file=cacheFile)
    result
  }
}

## Examples
#mydata <- quote(iris %>% filter(Species=="setosa")) %>% cache_it("tt")
#mydata <- quote(iris %>% filter(Species=="setosa")) %>% cache_it()
#mydata <- quote( { print("evaluate expr"); iris %>% filter(Species=="setosa") } ) %>% cache_it()

########################################################################################################################
#### File System

is.directory <- function(dirname) !is.na(file.info(dirname)$isdir)


mcdir <- function(dirname){
  if(!file.exists(dirname)){
    dir.create(dirname)
  }
  
  setwd(dirname)
}

locload <- function(fileName) local(get(load(fileName)))


rmerge <- function(LDF, by, ...){
  DF <- LDF[[1]]
  for (i in 2:length(LDF)) {
    DF <- merge(DF, LDF[[i]], by=by)
  }
  DF
}



## Deprecated: use trim_ext instead
trimEnd <- function(fileNames, ...) trim_ext(fileNames, ...)

trim_ext <-function(fileNames, exts=c()){
  for(fileExt in exts){
    fileNames <- str_replace(fileNames, paste(fileExt, "$",sep=""), "")
  }
  
  fileNames
}

write.delim <- function(df, file, header=TRUE,...){
  write.table(df, file, row.names=FALSE, col.names=header, sep="\t", ...)
}

rmSomeElements <- function(vec, toDel) vec[!(vec %in% toDel)]

rmLastElement <- function(vec) vec[-length(vec)]


########################################################################################################################
## Memory management


# improved list of objects
lsos <- function (pos = 1, pattern, order.by,
                  decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)/1000000
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  #        out
  out <- transform(out, var_name=rownames(out))
  rownames(out) <- NULL
  arrange(out, Size)
}

# shorthand that just shows top 1 results
lsosh <- function(..., n=10) {
  lsos(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

########################################################################################################################
### Statistics


## outlier handling
trim_outliers <- function(values, range=quantile(values, c(0.05, 0.95)))  pmax(range[1], pmin(range[2], values))

## use trim_outliers instead
#limit_range <- function(values, range)  pmax(range[1], pmin(range[2], values))

se<-function(x)	sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x)))
