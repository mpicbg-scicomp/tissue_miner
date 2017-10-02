#!/usr/bin/env Rscript --no-environ


argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: LineageGroupColoring.R  <movie_db_directory>")
}else{
    movieDir=normalizePath(argv[1])
    if(is.na(file.info(movieDir)$isdir)) stop(paste("movie directory does not exist"))
}


########################################################################################################################
### Setup environment

db_name=basename(movieDir)

scriptsDir=Sys.getenv("TM_HOME")

if(is.na(file.info(scriptsDir)$isdir)){
    stop(paste("TM_HOME  not correctly defined (",scriptsDir ,")"))
}

source(file.path(scriptsDir, "commons/TMCommons.R"))

db <- openMovieDb(movieDir)

mcdir(file.path(movieDir, "lg_color_optimization"))

# source("http://bioconductor.org/biocLite.R"); biocLite("RBGL")
sink(file=file("/dev/null", "w"), type="message")
  require.auto(RBGL)
  require.auto(graph)
sink(file=NULL, type="message")

########################################################################################################################
### Optimize colors for division groups
########################################################################################################################
echo('Optimize colors for division groups...')
dbonds <- dbGetQuery(db, "select cell_id, dbond_id, conj_dbond_id from directed_bonds")
#dbondsSlim <- with(dbonds, data.frame(frame, cell_id=cell_id, dbond_id))

## add division group
lgGroupInfo <- dbGetQuery(db, "select cell_id, lineage_group from cell_histories")
#lgGroupInfo <- dbGetQuery(db, "select * from cell_histories")

## combine bonds with lineage group info
dbondsWithLinInfo <- dt.merge(dbonds, lgGroupInfo, by="cell_id")

## remove background cell
dbondsWithLinInfo <- subset(dbondsWithLinInfo, cell_id != 10000)
neighbors <- dt.merge(dbondsWithLinInfo, plyr::rename(subset(dbondsWithLinInfo, select=-dbond_id), c(conj_dbond_id="dbond_id")), by=c("dbond_id"))

## reduce group interfaces per frame
lgInterfaces <- with(neighbors, data.frame(lineage_group.x, lineage_group.y)) %>% subset(!duplicated(paste(lineage_group.x, lineage_group.y)))

lgInterfcesNoDup <- subset(lgInterfaces, ac(lineage_group.x) > ac(lineage_group.y))


grpGraph <- ftM2graphNEL(as.matrix(lgInterfcesNoDup), edgemode="undirected")
graphColors <- sequential.vertex.coloring(grpGraph)
graphColorsDF <- data.frame(lineage_group=names(graphColors[[2]]), color=as.numeric(graphColors[[2]]))

gg <- ggplot(graphColorsDF, aes(as.factor(color))) + geom_bar() + ggtitle("color counts") + xlab("color") + ggtitle("lineage group coloring distribution")
ggsave2(gg)

## save for later (to be plugged in as color scheme
write.delim(graphColorsDF, file="lg_colors.txt")
# graphColorsDF <- read.delim("lg_colors.txt")

## gephi export
write.csv(data.frame(Id=names(graphColors[[2]]),Source=names(graphColors[[2]]), color=as.numeric(graphColors[[2]])), file=paste0(db_name, ".condesed_lineage_groups.gephi.nodes.csv"), row.names=F)
write.csv(with(lgInterfcesNoDup, data.frame(Source=lineage_group.x, Target=lineage_group.y)), file=paste0(db_name, ".condesed_lineage_groups.gephi.edges.csv"), row.names=F)

