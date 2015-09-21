#install.packages(c("plyr", "igraph", "ggplot2", "reshape", "stringr", "zoo"))
#source("http://bioconductor.org/biocLite.R")
#biocLite("graph")

#source("http://dl.dropbox.com/u/113630701/rlibs/base-commons.R")

require.auto(graph)
require.auto(zoo)

require.auto(igraph)

#######################################################################################################################
### Build the graph and determine division groups
#######################################################################################################################

print("infering cell division groups...")

cellDivs <- subset(cellinfo,  !is.na(left_daughter_cell_id))
cellDivsLong <- melt(with(cellDivs, data.frame(cell_id, left_daughter_cell_id, right_daughter_cell_id)), id.vars="cell_id", value.name="daughter_id")
cellDivsLong$variable <- NULL

cytoGraphData <-  with(cellDivsLong, data.frame(from=cell_id, to=daughter_id))

## clean up NAs because of fake daughters
with(cytoGraphData, as.data.frame(table(is.na(from), is.na(to))))
cytoGraphData <- subset(cytoGraphData, !is.na(to))

#### DEBUG todo some cells are missing: seems to divide but divided IDs don't show up anywhere
#graphCells <- unique(with(cytoGraphData, c(from, to)))
#allCells <- unique(cellinfo$cell_id)
#
#missingGraphCells <- graphCells[!(graphCells %in% allCells)]
#missingGraphCells <- allCells[!(allCells %in% graphCells)]
#
#
### compile table for matthias
#fakeDaughter <- subset(firstOccOrDiv, right_daughter_cell_id %in% missingGraphCells)
#fakeDaughter$row_in_cells_table=as.numeric(rownames(fakeDaughter))
#
### do they show up elsewhere
#subset(cells, fakeDaughter$daughter_id %in% tissue_analyzer_group_id)
### --> No
#
#write.delim(fakeDaughter, file="fakeDaughters.txt")
## fakeDaughter <- read.delim("fakeDaughter.txt")
### compile table for matthias END
#
#### DEBUG end

allCells <- unique(cellinfo$cell_id)

## unlead as it overloads some methods from graph package
#unloadNamespace("igraph")

cytoGraph <- graph.data.frame(cytoGraphData, directed=TRUE, vertices=data.frame(cell_id=allCells))

#if(F){ #### DEBUG new parser does not process debug movie
#longGraph <- melt(cytoGraphData)
#
#longGraphNotInCI <-subset(longGraph, !(value %in% unique(cellinfo$cell_id)))
#
#subset(cellinfo, left_daughter_cell_id==10386 | right_daughter_cell_id==10386)
#subset(cellinfo, cell_id==10386)
#subset(cellinfo, cell_id==10385)
#subset(cellinfo, cell_id==10384)
#
#subset(cells2, tissue_analyzer_group_id==3043859 & frame==101)
#
### check track group of missing daughter_id
#subset(cells2, tissue_analyzer_group_id==16386193)
#
#} #### DEBUG end

## split subgraphs and assign them to each division event
#g <- ftM2graphNEL(as.matrix(cytoGraphData), V=ac(allCells), edgemode="directed")
#divisionTrees <-  connComp(g)$membership
#names(divisionTrees)<- paste("lineage_group", 1:length(divisionTrees), sep="_")
#divTreesDF<-ldply(divisionTrees, function(x) as.data.frame(unlist(x)), .progress="text")
#names(divTreesDF)  <- c("lineage_group", "cell_id")

## http://igraph.sourceforge.net/doc/R/clusters.html
cytoGraphClusters <- clusters(cytoGraph)
divTreesDF <- with(cytoGraphClusters, data.frame(cell_id=as.numeric(V(cytoGraph)$name), lineage_group=paste("lg", membership, sep="_")))

#degree(ug)
## Visualize the distribution
lineage_groupCounts <- arrange(with(divTreesDF, as.data.frame(table(lineage_group))), -Freq)
ggplot(lineage_groupCounts, aes(factor(Freq))) + geom_bar() +xlab("number of division events in division group (graph sizes)")
ggsave2()
#subset(divTreesDF, lineage_group=="lineage_group_524")

save(divTreesDF, file="divTreesDF.RData")
# divTreesDF <- local(get(load("divTreesDF.RData")))


########################################################################################################################
### Calculate the generation number for each cell
########################################################################################################################

print("calculating cell generations...")

# caculate node degrees
inDegrees <- degree(cytoGraph, mode="in")
inDegreesDF <- data.frame(cell_id=names(inDegrees), in_degree=inDegrees)


## visualize the in-degree distribution (so which cell come from a division event)
ggplot(inDegreesDF, aes(factor(in_degree))) + geom_bar() + ggtitle("node in-degree distribution")
ggsave2()

## combine with division events table
cytoWithGrpDeg <- merge(divTreesDF, inDegreesDF, by="cell_id", all=T)


### count roots per division group
#rootCounts <- ddply(cytoWithGrpDeg, .(lineage_group), function(curSubTree){
#    # DEBUG curSubTree <- subset(cytoWithGrpDeg, lineage_group=="lineage_group_37")
#    c(root_count=nrow(subset(curSubTree, in_degree==0)))
#}, .progress="text")

rootCounts <- as.df(group_by(cytoWithGrpDeg, lineage_group) %>% filter(in_degree==0) %>% summarise(root_count=n()))
with(rootCounts, as.data.frame(table(root_count)))

ggsave2(ggplot(rootCounts, aes(as.factor(root_count))) + geom_bar() + ggtitle("root count distribution"))

## fix missing ids for no-further dividing cells
cytoWithGrpDeg  <- transform(cytoWithGrpDeg, in_degree=ifelse(is.na(in_degree), 0, in_degree))

procLG <<-1
numLGs <<-unlen(cytoWithGrpDeg$lineage_group)

### old implmementation
#if(F){
#
#calcGenerationForSubtree <- function(curSubTree){
#    if(procLG%%100==0) praste("percont done ", procLG/numLGs)
#    procLG <<- procLG+1
#
#    #    DEBUG curSubTree <- subset(cytoWithGrpDeg, lineage_group=="fakegroup_6")
#    if(nrow(curSubTree)==1) return(transform(curSubTree, generation=1)) ##  just speedup for single node lineage_groups
#
#    #if(nrow(subset(curSubTree, in_degree==0))>1) {stop(paste("division tree with multiple roots as in ", curSubTree$lineage_group[1], "should not happen")) }
#    ## todo remove this ugly hack by reenabling the original stop criterion
#    if(nrow(subset(curSubTree, in_degree==0))>1) return(merge(curSubTree, data.frame(cell_id=unique(curSubTree$cell_id), generation=NA),  all.x=T))
#
#    divtreeRoot <- ac(subset(curSubTree, in_degree==0)$cell_id)
#    nonRoots <- ac(curSubTree$cell_id)
##    divtreeRoot %in% V(cytoGraph)$name
##    nonRoots %in% V(cytoGraph)$name
#
#    short_paths <- get.shortest.paths(cytoGraph, from=divtreeRoot, to=nonRoots, mode="out")
##    browser()
#
#    names(short_paths) <- nonRoots
#    generations<-ldply(short_paths, length)
#    names(generations)  <- c("cell_id", "generation")
#
##     merge(curSubTree, generations, all.x=T)
#    merge(curSubTree, generations, all.x=T)
#}
#cytoWithGeneration <- ddply(cytoWithGrpDeg, .(lineage_group), calcGenerationForSubtree, .progress="text", .parallel=F)
##cytoWithGeneration <- ddply(head(cytoWithGrpDeg, 10000), .(lineage_group), calcGenerationForSubtree, .progress="text", .parallel=F)
#
#}

calcGenerationForSubtreeNEW <- function(cell_id_ST, in_degree_ST){
#    echo("processing", paste(cell_id_ST, collapse=","))
    if(procLG%%100==0) echo("percont done ", procLG/numLGs)
    procLG <<- procLG+1

    if(length(cell_id_ST)==1) return(0) ##  just speedup for single node lineage_groups

    #if(nrow(subset(curSubTree, in_degree==0))>1) {stop(paste("division tree with multiple roots as in ", curSubTree$lineage_group[1], "should not happen")) }
    ## todo remove this ugly hack by reenabling the original stop criterion
    if(sum(in_degree_ST>1)) return(as.numeric(NA))

    divtreeRoot <- ac(cell_id_ST[in_degree_ST==0])
    nonRoots <- ac(cell_id_ST)
#    divtreeRoot %in% V(cytoGraph)$name
#    nonRoots %in% V(cytoGraph)$name

    short_paths <- get.shortest.paths(cytoGraph, from=divtreeRoot, to=nonRoots, mode="out")$vpath

#browser()
    names(short_paths) <- nonRoots
    generations<-ldply(short_paths, length)
    names(generations)  <- c("cell_id", "generation")

    # todo necessary?
    genSorted <- generations$generation[with(generations, match(cell_id_ST, cell_id))]
#    print(genSorted)

    ## remove offset because from a geneticist point of view, cells that never divide correspond to the F0 parental line.
    genSorted <- genSorted - 1;

    return(as.numeric(genSorted))
}

cytoWithGeneration <- group_by(cytoWithGrpDeg, lineage_group) %>%
#    filter(n()==5) %>%
    mutate(generation=calcGenerationForSubtreeNEW(cell_id, in_degree)) %>%
    as.df()

#cytoWithGeneration <- ddply(fac2char(cytoWithGrpDeg), .(lineage_group), transform, generation=calcGenerationForSubtreeNEW(cell_id, in_degree), .progress="text", .parallel=T)


# TO SKIP GENERATION ANALYSIS: cytoWithGeneration <- cytoWithGrpDeg; cytoWithGeneration$generation <- -1

ggsave2(ggplot(cytoWithGeneration, aes(as.factor(generation))) + geom_bar())
cytoWithGeneration$in_degree <- NULL

save(cytoWithGeneration, file="cytoWithGeneration.RData")
# cytoWithGeneration <- local(get(load("cytoWithGeneration.RData")))


########################################################################################################################
### Graph Visualization
########################################################################################################################

if(F){ #### DEBUG

#library(igraph)
#cytoGraph <- graph.data.frame(cytoGraphData, directed=TRUE)
#plot(cytoGraph, vertex.label="", layout=layout.kamada.kawai(clusterGraph))
##plot(clusterGraph, vertex.label=V(clusterGraph)$name, edge.color=factor(get.edge.attribute(clusterGraph, "seqrev_weight")),
## plot(clusterGraph, vertex.label=V(clusterGraph)$name, edge.color=factor(get.edge.attribute(clusterGraph, "seqrev_weight")), layout=layout.kamada.kawai(clusterGraph))

plot(cytoGraph, vertex.label="", layout=layout.kamada.kawai(cytoGraph))


write.csv(with(subset(cytoWithGeneration, !is.na(rooted_parent_cell)), data.frame(Id=rooted_parent_cell,Source=rooted_parent_cell, generation, lineage_group )), file="cytoWithGenerationNodes.csv", row.names=F)
write.csv(with(subset(cytoWithGeneration, !is.na(rooted_parent_cell)), data.frame(Source=rooted_parent_cell, Target=left_progenitor)), file="cytoWithGenerationEdges.csv", row.names=F)

} #### DEBUG end
