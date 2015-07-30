#!/usr/bin/env Rscript

argv = commandArgs(TRUE)
if(length(argv) != 1){
    stop("Usage: FourWayVertices.R  <movie_db_directory>")
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

mcdir(file.path(movieDir, "4way_vertices"))

########################################################################################################################


multColors =c("2" = "yellow", "3"="lightblue", "4"="red", "5"="darkgreen", "6"="darkblue");
#multColors = create_palette(1:4)


dbonds <- dbGetQuery(db, "select frame,  frame from dbonds")


dbondVertices <- dbGetQuery(db, "select vertex_id from dbonds")
vertexMultiplicities <- with(dbondVertices, as.data.frame(table(vertex_id), responseName="multiplicity"))

#subset(vertexMultiplicities, multiplicity==6)

ggsave2(ggplot(vertexMultiplicities, aes(factor(multiplicity))) + geom_bar()+ ggtitle("vertex_multiplicity"))


## visualize them
vertices <- dbGetQuery(db, "select * from vertices")


## visualize proportions over time
framesPerVertex <- with(vertices, as.data.frame(table(frame), responseName="num_frame_vertices"))
vertexMultWithByFrame <- with(dt.merge(vertexMultiplicities, with(vertices, data.frame(vertex_id, frame))), as.data.frame(table(frame, multiplicity)))
relVertexMult <- transform(dt.merge(vertexMultWithByFrame, framesPerVertex), mult_proportion=Freq/num_frame_vertices)
gg <- ggplot(relVertexMult, aes(as.numeric(frame), mult_proportion, color=multiplicity)) + geom_line() + xlab("frame") + ggtitle("vertex_multiplicity_proportions_by_frame") + scale_color_manual(values=multColors)
ggsave2()


## filter for the bad guys
## todo also try <3 to see where the class 2 vertices are sitting
#twVertices <- subset(vertices, vertex_id %in% subset(vertexMultiplicities, Freq<3)$vertex_id)



fwVertices <- subset(vertices, vertex_id %in% subset(vertexMultiplicities, multiplicity>3)$vertex_id)
ggsave2(ggplot(fwVertices, aes(frame)) + geom_bar(binwidth=1)+ ggtitle("four_way_distribution_over_time"))


#frameOI=100
#render_source_image(frameOI) + geom_point(data=subset(fwVertices, frame==frameOI), aes(x_pos, y_pos),color="red", alpha=0.4, size=10)

#render_movie(fwVertices, "four_way_vertices.mp4", geom_point(aes(x_pos, y_pos),color="red", alpha=0.4, size=10))
#render_movie(twVertices, "two_way_vertices.mp4", geom_point(aes(x_pos, y_pos),color="red", alpha=0.4, size=10))

## do one movie for both types
not3Vertices <- subset(vertices, vertex_id %in% subset(vertexMultiplicities, multiplicity!=3)$vertex_id)
not3Vertices <- merge(not3Vertices, vertexMultiplicities)
not3Vertices <- transform(not3Vertices, multiplicity=as.factor(multiplicity))

#with(not3Vertices, as.data.frame(table(multiplicity)))
#multColors <- create_palette(not3Vertices$multiplicity)

render_movie(not3Vertices, "not3_vertices.mp4", list(
    geom_point(aes(x_pos, y_pos, color=multiplicity), alpha=0.4, size=10),
    scale_color_manual(values=multColors, drop = FALSE)
)) #sampleRate=50

#if(F){ #### DEBUG
#frameOI=3
#render_source_image(frameOI) + geom_point(data=subset(not3Vertices, frame==frameOI), aes(x_pos, y_pos, color=as.factor(multiplicity)), alpha=0.4, size=10) +   scale_color_manual(name="multiplicity", values=multColors, drop = FALSE)
#} #### DEBUG end

