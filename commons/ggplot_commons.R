require.auto(ggplot2)
require.auto(scales)
require.auto(grid)


scale_fill_redgreed <- function() scale_fill_manual(values = c("red","darkgreen"))

rot_x_lab <- function() theme(axis.text.x = element_text(angle = 90, hjust = 1))
rot_x_45 <- function() theme(axis.text.x = element_text(angle = 45, hjust = 1))
## DEPRACTED because of naming
rotXlab <- function() theme(axis.text.x = element_text(angle = 90, hjust = 1))



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, by.row=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


gg2Format="png"


## simplified save function for ggpltos
ggsave2 <- function(gplot=last_plot(), width=8, height=6, prefix="", saveData=FALSE, outputFormat=gg2Format, ...){
  title <- try(gplot$labels[["title"]])
  
  if(is.null(title)){
    varMapping <- gplot$labels
    varMapping <- varMapping[names(varMapping) %in% c("x", "y")]
    #        if(varMapping==NULL){
    #            varMapping <- gplot$mapping;
    #        }
    
    if(length(varMapping) == 1){
      title= paste("distribution of", varMapping)
    }else{
      title = try(paste(varMapping, collapse=" vs "))
      # stop("last plot had no title. Use ggsave() and give it a manual title")
    }
    
    rawFacetDesc <- format(gplot$facet)
    if(rawFacetDesc!="facet_null()"){
      title <- paste(title, "by", str_replace_all(str_match(rawFacetDesc, "facet_.*[(](.*)[)]")[,2], "~", "and"))
    }
  }
  
  
  fileBaseName <- ifelse(nchar(prefix)>0, paste0(prefix, " - ", title), title)
  
  ## clean up weired characters
  fileBaseName <- str_replace_all(fileBaseName, "[$%/?]", "_")
  
  fileName = paste0(fileBaseName, paste0(".", outputFormat))
  
  ## remove line-breaks and trim spaces
  fileName = str_replace_all(str_replace_all(fileName, "\\n", ""), "[ ]{2,}", " ")
  #    print(paste("saving plot to ", fileName))
  ggsave(fileName, width=width, height=height, ...)
  
  if(saveData){
    write.delim(gplot$data, file= paste0(fileBaseName, ".txt"))
  }
  
  return(fileName)
}


########################################################################################################################
### pca plots (http://largedata.blogspot.de/2011/07/plotting-pca-results-in-ggplot2.html)

makePcaPlot <- function(x = getData(), group = NA, items=rownames(x), title = "") {
  require(ggplot2)
  require(RColorBrewer)
  
  #  data <- x
  #  data <- t(apply(data, 1, scale))
  #  rownames(data) <- rownames(x)
  #  colnames(data) <- colnames(x)
  #  mydata <- t(data)
  mydata <-x
  mydata.pca <- prcomp(mydata, retx=TRUE, center=TRUE, scale.=TRUE)
  
  percent <- round((((mydata.pca$sdev)^2 / sum(mydata.pca$sdev^2))*100)[1:2])
  
  scores <- mydata.pca$x
  pc12 <- data.frame(PCA1=scores[,1], PCA2=scores[,2], group=group)
  
  #  ggplot(pc12, aes(PCA1, PCA2, colour = group)) + geom_point(size = 6, alpha = 3/4)
  ggplot(pc12, aes(PCA1, PCA2, colour = group, label=items)) +
    geom_point(size = 6, alpha = 3/4)   +
    geom_text(size = 6, alpha = 3/4)    +
    xlab(paste("PCA1 (", percent[2], "%)", sep = "")) +
    ylab(paste("PCA2 (", percent[2], "%)", sep = ""))
  
  qplot(PCA2, PCA1, geom="blank", main = title, xlab = paste("PCA2 (", percent[2], "%)", sep = ""), ylab = paste("PCA1 (", percent[1], "%)", sep = "")) +
    geom_point(aes(colour = group), size = 6, alpha = 3/4)
  #    theme(
  #      axis.text.x = element_text(size = base_size * 1.3 , lineheight = 0.9, colour = "grey50", hjust = 1, angle = 90),
  #      axis.text.y = element_text(size = base_size * 1.3, lineheight = 0.9, colour = "grey50", hjust = 1)
  #    )
}


## example
# makePcaPlot(getData(30,4,2,distort = 0.7))


########################################################################################################################
### Base-plot utils


plotPDF <- function(fileBaseName, expr, ...){ pdf(paste0(fileBaseName, ".pdf"), ...); expr; dev.off(); }
#plotPDF("test", plot(1:10))


## create a custom color palette for a fixed set of values
## scale_fill_manual(values = create_palette(unique(csWithTopoT1$t1_type)), drop = FALSE)
create_palette <- function(x, pal = 'Set1'){
  require.auto(RColorBrewer)
  
  ux <- sort(unique(x))
  n <-length(ux)
  
  if(n==0) return(c())
  
  setNames(brewer.pal(name = pal, n = n)[1:n], ux)
}

