## Make plot using the color scheme defined in misc/WingAlgnNcolorScheme.R

# mplot <- function(data,..., small=F){
#   algnData <- transform(merge(data, algnModel, by="movie"), dev_time=(time_sec+time_shift+54000)/3600) ## 54000 sec offset, namely 15hAPF
#   gg <- ggplot(algnData,...) + scale_color_manual(values=cols)+ facet_wrap(~roi)
#   if (small){
#     gg <- gg + scale_x_continuous(breaks=seq(16,40, 4),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")
#   }
#   else {gg <- gg + scale_x_continuous(breaks=seq(16,40, 2),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")}
#   return(gg)
# }
# 
# mplotcumsum <- function(data, ..., small=F){
#   algnData <- ddply(data, .(movie), transform, dev_time=(time_sec-min(time_sec)+refTime+54000)/3600) ## 54000 sec offset, namely 15hAPF
#   gg <- ggplot(algnData,...) + scale_color_manual(values=cols)+ facet_wrap(~roi) + xlab("time [hAPF]")
#   if (small){
#     gg <- gg + scale_x_continuous(breaks=seq(16,40, 4),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")
#   }
#   else {gg <- gg + scale_x_continuous(breaks=seq(16,40, 2),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")}
#   return(gg)
# }



### Prefix often used fuctions to get easily get them by autocompletion
mpf_plot <- function(data,..., small=F){
  algnData <- data #transform(merge(data, algnModel, by="movie"), dev_time=(time_sec+time_shift+54000)/3600) ## 54000 sec offset, namely 15hAPF
  gg <- ggplot(algnData,...) + scale_color_manual(values=cols)
  if (small){
    gg <- gg + scale_x_continuous(breaks=seq(16,40, 4),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")
  }
  else {gg <- gg + scale_x_continuous(breaks=seq(16,40, 2),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")}
  return(gg)
}

mpf_plotFacetByRoi <- function(data,..., small=F){
  algnData <- data #transform(merge(data, algnModel, by="movie"), dev_time=(time_sec+time_shift+54000)/3600) ## 54000 sec offset, namely 15hAPF
  gg <- ggplot(algnData,...) + scale_color_manual(values=cols)+ facet_wrap(~roi)
  if (small){
    gg <- gg + scale_x_continuous(breaks=seq(16,40, 4),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")
  }
  else {gg <- gg + scale_x_continuous(breaks=seq(16,40, 2),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")}
  return(gg)
}

mpf_plotFacetByMovie <- function(data,..., small=F){
  algnData <- data #transform(merge(data, algnModel, by="movie"), dev_time=(time_sec+time_shift+54000)/3600) ## 54000 sec offset, namely 15hAPF
  gg <- ggplot(algnData,...) + scale_color_manual(values=cols)+ facet_wrap(~movie)
  if (small){
    gg <- gg + scale_x_continuous(breaks=seq(16,40, 4),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")
  }
  else {gg <- gg + scale_x_continuous(breaks=seq(16,40, 2),limits=c(15,ceiling(max(algnData$dev_time)))) +xlab("time [hAPF]")}
  return(gg)
}

