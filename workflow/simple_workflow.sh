#/bin/sh

## An example workflow to analyze a cell-tracked movie

## First follow setup instructions on https://github.com/mpicbg-scicomp/tissue_miner/blob/master/README.md#how-to-run-locally

movieDbDir=.

## build the db
CreateMovieInfoDat.sh $movieDbDir
imageParser $(dirname {$movieDbDir}) $(basename {$movieDbDir}) %03d
CreateDbFromParser.R $movieDbDir

## roi tracking
UserRoiTracking.R $movieDbDir
RoiDeformation.R $movieDbDir

## analyze lineage
LineageGroupColoring.R $movieDbDir
LineageMovies.R $movieDbDir

## analyze cell topoloy
CountT1.R $movieDbDir
TopologyMovies.R $movieDbDir
UnbalanceT1Movie.R $movieDbDir
FourWayVertices.R $movieDbDir
PolygonClass.R $movieDbDir

## analyze shear contributions
CreateTriangles.R $movieDbDir
CategorizeTriangles.R $movieDbDir

## misc
AreaMovies.R $movieDbDir
DbElongationMovie.R $movieDbDir
DensityMovies.R $movieDbDir
StripeMovies.R $movieDbDir
