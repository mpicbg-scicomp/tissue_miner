#/bin/sh

## An example workflow to analyze a cell-tracked movie


movieDbDir=.

/sw/bin/xvfb-run imageParser $(dirname {$movieDbDir}) $(basename {$movieDbDir}) %03d
CreateDbFromParser.R $movieDbDir

LastFrameRoiBT.R $movieDbDir
RoiDeformation.R $movieDbDir

LineageGroupColoring.R $movieDbDir
LineageMovies.R $movieDbDir

CountT1.R $movieDbDir
TopologyMovies.R $movieDbDir
UnbalanceT1Movie.R $movieDbDir
FourWayVertices.R $movieDbDir
PolygonClass.R $movieDbDir

CreateTriangles.R $movieDbDir
CategorizeTriangles.R $movieDbDir

AreaMovies.R $movieDbDir
DbElongationMovie.R $movieDbDir
MakeMovies.R $movieDbDir
StripeMovies.R $movieDbDir
