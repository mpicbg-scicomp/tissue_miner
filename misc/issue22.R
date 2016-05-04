
"
export TM_HOME=/projects/project-raphael/scripts/tissue_miner
export PATH=$TM_HOME/db:$TM_HOME/shear:$TM_HOME/roi:$TM_HOME/misc:$TM_HOME/movies:$TM_HOME/shear_contributions:$TM_HOME/topology:$TM_HOME/triangles:$TM_HOME/lineage:$PATH
#CreateDbFromParser.R .

## buggy one
cd /projects/project-raphael/rapha_bug22



## wt for comparison
cd movie_dbs/db_tests/tissue_analyzer/WT_25deg_111102_ForTissueMiner/
"


movieDir="/Users/brandl/Desktop/WT_25deg_111102_ForTissueMiner"
db_name=basename(movieDir)

Sys.setenv(TM_HOME="/Volumes/projects/project-raphael/scripts/tissue_miner")
#Sys.setenv(TM_HOME="/projects/project-raphael/scripts/tissue_miner")
scriptsDir=Sys.getenv("TM_HOME")

source(file.path(scriptsDir, "commons/TMCommons.R"))

movieDb <- openMovieDb(movieDir)

cellinfo <- dbGetQuery(movieDb, "select * from cell_histories")