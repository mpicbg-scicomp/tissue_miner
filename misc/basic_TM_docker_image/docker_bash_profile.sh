## Define where the TissueMiner directory is located
export TM_HOME="/tissue_miner/"

export PATH=$TM_HOME/db:$TM_HOME/shear:$TM_HOME/roi:$TM_HOME/misc:$TM_HOME/movies:$TM_HOME/shear_contributions:$TM_HOME/topology:$TM_HOME/triangles:$TM_HOME/lineage:$TM_HOME/docs/quickstart/scripts:$PATH

## add quickstart R scripts to the path
export PATH=${TM_HOME}/docs/quickstart/scripts:$PATH

## Add gnu-parallel ot path (disabled since not strictly necessary
#export PATH=/parallel-20120522/src:$PATH

## define a small snakemake wrapper to ease usage
function sm() {
    snakemake --snakefile ${TM_HOME}/workflow/tm.snkmk --keep-going "$@"
    # use -j3 to run 3 rules in parallel
#    snakemake --snakefile ${TM_HOME}/workflow/tm.snkmk --keep-going -j3 "$@"
}
export -f sm


## add the image parser to the path
export PATH=${TM_HOME}/parser:$PATH

## start virtual x-server ## no longer needed
#Xvfb :1 -screen 0 1024x768x16 &> xvfb.log  &
#export DISPLAY=:1.0

#cd /movies
