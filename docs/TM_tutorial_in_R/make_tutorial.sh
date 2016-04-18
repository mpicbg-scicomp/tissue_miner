#export TM_HOME=/projects/project-raphael/scripts/tissue_miner
export TM_HOME=/Volumes/projects/project-raphael/scripts/tissue_miner

## define where the tutorial data is defined
export TUTORIAL_DATA=/projects/project-raphael/tutorial/MovieRepository_DB

echo '
require("rmarkdown")

render(file.path(Sys.getenv("TM_HOME"), "docs/TM_tutorial_in_R/TM_R-UserManual.Rmd"), "html_document")
' | R -q --vanilla


## checkout the gh-pages branch of the repository and commit the changes
cd $TM_HOME && cd ..
git clone https://github.com/mpicbg-scicomp/tissue_miner.git tm_static_html

cd tm_static_html
git checkout --orphan gh-pages

#### just needed for initial run:
##git rm -rf .
##mkdir tm_tutorial
##cp $TM_HOME/docs/TM_tutorial_in_R/R-tutorial.html tm_tutorial/
##git add -A tm_tutorial
##git commit -m 'started pages branch for static tutorial files'

git commit -m 'fixed some typos'
git push origin gh-pages


## check updated version https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html