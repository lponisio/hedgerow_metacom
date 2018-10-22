
## run from directory with git repository
mkdir -p  ../hedgerow_metacom_saved

## ************************************************************
## data prep, need access to original data, github should have all
## necessary files needed without this unless you are the authors
## ************************************************************
## 
Rscript -e 'install.packages("bipartite", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("fossil", repos="http://cran.r-project.org")'

## spatial data prep
Rscript -e 'install.packages("rgdal", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("maptools", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("sp", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("raster", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("spectralGP", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("dismo", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("viridis", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("rgeos", repos="http://cran.r-project.org")'


## ************************************************************
## packages needed for analysis
## ************************************************************
## spatial packages
Rscript -e 'install.packages("maptools", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("sp", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("spatstat", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("rgdal", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("raster", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("rgeos", repos="http://cran.r-project.org")'

## Bayesian analyses
Rscript -e 'install.packages("vegan", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("igraph", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("nimble", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("abind", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("RColorBrewer", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("viridis", repos="http://cran.r-project.org")'

## frequentist analyses
Rscript -e 'install.packages("lmerTest", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("lme4", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("bipartite", repos="http://cran.r-project.org")'
