
rm(list = ls())

# Install packages -----

# check if renv is installed
if(!("renv" %in% installed.packages()[,1])){
  install.packages('renv')
} # maybe a session restart will be needed 

renv::activate()                
renv::restore(prompt = FALSE) 

system('Rscript --version') # !
source('00_system.R')

dir.create('RData')

# Application ------

# Pre-processing

renv::run('01_preprocessing.R')

# Dense scenario

renv::run('02_a_App_Dense_gridsearch.R')
renv::run('02_b_App_Dense.R')

# Sparse scenario

renv::run('03_a_App_Sparse_gridsearch.R')
renv::run('03_b_App_Sparse.R')

# Figures

dir.create('Figures')
renv::run('04_Figures.R')
