
rm(list = ls())

# Install packages -----

# check if renv is installed
if(!("renv" %in% installed.packages()[,1])){
  install.packages('renv')
} # maybe a session restart will be needed 

renv::rebuild()
system('Rscript --version') # !
source('00_system.R')

dir.create('RData')

# Confidence intervals simulation study -----

renv::run('02_IC_Simulation.R')

# LRT simulation study -----

renv::run('03_LRT_Simulation.R')

# Application ------

# Pre-processing

renv::run('04_a_preprocessing.R')

# Dense scenario

renv::run('04_b_App_Dense.R')

# Sparse scenario

renv::run('04_b_App_Sparse.R')

# Figures

dir.create('Figures')
renv::run('05_Figures.R')
