#!/bin/bash
#####################################
##  Amber MDmix simulation
##  Project Variables
##  Adrià Cabello Blanque
####################################


# load project variables
source projectvars.sh


#Set the right environment variables:

module load amber12 python27 mdmix2

# Run simulation
for d in $MD/*; do cd $d/min; qsub min.q; cd -;done

# Run analysis

sbatch analisi.q
