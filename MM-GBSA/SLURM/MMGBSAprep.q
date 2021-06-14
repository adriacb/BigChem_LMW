#!/bin/bash

#SBATCH -e error.log
#SBATCH -o output.log
#SBATCH -J "MMGBSA_BigChem"
#SBATCH -p std
#SBATCH --time=72:00:00

bash prepare_MMGBSA.sh 40

