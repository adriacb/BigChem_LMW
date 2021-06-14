#!/bin/bash

#SBATCH -e error.log
#SBATCH -o output.log
#SBATCH -J "MMGBSA_BigChem"
#SBATCH -p std
#SBATCH --time=72:00:00


bash split_MMGBSA.sh 50 ligands_docking_H_2k.sdf
