#!/bin/bash
#SBATCH --job-name=EXCLUDE
#SBATCH -D .
#SBATCH --time=01:00:00
#SBATCH --output=test_cluster.o
#SBATCH --error=test_cluster.q.e
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5976

WD='/alexandria1/mminarro/BIGCHEM_FRAGMENTS/DOCKED/WHOLE_DOCKED_SET_EXCLUDE_H'
WD2='/alexandria1/mminarro/BIGCHEM_FRAGMENTS/DOCKED'

#perl $WD/filter_dist_exclude.pl area_exclude.txt $WD2/output/docked_${SLURM_ARRAY_TASK_ID}.sd > $WD/output/H_filter_${SLURM_
ARRAY_TASK_ID}.sd
perl $WD/filter_Hyd.pl area_exclude.txt $WD/output/H_filter_${SLURM_ARRAY_TASK_ID}.sd > $WD/1A_TOLERANCE/1A_TOLERANCE__${SLUR
M_ARRAY_TASK_ID}.sd
