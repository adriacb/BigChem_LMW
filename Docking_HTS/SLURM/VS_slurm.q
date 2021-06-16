#!/bin/bash

#SBATCH -D .
#SBATCH --output=log/rdock_2_filter.q.o
#SBATCH --error=log/rdock_2_filter.q.e
#SBATCH --job-name=rdock_2_filter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=602-700%30

module load rdock
#-------------------------------------------------------------------------
indir=/alexandria1/DB/BigChem_FRAGMENTS/SPLIT_DB/
outdir=VS_Results

num=$(( ${SLURM_ARRAY_TASK_ID} ))
curr_file=`ls $indir*.sd | awk -v line=$num '{if (NR == line) print $0}'`
file_name=`awk 'END{ var=FILENAME; n=split (var,a,/\//); print a[n]}' $curr_file`
#-------------------------------------------------------------------------

rbdock -i $curr_file -o $outdir/docked_filter_2_$file_name -r 4lr6.prm -p dock.prm -t VS_FILTER.txt > $outdir/log/rdock_2_${SLURM_ARRAY_TASK_ID}_filter.log


