#!/bin/bash

#SBATCH -D .
#SBATCH --output=log/mmgbsa.q.o
#SBATCH --error=log/mmgbsa.q.e
#SBATCH --job-name=mmgbsa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00

 
#-------------------------------------------------------------------------
# num=$(( ${SLURM_ARRAY_TASK_ID} ))

#-------------------------------------------------------------------------

initdir=$PWD

module load schrodinger/2018-2
module load tools/java/1.8.0_181
export LMUTIL=$SCHRODINGER/mmshare-v*/bin/Linux-x86*/lmutil
export prod_schrondinger=/prod/apps/schrodinger/2018-2

#-------------------------------------------------------------------------

#### check the license

lic_check=`$LMUTIL lmstat -a -c 1715@sdf1.cesca.cat | grep "Users of SUITE" | awk '{print $11}'`
echo "mm refs ${lic_check}" >> ${num}_mmgbsa_server.log
while [ $lic_check -gt "12" ]; do
	lic_check=`$LMUTIL lmstat -a -c 1715@sdf1.cesca.cat | grep "Users of SUITE" | awk '{print $11}'`
	wait
done;

#### launch mmgbsa

#cd ./tmp_${num}/

$prod_schrondinger/prime_mmgbsa -out_type LIGAND -job_type REAL_MIN -OVERWRITE -jobname mmgbsa_output -WAIT -LOCAL input.mae

#cd ..

#echo $initdir/tmp_${num}/ >> files.txt


