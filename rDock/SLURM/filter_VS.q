#!/bin/bash
​
#SBATCH --job-name=test_cluster
#SBATCH -D .
#SBATCH --time=01:00:00
#SBATCH --output=log/test_cluster-%a.q.o
#SBATCH --error=log/test_cluster-%a.q.e
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-196
​
module load rdock
​
​
OUTFILE='FILTERED_SCORE_INTER_LAST.sd'
​
for file in $PWD/VS_Results/*.sd
do
#FILTER DIFERENT MOLECULES
$RBT_ROOT/bin/sdsort -n -fSCORE -s $file |$RBT_ROOT/bin/sdfilter -f'$SCORE.RESTR < 0.7' |$RBT_ROOT/bin/sdfilter -f'$SCORE.INTER < -12.0'|$RBT_ROOT/bin/sdfilter -f'$_COUNT == 1'>> tmp.sd
done
​
#ORDER MOLECULES BY SCORE INTER
$RBT_ROOT/bin/sdsort -n -fSCORE.INTER tmp.sd |$RBT_ROOT/bin/sdfilter -f'$_COUNT == 1' >> $OUTFILE
​
#echo Compressing \${sdout}...
#rm tmp.sd
#gzip -9vf \${sdout}
