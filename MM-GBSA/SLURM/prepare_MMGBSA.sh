#!/bin/bash

module load schrodinger/2018-2

split=$1

#$SCHRODINGER/utilities/structconvert -imol2 receptor.mol2 -omae receptor.mae

for i in $(eval echo {1..$split});
do
if [[ -f "./tmp_${i}/tmp${i}.sd" ]]; then
	
	$SCHRODINGER/utilities/structconvert -isd ./tmp_${i}/tmp${i}.sd -omae ./tmp_${i}/ligands_${i}.mae;
	
	cat receptor.mae ./tmp_${i}/ligands_${i}.mae > ./tmp_${i}/input.mae
	
else
	break
fi
done;


