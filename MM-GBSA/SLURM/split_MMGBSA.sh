#!/bin/bash

module load rdock

split=$1
infile=$2


if [ $split ]  && [ $infile ]; then
	echo "File "$2 "splitted in $split parts.";
else 
	echo "split_MMGBSA.sh [number of splits] [sdf input file]";
	exit 0;
fi


cmd="sdsplit -${split} ${infile}"

if module load rdock; then
echo $cmd
$cmd
else
	echo "Cannot find rDock module";
	exit 0;
fi 




for i in $(eval echo {1..$split});
do
if [[ -f "tmp${i}.sd" ]]; then
	mkdir tmp_${i};
	mv tmp${i}.sd ./tmp_${i};
	
else
	break
fi
done;


