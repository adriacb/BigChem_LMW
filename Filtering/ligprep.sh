#!/bin/bash

module load Corina
module load rdock
module load ChemAxon
module load openbabel3

input_smi=$1
output_sdf=$2
#Filtering salts with corina + protonation states + tautomers + 3D conformers
corina -i t=smiles -o t=sdf -d rs ${input_smi} | \
obabel -isd -osd --unique | \
cxcalc microspeciesdistribution -H 7 -M True | sdfilter -f'$DISTR[pH=7] > 5' | \
cxcalc dominanttautomerdistribution --protectcharge true | sdfilter -f'$TAUTOMER_DISTRIBUTION > 5' > 2D_${output_sdf}

corina -i t=sdf -o t=sdf -t n -d stergen,axchir,msi=4,preserve,wh,rc,de=8,mc=5,sc,errorfile=Failed.sdf 2D_${output_sdf} ${output_sdf}
