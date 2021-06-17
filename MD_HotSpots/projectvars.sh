#!/bin/bash
#####################################
##  Amber MDmix simulation
##  Project Variables
##  Adrià Cabello Blanque
####################################


export WDR=$PWD; # IMPORTANT: If you provide the absolute path, make sure
                 #            that your path DOES NOT contains white-spaces
                 #            otherwise, you will get weird execution errors.
                 #            If you cannot fix the dir names containing such white-space
                 #            chars, you MUST set this var using the current folder '.'
                 #            instead of '$PWD', i.e:    export WDR=.;

mkdir md tleap

export MD=$WDR/md;
export TLEAP=$WDR/tleap;
export QueueName='SLURM'                     # Specify the name of the Queue system you are using (SGE to SLURM)
export PDB="your_file.pdb"                   # Specify PDB file name
export PName="project_name"                  # Name it as the struct you want to load the PDB



# 1. Prepare a PDB file with your sistem to simulate (we used MOE)

# 2.1 Create Object File Format
touch $TLEAP/leaprc.ff14SB

cat <<EOF > $TLEAP/leaprc.ff14SB
$PName = loadpdb $PDB
check $PName
saveoff $PName ${PName}.off
quit
EOF

export FORCEFIELD=$TLEAP/leaprc.ff14SB


# 2.2 Run Object file format

tleap -f $FORCEFIELD # or $AMBERHOME/exe/tleap -f leaprc.ff99SB if the program is not found in your PATH environmental variable

# > check $PName
# > saveOff $PName.off
# > quit



# 3. Create config file (CFG) Manually edit file “PTP_$PName.cfg”
cat <<EOF > $PWD/PTP_${PName}.cfg
# pyMDMix MD Settings Configuration File (MSCF)
# All non used options should be commented/removed
[SYSTEM]
# Name to identify the system to be loaded
NAME = $PName
# One of these two options should be given:
OFF = ${PName}.off
# PDB file
PDB = ${PName}.pdb

# Unit name containing the system inside the object file (default: first unit found)
UNAME = $PName

# Comma separated list of non-standard residues that should be considered as solute
# Used in automask detection and solute-solvent idetnfication (default: empty)
#EXTRARES =

# Forcefields or Forcefield modification files (frcmod) we should consider
# when parameterizing the solvated system
EXTRAFF = leaprc.ff14SB, frcmod.ionsjc_tip3p

[MDSETTINGS]
# Comma separated list of solvent box. MANDATORY
SOLVENTS = ETA, PYZ, WAT

# Number of replicas for each solvent. DEFAULT:1
NREPL = 3

# Number of nanoseconds to run. DEFAULT:20
NANOS = 50

# Temperature of each replica. DEFAULT:300K
TEMP = 300

# Restraining scheme (HA=heavyatoms , BB=back-bone atoms; FREE). DEFAULT: FREE.
RESTR= BB

# Restraining force when RESTR!=FREE. In kcal/mol.A^2. Default=0.0.
FORCE= 0.01

# Mask of residues where restraints should be applied (default: auto detect).
#RESTRMASK = 3-283,290-296
EOF

# 4. Create the project

mdmix create project -n $PName -f $PWD/PTP_${PName}.cfg

# mdmix info project

# Create MD files
cd $MD
mdmix queue write -n $QueueName
