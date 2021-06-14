# BigChem_LMW

## Intro

### 1. Filtering

Depending if the DB is compressed or not:

1. Not Compressed (SMILES format):
    * `python filterAtoms.py`
2. Compressed (ZIP format):
    * `python filterAtomsG.py`

The output is a file containing the filtered molecules in SMILES format .smi.

#### Data analysis

1. Molecular descriptors
    * `python PMI.py` Principal Moments of Inertia.
    * `Rscript HBA_HBD.r` Hydrogen Bond Donors and Acceptors.
    * `Rscript LogP_AMR.r` LogP and AMR.
    * `Rscript mweight.r` Molecular Weight.
    * `python rotBonds.py` Rotatable Bonds.
2. Histogram
    * `python FragmentsStats.py` Distributions of fragments.
3. Data visualization
    * `python DrawGrid.py` Visualization of molecules.

### 2. Structure and Library preparation

1. Structure preparation with MOE.

    * `MOE | File | Open | $MOE/4lr6.pdb.gz` Load 4LR6 into MOE.
    * `MOE | Window | Sequence Editor` open the editor.
    * `MOE | Compute | Prepare | Structure preparation` prepare structure.
    * `Structure Preparation | Popup | Clear Selection` and `Structure Preparation | Correct`.

2. Ligand preparation with Corina and ChemAxon.

    * `bash ligprep.sh`

### 3. pyMDMix

### 5. rDock

1. System Definition
    * `4lr6_pharma.prm` Containing the Cavity definition, cavity restraint penalty, and ph4 restraints.
    * `pharma.restr` Containing the ph4 restraints coordinates obtained from pyMDMix.
3. Cavity Generation
    * `rbcavity -was -d -r 4lr6_pharma.prm` rDock.
    * `pymol 4lr6_prep.mol2 4lr6_pharma_cav1.grd ligand_centered.sdf` Visualization with PyMol.
    * `isomesh cavity, 4lr6_pharma_cav1, 0.99` Visualize the cavity.
4. Docking
    1. We need a VS_FILTER file:
        * `VS_FILTER.txt` Containing the configuration for the HTVS protocol from rDock.
    2. We used SLURM.
        * `VS_slurm.q`, which basically contains: `rbdock -i $curr_file -o $outdir/docked_filter_2_$file_name -r 4lr6.prm -p dock.prm -t VS_FILTER.txt`.
    3. We can use an alternative for filtering after the Docking protocol.
        * `filter_VS.q`

### 6. MM-GBSA

We used SLURM to perform the computations.

1. Add hydrogens to the output sd file from rDock.

     * `obabel -isdf -----.sdf -osdf -h -O ligands_docking_H.sdf`

2. The MMGBSA requires a MAE format containing the RECEPTOR, and the LIGANDS in this specific order.

     * `$ SCHRODINGER/utilities/structconvert -imol2 receptor.mol2 -omae receptor.mae`
     * `$ SCHRODINGER/utilities/structconvert -isd ligands_docking_H.sdf -omae ligands_docking.mae`
     * `$ cat receptor.mae ligands_docking.mae > input.mae`

3. Run MM-GBSA Prime (Schrödinger)

     * `$ prod_schrondinger/prime_mmgbsa -out_type LIGAND -job_type REAL_MIN -OVERWRITE -jobname mmgbsa_output -WAIT -LOCAL $initdir/input.mae`

Two strategies:

* Perform all the computations together
   * `sbatch MMGBSA_slurm2.q`
* If the are enough licenses we can split the molecules.
