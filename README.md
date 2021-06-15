# BigChem_LMW

## Intro

### 1. Filtering
### 2. Structure and Library preparation
### 3. pyMDMix

1. Preparation
   * PyMol: clean structure (remove ligand and H20)
   * PDB Check: MOE (cappings, H, check aa, conforming names to AMBER)
   * TLEAP: prepare off (or PDB), input MDMix : check dpb, neutralize system, saveoff (or save pdb)
      * `tleap -f leaprc.ff99SB`
         >       4lr6 = loadpdb 4lr6.pdb
         >       check 4lr6
         >       saveoff 4lr6 4lr6.off
         >       quit

3. Create Project
   * cfg: solvents, ff14SB, name project, ns, replicas...
   * `mdmix create project -n 4lr6 -f 4lr6.cfg`
4. Run simulation
   * launch simulation (min.q) per replica
   * mdmix info project
5. Analysis (analisis.q)
   * align trajectories
   * calculate density grid
   * calculate energy
6. PyMol results
   * load reference structure
   * load energy grids
      * isosurface ETA_CT, ETA_CT_DGO, -1
      * color by probe
    
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
   * `MMGBSAsplit.q` splits the sdf file into desired number of molecules.
   * `MMGBSAprep.q` wich basically converts every split into MAE format.
   * `MMGBSA_slurm.q` runs the MMGBSA using a SLURM array, the number of jobs running at the same time cannot be greater than the number of CPUS.
