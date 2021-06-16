6. MM-GBSA

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
