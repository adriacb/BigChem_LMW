MD with solvent mixures | Identification of HotSpots | pyMDMix

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
    
