# 5. Docking

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

## Filtering H's and put Hyd
