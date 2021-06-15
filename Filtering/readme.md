# 1. Filtering
Depending if the DB is compressed or not:

1. Not Compressed (SMILES format):
    * `python filter_HA_R.py`
2. Compressed (ZIP format):
    * `python filter_HA_R_zip.py`

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
