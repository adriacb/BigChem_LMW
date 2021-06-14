import os
import re
import sys
from datetime import datetime
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


if __name__ == "__main__":

    start_time = datetime.now()
    # Get the names of the compressed files
    file_DB = "/NAS/DB/testUNIQUE/ERZ20_unique.smi"
    mols = list()
    c = 0
    with open(file_DB) as f:
        for line in f:
            if c < 50:
                l = line.decode('utf8').split()[0]
                mol = Chem.MolFromSmiles(l)
                mols.append(mol)
                c += 1



    img = Draw.MolsToGridImage(mols, molsPerRow=10, subImgSize=(250, 250))
    img.save('/home/acabello/Documents/ClustRdkit/totalDB.png')
    end_time = datetime.now()
    print("Total time: {}".format(end_time-start_time))