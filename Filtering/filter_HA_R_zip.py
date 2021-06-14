import os
import re
import sys
from datetime import datetime
import zipfile
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import multiprocessing as mp
from multiprocessing import Pool

## Use this script for compressed DB (EnamineReal)

directory_DB = ''
outdir = ''

try:
    # Change the current working Directory
    smiles = os.chdir(directory_DB)
    print("Directory changed: {}".format(os.getcwd()))
except OSError:
    print("Can't change the Current Working Directory")


def countHeavyAtomsSmiles(files):
    """It takes as input a list of path files
    returns the total amount of the DB and filtered
    and store the smiles and IDs in a file"""
    molsFilt = 0
    totalDB = 0
    c = 0

    for fn in files:
        with zipfile.ZipFile(fn) as f:
            for files in f.namelist():
                print(files)
                out = open(outdir+str(files), "a+")
                if files.endswith('.smiles'):
                    for line in f.open(files).readlines()[1:]:
                        l = line.decode("utf-8", "ignore").split()
                        smile = l[0]
                        ID = l[1]
                        mol = Chem.MolFromSmiles(smile)
                        rc = mol.GetRingInfo().NumRings()
                        HA = mol.GetNumHeavyAtoms()
                        totalDB += 1
                        if HA <= 14 and rc >= 1:
                            molsFilt += 1
                            if c <= 50:
                                Chem.Draw.MolToImageFile(mol, outdir+'1ring_{}.png'.format(c), size=(200,200))
                                c += 1
                                #write to file
                            out.write("{} {}\n".format(smile, ID)) 
                out.close()
        f.close()
    return molsFilt, totalDB



if __name__ == "__main__":

    start_time = datetime.now()
    # Get the names of the compressed files
    files = [f for f in os.listdir(os.getcwd()) if f.endswith(".zip")]
    totalFilt,totalDB = countHeavyAtomsSmiles(files)
    print("Total DB = {}\nTotal Filtered = {}".format(totalDB, totalFilt))
    end_time = datetime.now()
    print("Total time: {}".format(end_time-start_time))
