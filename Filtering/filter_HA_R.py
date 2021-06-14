import os
import re
import sys
from datetime import datetime
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

## Use this script for non-compressed DB



directory_DB = ''
outdir = ''

def getListFiles(folder, extension = '.smi'):
    """ Returns a list of files inside the <folder> directory """
    file_list = []
    for root, dirs, files in os.walk(folder):
        for dir in dirs:
            dirpath = os.path.join(root, dir)

        for filename in files:
            filepath = os.path.join(root, filename)
            # Check if it's file
            if os.path.isfile(filepath) and filepath.endswith(extension):
                file_list.append(filepath)

    return file_list

def countHeavyAtomsSmiles(files):
    """It takes as input a list of path files
    returns the total amount of the DB and filtered
    and store the smiles and IDs in a file"""
    molsFilt = 0
    totalDB = 0
    c = 0

    for fn in files:
        print(fn)
        n = fn.split('/')[-1]
        if os.path.getsize(fn) > 0: # If it is not empty
            with open(fn) as f:
                out = open(outdir+str(n), "a+")
                next(f)
                for line in f:
                    smile = line.split()[0]
		    ID = line.split()[1]
                    mol = Chem.MolFromSmiles(smile)
                    rc = mol.GetRingInfo().NumRings()
                    HA = mol.GetNumHeavyAtoms()
                    #print(HA, rc)
                    totalDB += 1
                    if HA <= 14 and rc >= 1:
                        #print(n)
                        molsFilt += 1
                        #if c <= 50: # Draw 51 first
                        #    Chem.Draw.MolToImageFile(mol, outdir+'1ring_{}_z20.png'.format(c), size=(200,200))
                        #    c += 1
                            #write to file
                        out.write(line+' '+ID+"\n") # write the smile and ID

            out.close()
            f.close()
    return molsFilt, totalDB



if __name__ == "__main__":

    start_time = datetime.now()
    # Get the names of the compressed files
    files = getListFiles(directory_DB)
    print(len(files))
    totalFilt,totalDB = countHeavyAtomsSmiles(files)
    print("Total DB = {}\nTotal Filtered = {}".format(totalDB, totalFilt))
    end_time = datetime.now()
    print("Total time: {}".format(end_time-start_time))
