import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Descriptors3D
from rdkit.Chem import AllChem

sns.set_style("whitegrid")

def getListFiles(folder, extension = '.sdf'):
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

def getMolsPMIs(files):
    PMIs = pd.DataFrame(columns=["SMILE","NPR1","NPR2"])
    failed = 0
    tot = 0
    for file in files:
        suppl = Chem.SDMolSupplier(file)
        for m in suppl:
            #mol = Chem.MolFromSmiles(m)
            try:
                NPR1, NPR2 = inertia(m)
                PMIs = PMIs.append({'NPR1':NPR1, 'NPR2':NPR2}, ignore_index=True)
                tot += 1
                print(tot)
            except:
                failed += 1
                print(failed)
                pass
    return PMIs, failed


def inertia( mol ):

    NPR1 = Descriptors3D.NPR1(mol, 0)
    NPR2 = Descriptors3D.NPR2(mol, 0)
    return (NPR1, NPR2)


if __name__ == "__main__":

    start_time = datetime.now()
    # Get the names of the compressed files
    file_DB_path = "" # MUST BE SPECIFIED
    sdf_files = getListFiles(file_DB_path)

    molsPMIs, failed = getMolsPMIs(sdf_files)
    print("Failed PMIs: {}".format(failed))
    fig,ax = plt.subplots()
    plane = plt.imread('/ImagesPMI/plane.png')
    ring = plt.imread('/ImagesPMI/ring.png')
    d3 = plt.imread('/ImagesPMI/d3.png')
    sns.regplot(data= molsPMIs, x = "NPR1", y = "NPR2", fit_reg=False, marker = "+")
    ax.imshow(d3, aspect='auto', extent=(1.00, 1.05, 1.00, 1.05))
    ax.imshow(plane, aspect='auto', extent=(-0.03, 0.02, 1.00, 1.05))
    ax.imshow(ring, aspect='auto', extent=(.466, 0.5, 0.5, .466))
    plt.ylim(bottom=0.466, top=1.05)
    plt.xlim(left = -0.03, right=1.10)
    plt.title("Principal Moments of Inertia")
    plt.savefig('PMIs.png', dpi=1200)
    #plt.show()


    end_time = datetime.now()
    print("Total time: {}".format(end_time-start_time))
