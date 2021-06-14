import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem import AllChem
from progress.bar import IncrementalBar
sns.set_style("white")

total = 4123967

def getMolsPMIs(pathfile):
    rtb = pd.DataFrame(columns=["RBonds"])
    c = 0
    bar = IncrementalBar('Progress:', max = total)
    with open(file_DB) as f:
    	for line in f:
        	bar.next()
        	smile = re.match(r"(.*)", line).group(0).replace("\t","")
        	mol = Chem.MolFromSmiles(smile)
        	try:
        		rbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        		rtb = rtb.append({'RBonds':rbonds}, ignore_index=True)
        	except:
        		print(smile)
        		pass
    bar.finish()
    return rtb




if __name__ == "__main__":

    start_time = datetime.now()
    # Get the names of the compressed files
    file_DB = "" # MUST BE CHANGED

    mols = getMolsPMIs(file_DB)
    sns.histplot(mols, x = "RBonds", color='#404080' ,bins= 6).set_title('Rotatable Bonds')
    plt.savefig('RotBonds.png', dpi=300)
    #plt.show()


    end_time = datetime.now()
    print("Total time: {}".format(end_time-start_time))
