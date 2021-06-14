#!/usr/bin/python

import os
import re
import sys
from rdkit import Chem, DataStructs
from rdkit.Chem import PandasTools, rdFMCS, MACCSkeys
from rdkit import RDConfig
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse



'''
https://www.rdkit.org/docs/GettingStartedInPython.html
'''

def sub_Matrix(patterns, ref):
	# Create M x N Dataframe to check substructure
    substructure = pd.DataFrame(columns=patterns['Name'], index=ref['Name'])

    # create vectors containing Mols, the first one we used the conversion from mol-->smarts--> mol
    patterns_smarts = [Chem.MolFromSmarts(Chem.MolToSmarts(x)) for x in patterns['ROMol']]
    # this second vector we used the direct mol object
    refs_smarts = [x for x in ref['ROMol']]



    for index, row in patterns.iterrows():  # iterating over the patterns

    	for index2, row2 in ref.iterrows():  # iterating over the fragments

    		substructure.loc[row2['ID'],row['ID']] = refs_smarts[index2].HasSubstructMatch(patterns_smarts[index])

    
    return substructure


def tanimoto_Matrix(patterns, ref):
	# Create M x N Dataframe to check substructure
    substructure = pd.DataFrame(columns=patterns['Name'], index=ref['Name'])

    # we use chemical similarity fingerprints (not topological)
    patterns_fps = [Chem.MACCSkeys.GenMACCSKeys(x) for x in patterns['ROMol']]
    refs_fps = [Chem.MACCSkeys.GenMACCSKeys(x) for x in ref['ROMol']]


    for index, row in patterns.iterrows():  # iterating over the patterns
    	for index2, row2 in ref.iterrows():  # iterating over the fragments

    		substructure.loc[row2['ID'],row['ID']] = DataStructs.TanimotoSimilarity(patterns_fps[index],refs_fps[index2])


    return substructure


def tanimoto_MCS(patterns, ref):
	# Create M x N Dataframe to check substructure
    substructure = pd.DataFrame(columns=patterns['Name'], index=ref['Name'])

    for index, row in patterns.iterrows():  # iterating over the columns of substructure
    	for index2, row2 in ref.iterrows():  # iterating over the indexes of substructure

    		substructure.loc[row2['ID'],row['ID']] = rdFMCS.FindMCS((row2['ROMol'],row['ROMol'])).numBonds
    		# .numAtoms, .numBonds, .canceled

    return substructure



def main():
    # Read patterns from sdf file and import into pandas dataframe
    patterns = PandasTools.LoadSDF(args.in_pat,molColName='ROMol', includeFingerprints=True, isomericSmiles=True, smilesName='smiles',embedProps=True)
    #print(patterns.info())
    
    # Read fragments sdf file and import into pandas dataframe
    reference = PandasTools.LoadSDF(args.in_ref,molColName='ROMol', includeFingerprints=True, isomericSmiles=True, smilesName='smiles',embedProps=True)
    #print(reference.info()) 

    # Create M x N Dataframe to check substructure
    subs = sub_Matrix(patterns, reference)

    # Create M x N Dataframe to check tanimoto similarity
    taani = tanimoto_Matrix(patterns, reference)

    # Create M x N Dataframe to check tanimoto MCS similarity
    # taani_MCS = tanimoto_MCS(patterns, reference)


    # Remove the ID that does not have coordinates and fill NA values
    subs = subs.drop(['5uvw'])
    subs = subs.fillna(99)
    taani = taani.drop(['5uvw'])
    taani = taani.fillna(99)
    #taani_MCS = taani_MCS.drop(['5uvw'])
    #taani_MCS = taani_MCS.fillna(99)
    #print(subs.head(20))
   
    with sns.axes_style("white"):
    	ax = sns.heatmap(subs, cmap="YlGnBu")
    	plt.show()




    #print(taani.head(20))
    with sns.axes_style("white"):
    	ax2 = sns.heatmap(taani, cmap="YlGnBu")
    	plt.show()



    print(max(taani.max()))
    # print(taani_MCS.head(20))
    # with sns.axes_style("white"):
    # 	ax3 = sns.heatmap(taani_MCS, cmap="YlGnBu")
    # 	plt.show()







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select molecules from .SDF file that exist in an .SMI file.')

    parser.add_argument('--ipat', dest='in_pat', help='Enter input .sdf file format containing patterns.')
    parser.add_argument('--iref', dest='in_ref',help='Enter input .sdf file format with fragments.')
    #parser.add_argument('--osdf', dest='out_sdf',help='Specify output .sdf file format.')

    args = parser.parse_args()

    if args.in_pat is None or args.in_ref is None: # or args.out_sdf is None:
        print("Type 'python select_ID_sdf_from_smi.py -h' to introduce the arguments.")
    else:
        main()
