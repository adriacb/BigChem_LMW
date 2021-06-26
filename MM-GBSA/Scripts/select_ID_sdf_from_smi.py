import os
import re
import sys
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
import argparse


def main():
    print(datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
    '''
    http://rdkit.org/docs/source/rdkit.Chem.PandasTools.html
    '''

    # Read smiles and IDs from -smi file
    columns = ['Smiles','ID(s)']
    smiles = pd.read_csv(args.in_smi, names=columns, sep='\t')
    print(smiles.info())
    # Read sdf file and import into pandas dataframe
    suppl = PandasTools.LoadSDF(args.in_sdf)
    print(suppl.info())
    # Check whether the Ids from the smile file exist in the sdf file and find the union of both datasets
    union = suppl.loc[suppl['ID'].isin(smiles['ID(s)'])]
    print(union.info())
    # Write the results into sdf file
    PandasTools.WriteSDF(union,args.out_sdf, properties=list(union.columns))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select molecules from .SDF file that exist in an .SMI file, and save it into an sdf file.')

    parser.add_argument('--isdf', dest='in_sdf', help='Enter input .sdf file format.')
    parser.add_argument('--ismi', dest='in_smi',help='Enter input .smi file format.')
    parser.add_argument('--osdf', dest='out_sdf',help='Specify output .sdf file format.')

    args = parser.parse_args()

    if args.in_sdf is None or args.in_smi is None or args.out_sdf is None:
        print("Type 'python select_ID_sdf_from_smi.py -h' to check the arguments.")
    else:
        main()
