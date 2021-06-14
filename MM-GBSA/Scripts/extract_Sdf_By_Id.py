import os
import re
import sys
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit import RDConfig
import pandas as pd
import argparse


def main():
    # Read smiles and IDs from -smi file
    columns = ['ID(s)']
    smiles = pd.read_csv(args.in_smi, names=columns)
    print(smiles.info())
    # Read sdf file and import into pandas dataframe
    suppl = PandasTools.LoadSDF(args.in_sdf)
    print(suppl.info())
    # Check whether the Ids from the smile file exist in the sdf file and merge both datasets
    merge = suppl.loc[suppl['ID'].isin(smiles['ID(s)'])]
    print(merge.info())
    # Write the results into sdf file
    PandasTools.WriteSDF(merge,args.out_sdf, properties=list(merge.columns))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select molecules from .SDF file that exist in an .SMI file.')

    parser.add_argument('--isdf', dest='in_sdf', help='Enter input .sdf file format.')
    parser.add_argument('--iids', dest='in_smi',help='Enter input .smi file format.')
    parser.add_argument('--osdf', dest='out_sdf',help='Specify output .sdf file format.')

    args = parser.parse_args()

    if args.in_sdf is None or args.in_smi is None or args.out_sdf is None:
        print("Type 'python select_ID_sdf_from_smi.py -h' to introduce the arguments.")
    else:
        main()
