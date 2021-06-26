import os
import re
import sys
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit import RDConfig
import pandas as pd
import matplotlib.pyplot as plt
import argparse


def plotHist(filtered):
    ax = filtered.hist(column='r_psp_MMGBSA_dG_Bind', grid=False, figsize=(12,8), color='#86bf91', zorder=2, rwidth=0.9)

    ax = ax[0]
    for x in ax:

        # Despine
        x.spines['right'].set_visible(False)
        x.spines['top'].set_visible(False)
        x.spines['left'].set_visible(False)

        # Switch off ticks
        x.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", left="off", right="off", labelleft="on")

        # Draw horizontal axis lines
        vals = x.get_yticks()
        for tick in vals:
            x.axhline(y=tick, linestyle='dashed', alpha=0.4, color='#eeeeee', zorder=1)

        # Remove title
        x.set_title("")

        # Set x-axis label
        x.set_xlabel("dG Bind", labelpad=20, weight='bold', size=12)

        # Set y-axis label
        x.set_ylabel("Counts", labelpad=20, weight='bold', size=12)
    plt.show()    

def main():
    # Read MMGBSA output
    #columns = ['ID(s)', 'r_psp_MMGBSA_dG_Bind', 'r_psp_MMGBSA_dG_Bind_Solv_GB','r_psp_MMGBSA_dG_Bind(NS)','r_psp_MMGBSA_dG_Bind(NS)_Solv_GB']
    mmgbsa = pd.read_csv(args.mmgbsa)
    print(mmgbsa.info())

    # We filter those dG_Bind values having a value less than [args.cutoff]
    filtered = mmgbsa.loc[mmgbsa['r_psp_MMGBSA_dG_Bind'] <= float(args.cutoff)]
    #print(filtered.info())
    print(filtered['r_psp_MMGBSA_dG_Bind'].describe())
    print()

    # Plot histogram dG bind
    plotHist(filtered)


    # Read sdf file and import into pandas dataframe
    suppl = PandasTools.LoadSDF(args.in_sdf)
    #print(suppl.info())

    # Check whether the Ids from the MMGBSA file exist in the sdf file and find the union between the datasets
    union = suppl.loc[suppl['ID'].isin(filtered['title'])]
    print(union.info())

    #Now merge results from MM/GBSA
    merged = pd.merge(union, filtered[['r_psp_MMGBSA_dG_Bind', 'r_psp_MMGBSA_dG_Bind_Solv_GB','r_psp_MMGBSA_dG_Bind(NS)','r_psp_MMGBSA_dG_Bind(NS)_Solv_GB']].astype(object), left_index=True, right_index=True)
    print(merged.info())

    # Write the results into sdf file adding the MM/GBSA information
    PandasTools.WriteSDF(merged, args.out_sdf, properties=list(merged.columns))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select molecules from .SDF file that exist in an .SMI file.')

    parser.add_argument('--isdf', dest='in_sdf', help='Enter input .sdf file format.')
    parser.add_argument('--cutoff', dest='cutoff', help='Enter input .sdf file format.')
    parser.add_argument('--imm', dest='mmgbsa',help='Enter input .smi file format.')
    parser.add_argument('--osdf', dest='out_sdf',help='Specify output .sdf file format.')

    args = parser.parse_args()

    if args.in_sdf is None or args.mmgbsa is None or args.out_sdf is None or args.cutoff is None:
        print("Type 'python select_ID_sdf_from_smi.py -h' to introduce the arguments.")
    else:
        main()
