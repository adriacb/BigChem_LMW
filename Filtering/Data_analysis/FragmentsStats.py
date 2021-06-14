#! /usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from rdkit import Chem
from rdkit import DataStructs
import multiprocessing as mp
from multiprocessing import Pool, Process, Queue, cpu_count

def parseF(file_path):
    """ Read first n lines from file """
    data = dict() 
    file_n = file_path.split("/")[-1]
    with open(file_path) as f:
        print(file_n)

        for line in f: # This step is High time consuming, maybe reading smiles directly and count length
            line = line.rstrip('\n')
            mol = Chem.MolFromSmiles(line)
            size = mol.GetNumHeavyAtoms() # heavyatomcount()
            data[size] = data.get(size, []) + [1]
    f.close()

    return data


def MpParsing(file_list):
    global cpu_cores
    #cpu_cores = cpu_count()
    q = Queue()
    procs = []

    # "Chunk" the input file list into sublists by CPU
    for i in range(0, cpu_cores):
        sub_list = [file_list[j] for j in range(0, len(file_list)) if j % cpu_cores == i]

        if len(sub_list) > 0:
            p = Process(target=procFiles, args=([sub_list, q]))
            p.start()
            procs.append(p)

    #all_results = []
    ALL = dict()
    for i in range(0, len(procs)):
        #all_results.append(q.get())
        currDict = q.get()[0]
        for key in currDict:
            if key not in ALL:
                ALL[key] = currDict[key]
            else:
                ALL[key] += currDict[key]

    return ALL


def procFiles(file_list, q):
    results = []
    try:
        for f in file_list:
            results.append(parseF(f))
    except:
        q.put([])
        raise
    # Put results into queue
    q.put(results)

def getListFiles(folder, format = ".smi"):
    """ Returns a list of files inside the <folder> directory """
    file_list = []
    for root, dirs, files in os.walk(folder):
        for dir in dirs:
            dirpath = os.path.join(root, dir)
            print(dirpath)
        for filename in files:
            if filename.endswith(format) and os.stat("/data1/acabello/BigChem_VS/fzinc20/"+filename).st_size != 0:
                filepath = os.path.join(root, filename)
                # Check if it's file
                if os.path.isfile(filepath):
                    file_list.append(filepath)

    return file_list

def renderPlots(dictionary):

    data = []
    for i, (key, val) in enumerate(dictionary.iteritems()):

        for v in val:
            data.append((key))


    counts, bins = np.histogram(data)
    plt.hist(bins[:-1], weights=counts, histtype='barstacked') # histtype='barstacked'
    plt.title("ZINC20 distribution of fragments")
    plt.xlabel("Number of atoms")
    plt.ylabel("Counts")
    mean = np.mean(data)
    plt.axvline( mean , color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(mean*0.9, max_ylim*(0.9), 'Mean: {:.2f}'.format( mean ))
    plt.show()

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='FragmentsStats')
    parser.add_argument('-c', action="store", dest="CP", help='Set number of CPUs.', type=int, required=True)
    args = parser.parse_args()
    path = "/data1/acabello/BigChem_VS/fzinc20/"
    cpu_cores = args.CP

    start_time = datetime.now()
    print("Start time: {}".format((start_time)))    
    file_list = getListFiles(path)
    Dataset = MpParsing(file_list) # Use a subset for testing purposes

    end_time = datetime.now()
    print("End time: {}".format((end_time)))
    print("Time elapsed: {}".format((end_time-start_time)))
    renderPlots(Dataset)
