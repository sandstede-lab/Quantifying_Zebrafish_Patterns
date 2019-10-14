'''
    get_barcodes.py

    Input: text file distance matrix
    Output: persistent homology barcodes from Ripser


    Author: Melissa R. McGuirl, Brown University, 2019.

'''

from ripser import ripser
import argparse
import numpy as np

def main():
    descriptor = "calls ripser to compute persistent homology"
    parser = argparse.ArgumentParser(description = descriptor)
    parser.add_argument('-i', '--indir',action = 'store',required = True,  help = '''provide path to directory containing distance matrix input''')
    parser.add_argument('-d', '--dim', action = 'store', required = True, help = '''provide maximum homology dimension for ripser computation''')
    parser.add_argument('-o', '--outdir',action = 'store',required = True,  help = '''provide path to directory for output''')


    args = parser.parse_args()
    inFile = args.indir
    outFile = args.outdir
    D = np.loadtxt(inFile, dtype='float', delimiter = ',')
    dim = int(args.dim)
    bars = ripser(D,  distance_matrix=True, maxdim=dim)
    for i in range(dim + 1): 
        PD = bars['dgms'][i]
        np.savetxt(outFile +  '_dim' +  str(i), PD)    
        
if __name__=="__main__":
    main()

