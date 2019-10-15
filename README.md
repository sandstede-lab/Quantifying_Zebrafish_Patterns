# Quantifying Zebrafish Patterns
Code for quantifying spots and stripes using topological data analysis and machine learning 

_Authors:_ Melissa R. McGuirl, Alexandria Volkening, Bjorn Sandstede 

For questions/comments please contact Melissa R. McGuirl at melissa_mcguirl@brown.edu.


## Description 

This software compute pattern statistics using topological data analysis and machine learning techniques. The input of this software is a collection of coordinate data from pigment cells.

This software is based upon the work presented in M.R. McGuirl, A. Volkening, and B. Sandstede, "Topological data analysis of zebrafish patterns" (2019). In preparation.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

#### Programs
*  Ripser.py v0.3.2 (https://ripser.scikit-tda.org/index.html), 
*  Python 
*  Matlab 

#### Python libraries
 * ripser
 * matplotlib
 * numpy

### Install the ripser program as follows: 
```
	pip install Cython
	pip install Ripser
```

### Source 
```
      cd 
      git clone https://github.com/MelissaMcguirl/Quantifying_Zebrafish_Patterns
      cd Quantifying_Zebrafish_Patterns
      pip install -r requirements.txt
```

### Input data format: 

```
	Input data is assumed to be a .mat file containing cell coordinate data over time for M, I^l, I^d, X^c (X^d), and X^s (X^l) cells under the naming scheme cellsM, cellsIl, cellsId, cellsXc, and cellsXsn, respectively. 
	The cell coordinate data is a 3-d matrix of sixe N X 2 X M, where N is the maxmimum number of cells present across all cells and all stages of development and M is the number of time points in your data. 
	The .mat file should also contain vectors of size M X 1 containing the number of cells of each cell type per time point, named numMel, numiril, numIrid, numXanc, and numXansn.  

```

### Data samples:    

```
      cd data/sample_inputs
      1) Out_WT_default_1.mat (simulation of wild-type pattern formation from agent-based model of A.V. and B.S. under the default parameter regime)
      2) Out_pfef_default_1.mat  (simulation of pfeffer pattern formation from agent-based model of A.V. and B.S. under the default parameter regime)
      3) Out_shady_default_1.mat  (simulation of shady pattern formation from agent-based model of A.V. and B.S. under the default parameter regime)
      4) Out_nacre_default_1.mat (simulation of nacre pattern formation from agent-based model of A.V. and B.S. under the default parameter regime)
```


### Matlab examples:    

```
      cd src/matlab
      1) test_WT.m (quantify wild-type stripes)
      2) test_pfeffer.m (quantify pfeffer spots)
      3) test_shady.m (quantify shady spots)
      4) test_nacre.m (quantify nacre spots)
```


## Notes

Code notation: X^C is equivalend to X^d in the paper, and X^S is equivalent to X^l. 

