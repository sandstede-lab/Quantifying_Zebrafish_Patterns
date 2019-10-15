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


### Data samples:    

MAT files from simulations of the agent-based model of A.V. and B.S. under the default parameter regime are provided for sample data. If you would like more data, please contact Melissa R. McGuirl at melissa_mcguirl@brown.edu.

```
      cd data/sample_inputs
```
1) Out_WT_default_1.mat 
2) Out_pfef_default_1.mat  
3) Out_shady_default_1.mat 
4) Out_nacre_default_1.mat 
      



### Matlab examples:    

Example scripts are provided to demonstrate how to generate input data and run the program from wild-type stripes and mutants. 

```
      cd src/matlab/examples    
```

1) test_WT.m (quantify wild-type stripes)
2) test_pfeffer.m (quantify pfeffer spots)
3) test_shady.m (quantify shady spots)
4) test_nacre.m (quantify nacre spots)



## Pipeline 
The two main files are quantify_spots.m and quantify_stripes.m for quantifying spots and stripes, respectively. 


### Generating input data
The following steps are illustrates in the test examples in src/matlab/examples
	
1) Load in cell-coordinate data to MATLAB
2) Extract cell-coordinate data at time point of interest 
3) Save cell-coordinates as cells_mel, cells_iriL, cells_xanD, cells_xanL
3) Generate distance matrices of cell-cell pairwise distances and save as text file
4) Run Ripser using get_barcodes.py in src/python to get persistent homology data into PD_dir
5) Compute boundaryX and boundaryY, the right and top boundaries of the input domain (assume domain starts at origin)
6) Specify persistence cut-off for betti number computation
7) For spots, specify cell-type used for quantifying spots
8) For stripes, get time series (1) of X^d cell locations (cellsXd_all), (2) number of X^d cells (numXand_all), and y-boundaries (boundayY_all) for identifying when new interstripes form. 


### Running the program

1) To quantify stripes: quantify_stripes(cells_mel, cells_iriL, cells_xanD, cells_xanL, mel1_dir, xanC1_dir, xanS1_dir, boundaryX, boundaryY, cellsXd_all, numXand_all, boundaryY_all)
2) To quantify spots: quantify_spots(cells_mel, cells_iriL, cells_xanD, cells_xanL, PD_dir, boundaryX, boundaryY, pers_cutoff, cell_type)




