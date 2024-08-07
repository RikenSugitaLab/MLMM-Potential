# MLMM-Potential Repository README

## This program is used to fit QM/MM potential to get a accurate and transferrable Machine Learning Force Field. 

**mlmm**: package of network.  
**server**: interface bewteen genesis and Machine learning program by socket.  
**my_genesis_multipole**: modified version of genesis used for simulation with ML potential.  
**scripts_for_prepare_data**: python program used to extract reference force, energy and charge from output of genesis  
**config.yml**: config file for training  
**run.inp**: genesis input sample for MLMM simulation  
**run.sh**: script for runing MLMM simulation  
**run_training.sh**: script for running traning program  
**mlmm_main_pl.py**: main program file for training  

# Required Library:
pytorch<2.0.0  
pytorch_lightning<2.0.0 (with CLI function)    
MDAnalysis  
mpi4py  
dgl  
dask  
netCDF4  
numpy  
h5py  
pmda  

Intall by pip install -r requirements.txt

# Usage
**training**: shown in run_traning.sh (more details can be found in website of pytorch lightning).  
**simulation**: copy server folder into working directory and put mlmm package into server folder. Then run simulation according to run.sh.  
Detailed explanation on options can be found in doc. 

## Training Data Preparation
For training, three files are needed
1. Topology file containing information of system such as psf, prmtop, pdb file (formats supported by MDAnalysis can be used)

2. Trajectory file containing coordinate and force such as dcd, nc file (formats supported by MDAnalysis can be used)  
   Recommend to use nc format, since force can be stored. For dcd file, force can not be stored.

3. Label File containing energy and qm charge (hdf5 file). The key in hdf5 file should be set as qm_energy and qm_charge.  
   If no force is stored in trajectory file, please store the force of qm/mm atoms in hdf5 file with the key qm_force and mm_force repectively.

### Default Unit: 
Energy: kcal/mol   
Force: `kcal/mol/angstrom`  
Coordinate: Angstrom  

# Notes
For the time being, the validate and predict function can't not be used. To generate prediction of any conformation, please store it into dcd format and use energy_analysis module in genesis to get the prediction.

# Citations
Yao-Kun Lei, Kiyoshi Yagi, Yuji Sugita; Learning QM/MM potential using equivariant multiscale model. J. Chem. Phys. 160 (21): 214109. https://doi.org/10.1063/5.0205123
