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
Detailed explanation on options can be found in doc (in preparation).

# Notes
For the time being, the validate and predict function can't not be used. To generate prediction of any conformation, please store it into dcd format and use energy_analysis module in genesis to get the prediction.