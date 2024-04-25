# MLMM-Potential Repository README

This program is used to fit QM/MM potential to get a accurate and transferrable Machine Learning Force Field. 

**mlmm**: package of network.  
server: interface bewteen genesis and Machine learning program by socket.  
my_genesis_multipole: modified version of genesis used for simulation with ML learning potential.  
scripts_for_prepare_data: python program used to extract reference force, energy and charge from output of genesis  
config.yml: config file for training  
run.inp: genesis input sample for MLMM simulation  
run.sh: script for runing MLMM simulation  
run_training.sh script for running traning program  
mlmm_main_pl.py: main program file for training  

# Required Library:
pytorch<=2.0.0  
pytorch_lightning<=2.0.0  
MDAnalysis  
mpi4py  
dgl  

Intall by pip install -r requirements.txt
