# Manual for ML/MM Simulation in GENESIS

To run ML/MM simulation
model file with .pt extension is needed 

## Related Option in GENESIS
[QMMM] Section  
qmtyp: set it to **ml**

max_order: The truncation order of Taylor Expansion of external potential. 1,2,3 is supported. Please use 3.

lock_id: During Initialization of interface, to synchronizea GENESIS with ML Server before socket connection, ML Server will wait for the generation of a file named by file.lock.{lock_id}. **Therefore if multiple tasks is ran on the same working directory, please set this option to different values. These files should be deleted before running new tasks.**

port: port for socket connection

## Related Option in ML Server
model_path: Path of model file ( .pt)  

n_atom: number of QM atoms

max_order: The truncation order of Taylor Expansion of external potential. 1,2,3 is supported. Please use 3.

precision: precision of numerical calculation ("float32" or "float64) 
model file should be stored in correpsponding precision

max_connections: number of requests for ML calculation. (n_steps + 1, n_step: number of integration steps)

port: The same as option port in GENESIS

lockindex: The same as option lock_id in GENESIS

cutoff: cut off distance for graph convolution (angstrom)

ngpu: number of avaiable gpus in each node. Default: 1

npergpu: number of models in each gpu. Default: 1

n_replicas: number of replicas in each node. Default: 1
This option should be set when Relica Exchange Sampling Methods is used, and ngpu as well as npergpu should be set since multiple processes (replica) may ran on one node.
The total number of models in each node is ngpu\*nergpu, this number should be equal to the number of replicas in one node.








