#$ -S /bin/bash
#$ -cwd
#$ -pe ompi 256
#$ -V
#$ -q gpu1.q

#mpiexec -n 32 --hostfile host2 python test.py
source /home/mdsoft/mpi-selector/data/ib-openmpi-4.0.3_intel-21.4.0_cuda-11.4_cent7.sh
bindir=/home/yklei/software/my_genesis_multipole/bin

source ~/.bashrc
#conda activate pytorch
conda activate old_gpu

rm file.lock
rm run_*.rst
rm run_*.rem
rm run_*.log
rm run_*.dcd
rm run.rst
rm -r  production.*

mpiexec -n 32 -npernode 2 -x OMP_NUM_THREADS=8 --report-bindings ${bindir}/atdyn run.inp > rp1.log : -n 32 -npernode 2 python ./server/main.py --cuda './model_more_window.pt' 6 'float32' 40001 --ngpu 32 --npergpu 1 --n_replica 32 --max_order 3 --fully_connected
echo "start"
