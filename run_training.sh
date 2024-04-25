#$ -S /bin/bash
#$ -cwd
#$ -N al
#$ -j y
#$ -pe impi 24
#$ -q all.q

rm /dev/shm/shared*
source ~/.bashrc
source /home/mdsoft/mpi-selector/data/ib-openmpi-4.0.3_intel-21.4.0_cuda-11.4_cent7.sh
bindir=/home/yklei/software/my_genesis/bin

conda activate pytorch_new
python mlmm_main_pl.py fit --config config.yml > log   #log_rank_loss2
df -h 
rm /dev/shm/shared*
