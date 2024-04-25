#$ -S /bin/bash
#$ -cwd
#$ -pe ompi 16
#$ -V
#$ -q all.q 

#mpiexec -n 32 --hostfile host2 python test.py
source /home/mdsoft/mpi-selector/data/ib-openmpi-4.0.3_intel-21.4.0_cuda-11.4_cent7.sh
bindir=/home/yklei/software/my_genesis_multipole/bin

source ~/.bashrc
conda activate pytorch

#export idx=`cat idx`

python read_label_data.py \
       --dir /swork/yklei/Tim/7.prod4_ob/prod5_ID.1 \
       --outfile Tim_test_ml_ID.h5 \
       --top ../step4_nvt_100.psf \
       --ref_model 0 \
       --trj prod5_ID.dcd \
       --logfile log_ID \
       --interval 1 \
       --drop_first \
       --store_nc \
       --linkbond 1.0 

