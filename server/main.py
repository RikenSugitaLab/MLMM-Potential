from qmmm_server import MLDriver
import socket
from sockets import Driver, InterfaceSocket, Status
import MDAnalysis as mda
import numpy as np
import time
import torch
import argparse
import logging
import time
from mpi4py import MPI
import os
import h5py

def write_data(filename,data):
    with h5py.File(filename,'w') as f:
        for i in data: 
           f[i]=data[i]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model_path', type=str, help='Path to ML model.')
    parser.add_argument('n_atom', type=int, help='number of QM atoms (including link atoms).')
    parser.add_argument('precision', type=str, help="precision of model, only support float32 or float64")
    parser.add_argument('max_connections', type=int, help="Maximum number of connections before server is shut down.")
    parser.add_argument('--port', type=int,default=31415, help="Port for socket.")
    parser.add_argument('--cuda', action="store_true", help="Device to run computation on.")
    parser.add_argument('--fully_connected', action="store_true", help="Device to run computation on.")
    parser.add_argument('--cutoff', type=float, help="Device to run computation on.")
    parser.add_argument('--ngpu', type=int,default=1, help="number of avaiable gpus in each node.")
    parser.add_argument('--npergpu', type=int,default=1, help="number of models in each gpu.")
    parser.add_argument('--n_replica',type=int,default=1, help="number of repica in each node")
    parser.add_argument('--log_level',type=str,default='info', help="log level (INFO, DEBUG).")
    parser.add_argument('--max_order',type=int,default=1,help="the max order of taylor expansion terms of electrostatic operator")
    parser.add_argument('--logfile',type=str,default='log_ml',help="the base name of log file")
    parser.add_argument('--lockindex',type=int,default=0,help="the name of lock file")
    args = parser.parse_args()

    color_program = 0
    max_order = args.max_order
    sh = (3,args.n_atom)
    sh2 = (args.n_atom)
    sh3 = (3,3,args.n_atom)
    sh4 = (3,3,3,args.n_atom)
    sh5 = (3,3,3,3,args.n_atom)
    A = np.zeros(sh, np.float64)
    A2 = np.zeros(sh2, np.float64)
    A3 = np.zeros(sh3, np.float64)
    A4 = np.zeros(sh4, np.float64)
    A5 = np.zeros(sh5, np.float64)
    Z = np.zeros(sh2, np.int32)
    accepted_connections = 0
    t_cal = 0
    t_all = 0

    try:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        comm2 = comm.Split(color=color_program, key = rank)
        size2 = comm2.Get_size()
        
        if size2 != size:
           comm = comm2
           rank = comm.Get_rank()
           size = size2
         
    except:
        rank = 0
        size = 1

    logger = logging.getLogger("main")
    if args.log_level == 'info':
        logger.setLevel(logging.INFO)
    elif args.log_level == 'debug':
        logger.setLevel(logging.DEBUG)

    dirname, filename = os.path.split(args.model_path)
    fh = logging.FileHandler(os.path.join(dirname,args.logfile+str(rank)),mode='w')
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(name)s: %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info('Starting FieldSchNet server...')

    n_model = args.npergpu * args.ngpu
    hostname = socket.gethostname()
    logger.info(f"node:{hostname}")
    if args.n_replica > 1:
        color1 = int(''.join(list(filter(str.isdigit, hostname))))
        comm_split1 = comm.Split(color=color1, key=rank)
        local_size1 = comm_split1.Get_size()
        local_rank1 = comm_split1.Get_rank()
        logger.info(f"local_rank1:{local_rank1}")

        color = int(local_rank1//(args.n_replica/n_model))
        comm_split = comm_split1.Split(color=color, key=local_rank1)
        local_size = comm_split.Get_size()
        gpu_id = color%args.ngpu
    else:
        local_rank1 = 0
        color = 0
        gpu_id = 0

    if args.precision == 'float32':
        dtype = np.float32
        torch.set_default_dtype(torch.float32)
    elif args.precision == 'float64':
        dtype = np.float64
        torch.set_default_dtype(torch.float64)

    # Set up the QM Driver
    try:
        rank_local = comm_split.Get_rank()
        local_size = comm_split.Get_size()
    except:
        rank_local = 0
        local_size = 1
    logger.info(f"local_rank:{rank_local}")

    if args.cuda:
        logger.info(f"device: gpu:{gpu_id}")
        device = torch.device("cuda:"+str(gpu_id))
    else:
        logger.info("device: cpu")
        device = torch.device("cpu")

    if rank_local == 0:
        logger.info(" set ml driver")
        mldriver = MLDriver(args.model_path, device=device, fully_connected=args.fully_connected,cutoff=args.cutoff)

    logger.info(f"open socket with port{args.port + local_rank1}")
    server = InterfaceSocket(port=args.port+local_rank1)
    server.open()
    comm.barrier()
    f = open('file.lock.'+str(args.lockindex),'w')
    #f = open('file.lock','w')
    f.close()

    i=0
    logger.info(f"accept connection")
    client, address = server.server.accept()
    driver = Driver(client)
    logger.info(" @SOCKET:   Client asked for connection from " + str(address) + ". Now hand-shaking.")
    has_data = False

    while(True):
        logger.debug(" get_statues")
        stat = driver.get_status()
        t0 = time.perf_counter()
        if stat == Status.Up | Status.NeedsInit:
            Z, qm_classic_charge = driver.initialise(Z,A2)
            Z[Z==0] = 88 #special for Tim
            if Z[2] == 17:
                Z[2] = Z[2] + 18 #special for SN2
            logger.info(qm_classic_charge)
            logger.info(Z)
            logger.info('Driver initialised.')

        elif stat == Status.Up | Status.HasData:
            logger.debug(" get data")
            coord = driver.get_data(A).T
            T = {}
            T0 = driver.get_data(A2,send_signal=False,need_dataready=False).T
            T1 = driver.get_data(A,send_signal=False,need_dataready=False).T
            if max_order >= 1:
                T2 = driver.get_data(A3,send_signal=False,need_dataready=False).transpose(2,1,0)
                logger.debug(" get T2")
            if max_order >= 2:
                T3 = driver.get_data(A4,send_signal=False,need_dataready=False).transpose(3,2,1,0)
                logger.debug(" get T3")
            if max_order >= 3:
                T4 = driver.get_data(A5,send_signal=False,need_dataready=False).transpose(4,3,2,1,0)
                logger.debug(" get T4")

            t2 = time.perf_counter()
            logger.debug(" perf_counter")

            if rank_local == 0:
                if local_size > 1:
                    logger.debug('gather')
                    coord = comm_split.gather(coord)
                    T0 = comm_split.gather(T0)
                    T1 = comm_split.gather(T1)
                    if max_order >= 1:
                        T2 = comm_split.gather(T2)
                    if max_order >= 2:
                        T3 = comm_split.gather(T3)
                    if max_order >= 3:
                        T4 = comm_split.gather(T4)

                logger.debug("create input ")
                T['T0'] = T0
                T['T1'] = T1
                if max_order >= 1:
                    T['T2'] = T2
                if max_order >= 2:
                    T['T3'] = T3
                if max_order >= 3:
                    T['T4'] = T4

                input_dict = {}
                input_dict = T
                input_dict['xyz'] = coord
                #write_data('input.h5',input_dict)

                mldriver._create_input(coord, Z, T, qm_classic_charge, dtype=dtype)
            
                logger.debug("calculation")
                results = mldriver._run_calculation(dtype)
                
                force_norm = np.linalg.norm(results['qm_force'].reshape([-1]))
                logger.info(f"Force Norm: {force_norm}")
                logger.debug("calculation done")
                t_cal = t_cal +  time.perf_counter() - t2
            else:
                results = {}
                results['qm_force'] = None
                results['qm_energy'] = None
                results['qm_charge'] = None
                results['T0_gradient'] = None
                results['T1_gradient'] = None
                if max_order > 1:
                    results['T2_gradient'] = None
                if max_order > 2:
                    results['T3_gradient'] = None
            
            if local_size > 1:
                logger.debug('scatter')
                results['qm_force'] = comm_split.scatter(results['qm_force'],root=0)
                results['qm_energy'] = comm_split.scatter(results['qm_energy'],root=0)
                results['qm_charge'] = comm_split.scatter(results['qm_charge'],root=0)
                results['T0_gradient'] = comm_split.scatter(results['T0_gradient'],root=0)
                results['T1_gradient'] = comm_split.scatter(results['T1_gradient'],root=0)
                if max_order > 1:
                    results['T2_gradient'] = comm_split.scatter(results['T2_gradient'],root=0)
                if max_order > 2:
                    results['T3_gradient'] = comm_split.scatter(results['T3_gradient'],root=0)

            #write_data('server_results.h5',results)

            accepted_connections += 1
            has_data = True

        elif stat == Status.Up | Status.Ready and has_data:
            logger.debug("send back results")
            driver.send_data(results['qm_energy'],change_status=False)
            driver.send_data(results['qm_force'].T,send_signal=False,change_status=False)
            driver.send_data(results['qm_charge'],send_signal=False,change_status=False)
            driver.send_data(results['T0_gradient'],send_signal=False,change_status=False)

            if max_order > 1:
                driver.send_data(results['T1_gradient'].T,send_signal=False,change_status=False)
                if max_order > 2:
                    driver.send_data(results['T2_gradient'].transpose(2,1,0),send_signal=False,change_status=False)
                else:
                    driver.send_data(results['T2_gradient'].transpose(2,1,0),send_signal=False)
            else:
                driver.send_data(results['T1_gradient'].T,send_signal=False)

            if max_order > 2:
                driver.send_data(results['T3_gradient'].transpose(3,2,1,0),send_signal=False)
            logger.debug("send back results done")
            has_data = False
            if accepted_connections == args.max_connections:
                break
        t_all = t_all +  time.perf_counter() - t0
    logger.info(f"time of ML calculation: {t_cal}")
    logger.info(f"time of data transfer: {t_all - t_cal}")
    logger.info(f"total time: {t_all}")
    server.close()
