import numpy as np
import adamod
import json
import sys
import os
from mlmm.trainer.trainer_alphachem import Trainer
import dgl
from mlmm.dataset.mydataset_alphachem import load_dataset
import torch
from torch import nn
from torch.utils.data import DataLoader
from mlmm.nn.activations import swish
from mlmm.nn.base import Dense
from mlmm.representation import SchNet, TDTNet, FieldNet
from mlmm.model.model_mlmm import model_mlmm
import argparse
from mlmm.utils.utils import set_seed, make_log_dir, moment_update, dict2namespace
from mlmm.opt.lr_scheduler import get_scheduler
from mlmm.opt.diffmod import DiffMod
import math
from torch.utils.tensorboard import SummaryWriter
from mlmm.utils.logger import setup_logger
import yaml
import shutil
#torch.set_default_dtype(torch.float64)

if "use_fp16" in os.environ:
    if os.environ['use_fp16'] == 'True':
        use_fp16=True
    elif os.environ['use_fp16'] == 'False':
        use_fp16=False
    else:
        raise AssertionError("wrong setting, use_fp16 can only be True or False")
else:
    use_fp16=False

def main(args1):
    args = dict2namespace(args1['Basic'])
    args_dataset = dict2namespace(args1['Dataset'])
    args_network = dict2namespace(args1['Network'])
    args_opt = dict2namespace(args1['Optimizer'])
    args_trainer = dict2namespace(args1['Trainer'])

    if args.use_fp16 != use_fp16:
        raise AssertionError("wrong setting, use_fp16 are not the same")

    if args.np > 0:
        gpu_id = int(math.fmod(args.local_rank,args.np))
    else:
        gpu_id = args.local_rank

    torch.cuda.set_device(gpu_id)
    set_seed(args)

    ### distributed training initiation ###
    if args.parallel:
        torch.distributed.init_process_group(
	    	'gloo',
            init_method='env://'
	    )
        args.local_rank = torch.distributed.get_rank()
    else:
        args.local_rank = 0

    ### log file ###
    log_dir,output_dir = make_log_dir(args)
    if args.local_rank>=0:
        if args.train:
            log_file = os.path.join(log_dir,"log_train.txt")
            logger = setup_logger(output=log_file,
                                  distributed_rank=args.local_rank, 
                                  name="rl",
                                  log_level=args.log_level)
        else:
            log_file = os.path.join(log_dir,"log_eval.txt")
            logger = setup_logger(output=log_file, 
                                  distributed_rank=args.local_rank, 
                                  name="rl",
                                  log_level=args.log_level)

    if gpu_id == 0:
       summary_writer = SummaryWriter(log_dir=os.path.join(output_dir,'summary'))
    else:
       summary_writer = None
    
    if not os.path.exists(os.path.join(output_dir,args.config)):
        shutil.copy(args.config,os.path.join(output_dir,'config.yml'))

    if args.precision == 'float32':
        torch.set_default_dtype(torch.float32)
        logger.info("set the precision: float32")
    elif args.precision == 'floa64':
        torch.set_default_dtype(torch.float64)
        logger.info("set the precision: float64")
    else:
        logger.info("defalt precision")

    logger.info("=> loading the data")

    ### loading the dataset, setting up the dataloader###
    train_loader= load_dataset(args_dataset)
    logger.info("=> done")

    ### model ###
    logger.info("=> building the model")
    logger.info(args.model_path)
    rep = FieldNet(
                   n_interactions=args_network.n_interactions,
                   n_atom_basis=args_network.n_atom_basis,
                   n_filters=args_network.n_filters,
                   n_gaussians=args_network.n_gaussians,
                   dipole_features=args_network.dipole_features,
                   cutoff=args_network.cutoff,
                   log_gaussian=args_network.log_gaussian,
                   graph_rep=args_network.graph_rep,
                   do_graph_norm=args_network.do_graph_norm
                   )

    model = model_mlmm(rep,n_hidden=args_network.n_hidden,
                       n_layers=args_network.n_layers,
                       mean_e=train_loader.dataset.mean_e,
                       std_e=train_loader.dataset.std_e)
    model.cuda(gpu_id)

    if args.parallel:
        logger.info('=>processing in parallel')
        model = torch.nn.parallel.DistributedDataParallel(model,device_ids=[gpu_id],find_unused_parameters=True)
    else:
        logger.info('=>processing in serial')

    logger.info("=> done")

    ### optimizer ###
    ### eps should be small or the training will bu unstable ###
    args_opt.lr = float(args_opt.lr)
    args_opt.eps = float(args_opt.eps)
    if args_opt.opt == "adamod":
        optimizer = adamod.AdaMod(model.parameters(),lr=args_opt.lr,beta3=0.999,eps=args_opt.eps)
    elif args_opt.opt == "adam":
        optimizer = torch.optim.Adam(model.parameters(),
                                     lr=args_opt.lr)#,
                                     #eps=args.eps,
                                     #amsgrad=True)
    elif args_opt.opt == "sgd":
        optimizer = torch.optim.SGD(model.parameters(),
                                    lr=args_opt.lr)#,
                                    #momentum=args.momentum)
    elif args_opt.opt == "diffmod":
        optimizer = DiffMod(model.parameters(), lr=args_opt.lr,eps=args_opt.eps) 


    #if args.lr_scheduler != "":
    #    scheduler = get_scheduler(optimizer,len(train_loader),args)
    #else:
    scheduler = None
    #scheduler = torch.optim.lr_scheduler.StepLR(optimizer,step_size=500,gamma=0.5)

    ### loss function ##

    trainer = Trainer( 
                        output_dir,
                        model,
                        optimizer,
                        train_loader,
                        args_trainer,
                        scheduler=scheduler,
                        summary_writer=summary_writer
                        )
    
    args_trainer.wq = float(args_trainer.wq)
    args_trainer.wf = float(args_trainer.wf)
    args_trainer.wl2 = float(args_trainer.wl2)
    if args.train:
        logger.info("=> start training")
        trainer.train(args_trainer)
    else:
        logger.info("=> start evaluating")
        trainer.eval(args_trainer)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    ### params ###
    parser.add_argument("--local_rank",type=int,default=0,
                        help="which GPU to use. If distributed training, this argument need not to be specified")
    parser.add_argument("--config",type=str,default="config_test.yml",help="the name of file containing commitor, reaction coordinate, tps when inferring")
    parser.add_argument("--datafile",type=str,default="/home/yklei/practice/mlmm_energy/md_test/test.dcd")

    args1 = parser.parse_args()
    f = open(args1.config)
    args2 = yaml.load(f)
    args2['Basic']['local_rank']=args1.local_rank
    args2['Basic']['config']=args1.config
    args2['Dataset']['datafile']=args1.datafile
    args2['Dataset']['cutoff']=args2['Network']['cutoff']
    args2['Dataset']['train']=args2['Basic']['train']
    args2['Dataset']['model_path']=args2['Basic']['model_path']
    args2['Dataset']['parallel']=args2['Basic']['parallel']
    args2['Trainer']['model_path']=args2['Basic']['model_path']
    args2['Trainer']['train']=args2['Basic']['train']
    args2['Trainer']['local_rank']=args1.local_rank
    args2['Trainer']['parallel']=args2['Basic']['parallel']
    args2['Trainer']['np']=args2['Basic']['np']
    args2['Trainer']['eval_size']=args2['Dataset']['eval_size']
    set_seed(dict2namespace(args2['Basic']))

    main(args2)

