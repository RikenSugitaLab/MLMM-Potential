from dgl.data.dgl_dataset import DGLDataset
from pytorch_lightning import LightningDataModule
from mlmm.dataset.mydataset_alphachem import  collate
from torch.utils.data import Dataset, DataLoader
from typing import Dict, Union, List, Optional
import dgl
__all__=['Molecule_DataModule']

class Molecule_DataModule(LightningDataModule):
    """molecular dataModule

    Args:
        Dataset(Dataset or DGLDataset) : Dataset Class storing all the data
        batch_size(int): size of mini-batch
        frac_list (list): fraction of samples of training, validating dataset
    """
    def __init__(self, 
                 Dataset: Dataset, #Union[Dataset, DGLDataset],
                 batch_size: int,
                 frac_list: list,
                 num_workers: int = 0,
                 pin_memory: bool = False,
                 persistent_workers: bool = False):
        super().__init__()
        self.dataset = Dataset
        self.batch_size = batch_size
        self.frac_list = frac_list
        self.num_workers = num_workers
        self.pin_memory = pin_memory
        self.persistent_workers = persistent_workers

    def setup(self,stage: Optional[str] = None):
        self.dataset.load_data()
        self.dataset.ckpt_path = self.trainer.log_dir
        if self.trainer.resume_from_checkpoint is None:
            self.trainer.model.std_e = self.dataset.std_e
            self.trainer.model.mean_e = self.dataset.mean_e
        # datalist = dgl.data.utils.split_dataset(self.dataset,
                                    # frac_list = self.frac_list,
                                    # shuffle = True)

        if stage == "fit" or stage is None: 
            self.dataset_train = self.dataset #datalist[0]

        elif stage == 'predict' or stage is None:
            self.dataset_predict = self.dataset

        elif stage == 'validate' or stage is None:
            self.dataset_validate = self.dataset
        return

    def train_dataloader(self):
        return DataLoader(self.dataset_train, batch_size=self.batch_size,collate_fn=collate,num_workers=self.num_workers,pin_memory=self.pin_memory,persistent_workers=self.persistent_workers,shuffle=True)
    
    def val_dataloader(self):
        return DataLoader(self.dataset, batch_size=self.batch_size,collate_fn=collate,num_workers=self.num_workers,pin_memory=self.pin_memory,persistent_workers=self.persistent_workers)

    def predict_dataloader(self):
        return DataLoader(self.dataset, batch_size=self.batch_size,collate_fn=collate,num_workers=self.num_workers,pin_memory=self.pin_memory,persistent_workers=self.persistent_workers)

