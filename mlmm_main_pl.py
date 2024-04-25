import pytorch_lightning as pl
from mlmm.model.model_pl import LitMLMM
from mlmm.dataset import Molecule_DataModule
from pytorch_lightning.utilities.cli import LightningCLI
import torch

class MyLightningCLI(LightningCLI):
    def add_arguments_to_parser(self, parser):
        parser.link_arguments("model.rep_model.init_args.n_atom_basis", "model.out_model.init_args.n_atom_basis")
        #parser.link_arguments("model.rep_model.init_args.cutoff", "data.Dataset.init_args.cutoff")
        parser.add_argument("--save_model_pt", default=False)

    def before_fit(self):
        if self.config['fit']['save_model_pt']:
            st = torch.load(self.trainer.resume_from_checkpoint)['state_dict'] 
            st['std_e'] = st['std_e'][None]
            st['mean_e'] = st['mean_e'][None]
            self.model.load_state_dict(st)
            torch.save(self.model,'model.pt')
            exit()
    
    def before_validate(self):
        self.model.std_e = self.datamodule.dataset.std_e #0.0334
        self.model.mean_e = self.datamodule.dataset.mean_e #-10.3949

if __name__=="__main__":
    cli = MyLightningCLI(LitMLMM, Molecule_DataModule)#, run=False)#, subclass_mode_model=True)
    # cli.trainer.fit(cli.model, cli.data)
