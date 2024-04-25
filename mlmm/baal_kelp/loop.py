# Copyright The PyTorch Lightning team.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from copy import deepcopy
from typing import Any, Dict, Optional
import os

import torch
from mlmm.baal.get_label import CONV_UNIT_ENE, CONV_UNIT_FORCE
from pytorch_lightning import LightningModule
from pytorch_lightning.loops import Loop
from pytorch_lightning.loops.dataloader import PredictionLoop
#from pytorch_lightning.loops.fit_loop import FitLoop
from mlmm.baal.fit_loop import FitLoop
from pytorch_lightning.trainer.progress import Progress
from pytorch_lightning.trainer.states import TrainerFn, TrainerStatus
from pytorch_lightning.utilities.exceptions import MisconfigurationException
from pytorch_lightning.utilities.model_helpers import is_overridden

import flash
from flash.core.data.utils import _STAGES_PREFIX
from flash.core.utilities.imports import _PL_GREATER_EQUAL_1_5_0, requires
from flash.core.utilities.stages import RunningStage
#from flash.image.classification.integrations.baal.data import ActiveLearningDataModule
from mlmm.baal.data import ActiveLearningDataModule
from mlmm.baal.dropout import InferenceMCDropoutTask
from mlmm.baal.sampler import OptimizerSampler, SamplerBase
import torch.utils.data as torchdata
# from flash.image.classification.integrations.baal.dropout import InferenceMCDropoutTask

if not _PL_GREATER_EQUAL_1_5_0:
    from pytorch_lightning.trainer.connectors.data_connector import _PatchDataLoader
else:
    from pytorch_lightning.trainer.connectors.data_connector import _DataLoaderSource


class ActiveLearningLoop(Loop):
    @requires("baal")
    def __init__(self, label_epoch_frequency: int, inference_iteration: int = 2, should_reset_weights: bool = True, mode: str = "database"):
        """The `ActiveLearning Loop` describes the following training procedure. This loop is connected with the
        `ActiveLearningTrainer`

        Example::

            while unlabelled data or budget critera not reached:

                if labelled data
                    trainer.fit(model, labelled data)

                if unlabelled data:
                    predictions = trainer.predict(model, unlabelled data)
                    uncertainties = heuristic(predictions)
                    request labellelisation for the sample with highest uncertainties under a given budget

        Args:
            label_epoch_frequency: Number of epoch to train on before requesting labellisation.
            inference_iteration: Number of inference to perform to compute uncertainty.
        """
        super().__init__()
        self.label_epoch_frequency = label_epoch_frequency
        self.inference_iteration = inference_iteration
        self.should_reset_weights = should_reset_weights
        self.fit_loop: Optional[FitLoop] = None
        self.progress = Progress()
        self._model_state_dict: Optional[Dict[str, torch.Tensor]] = None
        self._lightning_module: Optional[flash.Task] = None
        self.mode = mode

    @property
    def done(self) -> bool:
        return self.progress.current.completed >= self.max_epochs

    def connect(self, fit_loop: FitLoop):
        self.fit_loop = fit_loop
        self.max_epochs = self.fit_loop.max_epochs
        self.fit_loop.max_epochs = self.label_epoch_frequency #self.label_epoch_frequency 

    def on_run_start(self, *args: Any, **kwargs: Any) -> None:
        assert isinstance(self.trainer.datamodule, ActiveLearningDataModule)
        # self.trainer.lightning_module.std_e= self.trainer.datamodule._dataset._dataset.std_e
        if self.trainer.lightning_module.std_e == 1.0 and self.trainer.lightning_module.mean_e == 0.0:
            torch.nn.init.constant_(self.trainer.lightning_module.std_e, self.trainer.datamodule._dataset._dataset.std_e)
            torch.nn.init.constant_(self.trainer.lightning_module.mean_e , self.trainer.datamodule._dataset._dataset.mean_e)
        # self.trainer.lightning_module.mean_e = self.trainer.datamodule._dataset._dataset.mean_e
        self.trainer.predict_loop.energy = self.trainer.datamodule._dataset._dataset.label_data['qm_energy']
        self.trainer.predict_loop._return_predictions = True
        self._lightning_module = self.trainer.lightning_module
        self._model_state_dict = deepcopy(self._lightning_module.state_dict())
        self.inference_model = InferenceMCDropoutTask(self._lightning_module, self.inference_iteration)

    def reset(self) -> None:
        pass

    def on_advance_start(self, *args: Any, **kwargs: Any) -> None:
        if self.trainer.datamodule.has_labelled_data:
            self._reset_dataloader_for_stage(RunningStage.TRAINING)
            self._reset_dataloader_for_stage(RunningStage.VALIDATING)
            if self.trainer.datamodule.has_test:
                self._reset_dataloader_for_stage(RunningStage.TESTING)
        if self.trainer.datamodule.has_unlabelled_data or self.mode == 'generative':
            self._reset_dataloader_for_stage(RunningStage.PREDICTING)
        self.progress.increment_ready()

    def advance(self, *args: Any, **kwargs: Any) -> None:

        self.progress.increment_started()

        if self.trainer.datamodule.has_labelled_data:
            self.fit_loop.run()

        if self.trainer.datamodule.has_test:
            self._reset_testing()
            metrics = self.trainer.test_loop.run()
            if metrics:
                self.trainer.logger.log_metrics(metrics[0], step=self.trainer.global_step)

        if self.trainer.datamodule.has_unlabelled_data and self.mode == 'database':
            self._reset_predicting()
            probabilities = self.trainer.predict_loop.run()
            self.trainer.datamodule.label(probabilities=probabilities)
        elif self.mode == "generative":
            self._reset_predicting()
            samples = self.trainer.predict_loop.run()
            ref_value = self.trainer.datamodule.label(sample=samples)
            #self.trainer.predict_loop.energy = torch.cat([self.trainer.predict_loop.energy,ref_value['qm_energy']])
            #self.trainer.model.parent_module.std_e = torch.std(self.trainer.predict_loop.energy).to(samples.device)
            #self.trainer.model.parent_module.mean_e = torch.mean(self.trainer.predict_loop.energy).to(samples.device)

            predict_energy_error = abs(self.trainer.predict_loop.myresults['qm_energy'] - ref_value['qm_energy']).mean()*CONV_UNIT_ENE
            predict_force_error = abs(self.trainer.predict_loop.myresults['qm_force'] - ref_value['qm_force']).mean()*CONV_UNIT_FORCE
            print(f"predict energy error: {predict_energy_error}")
            print(f"predict force error: {predict_force_error}")
            a=3
        else:
            raise StopIteration

        self._reset_fitting()
        self.progress.increment_processed()

    def on_advance_end(self) -> None:
        # if self.trainer.datamodule.has_unlabelled_data and self.should_reset_weights:
            # reload the weights to retrain from scratch with the new labelled data.
            # self._lightning_module.load_state_dict(self._model_state_dict)
        filename = 'active_epoch='+str(self.progress.current.completed)
        filename = os.path.join(self.trainer.log_dir,filename)
        if self.mode == 'generative':
            self.trainer.datamodule._dataset._dataset.datasets[-1].save(path=filename)
        self.progress.increment_completed()
        return super().on_advance_end()

    def on_run_end(self):
        self._reset_fitting()
        self._teardown()
        return super().on_run_end()

    def on_save_checkpoint(self) -> Dict:
        return {"datamodule_state_dict": self.trainer.datamodule.state_dict()}

    def on_load_checkpoint(self, state_dict) -> None:
        self.trainer.datamodule.load_state_dict(state_dict.pop("datamodule_state_dict"))

    def __getattr__(self, key):
        if key not in self.__dict__:
            return getattr(self.fit_loop, key)
        return self.__dict__[key]

    def _connect(self, model: LightningModule):
        if _PL_GREATER_EQUAL_1_5_0:
            self.trainer.training_type_plugin.connect(model)
        else:
            self.trainer.accelerator.connect(model)

    def _reset_fitting(self):
        self.trainer.state.fn = TrainerFn.FITTING
        self.trainer.training = True
        self.trainer.lightning_module.on_train_dataloader()
        self._connect(self._lightning_module)
        self.fit_loop.epoch_progress.current.completed +=1
        self.fit_loop.max_epochs = self.fit_loop.epoch_progress.current.completed + self.label_epoch_frequency #self.label_epoch_frequency 
        # self.fit_loop.epoch_progress.current.completed = self.progress.current.ready
        self.global_step = self.fit_loop.epoch_loop.global_step
        # self.fit_loop.current_epoch = self.progress.current.ready
        # self.trainer.max_epochs = self.progress.current.ready + self.label_epoch_frequency
        # self.fit_loop.epoch_progress = Progress()
        # self.fit_loop.epoch_progress = self.progress #Progress()

    def _reset_predicting(self):
        self.trainer.state.fn = TrainerFn.PREDICTING
        self.trainer.predicting = True
        self.trainer.lightning_module.on_predict_dataloader()
        self._connect(self.inference_model)

    def _reset_testing(self):
        self.trainer.state.fn = TrainerFn.TESTING
        self.trainer.state.status = TrainerStatus.RUNNING
        self.trainer.testing = True
        self.trainer.lightning_module.on_test_dataloader()
        self._connect(self._lightning_module)

    def _reset_dataloader_for_stage(self, running_state: RunningStage):
        dataloader_name = f"{_STAGES_PREFIX[running_state]}_dataloader"
        # If the dataloader exists, we reset it.
        dataloader = (
            getattr(self.trainer.datamodule, dataloader_name)
            if is_overridden(dataloader_name, self.trainer.datamodule)
            else None
        )
        if dataloader:
            if _PL_GREATER_EQUAL_1_5_0:
                setattr(
                    self.trainer._data_connector,
                    f"_{dataloader_name}_source",
                    _DataLoaderSource(self.trainer.datamodule, dataloader_name),
                )
            else:
                setattr(
                    self.trainer.lightning_module,
                    dataloader_name,
                    _PatchDataLoader(dataloader(), running_state),
                )
            setattr(self.trainer, dataloader_name, None)
            # TODO: Resolve this within PyTorch Lightning.
            try:
                getattr(self.trainer, f"reset_{dataloader_name}")(self.trainer.lightning_module)
            except MisconfigurationException:
                pass

    def _teardown(self) -> None:
        self.trainer.train_dataloader = None
        self.trainer.val_dataloaders = None
        self.trainer.test_dataloaders = None
        self.trainer.predict_dataloaders = None
        # Hack
        self.trainer.lightning_module.train_dataloader = None
        self.trainer.lightning_module.val_dataloader = None
        self.trainer.lightning_module.test_dataloader = None
        self.trainer.lightning_module.predict_dataloader = None


class MyPredictLoop(PredictionLoop):
    def __init__(self, sampler:SamplerBase):
        super().__init__()
        self.min = 1e20
        self.sampler = sampler
        self.myresults = None
        self.energy = None
    
    def on_run_start(self) -> None:
        self.sampler.setup(self.trainer.model) 
        return super().on_run_start()

    def advance(self):
        dataloader = self.trainer.training_type_plugin.process_dataloader(self.current_dataloader)

        if hasattr(dataloader.dataset._dataset.dataset, 'datasets'):
            # minene = dataloader.dataset._dataset.dataset.datasets[-1].total_energy().min()
            minene = dataloader.dataset._dataset.dataset.datasets[-1].label_data['qm_energy'].min()
        else:
            # minene = self.trainer.datamodule._dataset._dataset.total_energy().min()
            minene = self.trainer.datamodule._dataset._dataset.label_data['qm_energy'].min()
        
        if minene<self.min:
            self.min = minene
        print(f"reference energy:{self.min}")

        dataloader_iter = iter(dataloader)
        batch = next(dataloader_iter)
        dl_max_batches = self.max_batches[self.current_dataloader_idx]
        with self.trainer.profiler.profile("predict_batch_to_device"):
            batch = self.trainer.accelerator.batch_to_device(batch, dataloader_idx=self.current_dataloader_idx)

        self.sampler.ref = self.min
        # self.sampler.alpha = self.alpha
        samples, myresults = self.sampler.run(batch)
        #self.energy = torch.cat([self.energy,myresults['qm_energy']])
        #self.trainer.model.parent_module.std_e = torch.std(self.energy).to(samples.device)
        #self.trainer.model.parent_module.mean_e = torch.mean(self.energy).to(samples.device)
        self.predictions.append(samples)
        self.myresults = myresults
        # self.alpha = self.sampler.alpha
        a=3

