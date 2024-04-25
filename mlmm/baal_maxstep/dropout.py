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
import torch
from pytorch_lightning.utilities.exceptions import MisconfigurationException

import flash
from flash.core.utilities.imports import _BAAL_AVAILABLE

from torch.autograd import grad
from collections import defaultdict


if _BAAL_AVAILABLE:
    from mlmm.baal.myheuristic import ForceVariance
    from baal.bayesian.dropout import _patch_dropout_layers

onlyforce = ['ForceVariance']

class InferenceMCDropoutTask(flash.Task):
    def __init__(self, module: flash.Task, inference_iteration: int):
        super().__init__()
        self.parent_module = module
        self.trainer = module.trainer
        changed = _patch_dropout_layers(self.parent_module)
        changed = True
        if not changed:
            raise MisconfigurationException("The model should contain at least 1 dropout layer.")
        self.inference_iteration = inference_iteration

    def predict_step(self, batch, batch_idx, dataloader_idx: int = 0):
        # with torch.no_grad():
        with torch.enable_grad():
            out = []
            feature = defaultdict(list)
            for _ in range(self.inference_iteration):
                results = self.parent_module.predict_step(batch, batch_idx, dataloader_idx=dataloader_idx)
                out2 = []
                for key in results.keys(): 
                    if key not in ['graph_feature','feature_vector']:
                        bz = results[key].shape[0]
                        out2.append(results[key].reshape(bz,-1))
                    else:
                        feature[key].append(results[key])
                out.append(torch.cat(out2,dim=-1))

            # out.append(self.parent_module.predict_step(batch, batch_idx, dataloader_idx=dataloader_idx)["preds"])
        feature['feature_vector'] = torch.stack(feature['feature_vector'],axis=-1).mean(-1)#.detach()
        if 'graph_feature' in results.keys():
            feature['graph_feature'] = torch.stack(feature['graph_feature'],axis=-1).mean(-1)
        else:
            feature['graph_feature'] = None

        # BaaL expects a shape [num_samples, num_classes, num_iterations]
        return torch.stack(out).permute((1, 2, 0)), feature
        # return torch.stack(out).permute((1, 2, 0)).detach().cpu()
