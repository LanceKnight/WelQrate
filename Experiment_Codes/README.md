# WelQrate: Defining the Gold Standard in Small Molecule Drug Discovery


## Installation
We provide the recommended environment, which were used for benchmarking in the original paper. Users can also build their own environment based on their own needs.
```
conda create -n welqrate python=3.9
```

```
conda install pytorch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 pytorch-cuda=12.1 -c pytorch -c nvidia
```

```
pip install torch_geometric==2.3.1
```

```
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
```


```
pip install -r requirements.txt
```

```
pip install -e .
```


## Load the Dataset
Users can download and preprocess the datasets by calling `WelQrateDataset` class. Available datasets include AID1798, AID435008, AID435034, AID1843, AID2258, AID463087, AID488997, AID2689, and AID485290. Please refer to our [website](https://www.welqrate.org/) for more details. Besides, users can choose between 2D and 3D molecular representations by setting `mol_repr` to `2dmol` or `3dmol`.

```python
from welqrate.dataset import WelQrateDataset
# Load the 2D dataset
AID1798_dataset_2d = WelQrateDataset(dataset_name = 'AID1798', root =f'./datasets', mol_repr ='2dmol')

# Load the 3D dataset 
AID1843_dataset_3d = WelQrateDataset(dataset_name = 'AID1843', root =f'./datasets', mol_repr ='3dmol')

# Load a split dictionary
split_dict = AID1798_dataset_2d.get_idx_split(split_scheme ='random_cv1') # or 'scaffold_seed1; we provide 1-5 for both random_cv and scaffold_seed

```

## Train a model
We can store hyperparameters related to model, training scheme, and dataset in a configuration file. Users can refer to configuration files in `./config/` for different models. Then, we can config the model and start training by calling `train` function. After training, results are automatically saved in `./results/` folder.

```python
from welqrate.train import train
from welqrate.models.gnn2d.GCN import GCN_Model
import configparser
import torch

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
config = configparser.ConfigParser()
config.read(['./config/train.cfg', './config/gcn.cfg'])

# Load the model configuration
hidden_channels = int(config['MODEL']['hidden_channels'])
num_layers = int(config['MODEL']['num_layers'])
model = GCN_Model(hidden_channels=hidden_channels, num_layers=num_layers).to(device)

# Train the model
train(model, AID1798_dataset_2d, config, device)
```
Please check out `example.py` for more examples.

## Citation
If you find our work helpful, please cite our paper:


```       
@article{dong2024welqrate,
  title={WelQrate: Defining the Gold Standard in Small Molecule Drug Discovery Benchmarking},
  author={Yunchao, Liu and Dong, Ha and Wang, Xin and Moretti, Rocco and Wang, Yu and Su, Zhaoqian and Gu, Jiawei and Bodenheimer, Bobby and Weaver, Charles David and Meiler, Jens and Derr, Tyler and others},
  journal={arXiv preprint arXiv:2411.09820},
  year={2024}
}

```