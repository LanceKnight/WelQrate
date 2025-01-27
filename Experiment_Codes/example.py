from welqrate.dataset import WelQrateDataset
from welqrate.models.gnn2d.GCN import GCN_Model 
from welqrate.models.gnn3d.SchNet import SchNet_Model
import torch
from welqrate.train import train
import configparser


AID1798_dataset_2d = WelQrateDataset(dataset_name = 'AID1798', root =f'./datasets', mol_repr ='2dmol')

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

config = configparser.ConfigParser()

config.read(['./config/train.cfg', './config/gcn.cfg'])

hidden_channels = int(config['MODEL']['hidden_channels'])
num_layers = int(config['MODEL']['num_layers'])

model = GCN_Model(hidden_channels=hidden_channels, 
                  num_layers=num_layers).to(device)

train(model, AID1798_dataset_2d, config, device)

# 3D graph representation and SchNet model
# AID1843_dataset_3d = WelQrateDataset(dataset_name = 'AID1798', root =f'./datasets', mol_repr ='3dmol')
# config = configparser.ConfigParser()
# config.read(['./config/train.cfg', './config/schnet.cfg'])

# model = SchNet_Model(hidden_channels=int(config['MODEL']['hidden_channels']),
#                     out_channels=int(config['MODEL']['out_channels']),
#                     num_layers=int(config['MODEL']['num_layers']),
#                     energy_and_force=bool(config['MODEL']['energy_and_force']),
#                     cutoff=float(config['MODEL']['cutoff']),
#                     num_filters=int(config['MODEL']['num_filters']),
#                     num_gaussians=int(config['MODEL']['num_gaussians'])).to(device)

# train(model, AID1843_dataset_3d, config, device)


