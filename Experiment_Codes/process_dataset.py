from dataset import WelQrateDataset
import torch
import random
import numpy as np
from loader import get_train_loader, get_valid_loader, get_test_loader
from torch_geometric.loader import DataLoader
from rdkit import Chem
from mol_utils.ChIRo_dataset import ChIRotDataset
from loader import get_train_loader, get_test_loader, get_valid_loader
from mol_utils.BCL_dataset import BCL_WelQrateDataset
from argparse import ArgumentParser
import configparser
import ast

from torch.utils.data import Dataset, DataLoader


seed = 0
torch.manual_seed(seed)
random.seed(seed)
np.random.seed(seed)

dataset_name = 'AID488997'
root = '../dataset'

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

### ----------------- Test the bcl_dataset ----------------- ###
bcl_data = BCL_WelQrateDataset(dataset_name, root = '../bcl_dataset')
print(len(bcl_data))
print(bcl_data[0].bcl.shape)

### ----------------- Test the argparser ----------------- ###
# parser = ArgumentParser()
# parser.add_argument('--config', type=str) 
# args = parser.parse_args()
# config_file = args.config
# config = configparser.ConfigParser()
# config.read(config_file)

# print(type(bool(config['DATA']['one_hot'])))
# print(type(ast.literal_eval(config['MODEL']['rnn_sizes'])[0]))
# print(type(ast.literal_eval(config['MODEL']['rnn_types'])))

### ----------------- Test the ChIRo Dataset ----------------- ###
# dataset_2d = WelQrateDataset(dataset_name, root, mol_repr='2dmol')

# dataset_chiro = ChIRotDataset(dataset_name, root)
# print(dataset_chiro[0])
# print(dataset_2d[0])

### ----------------- Test the Welqrate Dataset ----------------- ###

# dataset_2d = WelQrateDataset(dataset_name, root, mol_repr='2dmol')
# print(len(dataset_2d))


# dataset_3d = WelQrateDataset(dataset_name, root, mol_repr='3dmol')
# print(dataset_3d[100].pubchem_cid)
# print(dataset_3d[100])
# # 
# InChI=1S/C25H27N3O4S/c1-13-10-17(15(3)27(13)5)11-21-24(30)28-23(19-9-8-18(31-6)12-20(19)32-7)22(16(4)29)14(2)26-25(28)33-21/h8-12,23H,1-7H3/b21-11- 
# InChI=1S/C25H27N3O4S/c1-13-10-17(15(3)27(13)5)11-21-24(30)28-23(19-9-8-18(31-6)12-20(19)32-7)22(16(4)29)14(2)26-25(28)33-21/h8-12,23H,1-7H3/b21-11-

# get the pubchem_cid which y = 1
# well_active_cid = [data.pubchem_cid for data in dataset_2d if data.y == 1]
# well_cid_list = dataset_2d.pubchem_cid.tolist()

# print(well_cid_list[:10])

# split_scheme = 'AID1798_2d_random_cv1'
# split_dict = torch.load(f'../dataset/split/random/{split_scheme}.pt')

# max_seq = max([len(data.smiles) for data in dataset_2d])
# charset = set("".join(list(dataset_2d.smiles)))
# char_to_int = dict((c,i) for i,c in enumerate(charset))

# print(char_to_int)


# print(dataset_2d[0])
# max_edge_index = 0

# for data in dataset_2d:
#     data.validate(raise_on_error=True)

# print("Maximum edge_index in the dataset:", max_edge_index)

# train_loader = get_train_loader(dataset_2d[split_dict['train']], 1024, 4, 1)
# valid_loader = get_valid_loader(dataset_2d[split_dict['valid']], 1024, 4, 1)
# test_loader = get_test_loader(dataset_2d[split_dict['test']], 1024, 4, 1)
# train_loader = DataLoader(dataset_2d[split_dict['train']], batch_size=1024, shuffle=True, num_workers=4)


# from tqdm import tqdm
# for i in tqdm(range(50)):
#     for batch in train_loader:
#         batch = batch.to(device)
#         if batch.edge_index.max().item() > 500000:
#             print("Found a batch with edge_index greater than 500000")
#         if batch.edge_index.max().item() > batch.num_nodes:
#             print("Found a batch with edge_index greater than the maximum edge_index in the dataset")
#             break

### ---------------- Test the Poorly Curated Data ------------------###
# poor_dataset2d = WelQrateDataset('AID1798', '../poor_dataset', mol_repr='2dmol')
# poor_dataset3d = WelQrateDataset('AID1798', '../poor_dataset', mol_repr='3dmol')
# poor_active_cid = [data.pubchem_cid for data in poor_dataset2d if data.y == 1]

# false_positive = [cid.item() for cid in poor_active_cid if cid not in well_active_cid]

# false_positive = false_positive + [2882111]

# # check if the false positive is in the well curated dataset
# for cid in false_positive:
#     if cid in well_cid_list:
#         print(f"cid {cid} is in the well curated dataset")
#         break


# split_scheme = 'AID1798_2d_scaffold_seed5'

# split_dict = poor_dataset2d.get_idx_split(split_scheme)

# train_dataset = poor_dataset2d[split_dict['train']]
# valid_dataset = poor_dataset2d[split_dict['valid']]
# test_dataset = poor_dataset2d[split_dict['test']]

# num_train_actives = len([data for data in train_dataset if data.y == 1])
# num_train_inactives = len([data for data in train_dataset if data.y == 0])
# print("Number of actives in training set:", num_train_actives)
# print("Number of inactives in training set:", num_train_inactives)
# num_valid_actives = len([data for data in valid_dataset if data.y == 1])
# num_valid_inactives = len([data for data in valid_dataset if data.y == 0])
# print("Number of actives in validation set:", num_valid_actives)
# print("Number of inactives in validation set:", num_valid_inactives)
# num_test_actives = len([data for data in test_dataset if data.y == 1])
# num_test_inactives = len([data for data in test_dataset if data.y == 0])
# print("Number of actives in test set:", num_test_actives)
# print("Number of inactives in test set:", num_test_inactives)


### ----------------- Test the Regression Task ----------------- ###
# regression_dataset_2d = WelQrateDataset('AID488997', root, task_type='regression', mol_repr='2dmol')
# # regression_dataset_3d = WelQrateDataset('AID488997', root, task_type='regression', mol_repr='3dmol')

# print(regression_dataset_2d[0].activity_value.item())
# print(len(regression_dataset_2d))
# # print(len(regression_dataset_3d))

# # calculate the mean of the activity values below 999
# activity_values = [data.activity_value.item() for data in regression_dataset_2d if data.activity_value.item() < 999]
# print(np.mean(activity_values))

### ----------------- Test the Split ----------------- ###

# test the random split file of 3d dataset
# for i in range(1,6):
#     split_path = f'../dataset/split/random/AID1798_3d_random_cv{i}.pt'
#     split_dict = torch.load(split_path)
#     print(f'cv{i}: train: {len(split_dict["train"])}, valid: {len(split_dict["valid"])}, test: {len(split_dict["test"])}')

# # test the scaffold split file of 3d dataset
# for i in range(1,6):
#     split_path = f'../dataset/split/scaffold/AID1798_3d_scaffold_seed{i}.pt'
#     split_dict = torch.load(split_path)
#     print(f'cv{i}: train: {len(split_dict["train"])}, valid: {len(split_dict["valid"])}, test: {len(split_dict["test"])}')




## well-curated data -> well test index -> well test cid 
## poorly-curated data -> well test cid -> poorly test index -> poorly split




