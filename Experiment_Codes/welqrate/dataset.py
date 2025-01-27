from torch_geometric.data import InMemoryDataset
from torch_geometric.data import Dataset, Data
import torch
import os
import os.path as osp
from tqdm import tqdm
from rdkit import Chem, RDLogger
from argparse import ArgumentParser
import pandas as pd
from welqrate.mol_utils.preprocess import (smiles2graph, inchi2graph, sdffile2mol_conformer,
                                  mol_conformer2graph3d)
from glob import glob
import numpy as np
import random
import zipfile
import subprocess


class WelQrateDataset(InMemoryDataset):
    def __init__(self, dataset_name, root ='dataset', mol_repr ='2dmol', task_type='classification'):
        if task_type not in ['classification', 'regression']:
            raise ValueError("task_type must be either 'classification' or 'regression'")
        
        self.name = dataset_name
        self.root = osp.join(root, dataset_name)
        self.task_type = task_type
        self.mol_repr = mol_repr
        
        print(f'dataset stored in {self.root}')
        
        if self.task_type == 'regression':
            raise ValueError('Regression task is not supported yet')

        # if self.name not in ['AID2689', 'AID488997', 'AID435008'] and self.task_type == 'regression':
        #     raise ValueError(f'{self.name} is not a regression dataset')
        
        super(WelQrateDataset, self).__init__(self.root)
        
        self.data, self.slices = torch.load(self.processed_paths[0])

        y = self.y.squeeze()
        self.num_active = len([x for x in y if x == 1])
        self.num_inactive = len([x for x in y if x == 0])

        print(f"Dataset {self.name} loaded.")
        print(f"Number of active molecules: {self.num_active}")
        print(f"Number of inactive molecules: {self.num_inactive}")

    @property
    def dataset_info(self):
        return {
                "AID1798": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/cd2dpdinu8grvi8dye3gi9bzsdt659d1.zip",
                    "file_type": "zip",
                    "raw_files": ['AID1798_actives.csv', 'AID1798_inactives.csv', 'AID1798_actives.sdf', 'AID1798_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/68m9qigxd7kt0xtta3chrx270grd1l4p.zip",
                    "split_file_type": "zip"
                },
                "AID1843": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/4f9wl7t9pmkm5p6695owj6bsxin0hdo1.zip",
                    "file_type": "zip",
                    "raw_files": ['AID1843_actives.csv', 'AID1843_inactives.csv', 'AID1843_actives.sdf', 'AID1843_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/o3hmvwo0surg1vhettflhtxxgn2o6clw.zip",
                    "split_file_type": "zip"
                },
                "AID2258": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/b1cg619c4f35p9mkm42mgu8a9m86t4at.zip",
                    "file_type": "zip",
                    "raw_files": ['AID2258_actives.csv', 'AID2258_inactives.csv', 'AID2258_actives.sdf', 'AID2258_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/2hhz7o3vxwjq6p93yzxfrporeze5zzgf.zip",
                    "split_file_type": "zip"
                },
                "AID2689": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/5no4ut1vusxmxsbwxsmb68o6ix94mxgy.zip",
                    "file_type": "zip",
                    "raw_files": ['AID2689_actives.csv', 'AID2689_inactives.csv', 'AID2689_actives.sdf', 'AID2689_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/hml9qznmwhhb6rcddw18vp7yoh331wlb.zip",
                    "split_file_type": "zip"
                },
                "AID435008": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/mu6hfitlenanp11wgm3z7u5drae2ofsd.zip",
                    "file_type": "zip",
                    "raw_files": ['AID435008_actives.csv', 'AID435008_inactives.csv', 'AID435008_actives.sdf', 'AID435008_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/ybyorbb541gnefjvagbqg1p37jfae98z.zip",
                    "split_file_type": "zip"
                },
                "AID435034": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/ctk2vs70bsjmpznalqoj4uqbbdkdcbgj.zip",
                    "file_type": "zip",
                    "raw_files": ['AID435034_actives.csv', 'AID435034_inactives.csv', 'AID435034_actives.sdf', 'AID435034_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/onistu6vrizl4o2nu0xlfhowcf1q498a.zip",
                    "split_file_type": "zip"
                },
                "AID463087": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/y52n3pzf27ghqtxt0xz899gwuj3yw3m2.zip",
                    "file_type": "zip",
                    "raw_files": ['AID463087_actives.csv', 'AID463087_inactives.csv', 'AID463087_actives.sdf', 'AID463087_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/y9jkhhu4cljq4awlfkvlo4im6iz3kj5t.zip",
                    "split_file_type": "zip"
                },
                "AID485290": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/8w0njtoqmgs8g1c12qj0p52w0d3wwbop.zip",
                    "file_type": "zip",
                    "raw_files": ['AID485290_actives.csv', 'AID485290_inactives.csv', 'AID485290_actives.sdf', 'AID485290_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/2w5byzpf5pk4jqb9rknzdihzztb3u9e6.zip",
                    "split_file_type": "zip"
                },
                "AID488997": {
                    "raw_url": "https://vanderbilt.box.com/shared/static/rguy1gm36x7clq6822riconznoflc4xe.zip",
                    "file_type": "zip",
                    "raw_files": ['AID488997_actives.csv', 'AID488997_inactives.csv', 'AID488997_actives.sdf', 'AID488997_inactives.sdf'],
                    "split_url": "https://vanderbilt.box.com/shared/static/dceunwrlotxkz58usfgnxr7meziy3cdr.zip",
                    "split_file_type": "zip"
                }
                }   


    @property
    def raw_file_names(self):
        # List files in the raw directory that match the dataset pattern
        pattern_csv = os.path.join(self.raw_dir, f'raw_{self.name}_*.csv')
        pattern_sdf = os.path.join(self.raw_dir, f'raw_{self.name}_*.sdf')
        files_csv = glob(pattern_csv)
        files_sdf = glob(pattern_sdf)
        return [os.path.basename(f) for f in files_csv + files_sdf]
    
    # need to change inchi/smiles to 2dmol
    @property
    def processed_file_names(self):
        return [f'processed_{self.mol_repr}_{self.name}.pt'] # example: processed_2dmol_AID1798.pt, processed_3dmol_AID1798.pt
    
    @property
    def processed_dir(self):
        return osp.join(self.root, 'processed', self.mol_repr) # example: processed/2dmol

    def download(self):
        os.makedirs(osp.join(self.root, 'raw'), exist_ok=True)
        os.makedirs(osp.join(self.root, 'split'), exist_ok=True)
        
        dataset_info = self.dataset_info.get(self.name)
        if not dataset_info:
            print(f"Download information for {self.name} is not available.")
            return

        # Download and extract raw files
        if not all(osp.exists(osp.join(self.root, 'raw', f)) for f in dataset_info['raw_files']):
            print(f"Downloading {self.name} raw files...")
            raw_zip_path = osp.join(self.root, 'raw', f'{self.name}.zip')
            subprocess.run(['curl', '-L', dataset_info['raw_url'], '--output', raw_zip_path], check=True)
            
            print("Extracting raw files...")
            subprocess.run(['unzip', '-j', raw_zip_path, '-d', osp.join(self.root, 'raw')], check=True)
            
            print("Removing raw zip file...")
            os.remove(raw_zip_path)

        # Download and extract split files
        if not os.listdir(osp.join(self.root, 'split')):  # Check if split directory is empty
            print(f"Downloading {self.name} split files...")
            split_zip_path = osp.join(self.root, 'split', f'{self.name}_split.zip')
            subprocess.run(['curl', '-L', dataset_info['split_url'], '--output', split_zip_path], check=True)
            
            print("Extracting split files...")
            subprocess.run(['unzip', split_zip_path, '-d', osp.join(self.root, 'split')], check=True)
            
            print("Removing split zip file...")
            os.remove(split_zip_path)


    def process(self):

        print(f'molecule representation: {self.mol_repr}')
        if self.mol_repr == '2dmol':   # combine smiles and inchi to 2dmol later
            self.file_type = '.csv'
        elif self.mol_repr == '3dmol':
            self.file_type = '.sdf'
        print(f'processing dataset {self.name}')
        RDLogger.DisableLog('rdApp.*')
        
        data_list = []

        mol_id = 0
        
        for file_name, label in [(f'{self.name}_actives{self.file_type}', 1),
                                 (f'{self.name}_inactives{self.file_type}', 0)]:

            source_path = os.path.join(self.root, 'raw', file_name)
            print(f'loaded raw file from {source_path}')
            
            if self.mol_repr == '2dmol':
                inchi_list = pd.read_csv(source_path, sep=',')['InChI'].tolist()
                cid_list = pd.read_csv(source_path, sep=',')['CID'].tolist() 
                smiles_list = pd.read_csv(source_path, sep=',')['SMILES'].tolist()

                for i, mol in tqdm(enumerate(inchi_list), total = len(inchi_list)):
                    pyg_data = inchi2graph(mol)

                    pyg_data.y = torch.tensor([label], dtype=torch.int) 

                    pyg_data.pubchem_cid = torch.tensor([int(cid_list[i])], dtype=torch.int)
                    pyg_data.mol_id = torch.tensor([mol_id], dtype=torch.int)  # index of the molecule in the dataset
                    pyg_data.smiles = smiles_list[i]
                    pyg_data.mol_id = torch.tensor([mol_id], dtype=torch.int)
                    data_list.append(pyg_data)
                    mol_id += 1
                    
            elif self.mol_repr == '3dmol':
                
                df = pd.read_csv(f'{source_path[:-4]}.csv')
                smiles_dict = df.set_index('CID')['SMILES'].to_dict()
                inchi_dict = df.set_index('CID')['InChI'].to_dict()

                mol_conformer_list, cid_list = sdffile2mol_conformer(source_path)
                # get the smiles/inchi from the csv file with CID
                
                for i, (mol, conformer) in tqdm(enumerate(mol_conformer_list), total=len(mol_conformer_list)):
                    pyg_data = mol_conformer2graph3d(mol, conformer)

                    pyg_data.pubchem_cid = torch.tensor([int(cid_list[i])], dtype=torch.int)
                    pyg_data.y = torch.tensor([label], dtype=torch.int)
                    pyg_data.smiles = smiles_dict[int(cid_list[i])]
                    pyg_data.inchi = inchi_dict[int(cid_list[i])]
     
                    pyg_data.mol_id = torch.tensor([mol_id], dtype=torch.int)
                    data_list.append(pyg_data)
                    mol_id += 1
        
        
        data, slices = self.collate(data_list)
        processed_file_path = os.path.join(self.processed_dir, f'processed_{self.mol_repr}_{self.name}.pt')
        torch.save((data, slices), processed_file_path)

        
    def get_idx_split(self, split_scheme = 'random_cv1'):
        
        valid_schemes = [f'random_cv{i}' for i in range(1, 6)] + [f'scaffold_seed{i}' for i in range(1, 6)]
        if split_scheme not in valid_schemes:
            raise ValueError(f"Invalid split_scheme. Must be one of {valid_schemes}")
        
        num = int(split_scheme[-1])
        path = osp.join(self.root, 'split')
        try: 
            if 'random' in split_scheme:
                print(f'loading random split cv{num}')
                split_scheme = f'{self.name}_{self.mol_repr[:2]}_random_cv{num}'
                split_dict = torch.load(osp.join(path, 'random', f'{split_scheme}.pt'))
            elif 'scaffold' in split_scheme:
                print(f'loading scaffold split seed{num}')
                split_scheme = f'{self.name}_{self.mol_repr[:2]}_scaffold_seed{num}'
                split_dict = torch.load(osp.join(path, 'scaffold', f'{split_scheme}.pt'))

        except Exception as e:
            print(f'split file not found. Error msg: {e}')

        # count the number of active molecules in the split
        y = self.y.squeeze()
        num_active_train = len([x for x in y[split_dict['train']] if x == 1])
        print(f'train set: {num_active_train} actives and {len(split_dict["train"]) - num_active_train} inactives')
        num_active_valid = len([x for x in y[split_dict['valid']] if x == 1])
        print(f'valid set: {num_active_valid} actives and {len(split_dict["valid"]) - num_active_valid} inactives')
        num_active_test = len([x for x in y[split_dict['test']] if x == 1])
        print(f'test set: {num_active_test} actives and {len(split_dict["test"]) - num_active_test} inactives')


        
        return split_dict



