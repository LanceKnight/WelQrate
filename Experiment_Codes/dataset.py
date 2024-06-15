from torch_geometric.data import InMemoryDataset
from torch_geometric.data import Dataset, Data
import torch
import os
import os.path as osp
from tqdm import tqdm
from rdkit import Chem, RDLogger
from argparse import ArgumentParser
import pandas as pd
from mol_utils.preprocess import (smiles2graph, inchi2graph, sdffile2mol_conformer,
                                  mol_conformer2graph3d)
from glob import glob
import numpy as np
import random


class WelQrateDataset(InMemoryDataset):
    def __init__(self, dataset_name, root ='dataset', mol_repr ='2dmol', task_type='classification'):
        if task_type not in ['classification', 'regression']:
            raise ValueError("task_type must be either 'classification' or 'regression'")
        
        self.name = dataset_name
        self.root = root
        self.task_type = task_type
        self.mol_repr = mol_repr
        
        print(f'dataset stored in {self.root}')
        
        if self.name not in ['AID2689', 'AID488997', 'AID435008'] and self.task_type == 'regression':
            raise ValueError(f'{self.name} is not a regression dataset')
        
        
        super(WelQrateDataset, self).__init__(self.root)
        
        self.data, self.slices = torch.load(self.processed_paths[0])

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
        # wait for the website to be up
        pass

    def process(self):

        print(f'molecule representation:{self.mol_repr}.')
        if self.mol_repr == '2dmol':   # combine smiles and inchi to 2dmol later
            self.file_type = '.csv'
        elif self.mol_repr == '3dmol':
            self.file_type = '.sdf'
        print(f'processing dataset {self.name}')
        RDLogger.DisableLog('rdApp.*')
        
        data_list = []
        invalid_id_list = []
        mol_id = 0
        
        for file_name, label in [(f'{self.name}_actives{self.file_type}', 1),
                                 (f'{self.name}_inactives{self.file_type}', 0)]:

            source_path = os.path.join(self.root, 'raw', file_name)
            print(f'loaded raw file from {source_path}')
            
            if self.mol_repr == '2dmol':
                inchi_list = pd.read_csv(source_path, sep=',')['InChI'].tolist()
                cid_list = pd.read_csv(source_path, sep=',')['CID'].tolist() # will change for the last version of the dataset
                smiles_list = pd.read_csv(source_path, sep=',')['SMILES'].tolist()
                if self.task_type == 'regression':
                    activity_value_list = pd.read_csv(source_path, sep=',')['activity_value'].tolist()
                
                # extract the smiles/inchi column
                for i, mol in tqdm(enumerate(inchi_list), total = len(inchi_list)):

                    pyg_data = inchi2graph(mol)

                    if pyg_data.valid is False:
                        invalid_id_list.append([mol_id, mol])
                    
                    pyg_data.y = torch.tensor([label], dtype=torch.int) # why int
                    if self.task_type == 'regression':
                        pyg_data.activity_value = torch.tensor([activity_value_list[i]], dtype=torch.float)
                    pyg_data.pubchem_cid = torch.tensor([int(cid_list[i])], dtype=torch.int)
                    pyg_data.mol_id = torch.tensor([mol_id], dtype=torch.int)
                    pyg_data.smiles = smiles_list[i]
                    data_list.append(pyg_data)
                    mol_id += 1
                    
            elif self.mol_repr == '3dmol':
                

                df = pd.read_csv(f'{source_path[:-4]}.csv')
                smiles_dict = df.set_index('CID')['SMILES'].to_dict()
                inchi_dict = df.set_index('CID')['InChI'].to_dict()
                if self.task_type == 'regression':  
                    activity_dict = df.set_index('CID')['activity_value'].to_dict()

                
                mol_conformer_list, cid_list = sdffile2mol_conformer(source_path)
                # get the smiles/inchi from the csv file with CID

                
                for i, (mol, conformer) in tqdm(enumerate(mol_conformer_list), total=len(mol_conformer_list)):
                    pyg_data = mol_conformer2graph3d(mol, conformer)
                    if pyg_data.valid is False: # need to implement valid
                        invalid_id_list.append([mol_id, mol])
                    
                    # need to extract the smiles/inchi from the csv file with CID
                    # directly converting mol to smiles/inchi encounters some bugs
                    pyg_data.pubchem_cid = torch.tensor([int(cid_list[i])], dtype=torch.int)
                    pyg_data.y = torch.tensor([label], dtype=torch.int)
                    pyg_data.smiles = smiles_dict[int(cid_list[i])]
                    pyg_data.inchi = inchi_dict[int(cid_list[i])]
                    
                    if self.task_type == 'regression':
                        pyg_data.activity_value = torch.tensor([activity_dict[int(cid_list[i])]], dtype=torch.float)
                    pyg_data.mol_id = torch.tensor([mol_id], dtype=torch.int)
                    data_list.append(pyg_data)
                    mol_id += 1
        
        # save invalid_id_list
        pd.DataFrame(invalid_id_list).to_csv(
            os.path.join(self.processed_dir, f'{self.name}-{self.mol_repr}-invalid_id.csv')
            , header=None, index=None)
        

        data, slices = self.collate(data_list)
        processed_file_path = os.path.join(self.processed_dir, f'processed_{self.mol_repr}_{self.name}.pt')
        torch.save((data, slices), processed_file_path)

        
    def get_idx_split(self, split_scheme):
        print(f'loading {split_scheme} split')
        path = osp.join(self.root, 'split')
        active_smiles = pd.read_csv(osp.join(self.root, 'raw', f'{self.name}_actives.csv'), sep=',')['SMILES'].tolist()
        inactive_smiles = pd.read_csv(osp.join(self.root, 'raw', f'{self.name}_inactives.csv'), sep=',')['SMILES'].tolist()
        all_smiles = active_smiles + inactive_smiles
        activity_dict = {}

        num_active = len(active_smiles)
        num_inactive = len(all_smiles) - num_active
        for idx in range(len(all_smiles)):
            if idx < num_active:
                activity_dict[idx] = 1
            else:
                activity_dict[idx] = 0

        print(f'before split: Active: {num_active}, Inactive: {num_inactive}')

        try: 
            if 'random' in split_scheme:
                print(osp.join(path, 'random', f'{split_scheme}.pt'))
                # print the current working directory
                print(os.getcwd())
                split_dict = torch.load(osp.join(path, 'random', f'{split_scheme}.pt'))
            elif 'scaffold' in split_scheme:
                split_dict = torch.load(osp.join(path, 'scaffold', f'{split_scheme}.pt'))

        except Exception as e:
            print(f'split file not found. Error msg: {e}')

        try:
            invalid_id_list = pd.read_csv(os.path.join(self.processed_dir, 
                                                       f'{self.name}-{self.mol_repr}-invalid_id.csv')).values.tolist()
            if len(invalid_id_list) == 0:
                print(f'invalid_id_list is empty')
            else:
                for id, label in invalid_id_list:
                    print(f'checking invalid id {id}')
                    if label == 1:
                        print('====warning: a positive label is removed====')
                    if id in split_dict['train']:
                        split_dict['train'].remove(id)
                        print(f'found in train and removed')
                    if id in split_dict['test']:
                        split_dict['test'].remove(id)
                        print(f'found in test and removed')           
            
            print(f'using {split_scheme} split')
            split_dict = torch.load(f'data_split/{self.name}_{split_scheme}.pt')
        except Exception as e:
            print(f'Cannot open invalid mol file: {e}')
        


        return split_dict



