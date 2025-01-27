from torch_geometric.data import InMemoryDataset
from torch_geometric.data import Dataset, Data
import torch
import os
import os.path as osp
from tqdm import tqdm
from rdkit import Chem, RDLogger
from argparse import ArgumentParser
import pandas as pd
from glob import glob
import numpy as np
import random

class BCL_WelQrateDataset(InMemoryDataset):
    def __init__(self, dataset_name, root ='bcl_dataset'):

        self.name = dataset_name
        self.root = root
        
        print(f'dataset stored in {self.root}')
    
        
        super(BCL_WelQrateDataset, self).__init__(self.root)
        
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        # List files in the raw directory that match the dataset pattern
        pattern_actives = os.path.join(self.raw_dir, f'AID{self.name}_actives_bcl_feat.csv')
        pattern_inactives = os.path.join(self.raw_dir, f'AID{self.name}_inactives_bcl_feat.csv')

        files_actives = glob(pattern_actives)
        files_inactives = glob(pattern_inactives)
        return [os.path.basename(f) for f in files_actives + files_inactives]
    
    # need to change inchi/smiles to 2dmol
    @property
    def processed_file_names(self):
        return [f'processed_{self.name}_bcl.pt'] # example: processed_2dmol_AID1798.pt, processed_3dmol_AID1798.pt
    
    @property
    def processed_dir(self):
        return osp.join(self.root, 'processed') # example: processed/2dmol

    def download(self):
        # wait for the website to be up
        pass

    def process(self):

        print(f'processing dataset {self.name}')
        RDLogger.DisableLog('rdApp.*')
        
        data_list = []
        mol_id = 0
        
        for file_name, label in [(f'{self.name}_actives_bcl_feat.csv', 1),
                                 (f'{self.name}_inactives_bcl_feat.csv', 0)]:

            source_path = os.path.join(self.root, 'raw', file_name)
            print(f'loaded raw file from {source_path}')
            
            df = pd.read_csv(source_path, sep=',', header=None)
            for i in tqdm(range(len(df))):
                row = df.iloc[i]
                pyg_data = Data()
                
                if pd.isna(row.iloc[-1]):
                    print(f"Warning: NaN found at index {i} in the last column. Skipping this entry.")
                    continue
                
                pyg_data.y = torch.tensor([label], dtype=torch.int)
                pyg_data.bcl = torch.tensor(row[1:-2].values, dtype=torch.float32)
                pyg_data.mol_id = torch.tensor([mol_id], dtype=torch.int)
                pyg_data.pubchem_cid = torch.tensor([int(row.iloc[-1])], dtype=torch.int)
                data_list.append(pyg_data)
                mol_id += 1

 
        data, slices = self.collate(data_list)
        processed_file_path = os.path.join(self.processed_dir, f'processed_{self.name}_bcl.pt')
        torch.save((data, slices), processed_file_path)

        
    def get_idx_split(self, split_scheme):
        print(f'loading {split_scheme} split')
        path = osp.join(self.root, 'split')

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

        return split_dict