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
from mol_utils.embedding_functions import embedConformerWithAllPaths



class ChIRotDataset(InMemoryDataset):
    def __init__(self, dataset_name, root ='../dataset'):
        self.name = dataset_name
        self.root = root
        
        print(f'dataset stored in {self.root}')
        
        super(ChIRotDataset, self).__init__(self.root)
        
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
        return [f'processed_ChIRo_{self.name}.pt'] # example: processed_2dmol_AID1798.pt, processed_3dmol_AID1798.pt
    
    @property
    def processed_dir(self):
        return osp.join(self.root, 'processed', 'ChIRo') # example: processed/2dmol

    def download(self):
        # wait for the website to be up
        pass


    def chiro_process(self, mol):

        return_values = embedConformerWithAllPaths(mol, repeats=False)
        if return_values is not None:
            edge_index, edge_features, node_features, \
            bond_distances, bond_distance_index, bond_angles, \
            bond_angle_index, dihedral_angles, dihedral_angle_index = return_values
        else:
            return

        bond_angles = bond_angles % (2 * np.pi)
        dihedral_angles = dihedral_angles % (2 * np.pi)

        data = Data(
            x=torch.as_tensor(node_features),
            edge_index=torch.as_tensor(edge_index, dtype=torch.long),
            edge_attr=torch.as_tensor(edge_features))
        data.bond_distances = torch.as_tensor(bond_distances)
        data.bond_distance_index = torch.as_tensor(bond_distance_index,
                                                   dtype=torch.long).T
        data.bond_angles = torch.as_tensor(bond_angles)
        data.bond_angle_index = torch.as_tensor(bond_angle_index,
                                                dtype=torch.long).T
        data.dihedral_angles = torch.as_tensor(dihedral_angles)
        data.dihedral_angle_index = torch.as_tensor(
            dihedral_angle_index, dtype=torch.long).T

        return data


    def process(self):
        
        print(f'processing dataset {self.name}')
        RDLogger.DisableLog('rdApp.*')
        
        data_list = []
        invalid_id_list = []
        mol_id = 0
        
        for file_name, label in [(f'{self.name}_actives.sdf', 1),
                                 (f'{self.name}_inactives.sdf', 0)]:

            source_path = os.path.join(self.root, 'raw', file_name)
            print(f'loaded raw file from {source_path}')
            
            mol_conformer_list, cid_list = sdffile2mol_conformer(source_path)
            for i, (mol, conformer) in tqdm(enumerate(mol_conformer_list), total=len(mol_conformer_list)):
                pyg_data = self.chiro_process(mol)
 
                pyg_data.pubchem_cid = torch.tensor([int(cid_list[i])], dtype=torch.int)
                pyg_data.y = torch.tensor([label], dtype=torch.int)
                pyg_data.mol_id = torch.tensor([mol_id], dtype=torch.int)
                data_list.append(pyg_data)
                mol_id += 1

        # save invalid_id_list
        pd.DataFrame(invalid_id_list).to_csv(
            os.path.join(self.processed_dir, f'{self.name}-ChIRo-invalid_id.csv')
            , header=None, index=None)
        

        data, slices = self.collate(data_list)
        processed_file_path = os.path.join(self.processed_dir, f'processed_ChIRo_{self.name}.pt')
        torch.save((data, slices), processed_file_path)


if __name__ == '__main__':
    from torch_geometric.loader import DataLoader
    from models.gnn25d.ChIRoNet import ChIRoNet
    dataset = ChIRotDataset('AID9999', root='../dataset')
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    loader = DataLoader(dataset, batch_size=5, shuffle=True)

    model = ChIRoNet(
            F_z_list=[8,8,8],
            F_H=64,
            F_H_embed=28,
            F_E_embed=7, # original 14
            F_H_EConv=64,
            layers_dict={
                "EConv_mlp_hidden_sizes": [32, 32],
                "GAT_hidden_node_sizes": [64],
                "encoder_hidden_sizes_D": [64, 64],
                "encoder_hidden_sizes_phi": [64, 64],
                "encoder_hidden_sizes_c": [64, 64],
                "encoder_hidden_sizes_alpha": [64, 64],
                "encoder_hidden_sizes_sinusoidal_shift": [256, 256],
                "output_mlp_hidden_sizes": [128, 128]
            },
            activation_dict={
                "encoder_hidden_activation_D": torch.nn.LeakyReLU(negative_slope=0.01),
                "encoder_hidden_activation_phi": torch.nn.LeakyReLU(negative_slope=0.01),
                "encoder_hidden_activation_c": torch.nn.LeakyReLU(negative_slope=0.01),
                "encoder_hidden_activation_alpha": torch.nn.LeakyReLU(negative_slope=0.01),
                "encoder_hidden_activation_sinusoidal_shift": torch.nn.LeakyReLU(negative_slope=0.01),
                "encoder_output_activation_D": torch.nn.Identity(),
                "encoder_output_activation_phi": torch.nn.Identity(),
                "encoder_output_activation_c": torch.nn.Identity(),
                "encoder_output_activation_alpha": torch.nn.Identity(),
                "encoder_output_activation_sinusoidal_shift": torch.nn.Identity(),
                "EConv_mlp_hidden_activation": torch.nn.LeakyReLU(negative_slope=0.01),
                "EConv_mlp_output_activation": torch.nn.Identity(),
                "output_mlp_hidden_activation": torch.nn.LeakyReLU(negative_slope=0.01),
                "output_mlp_output_activation": torch.nn.Identity()
            },
            GAT_N_heads=4,
            chiral_message_passing=True,
            CMP_EConv_MLP_hidden_sizes=[256,256],
            CMP_GAT_N_layers=3,
            CMP_GAT_N_heads=2,
            c_coefficient_normalization="sigmoid",
            encoder_reduction='sum',
            output_concatenation_mode='molecule',
            EConv_bias=True,
            GAT_bias=True,
            encoder_biases=True,
            dropout=0.0
        ).to(device) 

    for batch in loader:
        batch = batch.to(device)
        out = model(batch)
        print(out.shape)

