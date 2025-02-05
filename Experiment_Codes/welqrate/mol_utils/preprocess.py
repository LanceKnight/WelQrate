
from torch_geometric.data import Data
from rdkit import Chem
import numpy as np
from .features import (allowable_features, atom_to_feature_vector, atom_to_one_hot_vector,
 bond_to_feature_vector, atom_feature_vector_to_dict, bond_feature_vector_to_dict)
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import rdkit.Chem.EState as EState
import rdkit.Chem.rdPartialCharges as rdPartialCharges
import torch
from torch_geometric.nn import radius_graph



def ReorderCanonicalRankAtoms(mol):
    order = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol))])))[1]
    mol_renum = Chem.RenumberAtoms(mol, order)
    return mol_renum, order


def atomized_mol_level_features(all_atom_features, mol):
    '''
    Get more atom features that cannot be calculated only with atom,
    but also with mol
    :param all_atom_features:
    :param mol:
    :return:
    '''
    # Crippen has two parts: first is logP, second is Molar Refactivity(MR)
    all_atom_crippen = rdMolDescriptors._CalcCrippenContribs(mol)
    all_atom_TPSA_contrib = rdMolDescriptors._CalcTPSAContribs(mol)
    all_atom_ASA_contrib = rdMolDescriptors._CalcLabuteASAContribs(mol)[0]
    all_atom_EState = EState.EStateIndices(mol)

    new_all_atom_features = []
    for atom_id, feature in enumerate(all_atom_features):
        crippen_logP = all_atom_crippen[atom_id][0]
        crippen_MR = all_atom_crippen[atom_id][1]
        atom_TPSA_contrib = all_atom_TPSA_contrib[atom_id]
        atom_ASA_contrib = all_atom_ASA_contrib[atom_id]
        atom_EState = all_atom_EState[atom_id]

        feature.append(crippen_logP)
        feature.append(crippen_MR)
        feature.append(atom_TPSA_contrib)
        feature.append(atom_ASA_contrib)
        feature.append(atom_EState)

        new_all_atom_features.append(feature)
    return new_all_atom_features


def smiles2graph(smiles_string, removeHs=True, reorder_atoms=False):
    """
    Converts SMILES string to 2D graph Data object
    :input: SMILES string (str)
    :return: graph object
    """
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        mol = mol if removeHs else Chem.AddHs(mol)
        if reorder_atoms:
            mol, _ = ReorderCanonicalRankAtoms(mol)

    except Exception as e:
        print(f'cannot generate mol, error: {e}, smiles: {smiles_string}')

    if mol is None:

        smiles = smiles_cleaner(smiles_string)
        try:
            mol = Chem.MolFromSmiles(smiles)

            mol = mol if removeHs else Chem.AddHs(mol)
            if reorder_atoms:
                mol, _ = ReorderCanonicalRankAtoms(mol)

        except Exception as e:
            print(f'cannot generate mol, error: {e}, smiles: {smiles_string}')
            mol = None

        if mol is None:
            raise ValueError(f'cannot generate molecule with smiles: {smiles_string}')
        # Create an empty data object
            # graph = Data(
            #     edge_index=np.empty((2, 0), dtype=np.int64),
            #     edge_attr=np.empty((0, 3), dtype=np.int64),
            #     x=np.empty((0, 39), dtype=np.int64),
            #     num_nodes=0,
            #     num_edges=0,
            #     smiles=smiles_string,
            #     valid=False
            # )
            # print('Cannot generate molecule with smiles:', smiles_string, 'Returning empty graph')
    else:
        # calculate Gasteiger charges
        rdPartialCharges.ComputeGasteigerCharges(mol)
        
        # atoms
        atom_features_list = []
        for atom in mol.GetAtoms():
            atom_features_list.append(atom_to_feature_vector(atom))

        atom_features_list = atomized_mol_level_features(atom_features_list, mol)
        x = torch.tensor(atom_features_list, dtype=torch.float32)

        # bonds
        num_bond_features = 3  # bond type, bond stereo, is_conjugated
        if len(mol.GetBonds()) > 0:  # mol has bonds
            edges_list = []
            edge_features_list = []
            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                
                bond_attr = []
                bond_attr += one_hot_vector(bond.GetBondTypeAsDouble(),
                                            [1.0, 1.5, 2.0, 3.0])
                is_aromatic = bond.GetIsAromatic()
                is_conjugated = bond.GetIsConjugated()
                is_in_ring = bond.IsInRing()
                bond_attr.append(is_aromatic)
                bond_attr.append(is_conjugated)
                bond_attr.append(is_in_ring)

                # add edges in both directions
                edges_list.append((i, j))
                edge_features_list.append(bond_attr)
                edges_list.append((j, i))
                edge_features_list.append(bond_attr)

            # data.edge_index: Graph connectivity in COO format with shape [2, num_edges]
            edge_index = torch.tensor(edges_list).t().contiguous()
            # data.edge_attr: Edge feature matrix with shape [num_edges, num_edge_features]
            edge_attr = torch.tensor(edge_features_list, dtype=torch.float32)

        else:  # mol has no bonds
            edge_index = torch.from_numpy(np.empty((2, 0)))
            edge_attr = torch.from_numpy(np.empty((0, num_bond_features)))
            print('Warning: molecule does not have bond:', smiles_string)

        graph = Data(
            edge_index=edge_index,
            edge_attr=edge_attr,
            x=x,
            num_nodes=torch.tensor([len(x)]),
            num_edges=len(edge_index[0]),
            smiles=smiles_string
                            )

    return graph


def inchi2graph(inchi_string, removeHs=True, reorder_atoms=False):
    """
    Converts Inchi string to 2D graph Data object
    :input: SMILES string (str)
    :return: graph object
    """
    try:
        mol = Chem.MolFromInchi(inchi_string)
        mol = Chem.RemoveHs(mol) if removeHs else mol
        if reorder_atoms:
            mol, _ = ReorderCanonicalRankAtoms(mol)

    except Exception as e:
        print(f'cannot generate mol, error: {e}, inchi: {inchi_string}')
        mol = None

    if mol is None:
        raise ValueError(f'cannot generate molecule with inchi: {inchi_string}')

    else:
        # calculate Gasteiger charges
        rdPartialCharges.ComputeGasteigerCharges(mol)
        
        # atoms
        atom_features_list = []
        atom_num_list = []
        one_hot_atom_list = []
        for atom in mol.GetAtoms():
            atom_num_list.append(atom.GetAtomicNum())
            atom_features_list.append(atom_to_feature_vector(atom))
            one_hot_atom_list.append(atom_to_one_hot_vector(atom))

        atom_features_list = atomized_mol_level_features(atom_features_list, mol)
        x = torch.tensor(atom_features_list, dtype=torch.float32)
        atom_num = torch.tensor(atom_num_list, dtype=torch.int64)
        one_hot_atom = torch.tensor(one_hot_atom_list, dtype=torch.int64)

        # bond features: bond type, bond stereo, is_conjugated
        if len(mol.GetBonds()) > 0:  # mol has bonds
            edges_list = []
            edge_features_list = []
            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                
                bond_attr = []
                bond_attr += one_hot_vector(bond.GetBondTypeAsDouble(),
                                            [1.0, 1.5, 2.0, 3.0])
                is_aromatic = bond.GetIsAromatic()
                is_conjugated = bond.GetIsConjugated()
                is_in_ring = bond.IsInRing()
                bond_attr.append(is_aromatic)
                bond_attr.append(is_conjugated)
                bond_attr.append(is_in_ring)

                # add edges in both directions
                edges_list.append((i, j))
                edge_features_list.append(bond_attr)
                edges_list.append((j, i))
                edge_features_list.append(bond_attr)

            edge_index = torch.tensor(edges_list, dtype = torch.long).t().contiguous()
            edge_attr = torch.tensor(edge_features_list, dtype=torch.float32)

        else: 
            print('molecule does not have bond:', inchi_string)
            edge_index = torch.empty((2, 0), dtype=torch.long)
            edge_attr = torch.empty((0, 7), dtype=torch.float32)

        graph = Data(
            edge_index=edge_index,
            edge_attr=edge_attr,
            x=x,
            one_hot_atom=one_hot_atom,
            atom_num=atom_num,
            num_nodes=torch.tensor([len(x)]),
            num_edges=len(edge_index[0]),
            inchi=inchi_string,
        )

    return graph
    

### Functions that convert sdf files to molecular graphs

def sdffile2graph3d_lst(sdffile):
    """convert SDF file into a list of 3D graph.

    Args:
      sdffile: SDF file

    Returns:
      graph3d_lst: a list of 3D graph.
            each graph has (i) idx2atom (dict); (ii) distance_adj_matrix (np.array); (iii) bondtype_adj_matrix (np.array)

    """
    mol_conformer_lst = sdffile2mol_conformer(sdffile)
    graph3d_lst = mol_conformer2graph3d(mol_conformer_lst)
    return graph3d_lst

def sdffile2mol_conformer(sdffile):
    """convert sdffile into a list of molecule conformers.
    Args:
      sdffile: str, file
    Returns:
      smiles_lst: a list of molecule conformers.
    """

    supplier = Chem.SDMolSupplier(sdffile, removeHs=True, sanitize=False)
    mol_lst = [mol for mol in supplier if mol is not None]
    conformer_lst = []
    cid_lst = []
    for mol in mol_lst:
        if mol is not None:
            conformer = mol.GetConformer(id=0)
            cid = mol.GetProp('_Name')
            conformer_lst.append(conformer)
            cid_lst.append(cid)
    mol_conformer_lst = list(zip(mol_lst, conformer_lst))
    return mol_conformer_lst, cid_lst


def mol_conformer2graph3d(mol, conformer, removeHs=True, reorder_atoms=False):
    """convert (molecule, conformer) into a 3D graph.
    Args:
      mol_conformer: tuple (molecule, conformer)

    Returns:
      graph3d: a 3D pyg graph.
    """
    if mol is None:
            # Create an empty data object
        # graph = Data(
        #     edge_index=torch.empty((2, 0), dtype=torch.long),
        #     edge_attr=torch.empty((0, 3), dtype=torch.float32),
        #     x=torch.empty((0, 28), dtype=torch.int64),
        #     num_nodes=0,
        #     num_edges=0,
        #     valid=False
        # )
        raise ValueError(f'cannot generate molecule with conformer: {conformer}')
    else:
        mol = Chem.RemoveHs(mol, sanitize=False) if removeHs else mol
        if reorder_atoms:
            mol, _ = ReorderCanonicalRankAtoms(mol)
        
        # calculate Gasteiger charges
        rdPartialCharges.ComputeGasteigerCharges(mol)
        
        # atoms
        atom_features_list = []
        atom_num_list = []
        one_hot_atom_list = []
        for atom in mol.GetAtoms():
            atom_num_list.append(atom.GetAtomicNum())
            atom_features_list.append(atom_to_feature_vector(atom))
            one_hot_atom_list.append(atom_to_one_hot_vector(atom))

        atom_features_list = atomized_mol_level_features(atom_features_list, mol)
        x = torch.tensor(atom_features_list, dtype=torch.float32)
        atom_num = torch.tensor(atom_num_list, dtype=torch.int64)
        one_hot_atom = torch.tensor(one_hot_atom_list, dtype=torch.int64)

        # get 3D features: pos + edge_index
        positions = []
        for i in range(len(atom_num)):
            pos = conformer.GetAtomPosition(i)
            coordinate = np.array([pos.x, pos.y, pos.z]).reshape(1, 3)
            positions.append(coordinate)
        positions = np.concatenate(positions, 0)
                      
        edge_index = radius_graph(torch.tensor(positions), r=6.0, max_num_neighbors=50).to(torch.long)
                    
        graph = Data(
            edge_index=edge_index,
            x=x,
            pos= torch.tensor(positions, dtype=torch.float32),
            one_hot_atom=one_hot_atom,
            atom_num=atom_num,
            num_edges=len(edge_index[0]),
            num_nodes = torch.tensor([len(x)]),
        )
    
    return graph


def distance3d(coordinate_1, coordinate_2):
    return np.sqrt(
        sum([(c1 - c2)**2 for c1, c2 in zip(coordinate_1, coordinate_2)]))

pattern_dict = {'[NH-]': '[N-]', '[OH2+]':'[O]'}


def smiles_cleaner(smiles):
    '''
    This function is to clean smiles for some known issues that makes
    rdkit:Chem.MolFromSmiles not working
    '''
    print('fixing smiles for rdkit...')
    new_smiles = smiles
    for pattern, replace_value in pattern_dict.items():
        if pattern in smiles:
            print('found pattern and fixed the smiles!')
            new_smiles = smiles.replace(pattern, replace_value)
    return new_smiles

def one_hot_vector(val, lst):
	'''
	Converts a value to a one-hot vector based on options in lst
	'''
	if val not in lst:
		val = lst[-1]
	return map(lambda x: x == val, lst)





### Test
if __name__ == '__main__':
    graph = smiles2graph('O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5')

    # load sdf file
    sdf_file = 'test.sdf'
    # sdf_supplier = Chem.SDMolSupplier(sdf_file)
    # print(len(sdf_supplier))
    mol_conformer_lst = sdffile2mol_conformer(sdf_file)
    for mol, conformer in mol_conformer_lst:
        
        data = mol_conformer2graph3d(mol, conformer)
        print(data)
        print(data.x[:5, :])
        print(data.pos[:5, :])
        print(data.edge_attr[:5, :]) 

        break
    

    
