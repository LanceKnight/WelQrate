import rdkit
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import TorsionFingerprints
import numpy as np
import networkx as nx
import random
from rdkit import Chem
import torch
import rdkit.Chem.EState as EState
import rdkit.Chem.rdPartialCharges as rdPartialCharges
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
from .features import (atom_to_feature_vector, atom_to_one_hot_vector)

atomTypes = ['H', 'C', 'B', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I']
formalCharge = [-1, -2, 1, 2, 0]
degree = [0, 1, 2, 3, 4, 5, 6]
num_Hs = [0, 1, 2, 3, 4]
local_chiral_tags = [0, 1, 2, 3] 
hybridization = [
    rdkit.Chem.rdchem.HybridizationType.S,
    rdkit.Chem.rdchem.HybridizationType.SP,
    rdkit.Chem.rdchem.HybridizationType.SP2,
    rdkit.Chem.rdchem.HybridizationType.SP3,
    rdkit.Chem.rdchem.HybridizationType.SP3D,
    rdkit.Chem.rdchem.HybridizationType.SP3D2,
    rdkit.Chem.rdchem.HybridizationType.UNSPECIFIED,
    ]
bondTypes = ['SINGLE', 'DOUBLE', 'TRIPLE', 'AROMATIC']

def one_hot_embedding(value, options):
    embedding = [0]*(len(options) + 1)
    index = options.index(value) if value in options else -1
    embedding[index] = 1
    return embedding

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



def one_hot_vector(val, lst):
	'''
	Converts a value to a one-hot vector based on options in lst
	'''
	if val not in lst:
		val = lst[-1]
	return map(lambda x: x == val, lst)

def adjacency_to_undirected_edge_index(adj):
    adj = np.triu(np.array(adj, dtype = int)) #keeping just upper triangular entries from sym matrix
    array_adj = np.array(np.nonzero(adj), dtype = int) #indices of non-zero values in adj matrix
    edge_index = np.zeros((2, 2*array_adj.shape[1]), dtype = int) #placeholder for undirected edge list
    edge_index[:, ::2] = array_adj
    edge_index[:, 1::2] = np.flipud(array_adj)
    return edge_index

def get_all_paths(G, N = 3):
    # adapted from: https://stackoverflow.com/questions/28095646/finding-all-paths-walks-of-given-length-in-a-networkx-graph
    
    def findPaths(G,u,n):
        if n==0:
            return [[u]]
        paths = [[u]+path for neighbor in G.neighbors(u) for path in findPaths(G,neighbor,n-1) if u not in path]
        return paths
    
    allpaths = []
    for node in G:
        allpaths.extend(findPaths(G,node,N))
    
    return allpaths

def getNodeFeatures(list_rdkit_atoms, owningMol):
    F_v = (len(atomTypes)+1) +\
        (len(degree)+1) + \
        (len(formalCharge)+1) +\
        (len(num_Hs)+1)+\
        (len(hybridization)+1) +\
        2 + 4 + 5 # 52
    
    global_tags = dict(rdkit.Chem.FindMolChiralCenters(owningMol, force=True, includeUnassigned=True, useLegacyImplementation=False))
    
    node_features = np.zeros((len(list_rdkit_atoms), F_v))
    for node_index, node in enumerate(list_rdkit_atoms):
        features = one_hot_embedding(node.GetSymbol(), atomTypes) # atom symbol, dim=12 + 1 
        features += one_hot_embedding(node.GetTotalDegree(), degree) # total number of bonds, H included, dim=7 + 1
        features += one_hot_embedding(node.GetFormalCharge(), formalCharge) # formal charge, dim=5+1 
        features += one_hot_embedding(node.GetTotalNumHs(), num_Hs) # total number of bonded hydrogens, dim=5 + 1
        features += one_hot_embedding(node.GetHybridization(), hybridization) # hybridization state, dim=7 + 1
        features += [int(node.GetIsAromatic())] # whether atom is part of aromatic system, dim = 1
        features += [node.GetMass()  * 0.01] # atomic mass / 100, dim=1
        
        ### chiral tags go last ###
        #global chiral tag
        idx = node.GetIdx()
        global_chiral_tag = 0
        if idx in global_tags:
            if global_tags[idx] == 'R':
                global_chiral_tag = 1
            elif global_tags[idx] == 'S':
                global_chiral_tag = 2
            else:
                global_chiral_tag = -1
        
        features += one_hot_embedding(global_chiral_tag, [0,1,2]) # chiral tag of atom, dim=3+1 (global chiral features)
        
        #local chiral tag
        features += one_hot_embedding(node.GetChiralTag(), local_chiral_tags) # chiral tag of atom, dim=4+1 (local chiral features)
        
        node_features[node_index,:] = features
        
    return np.array(node_features, dtype = np.float32)

def getEdgeFeatures(list_rdkit_bonds):
    F_e = (len(bondTypes)+1) + 2 + (6+1) # 14
    
    edge_features = np.zeros((len(list_rdkit_bonds)*2, F_e))
    for edge_index, edge in enumerate(list_rdkit_bonds):
        features = one_hot_embedding(str(edge.GetBondType()), bondTypes) # dim=4+1
        features += [int(edge.GetIsConjugated())] # dim=1
        features += [int(edge.IsInRing())] # dim=1
        features += one_hot_embedding(edge.GetStereo(), list(range(6))) #dim=6+1

        # Encode both directed edges to get undirected edge
        edge_features[2*edge_index: 2*edge_index+2, :] = features
        
    return np.array(edge_features, dtype = np.float32)

def getInternalCoordinatesFromAllPaths(mol, adj, repeats = False): 
    if isinstance(mol, rdkit.Chem.rdchem.Conformer):
        conformer = mol
    if isinstance(mol, rdkit.Chem.rdchem.Mol):
        conformer = mol.GetConformer()
        
    graph = nx.from_numpy_array(adj, parallel_edges=False, create_using=None)
    
    distance_paths, angle_paths, dihedral_paths = get_all_paths(graph, N = 1), get_all_paths(graph, N = 2), get_all_paths(graph, N = 3)
    
    if len(dihedral_paths) == 0:
        print('No Dihedral Angle Detected')
        return None
    
    bond_distance_indices = np.array(distance_paths, dtype = int)
    bond_angle_indices = np.array(angle_paths, dtype = int)
    dihedral_angle_indices = np.array(dihedral_paths, dtype = int)
    
    if not repeats: # only taking (0,1) vs. (1,0); (1,2,3) vs (3,2,1); (1,3,6,7) vs (7,6,3,1)
        bond_distance_indices = bond_distance_indices[bond_distance_indices[:, 0] < bond_distance_indices[:, 1]]
        bond_angle_indices = bond_angle_indices[bond_angle_indices[:, 0] < bond_angle_indices[:, 2]]
        dihedral_angle_indices = dihedral_angle_indices[dihedral_angle_indices[:, 1] < dihedral_angle_indices[:, 2]]

    bond_distances = np.array([rdMolTransforms.GetBondLength(conformer, int(index[0]), int(index[1])) for index in bond_distance_indices], dtype = np.float32)
    bond_angles = np.array([rdMolTransforms.GetAngleRad(conformer, int(index[0]), int(index[1]), int(index[2])) for index in bond_angle_indices], dtype = np.float32)
    dihedral_angles = np.array([rdMolTransforms.GetDihedralRad(conformer, int(index[0]), int(index[1]), int(index[2]), int(index[3])) for index in dihedral_angle_indices], dtype = np.float32)
   
    return bond_distances, bond_distance_indices, bond_angles, bond_angle_indices, dihedral_angles, dihedral_angle_indices

def embedConformerWithAllPaths(rdkit_mol3D, repeats = False, removeHs = True):

    if removeHs:
        mol = Chem.RemoveHs(rdkit_mol3D, sanitize=False)

    # Edge Index
    adj = rdkit.Chem.GetAdjacencyMatrix(mol)
    edge_index = adjacency_to_undirected_edge_index(adj)

    # Edge Features
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
            # edge_feature = bond_to_feature_vector(bond)

            # add edges in both directions
            edges_list.append((i, j))
            edge_features_list.append(bond_attr)
            edges_list.append((j, i))
            edge_features_list.append(bond_attr)
            
        edge_index = torch.tensor(edges_list).t().contiguous()
        edge_attr = torch.tensor(edge_features_list, dtype=torch.float32)
        
    else:  # mol has no bonds
        edge_index = torch.from_numpy(np.empty((2, 0)))
        edge_attr = torch.from_numpy(np.empty((0, num_bond_features), dtype=torch.float32))
        
    # Node Features 
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

    # chiro features
    return_values = getInternalCoordinatesFromAllPaths(mol, adj, repeats = repeats)
    if return_values is not None:
        bond_distances, bond_distance_indices, bond_angles, bond_angle_indices, dihedral_angles, dihedral_angle_indices =\
        return_values

        return edge_index, edge_attr, x, bond_distances, bond_distance_indices, bond_angles, bond_angle_indices, dihedral_angles, dihedral_angle_indices
    else:
        return None