from dataset import WelQrateDataset
import rdkit
import rdkit.Chem as Chem
from tqdm import tqdm
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from mol_utils.scaffold_split import generate_scaffolds_list



dataset_name='AID1798'
good_root='../dataset/'
bad_root='../poor_dataset'
mol_repr='3dmol'
split_scheme='AID1798_3d_scaffold_seed5'

good_dataset = WelQrateDataset(dataset_name, good_root, mol_repr=mol_repr)
bad_dataset = WelQrateDataset(dataset_name, bad_root, mol_repr=mol_repr)

good_split = good_dataset.get_idx_split(split_scheme)
bad_split = bad_dataset.get_idx_split(split_scheme)

bad_train = bad_dataset[bad_split['train']]
bad_valid = bad_dataset[bad_split['valid']]
good_test = good_dataset[good_split['test']]

good_test_smiles = good_test.smiles
bad_train_smiles = bad_train.smiles
bad_valid_smiles = bad_valid.smiles


good_test_scaffolds = generate_scaffolds_list(good_test_smiles)
bad_train_scaffolds = generate_scaffolds_list(bad_train_smiles)
bad_valid_scaffolds = generate_scaffolds_list(bad_valid_smiles)

print(f'good_test:{len(good_test_scaffolds)}')
print(f'bad_train:{len(bad_train_scaffolds)}')
print(f'bad_valid:{len(bad_valid_scaffolds)}')

leakge_count = 0
for scaffold in good_test_scaffolds:
    if scaffold in bad_train_scaffolds or scaffold in bad_valid_scaffolds:
        #print('scaffold leakage')
        leakge_count += 1
    
print(leakge_count)
