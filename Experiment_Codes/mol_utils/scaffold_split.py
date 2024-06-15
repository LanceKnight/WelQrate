import os
import json
import random
import hashlib
import torch
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles


def mol_to_smiles(sdf_path):
    """
    Read a SDF file and generate a list of SMILES strings.

    :param sdf_path: path to the SDF file.
    :return: a list of SMILES strings.
    """
    supplier = SDMolSupplier(sdf_path)
    smiles_list = []
    for mol in tqdm(supplier, desc=">Processing SMILES"):
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            smiles_list.append(smiles)
    return smiles_list


def generate_scaffolds(smiles_list, include_chirality=False):
    """
    Computes the Bemis-Murcko scaffold for a SMILES string list.

    :param mol: a list of SMILES strings.
    :param include_chirality: Whether to include chirality in the computed scaffold..
    :return: The Bemis-Murcko scaffolds for the molecules
    """
    scaffolds = {}
    not_parsed = []
    for idx, smiles, in tqdm(enumerate(smiles_list), desc=">Generating scaffolds"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            scaffold = MurckoScaffoldSmiles(mol=mol, includeChirality=include_chirality)
            if scaffold not in scaffolds:
                scaffolds[scaffold] = []
            scaffolds[scaffold].append(idx)
        else:
            not_parsed.append(idx)
    print(f'Number of scaffold generated: {len(scaffolds)}')
    print(f'Number of molecules not parsed: {len(not_parsed)}')
    return scaffolds, not_parsed

def generate_scaffolds_list(smiles_list, include_chirality=False):
    """
    Computes the unique Bemis-Murcko scaffold for a SMILES string list.

    :param smiles_list: a list of SMILES strings.
    :param include_chirality: Whether to include chirality in the computed scaffold.
    :return: A list of unique Bemis-Murcko scaffolds for the molecules
    """
    scaffolds = []
    for smiles in tqdm(smiles_list, desc=">Generating scaffolds"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            scaffold = MurckoScaffoldSmiles(mol=mol, includeChirality=include_chirality)
            scaffolds.append(scaffold)
    print(f'Number of unique scaffolds generated: {len(set(scaffolds))}')
    return scaffolds




import random

def ScaffoldSplitter(scaffolds, balanced, num_active, num_inactive, seed):
    """
    Split the indices of the molecules by scaffold so that no molecules sharing a scaffold are in different splits, 
    with ratios 0.6/0.2/0.2 for train/validation/test.

    :param scaffolds: a dictionary of scaffold: indices pairs.
    :param balanced: Whether to balance the splits by ensuring no single scaffold disproportionally affects any split size.
    :param num_active: number of active molecules.
    :param num_inactive: number of inactive molecules.
    :param seed: random seed.
    :return: train indices, validation indices, and test indices.
    """
    random.seed(seed)
    all_scaffold_indices = list(scaffolds.values())
    num_total = num_active + num_inactive
    num_train = int(num_total * 0.6)
    num_validation = int(num_total * 0.2)
    num_test = num_total - num_train - num_validation

    train_indices = []
    validation_indices = []
    test_indices = []

    if balanced:
        # Sort scaffolds into "big" and "small" by test size to ensure balance
        big_scaffold_indices = []
        small_scaffold_indices = []
        for scaffold_indices in all_scaffold_indices:
            if len(scaffold_indices) > num_test / 2:
                big_scaffold_indices.append(scaffold_indices)
            else:
                small_scaffold_indices.append(scaffold_indices)
        random.shuffle(small_scaffold_indices)
        random.shuffle(big_scaffold_indices)
        scaffold_indices_sets = big_scaffold_indices + small_scaffold_indices
    else:
        # Sort from largest to smallest scaffold sets:
        scaffold_indices_sets = sorted(all_scaffold_indices, key=len, reverse=True)

    # Distribute scaffolds into train, validation, and test sets
    for scaffold_indices in scaffold_indices_sets:
        if len(train_indices) + len(scaffold_indices) <= num_train:
            train_indices.extend(scaffold_indices)
        elif len(validation_indices) + len(scaffold_indices) <= num_validation:
            validation_indices.extend(scaffold_indices)
        else:
            test_indices.extend(scaffold_indices)

    # Reporting the distribution
    print(f'Total scaffolds: {len(all_scaffold_indices)}, '
          f'Train scaffolds: {len(train_indices)}, '
          f'Validation scaffolds: {len(validation_indices)}, '
          f'Test scaffolds: {len(test_indices)}')

    return train_indices, validation_indices, test_indices





import random

def PoorScaffoldSplitter(scaffolds, balanced, indice_list, seed):
    """
    Split the indices of the molecules by scaffold so that no molecules sharing a scaffold are in different splits,
    with a ratio of 3:1 for train/valid.

    :param scaffolds: a dictionary of scaffold: indices pairs.
    :param balanced: Whether to balance the train and validation sets by putting scaffolds that contain more than half of the validation size into train, the rest ordered randomly.
    :param indice_list: list of indices that need to be split.
    :param seed: random seed.
    :return: train indices and validation indices.
    """
    random.seed(seed)
    all_scaffold_indices = list(scaffolds.values())

    num_total = len(indice_list)
    num_train = round(num_total * 0.75)
    num_valid = num_total - num_train

    train_scaffold_count, valid_scaffold_count = 0, 0

    train_indices = []
    valid_indices = []

    if balanced:  # Put scaffolds that contain more than half of the validation size into train, the rest ordered randomly.
        big_scaffold_indices = []
        small_scaffold_indices = []
        for scaffold_indices in all_scaffold_indices:
            if len(scaffold_indices) > num_valid / 2:
                big_scaffold_indices.append(scaffold_indices)
            else:
                small_scaffold_indices.append(scaffold_indices)
        random.shuffle(small_scaffold_indices)
        random.shuffle(big_scaffold_indices)
        scaffold_indices_sets = big_scaffold_indices + small_scaffold_indices
    else:  # sort from largest to smallest scaffold sets
        scaffold_indices_sets = sorted(all_scaffold_indices, key=lambda x: len(x), reverse=True)

    for scaffold_indices_set in scaffold_indices_sets:
        if len(train_indices) + len(scaffold_indices_set) <= num_train:
            train_indices.extend(scaffold_indices_set)
            train_scaffold_count += 1
        else:
            valid_indices.extend(scaffold_indices_set)
            valid_scaffold_count += 1

    # Convert local scaffold indices to the actual indices in indice_list
    train_indices = [indice_list[i] for i in train_indices]
    valid_indices = [indice_list[i] for i in valid_indices]

    print(
        f'Total scaffolds: {len(all_scaffold_indices)}, train scaffolds: {train_scaffold_count}, validation scaffolds: {valid_scaffold_count}')

    return train_indices, valid_indices




def save_splits(train_indices, test_indices, dataset_name, seed, shrink=False):
    """
    Save the train and test indices to a file.

    :param train_indices: train indices.
    :param test_indices: test indices.
    :param dataset_name: dataset_name of the dataset.
    :param seed: random seed.
    """
    split_dict = {'train': train_indices, 'test': test_indices}
    directory = 'data_split'  # Changed to desired path
    if not os.path.exists(directory):
        os.makedirs(directory)
    if shrink:
        filename = os.path.join(directory, f'{dataset_name}-scaffold{seed}.pt')
    else:
        filename = os.path.join(directory, f'{dataset_name}_seed{seed}.pt')
    torch.save(split_dict, filename)
    data_md5 = hashlib.md5(json.dumps(split_dict, sort_keys=True).encode('utf-8')).hexdigest()
    # with open(f'{filename}.checksum', 'w') as checksum_file:
    #     checksum_file.write(data_md5)
    print(f'Completed: Dataset split for {dataset_name} with seed {seed} saved with checksum: {data_md5}')


def process_dataset(dataset_folder, dataset_name, seed_list, balanced=True, shrink=False):
    active_sdf_path = os.path.join(dataset_folder, f"{dataset_name}_actives_new.sdf")
    inactive_sdf_path = os.path.join(dataset_folder, f"{dataset_name}_inactives_new.sdf")

    smiles_active = mol_to_smiles(active_sdf_path)
    smiles_inactive = mol_to_smiles(inactive_sdf_path)

    all_smiles = smiles_active + smiles_inactive

    if shrink:
        all_smiles = all_smiles[:10000]
    else:
        all_smiles = all_smiles

    num_active = len(smiles_active)
    num_inactive = len(all_smiles) - num_active
    print(f'Before split: Active: {num_active}, Inactive: {num_inactive}')

    # make an index dictionary for activity (active/inacitive) for each in dex in all_smiles:
    activity_dict = {}
    for idx in range(len(all_smiles)):
        if idx < num_active:
            activity_dict[idx] = 1
        else:
            activity_dict[idx] = 0

    scaffolds, not_parsed = generate_scaffolds(all_smiles)

    for seed in seed_list:
        print(f'>Processing with seed {seed}')
        train_indices, test_indices = ScaffoldSplitter(scaffolds, balanced, num_active, num_inactive, seed)

        train_active = len([activity_dict[idx] for idx in train_indices if activity_dict[idx] == 1])
        train_inactive = len([activity_dict[idx] for idx in train_indices if activity_dict[idx] == 0])
        test_active = len([activity_dict[idx] for idx in test_indices if activity_dict[idx] == 1])
        test_inactive = len([activity_dict[idx] for idx in test_indices if activity_dict[idx] == 0])
        print(f'Train: Active: {train_active}, Inactive: {train_inactive}')
        print(f'Test: Active: {test_active}, Inactive: {test_inactive}')

        save_splits(train_indices, test_indices, dataset_name, seed, shrink)


if __name__ == '__main__':
    dataset_folder = '/home/dongh10/Documents/Datasets_cat/qsar/clean_sdf/raw'
    seed_list = [1, 2, 3, 4, 5]
    dataset_name_list = [ '1798', '435008', '435034', '1843', '2258', '463087', '488997', '2689', '485290']
    balanced_setting = True
    shrink = True

    # Update the loop to pass the 'shrink' parameter to the process_dataset function
    for dataset_name in tqdm(dataset_name_list, desc="Processing datasets"):
        print(f'Processing dataset {dataset_name}')
        process_dataset(dataset_folder, dataset_name, seed_list, balanced_setting, shrink)
        print(f'---------Finished processing dataset {dataset_name}---------')