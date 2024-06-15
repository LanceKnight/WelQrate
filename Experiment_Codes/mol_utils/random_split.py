from multiprocessing import Pool
import torch
import random
import numpy as np

from sklearn.model_selection import train_test_split

def create_random_train_valid_splits(active_idx, inactive_idx, train_size=0.75):
    # Ensure indices are sorted
    random.seed(0)
    np.random.seed(0)
    
    active_idx = sorted(active_idx)
    inactive_idx = sorted(inactive_idx)

    # Create labels for active and inactive indices
    active_labels = ['active'] * len(active_idx)
    inactive_labels = ['inactive'] * len(inactive_idx)

    # Combine indices and labels
    all_idx = active_idx + inactive_idx
    all_labels = active_labels + inactive_labels

    # Split the dataset using train_test_split while maintaining stratification
    train_idx, valid_idx, _, _ = train_test_split(
        all_idx, all_labels, train_size=train_size, stratify=all_labels, random_state=42
    )

    # Count active and inactive samples in the training set
    num_active_train = sum(idx in train_idx for idx in active_idx)
    num_inactive_train = sum(idx in train_idx for idx in inactive_idx)

    # Count active and inactive samples in the validation set
    num_active_valid = sum(idx in valid_idx for idx in active_idx)
    num_inactive_valid = sum(idx in valid_idx for idx in inactive_idx)

    # Print the counts for debugging
    print(f"Training Set: {num_active_train} active, {num_inactive_train} inactive")
    print(f"Validation Set: {num_active_valid} active, {num_inactive_valid} inactive")

    return {
        'train': sorted(train_idx),
        'valid': sorted(valid_idx)
    }








# excludeed seed from parameters and fix the seed to 0

def generate_random_split(num_active, num_inactive, dataset_name, inactive_active_ratio=None, k_cv = 5):

    seed = 0
    active_idx = list(range(num_active))
    inactive_idx = list(range(num_active, num_active + num_inactive))

    random.seed(seed)
    random.shuffle(active_idx)
    random.shuffle(inactive_idx)

    num_active_train = round(num_active * 0.8)
    # print(f'num_active_train:{num_active_train}')
    num_inactive_train = round(num_inactive * 0.8)
    # num_active_test = num_active - num_active_train
    # num_inactive_test = round(num_inactive - num_inactive_train)

    split_dict = {}

    for k in range(k_cv):
        if inactive_active_ratio is None:
            filename = f'{dataset_name}_random_cv{k}.pt'

            # split_dict['train'] = active_idx[:num_active_train] + inactive_idx[:num_inactive_train]
        else:
            filename = f'{dataset_name}_random_IAratio{inactive_active_ratio}_cv{k}.pt'
            num_inactive = num_active * inactive_active_ratio
            inactive_idx = inactive_idx[:num_active*inactive_active_ratio]
            # split_dict['train'] = active_idx[:num_active_train] + inactive_idx[:num_active_train * inactive_active_ratio]
        print(f'filename:{filename}')

        all_ids = active_idx + inactive_idx
        print(f'len(all_ids):{len(all_ids)}')
        active_trunk = int(len(active_idx) / k_cv)
        inactive_trunk = int(len(inactive_idx) / k_cv)
        split_dict['test'] = active_idx[k * active_trunk: (k + 1) * active_trunk] + \
                             inactive_idx[k * inactive_trunk: (k + 1) * inactive_trunk]
        split_dict['train'] = [id for id in all_ids if id not in split_dict['test']]
        # split_dict['test'] = active_idx[num_active_train : ] + inactive_idx[num_inactive_train:]
        # print(f'split_dict["train"]:{split_dict["train"]}\n\n split_dict["test"]:{split_dict["test"]}')

        num_train = len(split_dict['train'])
        num_test = len(split_dict['test'])
        print(f'num_active:{num_active}, num_inactive:{num_inactive}')
        print(f'num_train:{num_train}, num_test:{num_test}')
        print(f'train_active:{len([id for id in split_dict["train"] if id < num_active])}')
        print(f'train_inactive:{len([id for id in split_dict["train"] if id >= num_active])}')
        print(f'test_active:{len([id for id in split_dict["test"] if id < num_active])}')
        print(f'test_inactive:{len([id for id in split_dict["test"] if id >= num_active])}')
        print(f'first 10 train:{split_dict["train"][:10]}')
        print(f'first 10 test:{split_dict["test"][:10]}')
        print('\n')

        return split_dict

        # data_md5 = hashlib.md5(json.dumps(split_dict, sort_keys=True).encode('utf-8')).hexdigest()
        # print(f'data_md5_checksum:{data_md5}')
        # print(f'file saved at {filename}')

def create_random_ksplits(active_idx, inactive_idx, k, k_cv, inactive_active_ratio=None):
    if inactive_active_ratio is not None:
        num_active = len(active_idx)
        num_inactive = int(num_active * inactive_active_ratio)
        inactive_idx = inactive_idx[:num_inactive]

    # Calculate size of each segment per fold
    active_trunk = int(len(active_idx) / k_cv)
    inactive_trunk = int(len(inactive_idx) / k_cv)

    # Calculate start and end indices for test and validation segments
    start_test = k * active_trunk
    end_test = (k + 1) * active_trunk if k < k_cv - 1 else len(active_idx)
    start_valid = end_test % len(active_idx)  # Use modulo for wrapping
    end_valid = ((k + 2) * active_trunk) % len(active_idx)  # Use modulo for wrapping

    # Adjust start and end indices for inactive parts
    start_test_inactive = k * inactive_trunk
    end_test_inactive = (k + 1) * inactive_trunk if k < k_cv - 1 else len(inactive_idx)
    start_valid_inactive = end_test_inactive % len(inactive_idx)  # Use modulo for wrapping
    end_valid_inactive = ((k + 2) * inactive_trunk) % len(inactive_idx)  # Use modulo for wrapping

    # Create test and validation sets
    test_idx = active_idx[start_test:end_test] + inactive_idx[start_test_inactive:end_test_inactive]
    valid_idx = active_idx[start_valid:end_valid] + inactive_idx[start_valid_inactive:end_valid_inactive]

    # Create train set as all indices not in test or validation
    all_ids = active_idx + inactive_idx
    train_idx = [id for id in all_ids if id not in test_idx and id not in valid_idx]

    return {
        'train': train_idx,
        'valid': valid_idx,
        'test': test_idx
    }


if __name__ == '__main__':

    dataset_info = {# Minus the molecules that do not have BCL features (i.e., gives BCL processing error)
        '435008':{'num_active':233, 'num_inactive':217923-24},#{'num_active':233, 'num_inactive':217925},
        '1798':{'num_active':164, 'num_inactive':60542},#{'num_active':187, 'num_inactive':61645},
        '435034': {'num_active':362, 'num_inactive':61393-6},#{'num_active':362, 'num_inactive':61394},
        '1843': {'num_active':172, 'num_inactive':301318-30},#{'num_active':172, 'num_inactive':301321},
        '2258': {'num_active':213, 'num_inactive':302189-31},#{'num_active':213, 'num_inactive':302192},
        '463087': {'num_active':703, 'num_inactive':100171-17},#{'num_active':703, 'num_inactive':100172},
        '488997': {'num_active':252, 'num_inactive':302051-31},#{'num_active':252, 'num_inactive':302054},
        '2689': {'num_active':172, 'num_inactive':319617-29},#{'num_active':172, 'num_inactive':319620},
        '485290': {'num_active':278, 'num_inactive':341026-36},#{'num_active':281, 'num_inactive':341084},
        '9999':{'num_active':37, 'num_inactive':293-1},
    }


    dataset_name_list = ['435008', '1798', '435034', '1843', '2258', '463087', '488997','2689', '485290', '9999']
    # dataset_name_list = ['1798']
    inactive_active_ratio = 100
    for dataset_name in dataset_name_list:
        num_active = dataset_info[dataset_name]['num_active']
        num_inactive = dataset_info[dataset_name]['num_inactive']
        generate_random_split(num_active, num_inactive, dataset_name,
                              inactive_active_ratio=inactive_active_ratio)