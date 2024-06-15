import configparser
from argparse import ArgumentParser
import ast
from torch.optim import AdamW
from torch.nn import BCEWithLogitsLoss
from tqdm import tqdm
import torch
import random
import numpy as np
import os
from mol_utils.BCL_dataset import BCL_WelQrateDataset
from dataset import WelQrateDataset
from models.gnn2d.GCN import GCN_Model 
from models.gnn2d.GIN import GIN_Model
from models.gnn2d.GAT import GAT_Model
from models.gnn3d.DimeNet import DimeNet_Model, DimeNetplusplus_Model
from models.gnn3d.Spherenet import SphereNet_Model
from models.gnn3d.SchNet import SchNet
from models.MLP.mlp import MLP, bcl_MLP
from models.smiles.smiles2vec import Smiles2Vec
from models.smiles.textcnn import TextCNN
from loader import get_train_loader, get_test_loader, get_valid_loader, get_poor_train_loader
from scheduler import get_scheduler, get_lr
from utils.plot_loss import plot_epoch
import warnings
warnings.filterwarnings("ignore")




if __name__ == '__main__':
    
    # explore ways to refine the code
    parser = ArgumentParser()
    parser.add_argument('--config', type=str) 
    parser.add_argument('--no_train_eval', action='store_true')
    
    args = parser.parse_args()
    config_file = args.config
    config = configparser.ConfigParser()
    config.read(config_file)
    
    num_epochs = int(config['TRAIN']['num_epochs'])
    seed = int(config['GENERAL']['seed'])
    num_workers = int(config['GENERAL']['num_workers'])
    dataset_name = config['DATA']['dataset_name']
    batch_size = int(config['TRAIN']['batch_size'])
    model_type = config['MODEL']['model_type']
    split_scheme = config['DATA']['split_scheme']
    root = config['DATA']['root']  
    mol_repr = config['DATA']['mol_repr']
    one_hot = config['DATA']['one_hot']
    task_type = config['DATA']['task_type']
    weight_decay = float(config['TRAIN']['weight_decay'])

    
    if task_type == 'classification':
        from utils.train import train_class as train
        from utils.test import test_class as test
        loss_fn = BCEWithLogitsLoss()
        best_valid_logAUC = 0.0
    else:
        from utils.train import train_reg as train
        from utils.test import test_reg as test
        loss_fn = torch.nn.HuberLoss()
        best_valid_mae = 1e9
        

    if one_hot == 'True':
        one_hot = True
        in_channels = 12
    else:
        in_channels = 28
        one_hot = False


    print(f'====Model: {model_type}====')
    print(f'num_epochs: {num_epochs}')
    print(f'seed: {seed}')
    print(f'num_workers: {num_workers}')
    print(f'dataset_name: {dataset_name}')
    print(f'batch_size: {batch_size}')
    print(f'split_scheme: {split_scheme}')
    print(f'root: {root}')

    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if model_type in ['smiles2vec', 'textcnn']:
        assert mol_repr == '2dmol'
        
        if 'poor' in root:
            dataset = WelQrateDataset(dataset_name, root, mol_repr=mol_repr)
            good_root = root.replace('poor_dataset', 'dataset')
            good_dataset = WelQrateDataset(dataset_name, good_root, mol_repr=mol_repr)
            
            all_smiles = list(set(dataset.smiles + good_dataset.smiles))

            max_seq_len = max([len(data) for data in all_smiles])
            charset = set("".join(all_smiles))
            char_to_idx = dict((c,i) for i,c in enumerate(charset))            

            good_split_dict = good_dataset.get_idx_split(split_scheme)
            
            split_dict = dataset.get_idx_split(split_scheme)            
            
            train_loader = get_train_loader(dataset[split_dict['train']], batch_size, num_workers, seed)
            valid_loader = get_valid_loader(dataset[split_dict['valid']], batch_size, num_workers, seed)
            test_loader = get_test_loader(good_dataset[good_split_dict['test']], batch_size, num_workers, seed)

            
        else: 
            dataset = WelQrateDataset(dataset_name, root, mol_repr=mol_repr)
            max_seq_len = max([len(data.smiles) for data in dataset])
            charset = set("".join(list(dataset.smiles)))
            char_to_idx = dict((c,i) for i,c in enumerate(charset))
            
            split_dict = dataset.get_idx_split(split_scheme)
            train_loader = get_train_loader(dataset[split_dict['train']], batch_size, num_workers, seed)
            valid_loader = get_valid_loader(dataset[split_dict['valid']], batch_size, num_workers, seed)
            test_loader = get_test_loader(dataset[split_dict['test']], batch_size, num_workers, seed)    

    elif model_type in ['bcl_mlp']:
        dataset = BCL_WelQrateDataset(dataset_name, root)
        split_dict = dataset.get_idx_split(split_scheme)
        train_loader = get_train_loader(dataset[split_dict['train']], batch_size, num_workers, seed)
        valid_loader = get_valid_loader(dataset[split_dict['valid']], batch_size, num_workers, seed)
        test_loader = get_test_loader(dataset[split_dict['test']], batch_size, num_workers, seed)
    
    else:
        if 'poor' in root:
            dataset = WelQrateDataset(dataset_name, root, mol_repr=mol_repr)
            good_root = root.replace('poor_dataset', 'dataset')
            good_dataset = WelQrateDataset(dataset_name, good_root, mol_repr=mol_repr)
            good_split_dict = good_dataset.get_idx_split(split_scheme)
            split_dict = dataset.get_idx_split(split_scheme)
            
            train_loader = get_train_loader(dataset[split_dict['train']], batch_size, num_workers, seed)
            valid_loader = get_valid_loader(dataset[split_dict['valid']], batch_size, num_workers, seed)
            test_loader = get_test_loader(good_dataset[good_split_dict['test']], batch_size, num_workers, seed)
            
        else:
            dataset = WelQrateDataset(dataset_name, root, mol_repr=mol_repr)
            split_dict = dataset.get_idx_split(split_scheme)
            train_loader = get_train_loader(dataset[split_dict['train']], batch_size, num_workers, seed)
            valid_loader = get_valid_loader(dataset[split_dict['valid']], batch_size, num_workers, seed)
            test_loader = get_test_loader(dataset[split_dict['test']], batch_size, num_workers, seed)    

    
    print('Using One hot features:', one_hot, 'in_channels:', in_channels)
    
    if model_type == 'smiles2vec':
        model = Smiles2Vec(char_to_idx, 
                           max_seq_len, 
                           embedding_dim = int(config['MODEL']['embedding_dim']), 
                           use_bidir = bool(config['MODEL']['use_bidir']), 
                           use_conv = bool(config['MODEL']['use_conv']), 
                           rnn_sizes = ast.literal_eval(config['MODEL']['rnn_sizes']), 
                           rnn_types = ast.literal_eval(config['MODEL']['rnn_types']),).to(device)
        
    if model_type == 'textcnn':
        model = TextCNN(char_to_idx, 
                        max_seq_len, 
                        embedding_dim = int(config['MODEL']['embedding_dim']),
                        kernel_sizes = ast.literal_eval(config['MODEL']['kernel_sizes']),
                        num_filters = ast.literal_eval(config['MODEL']['num_filters']),
                        dropout = float(config['MODEL']['dropout'])).to(device)

    if model_type == 'gcn':
        # ensure mol_repr is 2dmol
        assert mol_repr == '2dmol'
        model = GCN_Model(in_channels= in_channels, 
                    hidden_channels=int(config['MODEL']['hidden_channels']),
                    num_layers=int(config['MODEL']['num_layers']),
                    one_hot = one_hot).to(device)
    
    elif model_type == 'gin':
        assert mol_repr == '2dmol'
        model = GIN_Model(in_channels=in_channels, 
                    hidden_channels=int(config['MODEL']['hidden_channels']),
                    num_layers=int(config['MODEL']['num_layers']),
                    one_hot = one_hot).to(device)
    
    elif model_type == 'gat':
        assert mol_repr == '2dmol'
        model = GAT_Model(in_channels=in_channels, 
                    hidden_channels=int(config['MODEL']['hidden_channels']), 
                    num_layers=int(config['MODEL']['num_layers']),
                    heads = int(config['MODEL']['heads']),
                    one_hot=one_hot).to(device)
    
    elif model_type == 'mlp':
        assert mol_repr == '2dmol'
        model = MLP(input_dim=in_channels, 
                    hidden_dim=int(config['MODEL']['hidden_channels']),
                    num_layers=int(config['MODEL']['num_layers']),
                    one_hot=one_hot).to(device)
    
    elif model_type == 'bcl_mlp':
        assert mol_repr == '2dmol'
        model = bcl_MLP(hidden_dim=int(config['MODEL']['hidden_channels']),
                        num_layers=int(config['MODEL']['num_layers']),).to(device)    

    
    elif model_type == 'schnet':
        assert mol_repr == '3dmol'
        model = SchNet(energy_and_force=False, 
                       cutoff=6.0, 
                       num_layers=int(config['MODEL']['num_layers']), 
                       in_channels= in_channels, 
                       hidden_channels=int(config['MODEL']['hidden_channels']), 
                       num_filters=int(config['MODEL']['num_filters']),
                       num_gaussians=int(config['MODEL']['num_gaussians']), 
                       out_channels=int(config['MODEL']['out_channels']),
                       one_hot=one_hot).to(device)
                       
    elif model_type == 'dimenet':
        model = DimeNet_Model(
                            hidden_channels = int(config['MODEL']['hidden_channels']), 
                            out_channels = int(config['MODEL']['out_channels']),
                            num_blocks = int(config['MODEL']['num_blocks']),
                            num_bilinear = int(config['MODEL']['num_bilinear']),
                            num_spherical = int(config['MODEL']['num_spherical']),
                            num_radial = int(config['MODEL']['num_radial']),
                            cutoff = float(config['MODEL']['cutoff']),
                            envelope_exponent = int(config['MODEL']['envelope_exponent']),
                            num_before_skip = int(config['MODEL']['num_before_skip']),
                            num_after_skip = int(config['MODEL']['num_after_skip']),
                            num_output_layers = int(config['MODEL']['num_output_layers']),
                            one_hot = one_hot).to(device)

    elif model_type == 'dimenet++':
        model = DimeNetplusplus_Model(
                            hidden_channels = int(config['MODEL']['hidden_channels']),
                            out_channels = int(config['MODEL']['out_channels']),
                            num_blocks = int(config['MODEL']['num_blocks']),
                            int_emb_size = int(config['MODEL']['int_emb_size']),
                            basis_emb_size = int(config['MODEL']['basis_emb_size']),
                            out_emb_channels = int(config['MODEL']['out_emb_channels']),
                            num_spherical = int(config['MODEL']['num_spherical']),
                            num_radial = int(config['MODEL']['num_radial']),
                            cutoff = float(config['MODEL']['cutoff']),
                            max_num_neighbors = int(config['MODEL']['max_num_neighbors']),
                            envelope_exponent = int(config['MODEL']['envelope_exponent']),
                            num_before_skip = int(config['MODEL']['num_before_skip']),
                            num_after_skip = int(config['MODEL']['num_after_skip']),
                            num_output_layers = int(config['MODEL']['num_output_layers']),
                            one_hot = one_hot).to(device)

    elif model_type == 'spherenet':
        model = SphereNet_Model(
            energy_and_force = config['MODEL']['energy_and_force'],
            cutoff = float(config['MODEL']['cutoff']),
            num_layers = int(config['MODEL']['num_layers']),
            hidden_channels = int(config['MODEL']['hidden_channels']),
            out_channels = int(config['MODEL']['out_channels']),
            int_emb_size = int(config['MODEL']['int_emb_size']),
            basis_emb_size_dist = int(config['MODEL']['basis_emb_size_dist']),
            basis_emb_size_angle = int(config['MODEL']['basis_emb_size_angle']),
            basis_emb_size_torsion = int(config['MODEL']['basis_emb_size_torsion']),
            out_emb_channels = int(config['MODEL']['out_emb_channels']),
            num_spherical = int(config['MODEL']['num_spherical']),
            num_radial = int(config['MODEL']['num_radial']),
            envelope_exponent = int(config['MODEL']['envelope_exponent']),
            num_before_skip = int(config['MODEL']['num_before_skip']),
            num_after_skip = int(config['MODEL']['num_after_skip']),
            num_output_layers = int(config['MODEL']['num_output_layers']),
            one_hot=one_hot).to(device)
    
    
    print(model)
    print(weight_decay)
    optimizer = AdamW(model.parameters(), weight_decay=weight_decay)
    scheduler = get_scheduler(optimizer, config, dataset[split_dict['train']])
    
    
    filename = f'saved_models/{dataset_name}_{model_type}_{split_scheme}_{task_type}.pt'
    
    best_epoch = 0
    
    print('\n' + '=' * 10 + ' STARTING TRAINING PROCESS... ' + '=' * 10 + '\n')
    os.makedirs(f'logs/{model_type}', exist_ok=True)
    os.makedirs('saved_models', exist_ok=True)
    with open(f'logs/{model_type}/loss_per_epoch.log', 'w+') as out_file:
        
        early_stopping_counter = 0
        early_stopping_limit = 30
        
        for epoch in range(num_epochs):
            loss = train(model, train_loader, optimizer, scheduler, device, loss_fn)
            lr = get_lr(optimizer)

            if not args.no_train_eval:
                if task_type == 'classification':
                    train_logAUC, train_EF, train_DCG, train_BEDROC = test(model=model, 
                                                                        loader=train_loader, 
                                                                        device=device, 
                                                                        type = 'train',
                                                                        model_name=model_type,
                                                                        save_result=True)
                else:
                    train_mae, train_rmse, train_r2 = test(model=model,
                                                            loader=train_loader,
                                                            device=device,
                                                            type = 'train',
                                                            model_name=model_type,
                                                            save_result=True)
                
                if epoch % 10 == 0:
                    if task_type == 'classification':
                        print(f'current_epoch={epoch} train_logAUC={train_logAUC} lr={lr}')
                        out_file.write(f'{epoch}\tloss={loss}\tlogAUC={train_logAUC}\tlr={lr}\t\n')
                    else:
                        print(f'current_epoch={epoch} train_mae={train_rmse} lr={lr}')
                        out_file.write(f'{epoch}\tloss={loss}\tmae={train_mae}\tlr={lr}\t\n')
                
            else:
                if epoch % 10 == 0:
                    print(f'current_epoch={epoch} loss={loss} lr={lr}')
                out_file.write(f'{epoch}\tloss={loss}\tlr={lr}\t\n')
            
            if task_type == 'classification':

                valid_logAUC, valid_EF, valid_DCG, valid_BEDROC = test(model=model, loader=valid_loader,
                                                                    device=device, type = 'valid',
                                                                    model_name=model_type,
                                                                    save_result=True)
                
                print(f'current_epoch={epoch} loss={loss} valid_logAUC={valid_logAUC} lr={lr}')
                out_file.write(f'\tvalid_logAUC={valid_logAUC}\tvalid_EF={valid_EF}\tvalid_DCG={valid_DCG}\tvalid_BEDROC={valid_BEDROC}\tlr={lr}\t\n')

                if valid_logAUC > best_valid_logAUC:
                    best_valid_logAUC = valid_logAUC
                    best_epoch = epoch
                    torch.save({'model': model.state_dict(),
                                'epoch': epoch}, filename)
                    early_stopping_counter = 0
                else:
                    early_stopping_counter += 1
                
            else:
                valid_mae, valid_rmse, valid_r2 = test(model=model, loader=valid_loader,
                                                                    device=device, type = 'valid',
                                                                    model_name=model_type,
                                                                    save_result=True)    
                print(f'current_epoch={epoch} loss={loss} valid_mae={valid_mae} lr={lr}')
                out_file.write(f'\tvalid_mae={valid_mae}\tvalid_rmse={valid_rmse}\tvalid_r2={valid_r2}\tlr={lr}\t\n')

                if valid_mae < best_valid_mae:
                    best_valid_mae = valid_mae
                    best_epoch = epoch
                    torch.save({'model': model.state_dict(),
                                'epoch': epoch}, filename)
                    
                    early_stopping_counter = 0
                else: 
                    early_stopping_counter += 1
                    
            if early_stopping_counter >= early_stopping_limit:
                print(f'Early stopping at epoch {epoch}')
                break

                
        print(f'finished training')
        
        
    # Testing

    if task_type == 'classification':    

        if os.path.exists(filename):
            model.load_state_dict(torch.load(filename)['model'])
        else:
            raise Exception(f'No trained model found. Please train the model first')
        print('testing...')
        logAUC, EF, DCG, BEDROC  = test(model, test_loader, device, type='test', model_name=model_type, save_result=True)
        print(f'{model_type.upper()} at epoch {best_epoch} Testing Results: logAUC={logAUC}\tEF={EF}\tDCG={DCG}\tBEDROC={BEDROC}\t')
        with open(f'result/{model_type}/test_result.txt', 'w+') as result_file: 
            result_file.write(f'logAUC={logAUC}\tEF={EF}\tDCG={DCG}\tBEDROC={BEDROC}')
            
    else:
        if os.path.exists(filename):
            model.load_state_dict(torch.load(filename)['model'])
        else:
            raise Exception(f'No trained model found. Please train the model first')
        print('testing...')
        mae, rmse, r2 = test(model, test_loader, device, type='test', model_name=model_type, save_result=True)
        print(f'{model_type.upper()} at epoch {best_epoch} Testing Results: MAE={mae}\tRMSE={rmse}\tR2={r2}\t')
        with open(f'result/{model_type}/test_result.txt', 'w+') as result_file: 
            result_file.write(f'MAE={mae}\tRMSE={rmse}\tR2={r2}')

