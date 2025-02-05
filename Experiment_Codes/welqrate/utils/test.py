import os
import sys
from tqdm import tqdm
import numpy as np
sys.path.append(os.path.abspath(os.getcwd()))
import torch
from torch_scatter import scatter_add
from utils.rank_prediction import rank_prediction
from utils.evaluation import calculate_logAUC, cal_EF, cal_DCG, cal_BEDROC_score, \
                       MAE, MSE, RMSE, R2

def test_class(model, loader, device, type, model_name, save_result=False):
    model.eval()
    filename = f'result/{model_name}/per_molecule_pred_of_{type}_set.txt'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    all_pred_y = []
    all_true_y = []

    for i, batch in enumerate(tqdm(loader)):
        batch.to(device)
        pred_y = model(batch).cpu().view(-1).detach().numpy()
        true_y = batch.y.view(-1).cpu().numpy()
        for j, _ in enumerate(pred_y):
            all_pred_y.append(pred_y[j])
            all_true_y.append(true_y[j])

    with open(filename, 'w') as out_file:
        if save_result:
            for k, _ in enumerate(all_pred_y):
                out_file.write(f'{all_pred_y[k]}\ty={all_true_y[k]}\n')
                
    if save_result:
        rank_prediction(type, model_name)
    all_pred_y = np.array(all_pred_y)
    all_true_y = np.array(all_true_y)
    logAUC = calculate_logAUC(all_true_y, all_pred_y)
    EF = cal_EF(all_true_y, all_pred_y, 100)
    DCG = cal_DCG(all_true_y, all_pred_y, 100)
    BEDROC = cal_BEDROC_score(all_true_y, all_pred_y)
    return logAUC, EF, DCG, BEDROC

# def test_class(model, loader, device, type, model_name, save_result=False):
#     model.eval()
#     filename = f'result/{model_name}/per_molecule_pred_of_{type}_set.txt'
#     os.makedirs(os.path.dirname(filename), exist_ok=True)
#     all_pred_y = []
#     all_true_y = []
#     for i, batch in enumerate(tqdm(loader)):
#         batch = batch.to(device)
        
#         # Check for isolated nodes
#         num_nodes = batch.num_nodes
#         connection_counts = torch.zeros(num_nodes, dtype=torch.int32, device=device)
#         ones = torch.ones_like(batch.edge_index[0], dtype=torch.int32)
#         connection_counts = scatter_add(ones, batch.edge_index[0], dim=0, dim_size=num_nodes)
#         connection_counts += scatter_add(ones, batch.edge_index[1], dim=0, dim_size=num_nodes)
        
#         # Add self-loops to isolated nodes
#         isolated_nodes = torch.where(connection_counts == 0)[0]
#         if len(isolated_nodes) > 0:
#             self_loops = torch.stack([isolated_nodes, isolated_nodes], dim=0)
#             batch.edge_index = torch.cat([batch.edge_index, self_loops], dim=1)
#             # Update node features if necessary (e.g., for attention mechanisms)
#             if hasattr(batch, 'edge_attr') and batch.edge_attr is not None:
#                 # Create self-loop edge attributes (you may need to adjust this based on your edge attributes)
#                 self_loop_attr = torch.zeros(len(isolated_nodes), batch.edge_attr.size(1), device=device)
#                 batch.edge_attr = torch.cat([batch.edge_attr, self_loop_attr], dim=0)
        
#         pred_y = model(batch).cpu().view(-1).detach().numpy()
#         true_y = batch.y.view(-1).cpu().numpy()
#         for j, _ in enumerate(pred_y):
#             all_pred_y.append(pred_y[j])
#             all_true_y.append(true_y[j])
    
#     with open(filename, 'w') as out_file:
#         if save_result:
#             for k, _ in enumerate(all_pred_y):
#                 out_file.write(f'{all_pred_y[k]}\ty={all_true_y[k]}\n')
                
#     if save_result:
#         rank_prediction(type, model_name, task_type='classification')
    
#     all_pred_y = np.array(all_pred_y)
#     all_true_y = np.array(all_true_y)
#     logAUC = calculate_logAUC(all_true_y, all_pred_y)
#     EF = cal_EF(all_true_y, all_pred_y, 100)
#     DCG = cal_DCG(all_true_y, all_pred_y, 100)
#     BEDROC = cal_BEDROC_score(all_true_y, all_pred_y)
#     return logAUC, EF, DCG, BEDROC


def test_reg(model, loader, device, type, model_name, save_result=False):
    model.eval()
    filename = f'result/{model_name}/per_molecule_pred_of_{type}_set.txt'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    all_pred_y = []
    all_true_y = []

    for i, batch in enumerate(tqdm(loader)):
        batch.to(device)
        pred_y = model(batch).cpu().view(-1).detach().numpy()
        true_y = batch.activity_value.view(-1).cpu().numpy()
        for j, _ in enumerate(pred_y):
            all_pred_y.append(pred_y[j])
            all_true_y.append(true_y[j])

    with open(filename, 'w') as out_file:
        if save_result:
            for k, _ in enumerate(all_pred_y):
                out_file.write(f'{all_pred_y[k]}\tactivity_value={all_true_y[k]}\n')
                
    if save_result:
        rank_prediction(type, model_name, task_type='regression')
    all_pred_y = np.array(all_pred_y)
    all_true_y = np.array(all_true_y)

    MAE_score = MAE(all_true_y, all_pred_y)
    # MSE_score = MSE(all_true_y, all_pred_y)
    
    RMSE_score = RMSE(all_true_y, all_pred_y)
    R2_score = R2(all_true_y, all_pred_y)
    return MAE_score, RMSE_score, R2_score


