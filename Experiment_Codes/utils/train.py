
from tqdm import tqdm
import numpy as np



def train_class(model, loader, optimizer, scheduler, device, loss_fn):
    
    model.train()
    loss_list = []

    for i, batch in enumerate(tqdm(loader, miniters=100)):
        batch.to(device)
        # assert batch.edge_index.max() < batch.x.size(0), f"Edge index {batch.edge_index.max()} exceeds number of nodes"
        y_pred = model(batch)
        
        loss= loss_fn(y_pred.view(-1), batch.y.view(-1).float())
            
        loss_list.append(loss.item())
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()

    loss = np.mean(loss_list)
    return loss


def train_reg(model, loader, optimizer, scheduler, device, loss_fn):
    
    model.train()
    loss_list = []

    for i, batch in enumerate(tqdm(loader, miniters=100)):
        batch.to(device)
        # assert batch.edge_index.max() < batch.x.size(0), f"Edge index {batch.edge_index.max()} exceeds number of nodes"
        y_pred = model(batch)
        
        loss = loss_fn(y_pred.view(-1), batch.activity_value.view(-1))
            
        loss_list.append(loss.item())
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()

    loss = np.mean(loss_list)
    return loss