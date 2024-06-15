import torch
import torch.nn.functional as F
from torch.nn import Linear, Module, ReLU, Dropout
from .GNNConv import GAT as Encoder
from torch_geometric.nn import global_add_pool
import torch


class GATv2_Model(torch.nn.Module):
    def __init__(
            self,
            in_channels,
            hidden_channels,
            num_layers,
            heads=1,
            dropout=0.0,
            act='relu',
            one_hot=False,
            act_first=False,
            act_kwargs=None,
            norm=None,
            norm_kwargs=None,
            jk=None
    ):
        super(GATv2_Model, self).__init__()

        # Initialize the encoder using the imported GAT class
        self.encoder = Encoder(
            in_channels=in_channels,
            hidden_channels=hidden_channels,
            num_layers=num_layers,
            heads=heads,
            dropout=dropout,
            act=act,
            act_first=act_first,
            act_kwargs=act_kwargs,
            norm=norm,
            norm_kwargs=norm_kwargs,
            jk=jk,
            v2=True
        )

        # Pooling operation
        self.pool = global_add_pool
        self.one_hot = one_hot
        # Feedforward network components
        # final_out_dim = hidden_channels * heads if not jk else in_channels
        self.ffn_dropout = Dropout(p=0.25)
        self.lin1 = Linear(hidden_channels, 64)
        self.lin2 = Linear(64, 1)
        self.activate_func = ReLU()

    def forward(self, batch_data):
        edge_index, edge_attr, batch = batch_data.edge_index, batch_data.edge_attr, batch_data.batch
        if self.one_hot:
            x = batch_data.one_hot_atom
            node_embedding = self.encoder(x, edge_index, edge_attr=edge_attr)
        else:
            x = batch_data.x
            node_embedding = self.encoder(x, edge_index, edge_attr=edge_attr)
        graph_embedding = self.pool(node_embedding, batch)
        graph_embedding = self.ffn_dropout(graph_embedding)
        # print(graph_embedding.shape)
        prediction = self.lin2(self.activate_func(self.lin1(graph_embedding)))
        return prediction
