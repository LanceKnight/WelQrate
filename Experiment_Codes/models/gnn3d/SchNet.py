import torch
import torch.nn.functional as F
from torch.nn import Embedding, Sequential, Linear
from torch_scatter import scatter
from torch_geometric.nn import radius_graph
from math import pi as PI
from torch_geometric.data import Data   
from torch.nn import Linear, Sequential, BatchNorm1d, ReLU
from torch.nn import Dropout, Linear


# need to create a new embedding function for the 28-dim node features which contain discrete and continuous features
class FeatureEmbedding(torch.nn.Module):
    def __init__(self, input_dim, output_dim):
        super(FeatureEmbedding, self).__init__()
        self.network = Sequential(
            Linear(input_dim, output_dim // 2),
            BatchNorm1d(output_dim // 2),
            ReLU(),
            Linear(output_dim // 2, output_dim),
            ReLU()
        )
    
    def forward(self, x):
        return self.network(x)

    def reset_parameters(self):
        for module in self.network:
            if hasattr(module, 'reset_parameters'):
                module.reset_parameters()

class update_e(torch.nn.Module):
    def __init__(self, hidden_channels, num_filters, num_gaussians, cutoff):
        super(update_e, self).__init__()
        self.cutoff = cutoff
        self.lin = Linear(hidden_channels, num_filters, bias=False)
        self.mlp = Sequential(
            Linear(num_gaussians, num_filters),
            ShiftedSoftplus(),
            Linear(num_filters, num_filters),
        )
        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.lin.weight)
        torch.nn.init.xavier_uniform_(self.mlp[0].weight)
        self.mlp[0].bias.data.fill_(0)
        torch.nn.init.xavier_uniform_(self.mlp[2].weight)
        self.mlp[0].bias.data.fill_(0)

    def forward(self, v, dist, dist_emb, edge_index):
        j, _ = edge_index
        C = 0.5 * (torch.cos(dist * PI / self.cutoff) + 1.0)
        W = self.mlp(dist_emb) * C.view(-1, 1)
        v = self.lin(v)
        e = v[j] * W
        return e

class update_v(torch.nn.Module):
    def __init__(self, hidden_channels, num_filters):
        super(update_v, self).__init__()
        self.act = ShiftedSoftplus()
        self.lin1 = Linear(num_filters, hidden_channels)
        self.lin2 = Linear(hidden_channels, hidden_channels)

        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.lin1.weight)
        self.lin1.bias.data.fill_(0)
        torch.nn.init.xavier_uniform_(self.lin2.weight)
        self.lin2.bias.data.fill_(0)

    def forward(self, v, e, edge_index):
        _, i = edge_index
        out = scatter(e, i, dim=0)
        out = self.lin1(out)
        out = self.act(out)
        out = self.lin2(out)
        return v + out

class update_u(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels):
        super(update_u, self).__init__()
        self.lin1 = Linear(hidden_channels, hidden_channels // 2)
        self.act = ShiftedSoftplus()
        self.lin2 = Linear(hidden_channels // 2, out_channels)

        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.lin1.weight)
        self.lin1.bias.data.fill_(0)
        torch.nn.init.xavier_uniform_(self.lin2.weight)
        self.lin2.bias.data.fill_(0)

    def forward(self, v, batch):
        v = self.lin1(v)
        v = self.act(v)
        v = self.lin2(v)
        u = scatter(v, batch, dim=0)
        return u

class emb(torch.nn.Module):
    def __init__(self, start=0.0, stop=5.0, num_gaussians=50):
        super(emb, self).__init__()
        offset = torch.linspace(start, stop, num_gaussians)
        self.coeff = -0.5 / (offset[1] - offset[0]).item() ** 2
        self.register_buffer('offset', offset)

    def forward(self, dist):
        dist = dist.view(-1, 1) - self.offset.view(1, -1)
        return torch.exp(self.coeff * torch.pow(dist, 2))

class ShiftedSoftplus(torch.nn.Module):
    def __init__(self):
        super(ShiftedSoftplus, self).__init__()
        self.shift = torch.log(torch.tensor(2.0)).item()

    def forward(self, x):
        return F.softplus(x) - self.shift

class SchnetLayer(torch.nn.Module):
    def __init__(self, energy_and_force=False, cutoff=6.0, num_layers=6, input_dim=28, hidden_channels=128, num_filters=128,
                 num_gaussians=50, out_channels=32, one_hot=False):
        # keep cutoff as 6.0 across all experiments
        
        super(SchnetLayer, self).__init__()
        self.energy_and_force = energy_and_force
        self.cutoff = cutoff
        self.num_layers = num_layers
        self.hidden_channels = hidden_channels
        self.num_filters = num_filters
        self.num_gaussians = num_gaussians
        self.one_hot = one_hot
        
        # Embedding for mixed feature types
        self.feature_embedding = FeatureEmbedding(input_dim, hidden_channels) 
        self.atom_embedding = Embedding(100, hidden_channels)
        
        self.dist_emb = emb(0.0, cutoff, num_gaussians)
        self.update_vs = torch.nn.ModuleList([update_v(hidden_channels, num_filters) for _ in range(num_layers)])
        self.update_es = torch.nn.ModuleList([update_e(hidden_channels, num_filters, num_gaussians, cutoff) for _ in range(num_layers)])
        self.update_u = update_u(hidden_channels, out_channels)
        self.reset_parameters()

    def reset_parameters(self):
        self.feature_embedding.network.apply(self.init_weights)
        self.atom_embedding.reset_parameters()
        for update_e in self.update_es:
            update_e.reset_parameters()
        for update_v in self.update_vs:
            update_v.reset_parameters()
        self.update_u.reset_parameters()

    @staticmethod
    def init_weights(m):
        if isinstance(m, torch.nn.Linear):
            torch.nn.init.xavier_uniform_(m.weight)
            torch.nn.init.constant_(m.bias, 0)

    def forward(self, batch_data):

        x, z, pos, batch, edge_index = batch_data.x , batch_data.one_hot_atom, batch_data.pos, batch_data.batch, batch_data.edge_index
        if self.energy_and_force:
            pos.requires_grad_()

        # edge_index = radius_graph(pos, r=self.cutoff, batch=batch)
        row, col = edge_index
        dist = (pos[row] - pos[col]).norm(dim=-1)
        
        dist_emb = self.dist_emb(dist)

        # one_hot: use 12-dimensional one-hot encoding for atom type; get the index and embed it
        # with torch.nn.Embedding; otherwise, use 28-dimensional continuous features and MLP to embed
        if self.one_hot:
            # print('using one-hot encoding')
            z_index = torch.argmax(z, dim=1)
            v = self.atom_embedding(z_index)
        else:
            # print('using all features')
            v = self.feature_embedding(x)
            
        for update_e, update_v in zip(self.update_es, self.update_vs):
            e = update_e(v, dist, dist_emb, edge_index)
            v = update_v(v, e, edge_index)

        u = self.update_u(v, batch)
        
        return u

### Schnet Model

class SchNet(torch.nn.Module):
    def __init__(
            self, energy_and_force=False, cutoff=6.0, num_layers=6, in_channels = 28, hidden_channels=128, num_filters=128,
            num_gaussians=50, out_channels=32, one_hot=False):

        super(SchNet, self).__init__()
        self.encoder = SchnetLayer(energy_and_force=energy_and_force,
                        cutoff=cutoff,
                        input_dim=in_channels,
                        num_layers=num_layers,
                        hidden_channels=hidden_channels,
                        num_filters=num_filters,
                        num_gaussians=num_gaussians,
                        out_channels=out_channels,
                        one_hot=one_hot
                    )
        self.one_hot = one_hot
        self.ffn_dropout = Dropout(p=0.25)
        self.lin1 = Linear(out_channels, 64)
        self.lin2 = Linear(64, 1)
        self.activate_func = ReLU()
        self.ff_dropout = Dropout(p=0.25)

    def forward(self, batch_data):
        # batch_data.z = batch_data.x.squeeze()
        graph_embedding = self.encoder(batch_data)

        graph_embedding = self.ffn_dropout(graph_embedding)
        prediction = self.lin2(self.activate_func(self.ff_dropout(self.lin1(graph_embedding))))
        
        return prediction



# import unittest
# from torch_geometric.data import Data

# class TestSchNet(unittest.TestCase):
#     def test_forward(self):
        
#         from torch_geometric.nn.pool import radius_graph

#         num_atoms = 5
#         feature_dim = 28  # Feature dimension per atom

#         features = torch.randn(2, num_atoms, feature_dim)
#         positions = torch.randn(2, num_atoms, 3)

#         batch_indices = torch.cat([torch.zeros(num_atoms, dtype=torch.long), 
#                                     torch.ones(num_atoms, dtype=torch.long)])

#         # Generate edge_index based on the Euclidean distance of nodes
#         edge_index = radius_graph(positions.view(-1, 3), r=10.0, batch=batch_indices)

#         batch_data = Data(x=features.view(-1, feature_dim), pos=positions.view(-1, 3), batch=batch_indices, edge_index=edge_index)
#         print('batch data x shape: ', batch_data.x.shape)
#         model = SchNet(energy_and_force=False, num_layers=3, hidden_channels=28, num_filters=128, num_gaussians=50, out_channels=32)

#         output = model(batch_data)
#         print("Output shape:", output.shape)
#         # self.assertEqual(output.shape, (2, 32))

# if __name__ == '__main__':
#     unittest.main()
