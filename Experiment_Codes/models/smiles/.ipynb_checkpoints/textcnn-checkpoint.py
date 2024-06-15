import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from typing import Callable, Union, List, Dict, Any

def get_activation(activation_fn: Union[Callable, str]) -> Callable:
    if isinstance(activation_fn, str):
        if activation_fn == 'relu':
            return F.relu
        elif activation_fn == 'sigmoid':
            return torch.sigmoid
        else:
            raise ValueError(f"Unknown activation function: {activation_fn}")
    return activation_fn

class HighwayLayer(nn.Module):
    """
    Highway layer from "Training Very Deep Networks" [1]

    y = H(x) * T(x) + x * C(x), where

    H(x): 1-layer neural network with non-linear activation
    T(x): 1-layer neural network with sigmoid activation
    C(X): 1 - T(X); As per the original paper

    The output will be of the same dimension as the input

    References
    ----------
    .. [1] Srivastava et al., "Training Very Deep Networks".https://arxiv.org/abs/1507.06228

    Examples
    --------
    >>> x = torch.randn(16, 20)
    >>> highway_layer = HighwayLayer(d_input=x.shape[1])
    >>> y = highway_layer(x)
    >>> x.shape
    torch.Size([16, 20])
    >>> y.shape
    torch.Size([16, 20])
    """

    def __init__(self,
                 d_input: int,
                 activation_fn: Union[Callable, str] = 'relu'):
        """
        Initializes the HighwayLayer.

        Parameters
        ----------
            d_input: int
                the dimension of the input layer
            activation_fn: str
                the activation function to use for H(x)
        """

        super(HighwayLayer, self).__init__()
        self.d_input = d_input
        self.activation_fn = get_activation(activation_fn)
        self.sigmoid_fn = get_activation('sigmoid')

        self.H = nn.Linear(d_input, d_input)
        self.T = nn.Linear(d_input, d_input)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass of the HighwayLayer.

        Parameters
        ----------
        x: torch.Tensor
            Input tensor of dimension (,input_dim).

        Returns
        -------
        output: torch.Tensor
            Output tensor of dimension (,input_dim)
        """

        H_out = self.activation_fn(self.H(x))
        T_out = self.sigmoid_fn(self.T(x))
        output = H_out * T_out + x * (1 - T_out)
        return output

class TextCNN(nn.Module):
    def __init__(self,
                 char_to_idx,
                 max_seq_len=270,
                 embedding_dim=75,
                 kernel_sizes=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20],
                 num_filters=[100, 200, 200, 200, 200, 100, 100, 100, 100, 100, 160, 160],
                 dropout=0.25,
                 ):
        super(TextCNN, self).__init__()

        self.char_to_idx = char_to_idx
        self.max_seq_len = max_seq_len
        self.embedding_dim = embedding_dim
        self.kernel_sizes = kernel_sizes
        self.num_filters = num_filters
        self.dropout = dropout
        # self.mode = mode

        self.embedding = nn.Embedding(num_embeddings=len(char_to_idx), embedding_dim=embedding_dim)
        
        self.conv_layers = nn.ModuleList()
        for kernel_size, num_filter in zip(kernel_sizes, num_filters):
            conv_layer = nn.Conv1d(in_channels=embedding_dim,
                                   out_channels=num_filter,
                                   kernel_size=kernel_size)
            self.conv_layers.append(conv_layer)
        
        self.fc1 = nn.Linear(sum(num_filters), 200)
        self.highway = HighwayLayer(200)
        

        self.fc2 = nn.Linear(200, 1)
        
        self.dropout_layer = nn.Dropout(dropout)
        self.relu = nn.ReLU()
    

    def forward(self, batch_data):
        smiles_list = batch_data.smiles
        smiles_seqs = []
        for smiles in smiles_list:
            seq = [self.char_to_idx[char] for char in smiles]
            seq += [0] * (self.max_seq_len - len(seq))
            smiles_seqs.append(seq)
        
        x = torch.tensor(smiles_seqs, dtype=torch.long, device=batch_data.x.device)
        
        x = self.embedding(x)
        x = x.permute(0, 2, 1)  # Conv1d expects (batch, channels, seq_len)
        
        conv_outputs = []
        for conv_layer in self.conv_layers:
            conv_out = F.relu(conv_layer(x))
            conv_out, _ = torch.max(conv_out, dim=2)
            conv_outputs.append(conv_out)
        
        x = torch.cat(conv_outputs, dim=1)
        x = self.dropout_layer(x)
        
        x = self.fc1(x)
        x = self.highway(x)
        x = self.relu(x)
    
        x = self.fc2(x)
        
        return x
