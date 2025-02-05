import torch
import random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

class Smiles2Vec(nn.Module):
    def __init__(self,
                 char_to_idx,
                 max_seq_len=270,
                 embedding_dim=50,
                 use_bidir=True,
                 use_conv=True,
                 filters=192,
                 kernel_size=3,
                 strides=1,
                 rnn_sizes=[300, 300,300],
                 rnn_types=["GRU", "GRU", "GRU"],
                 ):
        super(Smiles2Vec, self).__init__()
        
        self.char_to_idx = char_to_idx
        self.max_seq_len = max_seq_len
        self.embedding_dim = embedding_dim
        self.use_bidir = use_bidir
        self.use_conv = use_conv
        self.rnn_sizes = rnn_sizes
        self.rnn_types = rnn_types
        # self.mode = mode

        self.embedding = nn.Embedding(num_embeddings=len(char_to_idx), embedding_dim=embedding_dim)

        if use_conv:
            self.conv1d = nn.Conv1d(in_channels=embedding_dim,
                                    out_channels=filters,
                                    kernel_size=kernel_size,
                                    stride=strides)
        
        self.rnn_layers = nn.ModuleList()
        input_size = filters if use_conv else embedding_dim
        for idx, rnn_type in enumerate(rnn_types):
            rnn_class = nn.GRU if rnn_type == "GRU" else nn.LSTM
            rnn_layer = rnn_class(input_size, rnn_sizes[idx], batch_first=True, bidirectional=use_bidir)
            self.rnn_layers.append(rnn_layer)
            input_size = rnn_sizes[idx] * 2 if use_bidir else rnn_sizes[idx]

        self.fc = nn.Linear(input_size, 1)
        # if mode == "classification":
        #     self.sigmoid = nn.Sigmoid()
    
    def forward(self, batch_data):
        smiles_list = batch_data.smiles
        smiles_seqs = []
        for smiles in smiles_list:
            seq = [self.char_to_idx[char] for char in smiles]
            seq += [0] * (self.max_seq_len - len(seq))
            smiles_seqs.append(seq)
        
        x = torch.tensor(smiles_seqs, dtype=torch.long, device=batch_data.x.device)

        x = self.embedding(x)
        
        if self.use_conv:
            x = x.permute(0, 2, 1)  # Conv1d expects (batch, channels, seq_len)
            x = F.relu(self.conv1d(x))
            x = x.permute(0, 2, 1)  # Back to (batch, seq_len, channels)
        
        for rnn_layer in self.rnn_layers:
            x, _ = rnn_layer(x)
        
        output = self.fc(x[:, -1, :])
        # if self.mode == "classification":
        #     output = self.sigmoid(output)
        return output

# Test case
if __name__ == "__main__":
    # Sample character-to-index mapping
    # Set random seeds for reproducibility
    seed = 42
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    
    char_to_idx = {'[': 0, 'P': 1, 'e': 2, '@': 3, '1': 4, '#': 5, 
                   'n': 6, 'S': 7, '+': 8, '=': 9, ']': 10, 'c': 11, 
                   'o': 12, '2': 13, 'B': 14, '\\': 15, 'O': 16, 'H': 17, 
                   '5': 18, '-': 19, 'l': 20, '6': 21, 's': 22, '/': 23, 'C': 24, 
                   '4': 25, '3': 26, 'N': 27, 'I': 28, 'F': 29, '(': 30, ')': 31, 'r': 32}

    # Sample SMILES strings
    smiles_strings = ["Cc1ccc(C)c(OCC(O)CN2C(C)CCCC2C)c1", "CCOC(=O)N1CCN(C2=Nc3cccc4cccc2c34)CC1"]

    # Convert SMILES strings to integer sequences
    max_seq_len = 37  # Define the maximum sequence length
    smiles_seqs = []
    for smiles in smiles_strings:
        seq = [char_to_idx[char] for char in smiles]
        # Pad sequences to the maximum sequence length
        seq += [0] * (max_seq_len - len(seq))
        smiles_seqs.append(seq)

    smiles_seqs = torch.tensor(smiles_seqs, dtype=torch.long)

    # Instantiate the model
    model = Smiles2Vec(char_to_idx, max_seq_len=max_seq_len, mode="classification")

    # Pass the batch of SMILES sequences through the model
    outputs = model(smiles_seqs)

    # Print the outputs
    print("Model outputs:")
    print(outputs)
