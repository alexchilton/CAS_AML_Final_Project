import time
import torch
from torch import nn
import torch.nn.functional as F
import numpy as np
import inspect

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Conditional Prior VAE with deep conditional injection
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class FiLM(nn.Module):
    def __init__(self, label_size, hidden_size):
        super(FiLM, self).__init__()
        self.gamma = nn.Linear(label_size, hidden_size)
        self.beta = nn.Linear(label_size, hidden_size)

    def forward(self, h, y):
        gamma = self.gamma(y)
        beta = self.beta(y)
        return gamma * h + beta
        
class CPVAE(nn.Module):
    def __init__(self, input_size=12, label_size=10, label_hidden_size=10, hidden_size=32, code_size=2, dropout=0., epsilon="std"):
        super(CPVAE, self).__init__()
        
        self.input_size = input_size
        self.label_size = label_size
        self.label_hidden_size = label_hidden_size
        self.hidden_size = hidden_size
        self.code_size = code_size

        # Labels
        self.label_fc = nn.Linear(label_size, label_hidden_size)
        self.class_mu = nn.Parameter(torch.randn(label_size, code_size))
        self.class_logvar = nn.Parameter(torch.zeros(label_size, code_size))
        
        # Encoder
        self.encoder_fc1 = nn.Linear(input_size, 4 * hidden_size)  
        self.encoder_film1 = FiLM(label_hidden_size, 4 * hidden_size)
        self.encoder_fc2 = nn.Linear(4 * hidden_size, 2 * hidden_size)  
        self.encoder_film2 = FiLM(label_hidden_size, 2 * hidden_size)
        self.encoder_fc3 = nn.Linear(2 * hidden_size, hidden_size)  
        self.fc_mu = nn.Linear(hidden_size, code_size)
        self.fc_logvar = nn.Linear(hidden_size, code_size)

        self.encode_2d = nn.Linear(code_size, 2)
        
        # Decoder
        self.decoder_fc1 = nn.Linear(code_size, hidden_size) 
        self.decoder_film1 = FiLM(label_hidden_size, hidden_size)
        self.decoder_fc2 = nn.Linear(hidden_size, 2 * hidden_size) 
        self.decoder_film2 = FiLM(label_hidden_size, 2 * hidden_size)
        self.decoder_fc3 = nn.Linear(2 * hidden_size, 4 * hidden_size)  
        self.fc_out = nn.Linear(4 * hidden_size, input_size)   

    def encode(self, x, y):
        
        if x.dim() == 1:
            x = x.unsqueeze(0)
        if y.dim() == 1:
            y = y.unsqueeze(0)

        y_embed = self.label_fc(y)

        h = self.encoder_fc1(x)
        h = self.encoder_film1(h, y_embed)
        h = F.leaky_relu(h)
        
        h = self.encoder_fc2(h)
        h = self.encoder_film2(h, y_embed)
        
        h = self.encoder_fc3(h)
        h = F.leaky_relu(h, negative_slope=0.05) 

        mu = self.fc_mu(h)       
        logvar = self.fc_logvar(h) 
        logvar = torch.clamp(logvar, min=-10.0, max=10.0)
        return mu, logvar

    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std, requires_grad=True)
        return mu + eps * std

    def decode(self, z, y):
        
        if z.dim() == 1:
            z = z.unsqueeze(0)
        if y.dim() == 1:
            y = y.unsqueeze(0)

        y_embed = self.label_fc(y)

        h = self.decoder_fc1(z)
        h = self.decoder_film1(h, y_embed)
        
        h = self.decoder_fc2(h)
        h = self.decoder_film2(h, y_embed)
        
        h = self.decoder_fc3(h)
        h = F.leaky_relu(h, negative_slope=0.05) 
        h = self.fc_out(h)
        #h = F.sigmoid(h)
        return h
    
    def forward(self, d, verbose=False):
        x = d['x']
        y = d['y']
        
        class_indices = torch.argmax(y, dim=1)         # [B]
        mu_prior = self.class_mu[class_indices]        # [B, z_dim]
        logvar_prior = self.class_logvar[class_indices]

        mu, logvar = self.encode(x, y)
        z = self.reparameterize(mu, logvar)
        z_2d = self.encode_2d(z)
        
        recon_x = self.decode(z, y)  
        r = {'reconstructed_x': recon_x, 'z': z, 'z_2d': z_2d, 'mu': mu, 'logvar': logvar, 'mu_prior': mu_prior, 'logvar_prior': logvar_prior}
        return r

    @staticmethod
    def conditional_kl(mu_q, logvar_q, mu_p, logvar_p):
        return 0.5 * torch.sum(
            logvar_p - logvar_q +
            (torch.exp(logvar_q) + (mu_q - mu_p).pow(2)) / torch.exp(logvar_p) - 1
        )

    def weights_init(self, m):
        if isinstance(m, nn.Linear) or isinstance(m, nn.Conv1d):
            nn.init.xavier_normal_(m.weight)
            if m.bias is not None:
                nn.init.zeros_(m.bias)

def print_tensor_info(var):
    frame = inspect.currentframe().f_back 
    var_name = None
    
    for name, value in frame.f_locals.items():  
        if value is var:
            var_name = name
            break

    if var_name is None:
        var_name = "Unknown Variable"
    
    print(f"{var_name : <30} {str(var.size()) : <30} Dtype: {var.dtype}")
