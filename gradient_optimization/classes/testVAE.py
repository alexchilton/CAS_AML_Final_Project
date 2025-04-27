import time
import torch
from torch import nn
import torch.nn.functional as F

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simple VAE where input != output (XY-VAE)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class XYVAE(nn.Module):
    def __init__(self, input_size=12, hidden_size=32, code_size=2, channels=4, output_size=12, dropout=0., epsilon="std"):
        super(XYVAE, self).__init__()
        
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.code_size = code_size
        self.channels = channels
        
        # Encoder
        self.encoder_fc1 = nn.Linear(input_size, 4 * hidden_size)  
        self.encoder_fc2 = nn.Linear(4 * hidden_size, 2 * hidden_size)  
        self.encoder_fc3 = nn.Linear(2 * hidden_size, hidden_size)  
        self.fc_mu = nn.Linear(hidden_size, code_size)
        self.fc_logvar = nn.Linear(hidden_size, code_size)

        self.encode_2d = nn.Linear(code_size, 2)
        
        # Decoder
        self.decoder_fc1 = nn.Linear(code_size, hidden_size) 
        self.decoder_fc2 = nn.Linear(hidden_size, 2 * hidden_size)  
        self.decoder_fc3 = nn.Linear(2 * hidden_size, 4 * hidden_size)  
        self.fc_out = nn.Linear(4 * hidden_size, output_size)  
        
        self.decoder_activation = nn.Softmax(dim=-1)

    def encode(self, x):
        h = self.encoder_fc1(x)
        h = self.encoder_fc2(h)
        h = self.encoder_fc3(h)
        # h = F.leaky_relu(h, negative_slope=0.05) 
        h = F.sigmoid(h)
        mu = self.fc_mu(h)       
        logvar = self.fc_logvar(h) 
        logvar = torch.clamp(logvar, min=-10.0, max=10.0)
        return mu, logvar

    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std, requires_grad=True)
        return mu + eps * std

    def decode(self, z):
        h = self.decoder_fc1(z)
        h = self.decoder_fc2(h)
        h = self.decoder_fc3(h)
        #h = F.sigmoid(h)
        h = self.fc_out(h)
        h = h.view(h.shape[0], 50, 50, self.channels)
        h = self.decoder_activation(h)
        h = h.view(h.shape[0], -1)
        return h
    
    def forward(self, d, verbose=False):
        x = d['x']
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        z_2d = self.encode_2d(z)
        
        reconstructed = self.decode(z)  
        r = {'reconstructed': reconstructed, 'z': z, 'mu': mu, 'logvar': logvar, 'z_2d': z_2d}
        return r

    def kl_divergence(self, mu, log_var):
        return -0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())

    def weights_init(self, m):
        if isinstance(m, nn.Linear) or isinstance(m, nn.Conv1d):
            nn.init.xavier_normal_(m.weight)
            if m.bias is not None:
                nn.init.zeros_(m.bias)

def print_tensor_info(var):
    frame = inspect.currentframe().f_back  # Get the caller's frame
    var_name = None
    
    for name, value in frame.f_locals.items():  # Iterate over local variables
        if value is var:
            var_name = name
            break

    if var_name is None:
        var_name = "Unknown Variable"
    
    print(f"{var_name : <30} {str(var.size()) : <30} Dtype: {var.dtype}")
