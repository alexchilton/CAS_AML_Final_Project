import time
import numpy as np
import torch
from torch import nn
import torch.nn.functional as F

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Training function for XY-VAE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class Training():
    def __init__(self, 
                 name='Default', 
                 model=None, 
                 optimizer=None, 
                 with_kl=False,
                 device="cpu"):
        
        self.name = name
        self.model = model
        self.optimizer = optimizer
        self.with_kl = with_kl
        self.device = device
        self.history = {'train_kl_losses': [], 'train_recon_losses': [], 
                        'test_kl_losses': [], 'test_recon_losses': [], 'z_val': [], 'z_2d': []}

    def train(self, train_loader, test_loader, epochs=50, prints=10, 
              verbose=False, kl_sig=1., loss_type='MSE', loss_enlarge=1.):

        device = next(self.model.parameters()).device
        print(f"Training on device: {device}")
        
        t0 = time.time()
        print_mod = max(1, epochs//prints)  
        total_train_samples = len(train_loader.dataset)
        total_test_samples = len(test_loader.dataset)
        min_value, max_value = 0., kl_sig

        print (80*'~')
        print (f"Training model {self.name}")
        print (f"{total_train_samples} training samples, {total_test_samples} test samples, {epochs} epochs")
        if loss_enlarge != 1. :
            print (f"Loss per item multiplied {loss_enlarge} times.")
        print (80*'~')
        
        loss_fn = torch.nn.MSELoss(reduction='sum')
        
        for epoch in range(epochs):
            
            self.optimizer.zero_grad()
            train_kl_loss = torch.tensor([0.], device=device)
            test_kl_loss = torch.tensor([0.], device=device)
            train_recon_loss = torch.tensor([0.], device=device)
            test_recon_loss = torch.tensor([0.], device=device)
            train_loss = torch.tensor([0.], device=device)
            test_loss = torch.tensor([0.], device=device)
            Z = np.empty((0, self.model.code_size))
            Z_2d = np.empty((0, 2))
    
            ### train
            
            self.model.train()
            for data_x, data_y, labels in train_loader:

                data_x = data_x.to(device)
                data_y = data_y.to(device)
                labels = labels.to(device)
                
                d = {'x': data_x, 'labels': labels}
                out = self.model(d, verbose)

                Z = np.vstack([Z, out['z'].detach().cpu().numpy()])
                Z_2d = np.vstack([Z_2d, out['z_2d'].detach().cpu().numpy()])

                train_recon_loss += loss_fn(out['reconstructed'], data_y)
                
                if self.with_kl:
                    train_kl_loss += self.model.conditional_kl(out['mu'], out['logvar'], out['mu_prior'], out['logvar_prior'])

            train_loss = train_recon_loss
            
            if self.with_kl:
                β = min_value + 0.5 * (max_value - min_value) * (1 + torch.cos(torch.pi * torch.tensor(((epoch % print_mod) / print_mod), device=device)))
                train_loss += β * train_kl_loss
                
            train_loss.backward()
            self.optimizer.step()
            
            self.history['train_kl_losses'].append(train_kl_loss.item())
            self.history['train_recon_losses'].append(train_recon_loss.item())
            self.history['z_val'] = Z
            self.history['z_2d'] = Z_2d
    
            ### test
            
            self.model.eval()
            for data_x, data_y, labels in test_loader:

                data_x = data_x.to(device)
                data_y = data_y.to(device)
                labels = labels.to(device)
                
                d = {'x': data_x, 'labels': labels}
                out = self.model(d, verbose)

                test_recon_loss += loss_fn(out['reconstructed'], data_y)
                
                if self.with_kl:
                    test_kl_loss += self.model.conditional_kl(out['mu'], out['logvar'], out['mu_prior'], out['logvar_prior'])
                    
                    
            self.history['test_kl_losses'].append(test_kl_loss.item())
            self.history['test_recon_losses'].append(test_recon_loss.item())
            
            if epoch%print_mod == 0:
                s = f"Epoch [{epoch+1 : 6d}/{epochs}] Train KL-L: {loss_enlarge * train_kl_loss.item()/total_train_samples: 12.2f} | "
                s += f"Train RECON-L: {loss_enlarge * train_recon_loss.item()/total_train_samples: 12.2f} | "
                s += f"Test KL-L: {loss_enlarge * test_kl_loss.item()/total_test_samples: 12.2f} | "
                s += f"Test RECON-L: {loss_enlarge * test_recon_loss.item()/total_test_samples: 12.2f}"
                print(s)
                
        s = f"Epoch [{epoch+1 : 6d}/{epochs}] Train KL-L: {loss_enlarge * train_kl_loss.item()/total_train_samples: 12.2f} | "
        s += f"Train RECON-L: {loss_enlarge * train_recon_loss.item()/total_train_samples: 12.2f} | "
        s += f"Test KL-L: {loss_enlarge * test_kl_loss.item()/total_test_samples: 12.2f} | "
        s += f"Test RECON-L: {loss_enlarge * test_recon_loss.item()/total_test_samples: 12.2f}, final."
        print(s)
        
        t1 = time.time()
        print ('')
        print (f"Training time : {int(t1-t0)//60} minutes, {int(t1-t0)%60} seconds\n")