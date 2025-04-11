import time
import torch
import torch.nn.functional as F
import numpy as np

class Training():
    def __init__(self, name='Default', model=None, optimizer=None, model_class='VAE', with_kl=False):
        
        self.name = name
        self.model = model
        self.optimizer = optimizer
        self.model_class = model_class
        self.with_kl = with_kl
        self.history = {'train_kl_losses': [], 'train_recon_losses': [], 
                        'test_kl_losses': [], 'test_recon_losses': [], 'z_val': [], 'z_2d': []}

    def train(self, train_loader, test_loader, epochs=50, prints=10, 
              verbose=False, kl_sig=1., loss_type='MSE', loss_enlarge=1.):
        
        t0 = time.time()
        print_mod = epochs//prints
        total_train_samples = len(train_loader.dataset)
        total_test_samples = len(test_loader.dataset)
        min_value, max_value = 0., kl_sig
        
        print (f"{total_train_samples} training samples, {total_test_samples} test samples")
        print (f"Loss per item multiplied {loss_enlarge} times.\n")
    
        if loss_type == 'MSE':                    # default
            loss_fn = torch.nn.MSELoss(reduction='sum')
        if loss_type == 'CE':
            loss_fn = torch.nn.CrossEntropyLoss() # applies softmax internally, model outputs raw
        elif loss_type == 'BCE':
            loss_fn = torch.nn.BCEWithLogitsLoss()
        
        for epoch in range(epochs):
            self.optimizer.zero_grad()
            train_kl_loss, test_kl_loss = torch.tensor([0.]), torch.tensor([0.])
            train_recon_loss, test_recon_loss = torch.tensor([0.]), torch.tensor([0.])
            train_loss, test_loss = torch.tensor([0.]), torch.tensor([0.])
            Z = np.empty((0, self.model.code_size))
            Z_2d = np.empty((0, 2))
    
            ### train
            self.model.train()
            for data_x, data_y in train_loader:
                
                if self.model_class == 'LSTM': d = {'x': data[0]}
                elif self.model_class == 'VAE': d = {'x': data[0]}
                elif self.model_class in ['CVAE','CPVAE']: d = {'x': data_x, 'y': data_y}
                else:
                    x, edge_index = data.x, data.edge_index
                    batch = torch.zeros(x.shape[0], dtype=torch.long)     
                    d = {'x': x, 'edge_index': edge_index, 'batch': batch}
                    
                out = self.model(d, verbose)
                    
                Z = np.vstack([Z, out['z'].detach().numpy()])
                Z_2d = np.vstack([Z_2d, out['z_2d'].detach().numpy()])
    
                if loss_type == 'CE':
                    x_cat = torch.argmax(d['x'], dim=1)
                else:
                    x_cat = d['x']
    
                if self.with_kl:
                    train_recon_loss += loss_fn(out['reconstructed_x'], x_cat)
                    if self.model_class == 'CPVAE':
                        train_kl_loss += self.model.conditional_kl(mu_q=out['mu'], logvar_q=out['logvar'], 
                                                                   mu_p=out['mu_prior'], logvar_p=out['logvar_prior'])
                    else:
                        train_kl_loss += self.model.kl_divergence(out['mu'], out['logvar'])
                else:
                    train_recon_loss += loss_fn(out['reconstructed_x'], x_cat)
    
            β = min_value + 0.5 * (max_value - min_value) * (1 + torch.cos(torch.pi * torch.tensor(((epoch % print_mod) / print_mod))))
            if self.with_kl:
                train_loss = train_recon_loss + β * train_kl_loss
            else:
                train_loss = train_recon_loss
                
            train_loss.backward()
            self.optimizer.step()
            
            self.history['train_kl_losses'].append(train_kl_loss.item())
            self.history['train_recon_losses'].append(train_recon_loss.item())
            self.history['z_val'] = Z
            self.history['z_2d'] = Z_2d
    
            ### test
            self.model.eval()
            for data_x, data_y in test_loader:
                
                if self.model_class == 'LSTM': d = {'x': data[0]}
                elif self.model_class == 'VAE': d = {'x': data[0]}
                elif self.model_class in ['CVAE','CPVAE']: d = {'x': data_x, 'y': data_y}
                else:
                    x, edge_index = data.x, data.edge_index
                    batch = torch.zeros(x.shape[0], dtype=torch.long)     
                    d = {'x': x, 'edge_index': edge_index, 'batch': batch}
                    
                out = self.model(d, verbose)
    
                if loss_type == 'CE':
                    x_cat = torch.argmax(d['x'], dim=1)
                else:
                    x_cat = d['x']
    
                if self.with_kl:
                    test_recon_loss += loss_fn(out['reconstructed_x'], x_cat)
                    if self.model_class == 'CPVAE':
                        test_kl_loss += self.model.conditional_kl(mu_q=out['mu'], logvar_q=out['logvar'], 
                                                                   mu_p=out['mu_prior'], logvar_p=out['logvar_prior'])
                    else:
                        test_kl_loss += self.model.kl_divergence(out['mu'], out['logvar'])
                else:
                    test_recon_loss += loss_fn(out['reconstructed_x'], x_cat)
            
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
        print (f"Training time : {int(t1-t0)//60} minutes, {int(t1-t0)%60} seconds")
