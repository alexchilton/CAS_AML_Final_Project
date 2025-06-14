import numpy as np
import torch
from functions.proteins import *

def Predict(Trainer, sample_orig):
    x = sample_orig
    mu, logvar = Trainer.model.encode(x)
    z = Trainer.model.reparameterize(mu, logvar)
    out = Trainer.model.decode(z) 
    
    return out.detach()
    
def get_subsequences(dataset, protein):
    PC_SUB = Protein_Collection()
    aas, cas = protein.explode(dataset.chunk_length)
    for i in range(len(aas)):
        subprotein = Protein(aas=aas[i], cas=cas[i])
        for config in dataset.configs:
            subprotein.add_contact_map(config['lower'], config['upper'], True)
        PC_SUB.append(subprotein)
    
    return PC_SUB

def get_training_sets(unique_aa, data, idx=0):
    AA_train = torch.stack([torch.tensor(protein.get_oh(list(unique_aa)), dtype=torch.float32) for protein in data])
    CM_train = torch.stack([torch.tensor(np.stack(protein.contact_maps, axis=-1), dtype=torch.float32) for protein in data])
    TRAIN_x =    AA_train
    TRAIN_out =  CM_train
    TRAIN_x =    TRAIN_x.view(TRAIN_x.size(0),-1)
    TRAIN_out =  TRAIN_out.view(TRAIN_out.size(0),-1)
    
    return (TRAIN_x, TRAIN_out)

def get_protein_average_contact(Dist_m, protein, chunk_length):

    protein_DM = np.zeros((protein.len, protein.len))
    divider_DM =  np.zeros((protein.len, protein.len))
    hard_mask = np.zeros((protein.len, protein.len)).astype(bool)
    
    for i in range(chunk_length):
        
        Dist_matrix = Dist_m[i].view(chunk_length,chunk_length).detach().numpy()
        
        protein_DM[-chunk_length+i:,-chunk_length+i:] += Dist_matrix[:chunk_length-i,:chunk_length-i]
        divider_DM[-chunk_length+i:,-chunk_length+i:] += 1.
        hard_mask[-chunk_length+i:,-chunk_length+i:] = 1

        if i == 0:
            protein_DM[-chunk_length:,-chunk_length:] += Dist_matrix[-i:,-i:]
            divider_DM[-chunk_length:,-chunk_length:] += 1.
            hard_mask[-chunk_length:,-chunk_length:] = 1
            
        if i > 0:
            protein_DM[:i,:i] += Dist_matrix[-i:,-i:]
            divider_DM[:i,:i] += 1.
            hard_mask[:i,:i] = 1

            protein_DM[-chunk_length+i:,:i] += Dist_matrix[:chunk_length-i,-i:]
            divider_DM[-chunk_length+i:,:i] += 1.
            hard_mask[-chunk_length+i:,:i] = 1

            protein_DM[:i, -chunk_length+i:] += Dist_matrix[-i:,:chunk_length-i]
            divider_DM[:i, -chunk_length+i:] += 1.
            hard_mask[:i, -chunk_length+i:] = 1
        
    for i in range(chunk_length, len(Dist_m)):
        
        Dist_matrix = Dist_m[i].view(chunk_length,chunk_length).detach().numpy()
        protein_DM[i-chunk_length:i, i-chunk_length:i] += Dist_matrix
        divider_DM[i-chunk_length:i, i-chunk_length:i] += 1.
        hard_mask[i-chunk_length:i, i-chunk_length:i] = 1
    
    divider_DM[divider_DM == 0.] = 1.
    protein_DM /= divider_DM

    return (protein_DM, hard_mask)