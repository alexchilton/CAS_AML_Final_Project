from torch_geometric.data import Data
from torch_geometric.loader import DataLoader as DL
from torch.utils.data import TensorDataset

class CustomDataset(TensorDataset):
    def __init__(self, data_input, data_output, data_labels=None):
        self.input = data_input
        self.output = data_output
        self.labels = data_labels

    def __len__(self):
        return len(self.input)

    def __getitem__(self, idx):
        x = self.input[idx]
        y = self.output[idx]
        if self.labels is not None:
            z = self.labels[idx]
            return x, y, z
        else:
            return x, y

def split_train_test(contact_map, SUBPROTEINS):
    temp_train =  contact_map[:int(SUBPROTEINS.len*0.9)]
    temp_test =   contact_map[int(SUBPROTEINS.len*0.9):]
    temp_train =  temp_train.view(temp_train.size(0),-1)
    temp_test =   temp_test.view(temp_test.size(0),-1)
    return (temp_train, temp_test)

def get_dataloader(train_x, train_y, train_labels=None, batch_size=32):
    if train_labels is not None:
        temp_dataset =     CustomDataset(train_x, train_y, train_labels)
    else:
        temp_dataset =     CustomDataset(train_x, train_y)
    temp_dataloader =  DL(temp_dataset, batch_size=batch_size, shuffle=False)
    return temp_dataloader