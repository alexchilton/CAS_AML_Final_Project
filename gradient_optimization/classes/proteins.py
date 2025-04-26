import numpy as np
from scipy.spatial.distance import cdist
import random

class Protein:
    def __init__(self, atoms=[], aas=None, cas=None, previous_aa=['NUL']):
        self.atoms = atoms
        self.previous_aa = previous_aa
        self.contact_maps = []
        self.contact_maps_config = []

        if len(atoms) > 0:
            self.get_Ca()
            self.get_D()
        elif aas is not None and cas is not None:
            self.aa = aas
            self.ca = cas
            self.len = len(self.aa)
            self.get_D()
        else:
            print ('Data missing.')

    def get_Ca(self):
        coords, aa = [], []
        for atom in self.atoms:
            if atom[3] == 'CA':
                coords.append(list(atom[-1]))
                aa.append(atom[2])
        self.aa = aa
        self.ca = coords
        self.len = len(self.aa)

    def get_D(self):
        self.D = cdist(self.ca, self.ca, metric='euclidean')

    def add_contact_map(self, threshold_lower=0., threshold_upper=12., save=True):
        contact = self.D.copy()
        contact[self.D<threshold_lower] = 0.
        contact[self.D>threshold_upper] = 0.
        contact[(self.D>=threshold_lower) & (self.D<=threshold_upper)] = 1.
        if save:
            self.contact_maps.append(contact)
            self.contact_maps_config.append({'lower': threshold_lower, 'upper': threshold_upper})
        else:
            return contact

    def get_padded_D(self, max_len):
        arr = np.zeros((max_len, max_len))
        arr[:self.len,:self.len] = self.D
        return arr

    def get_padded_aa(self, max_len, aas):
        idx = np.array([aas.index(aa) for aa in self.aa])
        oneHot = np.eye(len(aas))[idx]
        padded_oneHot = np.zeros((max_len, len(aas)))
        padded_oneHot[:self.len] = oneHot
        return padded_oneHot
        
    def get_previous_aa(self, max_len, aas):
        idx = np.array([aas.index(aa) for aa in self.previous_aa])
        oneHot = np.eye(len(aas))[idx]
        padded_oneHot = np.zeros((1, len(aas)))
        padded_oneHot[0] = oneHot
        return padded_oneHot
        
    def explode(self, length):
        new_aas = []
        new_cas = []
        previous_aa = ['NUL']
        for i in range(length, self.len):
            new_aas.append(self.aa[i-length:i])
            new_cas.append(self.ca[i-length:i])
            if i > length:
                previous_aa.append(self.aa[i-length-1])
        return (new_aas, new_cas, previous_aa)

class Database:
    def __init__(self):
        self.proteins = []

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            new_db = Database()
            new_db.proteins = self.proteins[idx]
            return new_db
        return self.proteins[idx]

    def append(self, data):
        self.proteins.append(data)

    def extend(self, data):
        self.proteins.extend(data)

    @property
    def len(self):
        return len(self.proteins)

    @property
    def randomize(self):
        random.shuffle(self.proteins)
