import numpy as np
from scipy.spatial.distance import cdist
import random
from Bio.Align import substitution_matrices

class Protein:
    def __init__(self, atoms=[], aas=None, cas=None, blosum=False):
        self.atoms = atoms
        self.contact_maps = []
        self.contact_maps_config = []
        self.blosum = blosum

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

        if self.blosum:
            self.set_blosum()
        
    def get_oh(self, aas):
        idx = np.array([aas.index(aa) for aa in self.aa])
        oh = np.eye(len(aas))[idx]
        return oh

    def set_blosum(self):
        blosum = np.zeros(self.D.shape)
        for i, aa0 in enumerate(self.aa):
            for ii, aa1 in enumerate(self.aa):
                letter_1 = amino_acids[aa0]['one_letter'] if aa0 in amino_acids else 'X'
                letter_2 = amino_acids[aa1]['one_letter'] if aa1 in amino_acids else 'X'
                blosum[i,ii] = blosum62[(letter_1, letter_2)]
        self.blosum = blosum

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
        
    def explode(self, length):
        new_aas = []
        new_cas = []
        for i in range(length):
            new_aas.append(self.aa[-length+i:] + self.aa[:i])
            new_cas.append(self.ca[-length+i:] + self.ca[:i])
        for i in range(length, self.len):
            new_aas.append(self.aa[i-length:i])
            new_cas.append(self.ca[i-length:i])
        return (new_aas, new_cas)

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

blosum62 = substitution_matrices.load("BLOSUM62")

amino_acids = {
    "ALA": {"one_letter": "A", "name": "Alanine"},
    "ARG": {"one_letter": "R", "name": "Arginine"},
    "ASN": {"one_letter": "N", "name": "Asparagine"},
    "ASP": {"one_letter": "D", "name": "Aspartic Acid"},
    "CYS": {"one_letter": "C", "name": "Cysteine"},
    "GLN": {"one_letter": "Q", "name": "Glutamine"},
    "GLU": {"one_letter": "E", "name": "Glutamic Acid"},
    "GLY": {"one_letter": "G", "name": "Glycine"},
    "HIS": {"one_letter": "H", "name": "Histidine"},
    "ILE": {"one_letter": "I", "name": "Isoleucine"},
    "LEU": {"one_letter": "L", "name": "Leucine"},
    "LYS": {"one_letter": "K", "name": "Lysine"},
    "MET": {"one_letter": "M", "name": "Methionine"},
    "PHE": {"one_letter": "F", "name": "Phenylalanine"},
    "PRO": {"one_letter": "P", "name": "Proline"},
    "SER": {"one_letter": "S", "name": "Serine"},
    "THR": {"one_letter": "T", "name": "Threonine"},
    "TRP": {"one_letter": "W", "name": "Tryptophan"},
    "TYR": {"one_letter": "Y", "name": "Tyrosine"},
    "VAL": {"one_letter": "V", "name": "Valine"},
    "ASX": {"one_letter": "B", "name": "Aspartic Acid or Asparagine"},
    "GLX": {"one_letter": "Z", "name": "Glutamic Acid or Glutamine"},
    "UNK": {"one_letter": "X", "name": "Unknown or Other"}
}
