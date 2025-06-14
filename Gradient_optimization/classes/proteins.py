import numpy as np
from scipy.spatial.distance import cdist
import random
from Bio.Align import substitution_matrices
from collections import Counter
import math
from Bio.SeqUtils import IsoelectricPoint
from Bio.Seq import Seq

class Protein:
    def __init__(self, atoms=[], aas=None, cas=None, blosum=False, source_idx=None):
        self.atoms = atoms
        self.contact_maps = []
        self.contact_maps_config = []
        self.blosum = blosum
        self.source_idx = source_idx

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

        self.fasta = ''.join([amino_acids[aa]['one_letter'] if aa in amino_acids else 'X' for aa in self.aa])

    @property
    def hydrophobicity_score(self):
        values = [hydropathy[aa] for aa in self.fasta if aa in hydropathy]
        return sum(values) / len(values) if values else 0.0

    @property
    def charge(self, pH=7.0):
        charge = net_charge(self.fasta, pH=pH)
        return charge

    @property
    def cysteine_fraction(self):
        if len(self.fasta) == 0:
            return 0.0
        return self.fasta.count('C') / len(self.fasta)

    @property
    def isoelectric_point(self):
        seq_obj = Seq(self.fasta)
        ip = IsoelectricPoint.IsoelectricPoint(seq_obj)
        return ip.pi()
        
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
        coords = np.array(self.ca)
        assert coords.ndim == 2 and coords.shape[1] == 3
        self.D = cdist(coords, coords, metric='euclidean')

    def add_contact_map(self, threshold_lower=0., threshold_upper=12., save=True):
        contact = np.zeros(self.D.shape)
        contact[(self.D>=threshold_lower) & (self.D<=threshold_upper)] = 1.
        if save:
            self.contact_maps.append(contact)
            self.contact_maps_config.append({'lower': threshold_lower, 'upper': threshold_upper})
        else:
            return contact

    def fuzzy_band(self, lower=0.0, upper=12.0, sigma=0.5):
        return 1 / (1 + np.exp(-(D - lower)/sigma)) - 1 / (1 + np.exp(-(D - upper)/sigma))

    def soft_binning(self, bins, sigma=0.5):
        """
        Transforme une matrice de distances en une pile de cartes de contact floues.
        
        bins : liste de (lower, upper)
        Retourne : shape (N, N, len(bins))
        """
        N = self.D.shape[0]
        contact_stack = np.zeros((N, N, len(bins)))
    
        for idx, (low, high) in enumerate(bins):
            contact_stack[:, :, idx] = fuzzy_band(D, low, high, sigma)
        
        # Normalise les scores sur l'axe des plages
        contact_stack /= np.sum(contact_stack, axis=-1, keepdims=True) + 1e-8
        self.soft_contact_maps = contact_stack 
        return contact_stack

    def add_soft_contact_map(self, threshold_lower=0.0, threshold_upper=12.0, sigma=0.5, save=True):
        """
        Crée une carte de contact floue où la valeur est proche de 1 si la distance D_ij est dans l'intervalle,
        et décroît en dehors selon une fonction sigmoïde.
        
        Args:
            threshold_lower: borne inférieure de la plage de contact
            threshold_upper: borne supérieure de la plage de contact
            sigma: paramètre de flou, plus petit = transition plus abrupte
            save: stocke la carte dans self.contact_maps si True
        
        Returns:
            soft_contact: matrice NxN avec valeurs entre 0 et 1
        """
        D = self.D.copy()
    
        # Sigmoïde inversée pour modéliser la plage floue
        lower_sigmoid = 1.0 / (1.0 + np.exp(-(D - threshold_lower) / sigma))
        upper_sigmoid = 1.0 / (1.0 + np.exp(-(D - threshold_upper) / sigma))
    
        # La "fenêtre" entre les deux sigmoïdes
        soft_contact = lower_sigmoid - upper_sigmoid
        print (soft_contact)
    
        if save:
            self.contact_maps.append(soft_contact)
            self.contact_maps_config.append({
                'lower': threshold_lower,
                'upper': threshold_upper,
                'sigma': sigma
            })
        else:
            return soft_contact
        
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

class Protein_Collection:
    def __init__(self):
        self.proteins = []

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            new_db = Protein_Collection()
            new_db.proteins = self.proteins[idx]
            return new_db
        elif isinstance(idx, list):
            new_db = Protein_Collection()
            [new_db.append(self.proteins[item]) for item in idx]
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
# Kyte-Doolittle hydropathy index
hydropathy = {
    'A': 1.8,  'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8,  'K': -3.9, 'M': 1.9,  'F': 2.8,  'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
pKa = {
    'K': 10.5, 'R': 12.5, 'H': 6.0,    # Basic side chains
    'D': 3.9,  'E': 4.1,               # Acidic side chains
    'C_term': 3.1, 'N_term': 8.0      # Termini
}

def charge_at_pH(pH, pKa_val, is_acid=False):
    if is_acid:
        return -1 / (1 + 10**(pKa_val - pH))
    else:
        return 1 / (1 + 10**(pH - pKa_val))

def net_charge(sequence, pH=7.0):
    counts = Counter(sequence)
    # Side chains
    charge = 0.0
    charge += counts.get('K', 0) * charge_at_pH(pH, pKa['K'], is_acid=False)
    charge += counts.get('R', 0) * charge_at_pH(pH, pKa['R'], is_acid=False)
    charge += counts.get('H', 0) * charge_at_pH(pH, pKa['H'], is_acid=False)
    charge += counts.get('D', 0) * charge_at_pH(pH, pKa['D'], is_acid=True)
    charge += counts.get('E', 0) * charge_at_pH(pH, pKa['E'], is_acid=True)
    # Termini
    charge += charge_at_pH(pH, pKa['N_term'], is_acid=False)
    charge += charge_at_pH(pH, pKa['C_term'], is_acid=True)
    return charge