"""
Protein ML Library - Tools for protein structure generation and analysis
"""

__version__ = '0.1.0'
__author__ = ''

# Import key classes for easy access
from .preprocessing import ProteinPreprocessor, ProteinVisualizer
from .model import Encoder, Decoder, TrainingManager
from .generation import ProteinGenerator
from .constraints import ConstraintManager