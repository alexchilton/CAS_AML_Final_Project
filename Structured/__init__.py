"""
Protein ML Library - Tools for protein structure generation and analysis
"""

__version__ = '0.1.0'
__author__ = ''

# Import key classes and modules for easy access
from .preprocessing.ProteinPreprocessor import ProteinPreprocessor
from .preprocessing.ProteinVisualizer import ProteinVisualizer
from .preprocessing.GraphBuilder import GraphBuilder
from .preprocessing.FeatureExtractor import FeatureExtractor

from .model.Decoder import Decoder
from .model.Encoder import Encoder
from .model.LatentSpaceManager  import LatentSpaceManager
from .model.TrainingManager import TrainingManager
from .model.LossVisualizer import LossVisualizer

from .generation.GenerationVisualizer import GenerationVisualizer
from .generation.ProteinGenerator import ProteinGenerator
from .generation.GenerationOptimizer import GenerationOptimizer

from .utils.io_utils import io_utils
from .utils.visualization_utils import visualization_utils
from .utils.metrics import metrics

# Define __all__ for import *
__all__ = [
    "ProteinPreprocessor",
    "ProteinVisualizer",
    "GraphBuilder",
    "FeatureExtractor",
    "Decoder",
    "Encoder",
    "LatentSpaceManager",
    "TrainingManager",
    "LossVisualizer",
    "GenerationVisualizer",
    "ProteinGenerator",
    "GenerationOptimizer",
    "io_utils",
    "visualization_utils",
    "metrics",
]