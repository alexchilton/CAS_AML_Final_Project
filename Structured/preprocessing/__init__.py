# Import main classes to make them available directly from the package
from .protein_preprocessor import ProteinPreprocessor
from .protein_visualizer import ProteinVisualizer
from .graph_builder import GraphBuilder
from .feature_extractor import FeatureExtractor

# Define what's exported when using "from preprocessing import *"
__all__ = [
    'ProteinPreprocessor',
    'ProteinVisualizer',
    'GraphBuilder',
    'FeatureExtractor'
]

# You can add package-level variables if needed
__version__ = '0.1.0'