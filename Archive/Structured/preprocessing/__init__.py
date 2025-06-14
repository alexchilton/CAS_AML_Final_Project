# Import main classes to make them available directly from the package
from .ProteinPreprocessor import ProteinPreprocessor
from .ProteinVisualizer import ProteinVisualizer
from .GraphBuilder import GraphBuilder
from .FeatureExtractor import FeatureExtractor
from .ProteinFeatureProcessor import ProteinFeatureProcessor
from .ProteinGraphDataset import ProteinGraphDataset


# Define what's exported when using "from preprocessing import *"
__all__ = [
    'ProteinPreprocessor',
    'ProteinVisualizer',
    'GraphBuilder',
    'FeatureExtractor',
    'ProteinFeatureProcessor',
    'ProteinGraphDataset'
]

# You can add package-level variables if needed
__version__ = '0.1.1'