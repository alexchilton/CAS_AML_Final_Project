import os
from setuptools import setup, find_packages

setup(
    name="Structured",
    version="0.1.1",
    packages=find_packages(),  # This will automatically find all packages in your structure
    install_requires=[

    ],
    description="Protein structure analysis and generation tools",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/protein_analysis",
    python_requires=">=3.9",
)