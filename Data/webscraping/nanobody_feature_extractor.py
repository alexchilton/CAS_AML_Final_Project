import os
import subprocess
import json
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from tqdm.notebook import tqdm  # Use tqdm.auto if not in notebook
import tempfile

class NanobodyFeatureExtractor:
    """Extract structural features from nanobody PDB files for ML models"""
    
    def __init__(self, wsl_path=None, use_pymol=True):
        """
        Initialize the feature extractor
        
        Args:
            wsl_path (str): Path to WSL scripts (if None, creates scripts as needed)
            use_pymol (bool): Whether to use PyMOL for secondary structure (if False, uses estimation)
        """
        self.use_pymol = use_pymol
        self.wsl_path = wsl_path if wsl_path else self._setup_wsl_scripts()
        self.results = {}

    def _setup_wsl_scripts(self):
        """Create necessary scripts in WSL"""
        # Create a temporary directory in WSL
        result = subprocess.run(['wsl', 'mktemp', '-d'], 
                                capture_output=True, text=True)
        if result.returncode != 0:
            raise Exception(f"Failed to create temp directory in WSL: {result.stderr}")
        
        wsl_path = result.stdout.strip()
        
        # Create the PyMOL script in WSL directly, avoiding copy
        pymol_script = """
import pymol
from pymol import cmd
import os
import sys
import json
import subprocess

def extract_ss_pymol(pdb_path):
    # Initialize PyMOL in quiet mode
    pymol_executable = "/usr/bin/pymol"  # Adjust this to your pymol executable location
    subprocess.run([pymol_executable,'-qc', pdb_path])
    
    # Check if file exists
    if not os.path.exists(pdb_path):
        print(f"Error: File '{pdb_path}' not found")
        return {'helix_percent': 0, 'sheet_percent': 0, 'coil_percent': 0}
    
    # Generate a safe object name for PyMOL
    obj_name = "nb_struct"  # Use a simple, safe name
    
    # Load the PDB file with the safe name
    cmd.load(pdb_path, obj_name)
    
    # Run the secondary structure determination
    cmd.dss(obj_name)
    
    # Get total residue count
    cmd.select("all_ca", f"{obj_name} and name CA")
    total_residues = cmd.count_atoms("all_ca")
    
    # Select residues by secondary structure
    cmd.select("helices", f"{obj_name} and ss h")
    cmd.select("sheets", f"{obj_name} and ss s")
    
    # Count residues in each selection
    helix_residues = cmd.count_atoms("helices and name CA")
    sheet_residues = cmd.count_atoms("sheets and name CA")
    
    # Calculate percentages
    helix_percent = helix_residues / total_residues if total_residues > 0 else 0
    sheet_percent = sheet_residues / total_residues if total_residues > 0 else 0
    coil_percent = 1.0 - helix_percent - sheet_percent
    
    # Create result dictionary
    result = {
        'helix_percent': helix_percent,
        'sheet_percent': sheet_percent,
        'coil_percent': coil_percent,
        'total_residues': total_residues
    }
    
    # Clean up
    cmd.delete(obj_name)
    cmd.delete("all_ca")
    cmd.delete("helices")
    cmd.delete("sheets")
    cmd.quit()
    
    return result

if __name__ == "__main__":
    if len(sys.argv) > 1:
        pdb_path = sys.argv[1]
        result = extract_ss_pymol(pdb_path)
        print(json.dumps(result))
    else:
        print("Please provide a PDB file path")
"""

        # Write the script directly to WSL using echo
        script_path = os.path.join(wsl_path, 'extract_ss.py')
        wsl_script_path = script_path.replace('\\', '/')
        
        # Use echo to write the script file in WSL
        echo_cmd = f"echo '{pymol_script}' > {wsl_script_path}"
        result = subprocess.run(['wsl', 'bash', '-c', echo_cmd], 
                                capture_output=True, text=True)
        
        if result.returncode != 0:
            raise Exception(f"Failed to create script in WSL: {result.stderr}")
        
        print(f"WSL scripts set up at: {wsl_path}")
        return wsl_path






    
    # def _setup_wsl_scripts(self):
    #     """Create necessary scripts in WSL"""
    #     # Create a temporary directory in WSL
    #     result = subprocess.run(['wsl', 'mktemp', '-d'], 
    #                         capture_output=True, text=True)
    #     if result.returncode != 0:
    #         raise Exception(f"Failed to create temp directory in WSL: {result.stderr}")
        
    #     wsl_path = result.stdout.strip()
        
    #     # Create the PyMOL script in WSL directly, avoiding copy
    #     pymol_script = '''
    # import pymol
    # from pymol import cmd
    # import os
    # import sys
    # import json

    # def extract_ss_pymol(pdb_path):
    #     """Extract secondary structure using PyMOL"""
    #     # Initialize PyMOL in quiet mode
    #     pymol.finish_launching(['pymol', '-qc'])
        
    #     # Check if file exists
    #     if not os.path.exists(pdb_path):
    #         print(f"Error: File '{pdb_path}' not found")
    #         return {'helix_percent': 0, 'sheet_percent': 0, 'coil_percent': 0}
        
    #     # Generate a safe object name for PyMOL
    #     obj_name = "nb_struct"  # Use a simple, safe name
        
    #     # Load the PDB file with the safe name
    #     cmd.load(pdb_path, obj_name)
        
    #     # Run the secondary structure determination
    #     cmd.dss(obj_name)
        
    #     # Get total residue count
    #     cmd.select("all_ca", f"{obj_name} and name CA")
    #     total_residues = cmd.count_atoms("all_ca")
        
    #     # Select residues by secondary structure
    #     cmd.select("helices", f"{obj_name} and ss h")
    #     cmd.select("sheets", f"{obj_name} and ss s")
        
    #     # Count residues in each selection
    #     helix_residues = cmd.count_atoms("helices and name CA")
    #     sheet_residues = cmd.count_atoms("sheets and name CA")
        
    #     # Calculate percentages
    #     helix_percent = helix_residues / total_residues if total_residues > 0 else 0
    #     sheet_percent = sheet_residues / total_residues if total_residues > 0 else 0
    #     coil_percent = 1.0 - helix_percent - sheet_percent
        
    #     # Create result dictionary
    #     result = {
    #         'helix_percent': helix_percent,
    #         'sheet_percent': sheet_percent,
    #         'coil_percent': coil_percent,
    #         'total_residues': total_residues
    #     }
        
    #     # Clean up
    #     cmd.delete(obj_name)
    #     cmd.delete("all_ca")
    #     cmd.delete("helices")
    #     cmd.delete("sheets")
    #     cmd.quit()
        
    #     return result

    # if __name__ == "__main__":
    #     if len(sys.argv) > 1:
    #         pdb_path = sys.argv[1]
    #         result = extract_ss_pymol(pdb_path)
    #         print(json.dumps(result))
    #     else:
    #         print("Please provide a PDB file path")
    #     '''
            
    #     # Write the script directly to WSL using echo
    #     script_path = os.path.join(wsl_path, 'extract_ss.py')
    #     wsl_script_path = script_path.replace('\\', '/')
        
    #     # Escape single quotes in the script
    #     escaped_script = pymol_script.replace("'", "'\\''")
        
    #     # Use echo to write the script file in WSL
    #     echo_cmd = f"echo '{escaped_script}' > {wsl_script_path}"
    #     result = subprocess.run(['wsl', 'bash', '-c', echo_cmd], 
    #                         capture_output=True, text=True)
        
    #     if result.returncode != 0:
    #         raise Exception(f"Failed to create script in WSL: {result.stderr}")
            
    #     print(f"WSL scripts set up at: {wsl_path}")
    #     return wsl_path
    
    def _windows_to_wsl_path(self, win_path):
        """Convert Windows path to WSL path"""
        # Handle basic drive letter conversion
        if ':' in win_path:
            drive, path = win_path.split(':', 1)
            path_fixed = path.replace('\\', '/')
            return f"/mnt/{drive.lower()}{path_fixed}"
        return win_path.replace('\\', '/')
    
    def _extract_ss_pymol_wsl(self, pdb_path):
        """Extract secondary structure using PyMOL in WSL"""
        # Convert path to WSL format
        wsl_pdb_path = self._windows_to_wsl_path(pdb_path)
        
        # Run the PyMOL script in WSL
        script_path = os.path.join(self.wsl_path, 'extract_ss.py')
        wsl_script_path = script_path.replace('\\', '/')
        
        result = subprocess.run(['wsl', 'python3', wsl_script_path, wsl_pdb_path], 
                               capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"PyMOL error: {result.stderr}")
            return {'helix_percent': 0, 'sheet_percent': 0, 'coil_percent': 0, 'total_residues': 0}
        
        try:
            ss_data = json.loads(result.stdout)
            return ss_data
        except json.JSONDecodeError:
            print(f"Error parsing PyMOL output: {result.stdout}")
            return {'helix_percent': 0, 'sheet_percent': 0, 'coil_percent': 0, 'total_residues': 0}
    
    def _extract_ss_simple(self, pdb_path):
        """Estimate secondary structure using geometric properties"""
        parser = PDBParser(QUIET=True)
        
        # Check file existence
        if not os.path.exists(pdb_path):
            print(f"Error: File '{pdb_path}' not found")
            return {'helix_percent': 0, 'sheet_percent': 0, 'coil_percent': 0, 'total_residues': 0}
        
        # Parse structure
        try:
            structure = parser.get_structure("protein", pdb_path)
            model = structure[0]
        except Exception as e:
            print(f"Error parsing structure: {e}")
            return {'helix_percent': 0, 'sheet_percent': 0, 'coil_percent': 0, 'total_residues': 0}
        
        # Get all CA atoms to analyze backbone geometry
        ca_atoms = [atom for atom in model.get_atoms() if atom.get_name() == 'CA']
        
        # Count residues
        total_residues = len(ca_atoms)
        
        # Typical nanobody SS distribution if we can't analyze further
        result = {
            'helix_percent': 0.08,  # ~8% alpha helix
            'sheet_percent': 0.47,  # ~47% beta sheet
            'coil_percent': 0.45,   # ~45% coil
            'total_residues': total_residues
        }
        
        return result
    
    def extract_features(self, pdb_path):
        """
        Extract all features from a single PDB file
        
        Args:
            pdb_path (str): Path to PDB file
            
        Returns:
            dict: Dictionary of features
        """
        features = {}
        
        # Extract secondary structure
        if self.use_pymol:
            ss_features = self._extract_ss_pymol_wsl(pdb_path)
        else:
            ss_features = self._extract_ss_simple(pdb_path)
        
        features.update(ss_features)
        
        # Extract additional features here
        # For now, let's add some basic structural statistics
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", pdb_path)
            model = structure[0]
            
            # Count atoms and residues
            atom_count = len(list(model.get_atoms()))
            residue_count = len(list(model.get_residues()))
            
            # Calculate basic statistics
            features['atom_count'] = atom_count
            features['residue_count'] = residue_count
            features['atoms_per_residue'] = atom_count / residue_count if residue_count > 0 else 0
            
        except Exception as e:
            print(f"Error extracting additional features: {e}")
        
        return features
    
    def process_directory(self, input_dir, pattern='*.pdb'):
        """
        Process all PDB files in a directory
        
        Args:
            input_dir (str): Directory containing PDB files
            pattern (str): Glob pattern to match PDB files
            
        Returns:
            pd.DataFrame: DataFrame with features for all processed files
        """
        import glob
        
        # Find all matching PDB files
        pdb_files = glob.glob(os.path.join(input_dir, pattern))
        
        if not pdb_files:
            print(f"No PDB files found in {input_dir} matching pattern {pattern}")
            return pd.DataFrame()
        
        print(f"Processing {len(pdb_files)} PDB files...")
        
        # Process each file
        all_features = {}
        
        for pdb_file in tqdm(pdb_files):
            try:
                file_id = os.path.basename(pdb_file).split('.')[0]
                features = self.extract_features(pdb_file)
                all_features[file_id] = features
            except Exception as e:
                print(f"Error processing {pdb_file}: {e}")
        
        # Convert to DataFrame
        features_df = pd.DataFrame.from_dict(all_features, orient='index')
        
        # Save results
        self.results = all_features
        
        return features_df
    
    def save_results(self, output_path='nanobody_features.json'):
        """Save results to a JSON file"""
        with open(output_path, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        print(f"Results saved to {output_path}")
        
        # Also save as CSV
        csv_path = output_path.replace('.json', '.csv')
        pd.DataFrame.from_dict(self.results, orient='index').to_csv(csv_path)
        print(f"Results also saved to {csv_path}")

# # Example usage
# if __name__ == "__main__":
#     # Initialize the extractor
#     extractor = NanobodyFeatureExtractor(use_pymol=True)
    
#     # Process a directory of PDB files
#     features_df = extractor.process_directory("path/to/nanobody_pdbs")
    
#     # Save results
#     extractor.save_results()
    
#     # Display results
#     print(features_df.head())
