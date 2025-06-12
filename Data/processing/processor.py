import os
import numpy as np
import pickle
import json
import traceback
import gc
import geometry, features, parsing, io_utils, graph

class ProteinProcessor:
    def __init__(self, folder_path: str, output_path: str, max_file_size_mb: float = 25):
        self.folder_path = folder_path
        self.output_path = output_path
        self.max_file_size_mb = max_file_size_mb

    def process_pdb_full_features_memory_efficient(self):
        """
        Process PDB files with all features while being memory-efficient.
        
        Parameters:
        -----------
        folder_path : str
            Path to folder with PDB files
        output_path : str, optional
            Output directory
        max_file_size_mb : float
            Skip files larger than this
            
        Returns:
        --------
        dict
            Processing summary
        """
        if output_path is None:
            output_path = os.path.join(self.folder_path, "processed_data")
        
        os.makedirs(output_path, exist_ok=True)
        
        # Get PDB files
        pdb_files = [f for f in os.listdir(self.folder_path) if f.endswith('.pdb')]
        print(f"Found {len(pdb_files)} PDB files")
        
        summary = {}
        
        for pdb_file in pdb_files:
            pdb_path = os.path.join(self.folder_path, pdb_file)
            pdb_id = os.path.splitext(pdb_file)[0]
            
            # Create PDB-specific output directory
            pdb_output = os.path.join(output_path, pdb_id)
            os.makedirs(pdb_output, exist_ok=True)
            
            try:
                # Check file size
                file_size_mb = os.path.getsize(pdb_path) / (1024 * 1024)
                if file_size_mb > self.max_file_size_mb:
                    print(f"Skipping {pdb_file} - too large ({file_size_mb:.2f} MB)")
                    summary[pdb_id] = {"status": "skipped", "size_mb": file_size_mb}
                    continue
                
                print(f"Processing {pdb_file} ({file_size_mb:.2f} MB)")
                
                # Step 1: Parse basic structure and create graph (low memory impact)
                structure, residues_by_chain = parsing.parse_basic_structure(pdb_path)
                
                if not structure or not residues_by_chain:
                    print(f"Failed to parse {pdb_file}")
                    summary[pdb_id] = {"status": "error", "error": "Failed to parse structure"}
                    continue
                
                print(f"  Found {len(residues_by_chain)} chains")
                
                # Step 2: Calculate secondary structure once (low memory impact)
                ss_data = geometry.calculate_secondary_structure(pdb_path)
                print(f"  Calculated secondary structure for {len(ss_data)} residues")
                
                # Step 3: Create protein graph (moderate memory impact)
                nx_graph = graph.create_protein_graph(pdb_path)
                nx_graph = graph.add_structure_to_graph(nx_graph, ss_data)
                print(f"  Created graph with {nx_graph.number_of_nodes()} nodes and {nx_graph.number_of_edges()} edges")
                
                # Save graph immediately to free memory
                graph_path = os.path.join(pdb_output, f"{pdb_id}_graph.pkl")
                with open(graph_path, 'wb') as f:
                    pickle.dump(nx_graph, f)
                print(f"  Saved graph to {graph_path}")
                
                # Create PDB data dictionary (without the graph)
                pdb_data = {
                    'pdb_id': pdb_id,
                    'residues_by_chain': residues_by_chain,
                    'secondary_structure': {},
                    'backbone_atoms': {},
                    'edge_indices': {},
                    'bond_lengths': {},
                    'angles': {},
                    'torsions': {},
                    'charges': {},
                    'hydrophobic_residues': {}
                }
                
                # Process each chain separately with incremental save
                for chain_id in residues_by_chain:
                    print(f"  Processing chain {chain_id}...")
                    
                    # Create output path for this chain
                    chain_output = os.path.join(pdb_output, f"{pdb_id}_{chain_id}")
                    
                    # Step 4: Process one chain at a time with memory-efficient functions
                    #         Save results incrementally to disk
                    self.process_single_chain_full_features(
                        pdb_path, 
                        chain_id, 
                        chain_output, 
                        ss_data,
                        incremental_save=True
                    )
                    
                    # After processing is complete, load the chain results
                    chain_data = io_utils.load_chain_data(chain_output)
                    
                    # Add to main data dictionary
                    if chain_data:
                        for key in ['secondary_structure', 'backbone_atoms', 'edge_indices', 
                                'bond_lengths', 'angles', 'torsions', 'charges', 
                                'hydrophobic_residues']:
                            if key in chain_data:
                                pdb_data[key][chain_id] = chain_data[key]
                    
                    # Force garbage collection
                    gc.collect()
                
                # Save the main data dictionary
                data_path = os.path.join(pdb_output, f"{pdb_id}_data.pkl")
                with open(data_path, 'wb') as f:
                    pickle.dump(pdb_data, f)
                print(f"  Saved data to {data_path}")
                
                # Update summary
                summary[pdb_id] = {
                    "status": "success",
                    "chains": list(residues_by_chain.keys()),
                    "size_mb": file_size_mb,
                    "num_nodes": nx_graph.number_of_nodes(),
                    "num_edges": nx_graph.number_of_edges()
                }
                
                # Clear memory for next PDB
                del nx_graph
                del pdb_data
                gc.collect()
                
            except Exception as e:
                print(f"Error processing {pdb_file}: {str(e)}")
                traceback.print_exc()
                summary[pdb_id] = {"status": "error", "error": str(e)}
        
        # Save summary
        summary_path = os.path.join(output_path, "processing_summary.json")
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        return summary

    def process_single_chain_full_features(pdb_path, chain_id, output_prefix, ss_data, incremental_save=True):
        """
        Process a single chain with all features but in a memory-efficient way.
        
        Parameters:
        -----------
        pdb_path : str
            Path to PDB file
        chain_id : str
            Chain ID to process
        output_prefix : str
            Prefix for output files
        ss_data : dict
            Secondary structure data
        incremental_save : bool
            Whether to save results incrementally
            
        Returns:
        --------
        dict
            Chain data if not saved incrementally
        """
        try:
            # Step 1: Extract atoms for this chain only
            atoms = parsing.extract_chain_atoms(pdb_path, chain_id)
            
            if not atoms:
                print(f"  No atoms found for chain {chain_id}")
                return {} if not incremental_save else None
            
            print(f"  Extracted {len(atoms)} atoms for chain {chain_id}")
            
            # Save atoms if incremental
            if incremental_save:
                with open(f"{output_prefix}_atoms.pkl", 'wb') as f:
                    pickle.dump(atoms, f)
            
            # Step 2: Get chain residues
            chain_residues = []
            res_ids_seen = set()
            
            for atom in atoms:
                _, res_id, res_name, _, _, _ = atom
                if res_id not in res_ids_seen:
                    chain_residues.append((res_id, res_name))
                    res_ids_seen.add(res_id)
            
            # Step 3: Extract backbone atoms (low memory impact)
            backbone = geometry.extract_backbone_atoms(atoms)
            print(f"  Found backbone atoms for {len(backbone)} residues")
            
            if incremental_save:
                with open(f"{output_prefix}_backbone.pkl", 'wb') as f:
                    pickle.dump(backbone, f)
            
            # Step 4: Get secondary structure for this chain
            chain_ss = {}
            for res_id, _ in chain_residues:
                chain_ss[res_id] = ss_data.get((chain_id, res_id), '?')
            
            if incremental_save:
                with open(f"{output_prefix}_ss.pkl", 'wb') as f:
                    pickle.dump(chain_ss, f)
            
            # Step 5: Extract atomic positions (needed for multiple calculations)
            positions = np.array([atom[5] for atom in atoms])
            
            # Step 6: Create edge indices for different calculations
            print("  Creating edge indices...")
            
            # Pairs for bond lengths (with efficient memory)
            edge_index_pairs = geometry.create_edge_index_memory_efficient(atoms, positions, mode='pairs')
            
            if incremental_save:
                with open(f"{output_prefix}_edge_pairs.pkl", 'wb') as f:
                    pickle.dump(edge_index_pairs, f)
            
            # Step 7: Calculate bond lengths (moderate memory impact)
            print("  Calculating bond lengths...")
            bond_lengths = geometry.calculate_bond_lengths_efficient(atoms, edge_index_pairs)
            
            if incremental_save:
                with open(f"{output_prefix}_bond_lengths.pkl", 'wb') as f:
                    pickle.dump(bond_lengths, f)
                
                # Clear memory
                del edge_index_pairs
                gc.collect()
            
            # Step 8: Calculate angles (higher memory impact)
            print("  Calculating angles...")
            angles = geometry.calculate_angles_memory_efficient(atoms, positions, backbone_only=True)
            
            if incremental_save:
                with open(f"{output_prefix}_angles.pkl", 'wb') as f:
                    pickle.dump(angles, f)
                
                # Clear memory
                gc.collect()
            
            # Step 9: Calculate torsions (higher memory impact)
            print("  Calculating torsions...")
            torsions = geometry.calculate_torsions_memory_efficient(atoms, positions)
            
            if incremental_save:
                with open(f"{output_prefix}_torsions.pkl", 'wb') as f:
                    pickle.dump(torsions, f)
                
                # Clear memory
                gc.collect()
            
            # Step 10: Calculate charges and hydrophobicity (low memory impact)
            print("  Assigning charges and identifying hydrophobic residues...")
            charges = features.assign_charges(atoms)
            hydrophobic = features.identify_hydrophobic_residues(chain_residues)
            
            if incremental_save:
                with open(f"{output_prefix}_charges.pkl", 'wb') as f:
                    pickle.dump(charges, f)
                
                with open(f"{output_prefix}_hydrophobic.pkl", 'wb') as f:
                    pickle.dump(hydrophobic, f)
                
                # Create index file
                index = {
                    "atoms": f"{output_prefix}_atoms.pkl",
                    "backbone": f"{output_prefix}_backbone.pkl",
                    "secondary_structure": f"{output_prefix}_ss.pkl",
                    "edge_pairs": f"{output_prefix}_edge_pairs.pkl",
                    "bond_lengths": f"{output_prefix}_bond_lengths.pkl",
                    "angles": f"{output_prefix}_angles.pkl",
                    "torsions": f"{output_prefix}_torsions.pkl",
                    "charges": f"{output_prefix}_charges.pkl",
                    "hydrophobic": f"{output_prefix}_hydrophobic.pkl"
                }
                
                with open(f"{output_prefix}_index.json", 'w') as f:
                    json.dump(index, f)
                
                return None
            
            # If not incremental save, return all data
            return {
                'backbone_atoms': backbone,
                'secondary_structure': chain_ss,
                'edge_indices': edge_index_pairs,
                'bond_lengths': bond_lengths,
                'angles': angles,
                'torsions': torsions,
                'charges': charges,
                'hydrophobic_residues': hydrophobic
            }
                
        except Exception as e:
            print(f"  Error processing chain {chain_id}: {str(e)}")
            traceback.print_exc()
            return {} if not incremental_save else None

    def run(self):
        self.process_pdb_full_features_memory_efficient()
