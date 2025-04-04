from graphein.ml import InMemoryProteinGraphDataset
from graphein.ml.conversion import GraphFormatConvertor
import graphein.protein as gp
import os
from functools import partial
import networkx as nx
import traceback

# Your existing configuration
dist_edge_func = {"edge_construction_functions": [partial(gp.add_distance_threshold, threshold=5, long_interaction_threshold=0)]}
one_hot = {"node_metadata_functions" : [gp.amino_acid_one_hot, gp.meiler_embedding,
                                        partial(gp.expasy_protein_scale, add_separate=True)]}
config_1A = gp.ProteinGraphConfig(**{**dist_edge_func, **one_hot})

convertor = GraphFormatConvertor(src_format="nx", dst_format="pyg", columns=["coords", "edge_index",
                                                                             "amino_acid_one_hot", "bulkiness",
                                                                             "meiler", "rsa", "pka_rgroup",
                                                                             "isoelectric_points",
                                                                             "polaritygrantham",
                                                                             "hphob_black", "transmembranetendency"])

# Get paths to all your PDB files
pdb_dir = os.path.expanduser("~/Downloads/nanobody_extracted")  # Update this path
pdb_paths = [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f.endswith('.pdb')]

# Create label map
label_map = {os.path.splitext(os.path.basename(path))[0]: 1 if "nanobody" in path else 0
             for path in pdb_paths}

# Modified approach: Process graphs one by one, but during dataset creation
def create_dataset_with_error_tracking():
    problematic_files = []
    working_files = []

    # Create a dummy empty dataset first to set up directories
    try:
        # Create the empty dataset object
        dataset = InMemoryProteinGraphDataset(
            root="data/",
            name="train",
            paths=[],  # Start with no paths
            graph_label_map={},
            graphein_config=config_1A,
            graph_format_convertor=convertor,
            graph_transformation_funcs=[]
        )

        # Now process each file individually and add to the dataset
        for pdb_path in pdb_paths:
            filename = os.path.basename(pdb_path)
            print(f"Processing {filename}...")

            try:
                # Construct the graph directly using graphein
                graph = gp.construct_graph(config=config_1A, pdb_path=pdb_path)

                # Check for pka_rgroup in the graph nodes
                has_pka_rgroup = True
                for node in graph.nodes:
                    if 'pka_rgroup' not in graph.nodes[node]:
                        has_pka_rgroup = False
                        break

                if not has_pka_rgroup:
                    raise ValueError(f"pka_rgroup missing from some nodes in {filename}")

                # Convert the graph using your convertor
                converted_graph = convertor(graph)

                # This file worked
                working_files.append(filename)
                print(f"✓ Successfully processed {filename}")

            except Exception as e:
                # Track the problematic file and error
                problematic_files.append((filename, str(e)))
                print(f"✗ Error processing {filename}: {str(e)}")

        # Report results
        print(f"\nResults: {len(working_files)} working files, {len(problematic_files)} problematic files")

        # Save lists of files to disk
        with open("working_pdb_files.txt", "w") as f:
            for name in working_files:
                f.write(f"{name}\n")

        with open("problematic_pdb_files.txt", "w") as f:
            f.write("Problematic PDB files:\n")
            f.write("="*50 + "\n")
            for name, error in problematic_files:
                f.write(f"{name}: {error}\n")

        print("Lists saved to working_pdb_files.txt and problematic_pdb_files.txt")

        # If we have working files, create a dataset with them
        if working_files:
            # Filter paths to only working files
            working_paths = [p for p in pdb_paths if os.path.basename(p) in working_files]

            # Create the final dataset
            final_dataset = InMemoryProteinGraphDataset(
                root="data_final/",
                name="train",
                paths=working_paths,
                graph_label_map={k: label_map[k] for k in label_map if k in [os.path.splitext(os.path.basename(p))[0] for p in working_paths]},
                graphein_config=config_1A,
                graph_format_convertor=convertor,
                graph_transformation_funcs=[],
            )

            print(f"Final dataset created with {len(working_paths)} files.")
            return final_dataset
        else:
            print("No working files found, cannot create dataset.")
            return None

    except Exception as e:
        print(f"Error setting up dataset: {str(e)}")
        traceback.print_exc()

        # Let's try a different approach if the first one fails
        print("\nAttempting alternate approach to identify problematic files...")
        problematic_files = []
        working_files = []

        for pdb_path in pdb_paths:
            filename = os.path.basename(pdb_path)
            print(f"Processing {filename}...")

            try:
                # Just try to construct the graph directly
                graph = gp.construct_graph(config=config_1A, path=pdb_path)

                # Check each node
                missing_pka = False
                for node in graph.nodes:
                    if 'pka_rgroup' not in graph.nodes[node]:
                        missing_pka = True
                        break

                if missing_pka:
                    problematic_files.append((filename, "Missing pka_rgroup property"))
                    print(f"✗ {filename}: Missing pka_rgroup property")
                else:
                    working_files.append(filename)
                    print(f"✓ {filename}: All properties present")

            except Exception as e:
                problematic_files.append((filename, str(e)))
                print(f"✗ {filename}: Error: {str(e)}")

        # Save results
        with open("working_pdb_files.txt", "w") as f:
            for name in working_files:
                f.write(f"{name}\n")

        with open("problematic_pdb_files.txt", "w") as f:
            f.write("Problematic PDB files:\n")
            f.write("="*50 + "\n")
            for name, error in problematic_files:
                f.write(f"{name}: {error}\n")

        print(f"\nResults: {len(working_files)} working files, {len(problematic_files)} problematic files")
        print("Lists saved to working_pdb_files.txt and problematic_pdb_files.txt")

        return None

# Run the modified approach
dataset = create_dataset_with_error_tracking()
if dataset is not None:
    print("Dataset creation successful!")
else:
    print("Dataset creation failed.")