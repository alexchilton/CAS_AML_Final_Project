import os
import pytest
import numpy as np
import networkx as nx
import pickle
import json
import tempfile
import shutil
from unittest.mock import patch, MagicMock

# Import the functions to test
from preprocess_mem_eff import (
    parse_pdb,
    parse_basic_structure,
    calculate_secondary_structure,
    create_protein_graph,
    add_structure_to_graph,
    extract_backbone_atoms,
    extract_chain_atoms,
    calculate_bond_lengths_efficient,
    calculate_angles_memory_efficient,
    calculate_torsions_memory_efficient,
    calculate_single_torsion,
    assign_charges,
    identify_hydrophobic_residues,
    create_edge_index_memory_efficient,
    add_backbone_from_pdb,
    process_single_chain_full_features,
    load_chain_data,
    process_pdb_full_features_memory_efficient
)

# Fixtures for testing

@pytest.fixture
def test_pdb_file():
    """Create a minimal PDB file for testing"""
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        # Write a minimal PDB content
        tmp.write(b"""ATOM      1  N   VAL A   1      11.804  37.253  30.159  1.00  0.00           N
ATOM      2  CA  VAL A   1      11.870  36.283  31.232  1.00  0.00           C
ATOM      3  C   VAL A   1      11.031  36.741  32.425  1.00  0.00           C
ATOM      4  O   VAL A   1      11.211  37.862  32.905  1.00  0.00           O
ATOM      5  CB  VAL A   1      13.330  36.156  31.705  1.00  0.00           C
ATOM      6  N   PRO A   2      10.111  35.876  32.897  1.00  0.00           N
ATOM      7  CA  PRO A   2       9.303  36.104  34.094  1.00  0.00           C
ATOM      8  C   PRO A   2      10.129  36.067  35.369  1.00  0.00           C
ATOM      9  O   PRO A   2      10.988  35.210  35.489  1.00  0.00           O
ATOM     10  N   GLY B   1      37.068  26.613  37.468  1.00  0.00           N
ATOM     11  CA  GLY B   1      37.843  26.110  36.344  1.00  0.00           C
ATOM     12  C   GLY B   1      39.143  26.876  36.150  1.00  0.00           C
ATOM     13  O   GLY B   1      39.688  27.436  37.111  1.00  0.00           O
END
""")
        tmp_name = tmp.name

    yield tmp_name

    # Clean up
    if os.path.exists(tmp_name):
        os.unlink(tmp_name)

@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for output files"""
    output_dir = tempfile.mkdtemp()
    yield output_dir
    # Clean up
    shutil.rmtree(output_dir)

@pytest.fixture
def mock_graph():
    """Create a mock NetworkX graph"""
    G = nx.Graph()

    # Add some nodes
    G.add_node("A:VAL:1", chain_id="A", residue_number=1, residue_name="VAL",
               coords=np.array([1.0, 2.0, 3.0]), CA=np.array([1.0, 2.0, 3.0]))
    G.add_node("A:PRO:2", chain_id="A", residue_number=2, residue_name="PRO",
               coords=np.array([2.0, 3.0, 4.0]), CA=np.array([2.0, 3.0, 4.0]))
    G.add_node("B:GLY:1", chain_id="B", residue_number=1, residue_name="GLY",
               coords=np.array([5.0, 6.0, 7.0]), CA=np.array([5.0, 6.0, 7.0]))

    # Add some edges
    G.add_edge("A:VAL:1", "A:PRO:2", kind={"peptide_bond"})
    G.add_edge("A:VAL:1", "B:GLY:1", kind={"contact"}, distance=6.5)

    return G

# Tests for parsing functions

def test_parse_pdb(test_pdb_file):
    """Test the parse_pdb function with a minimal PDB file"""
    atoms, residues_by_chain = parse_pdb(test_pdb_file)

    # Check we parsed the expected atoms
    assert len(atoms) == 13

    # Check the chain information
    assert set(residues_by_chain.keys()) == {"A", "B"}
    assert len(residues_by_chain["A"]) == 2  # VAL1, PRO2
    assert len(residues_by_chain["B"]) == 1  # GLY1

    # Check first atom structure
    first_atom = atoms[0]
    assert first_atom[0] == "A"  # chain_id
    assert first_atom[1] == 1    # residue_id
    assert first_atom[2] == "VAL" # residue_name
    assert first_atom[3] == "N"   # atom_name
    assert first_atom[4] == "N"   # element
    assert isinstance(first_atom[5], np.ndarray) # position

def test_parse_basic_structure(test_pdb_file):
    """Test the parse_basic_structure function with a minimal PDB file"""
    structure, residues_by_chain = parse_basic_structure(test_pdb_file)

    # Check we got a structure object
    assert structure is not None

    # Check the residue information
    assert set(residues_by_chain.keys()) == {"A", "B"}
    assert len(residues_by_chain["A"]) == 2  # VAL1, PRO2
    assert len(residues_by_chain["B"]) == 1  # GLY1

    # Check first residue in chain A
    assert residues_by_chain["A"][0] == (1, "VAL")

def test_extract_chain_atoms(test_pdb_file):
    """Test the extract_chain_atoms function for a specific chain"""
    chain_atoms = extract_chain_atoms(test_pdb_file, "A")

    # Check we got only chain A atoms
    assert len(chain_atoms) == 9  # 5 from VAL, 4 from PRO

    # All atoms should be from chain A
    for atom in chain_atoms:
        assert atom[0] == "A"

    # Check specific residues
    res_ids = set(atom[1] for atom in chain_atoms)
    assert res_ids == {1, 2}  # VAL1, PRO2

# Tests for structure and graph functions

def test_create_protein_graph(test_pdb_file):
    """Test the create_protein_graph function"""
    # Mock the graphein import to ensure we use the fallback
    with patch('preprocess_mem_eff.gp', side_effect=ImportError("No graphein")):
        graph = create_protein_graph(test_pdb_file)

    # Check the graph structure
    assert isinstance(graph, nx.Graph)
    assert graph.number_of_nodes() == 3  # 3 residues

    # Check node attributes
    for node in graph.nodes():
        data = graph.nodes[node]
        assert 'chain_id' in data
        assert 'residue_number' in data
        assert 'residue_name' in data
        assert 'has_backbone' in data

    # Check we have the expected nodes
    node_ids = set(graph.nodes())
    expected_nodes = {"A:VAL:1", "A:PRO:2", "B:GLY:1"}
    assert node_ids == expected_nodes

@patch('preprocess_mem_eff.DSSP')
def test_calculate_secondary_structure(mock_dssp, test_pdb_file):
    """Test the calculate_secondary_structure function with a mocked DSSP"""
    # Create a mock DSSP result
    mock_dssp_instance = MagicMock()
    mock_dssp_instance.keys.return_value = [('A', ('_', 1, '_')), ('A', ('_', 2, '_')), ('B', ('_', 1, '_'))]
    mock_dssp_instance.__getitem__.side_effect = lambda key: {'secstruct': 'H'} if key[0] == 'A' else {'secstruct': 'E'}
    mock_dssp.return_value = mock_dssp_instance

    ss_data = calculate_secondary_structure(test_pdb_file)

    # Check we have SS data for all residues
    assert len(ss_data) == 3
    assert ss_data[('A', 1)] == 'H'
    assert ss_data[('A', 2)] == 'H'
    assert ss_data[('B', 1)] == 'E'

def test_add_structure_to_graph(mock_graph):
    """Test the add_structure_to_graph function"""
    # Create mock secondary structure data
    ss_data = {
        ('A', 1): 'H',
        ('A', 2): 'H',
        ('B', 1): 'E'
    }

    # Add SS info to the graph
    graph = add_structure_to_graph(mock_graph, ss_data)

    # Check the SS attributes were added
    assert graph.nodes["A:VAL:1"]["ss"] == 'H'
    assert graph.nodes["A:PRO:2"]["ss"] == 'H'
    assert graph.nodes["B:GLY:1"]["ss"] == 'E'

    # Check backbone flag is correctly set
    for node in graph.nodes():
        assert 'has_backbone' in graph.nodes[node]

def test_add_backbone_from_pdb(mock_graph, test_pdb_file):
    """Test the add_backbone_from_pdb function"""
    # Start with a graph without backbone info
    for node in mock_graph.nodes():
        mock_graph.nodes[node]['has_backbone'] = False
        for key in ['CA', 'N', 'C', 'O']:
            if key in mock_graph.nodes[node]:
                del mock_graph.nodes[node][key]

    # Add backbone info from PDB
    graph = add_backbone_from_pdb(mock_graph, test_pdb_file)

    # Check backbone info was added for at least some nodes
    backbone_present = sum(1 for node, data in graph.nodes(data=True) if data.get('has_backbone', False))
    assert backbone_present > 0

# Tests for calculation functions

def test_extract_backbone_atoms():
    """Test the extract_backbone_atoms function"""
    # Create mock chain atoms
    chain_atoms = [
        ("A", 1, "VAL", "N", "N", np.array([1.0, 2.0, 3.0])),
        ("A", 1, "VAL", "CA", "C", np.array([2.0, 3.0, 4.0])),
        ("A", 1, "VAL", "C", "C", np.array([3.0, 4.0, 5.0])),
        ("A", 1, "VAL", "O", "O", np.array([4.0, 5.0, 6.0])),
        ("A", 1, "VAL", "CB", "C", np.array([5.0, 6.0, 7.0])),
        ("A", 2, "PRO", "N", "N", np.array([6.0, 7.0, 8.0])),
        ("A", 2, "PRO", "CA", "C", np.array([7.0, 8.0, 9.0]))
    ]

    backbone = extract_backbone_atoms(chain_atoms)

    # Check we have backbone atoms for both residues
    assert set(backbone.keys()) == {1, 2}

    # Check first residue has all backbone atoms
    assert set(backbone[1].keys()) == {"N", "CA", "C", "O"}

    # Check second residue has partial backbone
    assert set(backbone[2].keys()) == {"N", "CA"}

def test_create_edge_index_memory_efficient():
    """Test the create_edge_index_memory_efficient function"""
    # Create mock chain atoms
    chain_atoms = [
        ("A", 1, "VAL", "N", "N", np.array([1.0, 2.0, 3.0])),
        ("A", 1, "VAL", "CA", "C", np.array([2.0, 3.0, 4.0])),
        ("A", 2, "PRO", "N", "N", np.array([10.0, 11.0, 12.0])),
        ("A", 2, "PRO", "CA", "C", np.array([11.0, 12.0, 13.0])),
        ("A", 3, "GLY", "N", "N", np.array([20.0, 21.0, 22.0]))
    ]

    positions = np.array([atom[5] for atom in chain_atoms])

    # Test with different modes
    edge_index = create_edge_index_memory_efficient(chain_atoms, positions, mode='pairs', distance_cutoff=10.0)

    # Check we got the expected edges
    assert isinstance(edge_index, np.ndarray)
    assert edge_index.shape[0] == 2  # [source, target] format

    # There should be connections within same residue and between residues 1-2
    # But not between residues 1-3 or 2-3 due to distance cutoff
    edge_count = edge_index.shape[1]
    assert edge_count > 0

def test_calculate_bond_lengths_efficient():
    """Test the calculate_bond_lengths_efficient function"""
    # Create mock chain atoms
    chain_atoms = [
        ("A", 1, "VAL", "N", "N", np.array([1.0, 2.0, 3.0])),
        ("A", 1, "VAL", "CA", "C", np.array([2.0, 3.0, 4.0])),
        ("A", 1, "VAL", "C", "C", np.array([3.0, 4.0, 5.0]))
    ]

    # Create edge index
    edge_index = np.array([[0, 1, 0, 2, 1, 2], [1, 0, 2, 0, 2, 1]])

    # Calculate bond lengths
    bond_lengths = calculate_bond_lengths_efficient(chain_atoms, edge_index)

    # Check we got the expected bond lengths
    assert len(bond_lengths) == 3  # Three unique pairs

    # Check one specific bond length (N-CA)
    assert ((0, 1) in bond_lengths) or ((1, 0) in bond_lengths)
    if (0, 1) in bond_lengths:
        n_ca_length = bond_lengths[(0, 1)]
    else:
        n_ca_length = bond_lengths[(1, 0)]

    # Should be approximately sqrt(3) since positions differ by 1 in each dimension
    assert np.isclose(n_ca_length, np.sqrt(3))

def test_calculate_angles_memory_efficient():
    """Test the calculate_angles_memory_efficient function"""
    # Create mock chain atoms with positions forming a right angle
    chain_atoms = [
        ("A", 1, "VAL", "N", "N", np.array([0.0, 0.0, 0.0])),
        ("A", 1, "VAL", "CA", "C", np.array([1.0, 0.0, 0.0])),
        ("A", 1, "VAL", "C", "C", np.array([1.0, 1.0, 0.0]))
    ]

    positions = np.array([atom[5] for atom in chain_atoms])

    # Calculate angles
    angles = calculate_angles_memory_efficient(chain_atoms, positions, distance_cutoff=2.0)

    # Check we got some angles
    assert len(angles) > 0

    # Check if we have the right angle (N-CA-C = 90 degrees)
    angle_key = (0, 1, 2)
    if angle_key in angles:
        assert np.isclose(angles[angle_key], 90.0, atol=1e-6)

def test_calculate_single_torsion():
    """Test the calculate_single_torsion function"""
    # Create positions that form a known torsion angle
    # Trans configuration = 180 degrees
    pos_i = np.array([0.0, 0.0, 0.0])
    pos_j = np.array([1.0, 0.0, 0.0])
    pos_k = np.array([2.0, 0.0, 0.0])
    pos_l = np.array([3.0, 0.0, 0.0])

    # Calculate torsion
    torsion = calculate_single_torsion(pos_i, pos_j, pos_k, pos_l)

    # Should be 180 or -180 degrees for trans
    assert np.isclose(abs(torsion), 180.0, atol=1e-6)

    # Test with a cis configuration = 0 degrees
    pos_l_cis = np.array([2.0, 1.0, 0.0])
    torsion_cis = calculate_single_torsion(pos_i, pos_j, pos_k, pos_l_cis)

    # Should be 0 degrees for cis
    assert np.isclose(torsion_cis, 0.0, atol=1e-6)

def test_calculate_torsions_memory_efficient():
    """Test the calculate_torsions_memory_efficient function"""
    # Create mock chain atoms forming a simple backbone-like structure
    chain_atoms = [
        ("A", 1, "VAL", "N", "N", np.array([0.0, 0.0, 0.0])),
        ("A", 1, "VAL", "CA", "C", np.array([1.0, 0.0, 0.0])),
        ("A", 1, "VAL", "C", "C", np.array([2.0, 0.0, 0.0])),
        ("A", 2, "PRO", "N", "N", np.array([3.0, 0.0, 0.0])),
        ("A", 2, "PRO", "CA", "C", np.array([4.0, 0.0, 0.0]))
    ]

    positions = np.array([atom[5] for atom in chain_atoms])

    # Calculate torsions with backbone focus
    torsions = calculate_torsions_memory_efficient(chain_atoms, positions,
                                                   distance_cutoff=2.0, backbone_only=True)

    # Check we got some torsions
    assert len(torsions) > 0

def test_assign_charges():
    """Test the assign_charges function"""
    # Create mock chain atoms
    chain_atoms = [
        ("A", 1, "VAL", "N", "N", np.array([1.0, 2.0, 3.0])),
        ("A", 1, "VAL", "CA", "C", np.array([2.0, 3.0, 4.0])),
        ("A", 1, "VAL", "O", "O", np.array([3.0, 4.0, 5.0])),
        ("A", 1, "VAL", "CB", "C", np.array([4.0, 5.0, 6.0]))
    ]

    charges = assign_charges(chain_atoms)

    # Check we have charges for all atoms
    assert len(charges) == 4

    # Check specific charge assignments from charge_map
    assert charges[0] == -0.5  # N
    assert charges[1] == 0.5   # C (CA)
    assert charges[2] == -0.5  # O
    assert charges[3] == 0.5   # C (CB)

def test_identify_hydrophobic_residues():
    """Test the identify_hydrophobic_residues function"""
    # Create mock chain residues
    chain_residues = [
        (1, "ALA"),  # Hydrophobic
        (2, "ASP"),  # Hydrophilic
        (3, "VAL"),  # Hydrophobic
        (4, "LYS")   # Hydrophilic
    ]

    hydrophobic = identify_hydrophobic_residues(chain_residues)

    # Check we identified the correct hydrophobic residues
    assert hydrophobic == {1, 3}

# Tests for processing functions

@patch('preprocess_mem_eff.calculate_secondary_structure')
@patch('preprocess_mem_eff.extract_chain_atoms')
def test_process_single_chain_full_features(mock_extract_chain, mock_calc_ss, temp_output_dir):
    """Test the process_single_chain_full_features function"""
    # Create mock data
    mock_extract_chain.return_value = [
        ("A", 1, "VAL", "N", "N", np.array([1.0, 2.0, 3.0])),
        ("A", 1, "VAL", "CA", "C", np.array([2.0, 3.0, 4.0])),
        ("A", 2, "PRO", "N", "N", np.array([6.0, 7.0, 8.0])),
        ("A", 2, "PRO", "CA", "C", np.array([7.0, 8.0, 9.0]))
    ]

    mock_calc_ss.return_value = {
        ('A', 1): 'H',
        ('A', 2): 'H'
    }

    # Process a single chain
    output_prefix = os.path.join(temp_output_dir, "test_A")
    process_single_chain_full_features(
        "dummy.pdb", "A", output_prefix, mock_calc_ss.return_value, incremental_save=True
    )

    # Check output files were created
    expected_files = [
        "test_A_atoms.pkl",
        "test_A_backbone.pkl",
        "test_A_ss.pkl",
        "test_A_edge_pairs.pkl",
        "test_A_bond_lengths.pkl",
        "test_A_angles.pkl",
        "test_A_torsions.pkl",
        "test_A_charges.pkl",
        "test_A_hydrophobic.pkl",
        "test_A_index.json"
    ]

    for filename in expected_files:
        full_path = os.path.join(temp_output_dir, filename)
        assert os.path.exists(full_path)

def test_load_chain_data(temp_output_dir):
    """Test the load_chain_data function"""
    # Create mock chain data files
    output_prefix = os.path.join(temp_output_dir, "test_chain")

    # Create a mock index file
    index = {
        "backbone": f"{output_prefix}_backbone.pkl",
        "secondary_structure": f"{output_prefix}_ss.pkl",
        "charges": f"{output_prefix}_charges.pkl"
    }

    with open(f"{output_prefix}_index.json", 'w') as f:
        json.dump(index, f)

    # Create mock data files
    backbone_data = {1: {"CA": np.array([1.0, 2.0, 3.0])}}
    with open(f"{output_prefix}_backbone.pkl", 'wb') as f:
        pickle.dump(backbone_data, f)

    ss_data = {1: "H", 2: "E"}
    with open(f"{output_prefix}_ss.pkl", 'wb') as f:
        pickle.dump(ss_data, f)

    charges_data = {0: -0.5, 1: 0.5}
    with open(f"{output_prefix}_charges.pkl", 'wb') as f:
        pickle.dump(charges_data, f)

    # Load the chain data
    chain_data = load_chain_data(output_prefix)

    # Check we loaded the expected data
    assert 'backbone_atoms' in chain_data
    assert 'secondary_structure' in chain_data
    assert 'charges' in chain_data

    # Check the content matches what we saved
    assert chain_data['backbone_atoms'] == backbone_data
    assert chain_data['secondary_structure'] == ss_data
    assert chain_data['charges'] == charges_data

# Integration test for the main processing function

@patch('preprocess_mem_eff.parse_basic_structure')
@patch('preprocess_mem_eff.calculate_secondary_structure')
@patch('preprocess_mem_eff.create_protein_graph')
def test_process_pdb_full_features_memory_efficient(mock_create_graph, mock_calc_ss, mock_parse_basic,
                                                    mock_graph, test_pdb_file, temp_output_dir):
    """Test the process_pdb_full_features_memory_efficient function"""
    # Set up mocks
    mock_parse_basic.return_value = (None, {"A": [(1, "VAL"), (2, "PRO")], "B": [(1, "GLY")]})
    mock_calc_ss.return_value = {('A', 1): 'H', ('A', 2): 'H', ('B', 1): 'E'}
    mock_create_graph.return_value = mock_graph

    # Create a test folder with a PDB file
    test_folder = os.path.join(temp_output_dir, "test_folder")
    os.makedirs(test_folder, exist_ok=True)

    test_pdb_path = os.path.join(test_folder, "test.pdb")
    shutil.copy(test_pdb_file, test_pdb_path)

    # Process the PDB
    with patch('preprocess_mem_eff.process_single_chain_full_features',
               return_value=None) as mock_process_chain:
        summary = process_pdb_full_features_memory_efficient(
            test_folder,
            output_path=os.path.join(temp_output_dir, "output"),
            max_file_size_mb=100
        )

    # Check the output
    assert "test" in summary
    assert summary["test"]["status"] == "success"

    # Check that process_single_chain_full_features was called for each chain
    assert mock_process_chain.call_count >= 2  # Once for chain A, once for chain B

# Optional: More specialized tests to test each feature more intensively

def test_memory_efficiency_with_large_data():
    """Test the memory efficiency aspects with larger data"""
    # Create larger test data
    n_atoms = 1000
    chain_atoms = []

    for i in range(n_atoms):
        res_id = i // 10 + 1  # 10 atoms per residue
        chain_atoms.append((
            "A", res_id, "ALA", f"ATOM{i}", "C",
            np.array([float(i), float(i+1), float(i+2)])
        ))

    positions = np.array([atom[5] for atom in chain_atoms])

    # Track memory usage during operation
    import tracemalloc

    tracemalloc.start()

    # Test memory-intensive operations
    edge_index = create_edge_index_memory_efficient(chain_atoms, positions,
                                                    mode='pairs', distance_cutoff=10.0)

    current, peak = tracemalloc.get_traced_memory()

    tracemalloc.stop()

    # Simple assertion that we're not using an excessive amount of memory
    # The exact values will depend on system, but we're checking the test completes
    assert edge_index.shape[0] == 2

    # We're mainly checking that the test completes without OOM errors
    print(f"Memory usage for {n_atoms} atoms - Current: {current / 1024**2:.2f}MB, Peak: {peak / 1024**2:.2f}MB")