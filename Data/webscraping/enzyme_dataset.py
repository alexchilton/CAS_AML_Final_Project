#!/usr/bin/env python3
"""
Script to collect enzyme structures with pH data, organized by category
Saves PDB files in separate folders based on enzyme families
"""

import requests
import pandas as pd
import time
import os
from Bio.PDB import PDBList

# Enzyme database organized by categories
ENZYME_CATEGORIES = {
    "acid_proteases": {
        "pepsin_human": {"ec": "3.4.23.1", "ph_opt": 2.0, "pdb": "1PSO"},
        "pepsin_porcine": {"ec": "3.4.23.1", "ph_opt": 2.0, "pdb": "4PEP"},
        "pepsin_a": {"ec": "3.4.23.1", "ph_opt": 1.8, "pdb": "5PEP"},
        "renin": {"ec": "3.4.23.15", "ph_opt": 3.5, "pdb": "1BIL"},
        "cathepsin_d": {"ec": "3.4.23.5", "ph_opt": 3.5, "pdb": "1LYA"},
        "plasmepsin": {"ec": "3.4.23.39", "ph_opt": 4.5, "pdb": "1SME"},
        "aspartic_protease": {"ec": "3.4.23.-", "ph_opt": 3.0, "pdb": "1APT"},
        "hiv_protease": {"ec": "3.4.23.16", "ph_opt": 4.5, "pdb": "1A30"},
    },
    
    "acid_phosphatases": {
        "acid_phosphatase_human": {"ec": "3.1.3.2", "ph_opt": 5.0, "pdb": "1RPT"},
        "acid_phosphatase_yeast": {"ec": "3.1.3.2", "ph_opt": 5.5, "pdb": "1CVI"},
        "purple_acid_phosphatase": {"ec": "3.1.3.2", "ph_opt": 5.0, "pdb": "1UTE"},
    },
    
    "glycosidases": {
        "glucoamylase": {"ec": "3.2.1.3", "ph_opt": 4.5, "pdb": "1AYX"},
        "beta_glucosidase": {"ec": "3.2.1.21", "ph_opt": 5.0, "pdb": "1GNX"},
        "alpha_glucosidase": {"ec": "3.2.1.20", "ph_opt": 5.5, "pdb": "1OBB"},
        "cellulase_endo": {"ec": "3.2.1.4", "ph_opt": 5.0, "pdb": "1CEM"},
        "cellulase_exo": {"ec": "3.2.1.91", "ph_opt": 5.0, "pdb": "1CEL"},
        "xylanase_family10": {"ec": "3.2.1.8", "ph_opt": 5.5, "pdb": "1XYN"},
        "xylanase_family11": {"ec": "3.2.1.8", "ph_opt": 5.0, "pdb": "1XYF"},
        "xylanase_thermophilic": {"ec": "3.2.1.8", "ph_opt": 6.0, "pdb": "1BVV"},
    },
    
    "lysozymes": {
        "lysozyme_chicken": {"ec": "3.2.1.17", "ph_opt": 5.2, "pdb": "1DPX"},
        "lysozyme_human": {"ec": "3.2.1.17", "ph_opt": 5.0, "pdb": "1REX"},
        "lysozyme_t4": {"ec": "3.2.1.17", "ph_opt": 5.5, "pdb": "3LZM"},
        "muramidase": {"ec": "3.2.1.17", "ph_opt": 5.5, "pdb": "1LYS"},
    },
    
    "lipases": {
        "gastric_lipase": {"ec": "3.1.1.3", "ph_opt": 4.0, "pdb": "1HLG"},
        "pancreatic_lipase": {"ec": "3.1.1.3", "ph_opt": 8.0, "pdb": "1LPA"},
        "lipase_fungal": {"ec": "3.1.1.3", "ph_opt": 6.0, "pdb": "1TIC"},
        "lipase_bacterial": {"ec": "3.1.1.3", "ph_opt": 7.0, "pdb": "1CVL"},
    },
    
    "oxidoreductases": {
        "glucose_oxidase": {"ec": "1.1.3.4", "ph_opt": 5.5, "pdb": "1CF3"},
        "lactate_oxidase": {"ec": "1.1.3.15", "ph_opt": 6.5, "pdb": "2E77"},
        "alcohol_oxidase": {"ec": "1.1.3.13", "ph_opt": 7.5, "pdb": "1AHZ"},
        "catalase_human": {"ec": "1.11.1.6", "ph_opt": 7.0, "pdb": "1DGF"},
        "catalase_bacterial": {"ec": "1.11.1.6", "ph_opt": 7.0, "pdb": "1IPH"},
        "catalase_fungal": {"ec": "1.11.1.6", "ph_opt": 7.0, "pdb": "1A4E"},
        "peroxidase_horseradish": {"ec": "1.11.1.7", "ph_opt": 6.8, "pdb": "1ATJ"},
        "peroxidase_cytochrome": {"ec": "1.11.1.5", "ph_opt": 5.5, "pdb": "1CYF"},
        "myeloperoxidase": {"ec": "1.11.1.7", "ph_opt": 7.0, "pdb": "1DNU"},
    },
    
    "amylases": {
        "alpha_amylase_human": {"ec": "3.2.1.1", "ph_opt": 6.9, "pdb": "1BLI"},
        "alpha_amylase_bacterial": {"ec": "3.2.1.1", "ph_opt": 6.5, "pdb": "1BAG"},
        "alpha_amylase_fungal": {"ec": "3.2.1.1", "ph_opt": 5.5, "pdb": "2GUY"},
        "beta_amylase": {"ec": "3.2.1.2", "ph_opt": 5.5, "pdb": "1BYA"},
    },
    
    "dehydrogenases": {
        "lactate_dehydrogenase": {"ec": "1.1.1.27", "ph_opt": 7.5, "pdb": "1LDM"},
        "malate_dehydrogenase": {"ec": "1.1.1.37", "ph_opt": 7.5, "pdb": "1MLD"},
        "alcohol_dehydrogenase": {"ec": "1.1.1.1", "ph_opt": 8.5, "pdb": "1CDO"},
        "glucose_dehydrogenase": {"ec": "1.1.1.47", "ph_opt": 7.5, "pdb": "1GCO"},
    },
    
    "kinases": {
        "hexokinase": {"ec": "2.7.1.1", "ph_opt": 7.5, "pdb": "1HKG"},
        "pyruvate_kinase": {"ec": "2.7.1.40", "ph_opt": 7.5, "pdb": "1PKM"},
        "creatine_kinase": {"ec": "2.7.3.2", "ph_opt": 7.0, "pdb": "1CRK"},
    },
    
    "alkaline_proteases": {
        "trypsin_bovine": {"ec": "3.4.21.4", "ph_opt": 8.0, "pdb": "1S0Q"},
        "trypsin_human": {"ec": "3.4.21.4", "ph_opt": 8.0, "pdb": "1TRN"},
        "chymotrypsin": {"ec": "3.4.21.1", "ph_opt": 8.0, "pdb": "1CHG"},
        "elastase": {"ec": "3.4.21.36", "ph_opt": 8.5, "pdb": "1ELA"},
        "thrombin": {"ec": "3.4.21.5", "ph_opt": 8.0, "pdb": "1HRT"},
        "subtilisin_bpn": {"ec": "3.4.21.62", "ph_opt": 10.5, "pdb": "1CSE"},
        "subtilisin_carlsberg": {"ec": "3.4.21.62", "ph_opt": 10.0, "pdb": "1SBC"},
        "subtilisin_e": {"ec": "3.4.21.62", "ph_opt": 10.5, "pdb": "1SCB"},
    },
    
    "alkaline_phosphatases": {
        "alkaline_phosphatase_ecoli": {"ec": "3.1.3.1", "ph_opt": 9.0, "pdb": "1ALK"},
        "alkaline_phosphatase_calf": {"ec": "3.1.3.1", "ph_opt": 9.5, "pdb": "1AJA"},
        "alkaline_phosphatase_shrimp": {"ec": "3.1.3.1", "ph_opt": 8.5, "pdb": "1SHN"},
    },
    
    "carbonic_anhydrases": {
        "carbonic_anhydrase_i": {"ec": "4.2.1.1", "ph_opt": 7.0, "pdb": "1AZM"},
        "carbonic_anhydrase_ii": {"ec": "4.2.1.1", "ph_opt": 7.5, "pdb": "1CA2"},
        "carbonic_anhydrase_iii": {"ec": "4.2.1.1", "ph_opt": 8.0, "pdb": "1FLJ"},
    },
    
    "extremophiles": {
        "thermolysin": {"ec": "3.4.24.27", "ph_opt": 8.0, "pdb": "1LNF"},
        "psychrophilic_protease": {"ec": "3.4.21.112", "ph_opt": 8.0, "pdb": "1SH7"},
        "halophilic_protease": {"ec": "3.4.21.61", "ph_opt": 8.5, "pdb": "1NLN"},
        "thermophilic_xylanase": {"ec": "3.2.1.8", "ph_opt": 6.0, "pdb": "1BVV"},
        "acidophilic_protease": {"ec": "3.4.23.34", "ph_opt": 2.5, "pdb": "1PSO"},
    },
    
    "miscellaneous": {
        "papain": {"ec": "3.4.22.2", "ph_opt": 6.5, "pdb": "1PPN"},
        "bromelain": {"ec": "3.4.22.32", "ph_opt": 7.0, "pdb": "1W0Q"},
        "ficin": {"ec": "3.4.22.3", "ph_opt": 7.0, "pdb": "4YYW"},
        "invertase": {"ec": "3.2.1.26", "ph_opt": 4.5, "pdb": "4EQV"},
        "urease": {"ec": "3.5.1.5", "ph_opt": 7.0, "pdb": "1UBP"},
        "arginase": {"ec": "3.5.3.1", "ph_opt": 9.5, "pdb": "1RLA"},
        "asparaginase": {"ec": "3.5.1.1", "ph_opt": 8.5, "pdb": "1NNS"},
    }
}

def download_pdb(pdb_id: str, output_dir: str):
    """Download PDB file"""
    pdb_list = PDBList()
    try:
        file_path = pdb_list.retrieve_pdb_file(
            pdb_id, 
            pdir=output_dir, 
            file_format='pdb',
            overwrite=True
        )
        if file_path and os.path.exists(file_path):
            new_name = os.path.join(output_dir, f"{pdb_id.upper()}.pdb")
            if os.path.exists(new_name):
                os.remove(new_name)
            os.rename(file_path, new_name)
            return True
    except Exception as e:
        print(f"    Download error for {pdb_id}: {e}")
        return False

def main():
    # Create main output directory
    base_dir = "enzyme_ph_structures"
    os.makedirs(base_dir, exist_ok=True)
    
    # Collect data
    dataset = []
    downloaded_count = 0
    total_enzymes = sum(len(enzymes) for enzymes in ENZYME_CATEGORIES.values())
    
    print(f"Collecting structures for {total_enzymes} enzymes in {len(ENZYME_CATEGORIES)} categories...\n")
    
    for category, enzymes in ENZYME_CATEGORIES.items():
        # Create category folder
        category_dir = os.path.join(base_dir, category)
        os.makedirs(category_dir, exist_ok=True)
        
        print(f"\nCategory: {category.replace('_', ' ').title()} ({len(enzymes)} enzymes)")
        print("-" * 50)
        
        for enzyme_name, info in enzymes.items():
            print(f"  {enzyme_name} (pH {info['ph_opt']})")
            
            pdb_id = info["pdb"]
            
            # Add to dataset
            dataset.append({
                'category': category,
                'enzyme_name': enzyme_name,
                'ec_number': info["ec"],
                'pdb_id': pdb_id,
                'ph_optimum': info["ph_opt"]
            })
            
            # Download structure to category folder
            if download_pdb(pdb_id, category_dir):
                downloaded_count += 1
                print(f"    ✓ Downloaded {pdb_id}")
            else:
                print(f"    ✗ Failed to download {pdb_id}")
            
            time.sleep(0.5)
    
    # Save dataset
    df = pd.DataFrame(dataset)
    df.to_csv("enzyme_ph_dataset.csv", index=False)
    
    # Create summary by category
    category_summary = df.groupby('category').agg({
        'enzyme_name': 'count',
        'ph_optimum': ['min', 'max', 'mean']
    })
    category_summary.columns = ['enzyme_count', 'ph_min', 'ph_max', 'ph_mean']
    category_summary.to_csv("category_summary.csv")
    
    # Create pH distribution
    ph_counts = df['ph_optimum'].value_counts().sort_index()
    ph_distribution = pd.DataFrame({
        'pH': ph_counts.index,
        'count': ph_counts.values
    })
    ph_distribution.to_csv("ph_distribution.csv", index=False)
    
    print(f"\n{'='*50}")
    print(f"Dataset complete!")
    print(f"Total enzymes: {len(dataset)}")
    print(f"Structures downloaded: {downloaded_count}")
    print(f"pH range: {df['ph_optimum'].min()} - {df['ph_optimum'].max()}")
    print(f"Unique pH values: {df['ph_optimum'].nunique()}")
    print(f"\nFiles created:")
    print(f"  - enzyme_ph_dataset.csv (main dataset)")
    print(f"  - category_summary.csv (summary by category)")
    print(f"  - ph_distribution.csv (pH distribution)")
    print(f"  - {base_dir}/ (main folder)")
    for category in ENZYME_CATEGORIES.keys():
        print(f"    - {category}/ (subfolder)")
    print(f"{'='*50}")
    
    # Print category summary
    print("\nCategory Summary:")
    for category in ENZYME_CATEGORIES.keys():
        cat_data = df[df['category'] == category]
        print(f"  {category.replace('_', ' ').title()}: {len(cat_data)} enzymes, pH {cat_data['ph_optimum'].min()}-{cat_data['ph_optimum'].max()}")

if __name__ == "__main__":
    main()