#!/usr/bin/env python3
"""
Script to collect ONE structure per enzyme with known pH optimum
Fixed version with better search parameters
"""

import requests
import pandas as pd
import time
import os
from Bio.PDB import PDBList

# Diverse enzymes with specific search terms and known PDB examples
ENZYME_PH_DATA = {
    # Acidic pH enzymes (pH 1-4)
    "pepsin": {"ec": "3.4.23.1", "ph_opt": 2.0, "search_term": "pepsin", "example_pdb": "1PSO"},
    "acid_phosphatase": {"ec": "3.1.3.2", "ph_opt": 5.0, "search_term": "acid phosphatase", "example_pdb": "1RPT"},
    "glucoamylase": {"ec": "3.2.1.3", "ph_opt": 4.5, "search_term": "glucoamylase", "example_pdb": "1AYX"},
    
    # Slightly acidic (pH 4-6)
    "lysozyme": {"ec": "3.2.1.17", "ph_opt": 5.2, "search_term": "lysozyme", "example_pdb": "1DPX"},
    "glucose_oxidase": {"ec": "1.1.3.4", "ph_opt": 5.5, "search_term": "glucose oxidase", "example_pdb": "1CF3"},
    
    # Neutral pH enzymes (pH 6-8)
    "alpha_amylase": {"ec": "3.2.1.1", "ph_opt": 6.9, "search_term": "alpha amylase", "example_pdb": "1BLI"},
    "catalase": {"ec": "1.11.1.6", "ph_opt": 7.0, "search_term": "catalase", "example_pdb": "1DGF"},
    "carbonic_anhydrase": {"ec": "4.2.1.1", "ph_opt": 7.5, "search_term": "carbonic anhydrase", "example_pdb": "1CA2"},
    
    # Alkaline pH enzymes (pH 8-10)
    "trypsin": {"ec": "3.4.21.4", "ph_opt": 8.0, "search_term": "trypsin", "example_pdb": "1S0Q"},
    "alkaline_phosphatase": {"ec": "3.1.3.1", "ph_opt": 9.0, "search_term": "alkaline phosphatase", "example_pdb": "1ALK"},
    
    # High alkaline enzymes (pH 10+)
    "subtilisin": {"ec": "3.4.21.62", "ph_opt": 10.5, "search_term": "subtilisin", "example_pdb": "1CSE"},
}

def search_pdb_by_text(search_term):
    """Search PDB using full text search"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": search_term
            }
        },
        "request_options": {
            "return_all_hits": False,
            "results_content_type": ["compact"],
            "rows": 5
        },
        "return_type": "entry"
    }
    
    try:
        response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
        if response.status_code == 200:
            data = response.json()
            return [result["identifier"] for result in data.get("result_set", [])]
        else:
            print(f"  Search failed with status {response.status_code}")
    except Exception as e:
        print(f"  Error searching for {search_term}: {e}")
    
    return []

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
        print(f"  Download error: {e}")
        return False

def main():
    # Create output directory
    output_dir = "enzyme_ph_structures"
    os.makedirs(output_dir, exist_ok=True)
    
    # Collect data
    dataset = []
    downloaded_count = 0
    
    print(f"Collecting structures for {len(ENZYME_PH_DATA)} enzymes...\n")
    
    for enzyme_name, info in ENZYME_PH_DATA.items():
        print(f"Processing {enzyme_name} (pH {info['ph_opt']})...")
        
        # Try the example PDB first
        pdb_id = info.get("example_pdb", None)
        
        if pdb_id:
            # Use the known example
            print(f"  Using known structure: {pdb_id}")
        else:
            # Search for structures
            pdb_ids = search_pdb_by_text(info["search_term"])
            if pdb_ids:
                pdb_id = pdb_ids[0]
                print(f"  Found structure: {pdb_id}")
            else:
                print(f"  No structures found")
                continue
        
        # Add to dataset
        dataset.append({
            'enzyme_name': enzyme_name,
            'ec_number': info["ec"],
            'pdb_id': pdb_id,
            'ph_optimum': info["ph_opt"]
        })
        
        # Download structure
        if download_pdb(pdb_id, output_dir):
            downloaded_count += 1
            print(f"  Successfully downloaded {pdb_id}")
        else:
            print(f"  Failed to download {pdb_id}")
        
        time.sleep(1)  # Be nice to PDB server
    
    # Save dataset
    if dataset:
        df = pd.DataFrame(dataset)
        df.to_csv("enzyme_ph_dataset.csv", index=False)
        
        print(f"\n{'='*50}")
        print(f"Dataset complete!")
        print(f"Total enzymes processed: {len(ENZYME_PH_DATA)}")
        print(f"Structures found: {len(dataset)}")
        print(f"Structures downloaded: {downloaded_count}")
        print(f"pH range: {df['ph_optimum'].min()} - {df['ph_optimum'].max()}")
        print(f"Files saved in: {output_dir}/")
        print(f"Dataset saved as: enzyme_ph_dataset.csv")
        print(f"{'='*50}")
    else:
        print("\nNo structures found. Dataset empty.")

if __name__ == "__main__":
    main()