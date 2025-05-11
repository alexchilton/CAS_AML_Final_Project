#!/usr/bin/env python3
"""
Working script to find enzyme structures from PDB
"""

import requests
import pandas as pd
import time
import os
from Bio.PDB import PDBList
import json

def search_pdb_enzymes():
    """Search PDB for enzyme structures using working queries"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # List of EC classes to search for
    ec_classes = [
        # Oxidoreductases
        "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "1.10", "1.11", "1.14",
        # Transferases  
        "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8",
        # Hydrolases
        "3.1", "3.2", "3.3", "3.4", "3.5", "3.6",
        # Lyases
        "4.1", "4.2", "4.3",
        # Isomerases
        "5.1", "5.2", "5.3", "5.4",
        # Ligases
        "6.1", "6.2", "6.3", "6.4"
    ]
    
    all_pdbs = []
    
    for ec_class in ec_classes:
        # Search by EC lineage
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity.rcsb_ec_lineage.id",
                    "operator": "exact_match",
                    "value": ec_class
                }
            },
            "return_type": "entry"
        }
        
        try:
            response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
            if response.status_code == 200:
                data = response.json()
                results = data.get("result_set", [])
                pdb_ids = [r["identifier"] for r in results]
                all_pdbs.extend(pdb_ids)
                print(f"EC {ec_class}: Found {len(pdb_ids)} structures")
            time.sleep(0.5)
        except:
            continue
    
    # Also search by enzyme names
    enzyme_keywords = [
        "dehydrogenase", "kinase", "phosphatase", "protease", "lipase",
        "amylase", "oxidase", "reductase", "transferase", "hydrolase",
        "lyase", "isomerase", "ligase", "synthase", "peptidase", "catalase"
    ]
    
    for keyword in enzyme_keywords:
        query = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": keyword
                }
            },
            "return_type": "entry"
        }
        
        try:
            response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
            if response.status_code == 200:
                data = response.json()
                results = data.get("result_set", [])[:50]  # Limit to 50 per keyword
                pdb_ids = [r["identifier"] for r in results]
                all_pdbs.extend(pdb_ids)
                print(f"{keyword}: Found {len(pdb_ids)} structures")
            time.sleep(0.5)
        except:
            continue
    
    # Remove duplicates
    return list(set(all_pdbs))

def get_enzyme_details(pdb_id):
    """Get enzyme details including full EC number"""
    try:
        # Get polymer entity data
        polymer_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
        polymer_response = requests.get(polymer_url)
        
        if polymer_response.status_code == 200:
            polymer_data = polymer_response.json()
            
            # Get full EC number
            ec_lineage = polymer_data.get("rcsb_polymer_entity", {}).get("rcsb_ec_lineage", [])
            ec_numbers = []
            for ec_info in ec_lineage:
                ec_id = ec_info.get("id", "")
                if ec_id and "." in ec_id:  # Full EC number
                    ec_numbers.append(ec_id)
            
            # Get entry data
            entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            entry_response = requests.get(entry_url)
            
            if entry_response.status_code == 200:
                entry_data = entry_response.json()
                
                info = {
                    'pdb_id': pdb_id,
                    'title': entry_data.get("struct", {}).get("title", ""),
                    'ec_number': ec_numbers[0] if ec_numbers else None,
                    'all_ec_numbers': ec_numbers,
                    'organism': "",
                    'resolution': None,
                    'method': entry_data.get("exptl", [{}])[0].get("method", "") if entry_data.get("exptl") else ""
                }
                
                # Get organism
                source_organisms = entry_data.get("rcsb_entity_source_organism", [])
                if source_organisms:
                    info['organism'] = source_organisms[0].get("scientific_name", "")
                
                # Get resolution
                if "rcsb_entry_info" in entry_data:
                    res_list = entry_data["rcsb_entry_info"].get("resolution_combined", [])
                    if res_list:
                        info['resolution'] = res_list[0]
                
                return info
    except:
        pass
    
    return None

def categorize_by_ec(ec_number):
    """Categorize enzyme by EC number"""
    if not ec_number:
        return "unclassified"
    
    ec_parts = ec_number.split(".")
    if not ec_parts:
        return "unclassified"
    
    major_class = ec_parts[0]
    categories = {
        "1": "oxidoreductases",
        "2": "transferases", 
        "3": "hydrolases",
        "4": "lyases",
        "5": "isomerases",
        "6": "ligases"
    }
    
    return categories.get(major_class, "unclassified")

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
    except:
        return False

def main():
    """Main function"""
    print("Finding Enzyme Structures from PDB")
    print("=" * 50)
    
    # Load existing dataset
    existing_pdbs = set()
    if os.path.exists("enzyme_ph_dataset.csv"):
        df = pd.read_csv("enzyme_ph_dataset.csv")
        existing_pdbs = set(df['pdb_id'].values)
        print(f"Loaded {len(existing_pdbs)} existing PDB IDs")
    
    # Search for enzymes
    print("\nSearching for enzyme structures...")
    all_pdbs = search_pdb_enzymes()
    print(f"\nFound {len(all_pdbs)} total structures")
    
    # Filter new ones
    new_pdbs = [pdb for pdb in all_pdbs if pdb not in existing_pdbs]
    print(f"{len(new_pdbs)} are new structures")
    
    # Get details for new structures
    enzyme_data = []
    print("\nGetting enzyme details...")
    
    for i, pdb_id in enumerate(new_pdbs[:300]):  # Limit to 300
        if i % 20 == 0:
            print(f"Processing {i}/{min(len(new_pdbs), 300)}...")
        
        details = get_enzyme_details(pdb_id)
        if details and details['ec_number']:
            details['category'] = categorize_by_ec(details['ec_number'])
            
            # Extract enzyme name from title
            title_words = details['title'].split()
            enzyme_name = "enzyme"
            for word in title_words:
                if any(suffix in word.lower() for suffix in ["ase", "in", "gen"]):
                    enzyme_name = word.lower().strip(",.;:")
                    break
            details['enzyme_name'] = enzyme_name
            
            # Add placeholder pH
            details['ph_optimum'] = 7.0
            
            enzyme_data.append(details)
        
        time.sleep(0.1)
    
    if enzyme_data:
        # Create DataFrame
        df = pd.DataFrame(enzyme_data)
        
        # Save dataset
        columns = ['category', 'enzyme_name', 'ec_number', 'pdb_id', 'ph_optimum', 'title', 'organism', 'resolution', 'method']
        save_df = df[[col for col in columns if col in df.columns]]
        save_df.to_csv("expanded_enzyme_dataset.csv", index=False)
        
        # Create directory structure
        base_dir = "expanded_enzyme_structures"
        os.makedirs(base_dir, exist_ok=True)
        
        # Download PDB files
        print("\nDownloading PDB files...")
        downloaded = 0
        
        for idx, row in save_df.iterrows():
            category_dir = os.path.join(base_dir, row['category'])
            os.makedirs(category_dir, exist_ok=True)
            
            if download_pdb(row['pdb_id'], category_dir):
                downloaded += 1
                if downloaded % 10 == 0:
                    print(f"Downloaded {downloaded}/{len(save_df)}")
            
            time.sleep(0.5)
        
        # Summary
        print(f"\n{'='*50}")
        print(f"Expansion complete!")
        print(f"New structures with EC numbers: {len(save_df)}")
        print(f"Downloaded: {downloaded}")
        print(f"Categories: {', '.join(save_df['category'].unique())}")
        print(f"\nFiles created:")
        print(f"  - expanded_enzyme_dataset.csv")
        print(f"  - {base_dir}/ (organized PDB files)")
        print(f"\nNext: Run BRENDA matcher on expanded_enzyme_dataset.csv")
        print(f"{'='*50}")
    else:
        print("No new enzyme structures with EC numbers found")

if __name__ == "__main__":
    main()