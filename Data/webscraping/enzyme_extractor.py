#!/usr/bin/env python3
"""
Fixed script to find enzyme structures from PDB with proper search queries
"""

import requests
import pandas as pd
import time
import os
from Bio.PDB import PDBList
import json

def search_pdb_by_ec_class(ec_class):
    """Search PDB for structures with specific EC class using proper query format"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Search specifically for EC number
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
        "request_options": {
            "return_all_hits": False,
            "results_content_type": ["compact"],
            "rows": 50
        },
        "return_type": "entry"
    }
    
    try:
        response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
        if response.status_code == 200:
            data = response.json()
            return [result["identifier"] for result in data.get("result_set", [])]
    except Exception as e:
        print(f"Error searching for EC class {ec_class}: {e}")
    
    return []

def search_pdb_alternative(ec_class):
    """Alternative search method using full-text search"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": f"EC:{ec_class}"
            }
        },
        "request_options": {
            "return_all_hits": False,
            "results_content_type": ["compact"],
            "rows": 50
        },
        "return_type": "entry"
    }
    
    try:
        response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
        if response.status_code == 200:
            data = response.json()
            return [result["identifier"] for result in data.get("result_set", [])]
    except:
        return []

def get_enzyme_info_simple(pdb_id):
    """Get basic enzyme information including EC number"""
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            
            # Get EC number
            ec_number = None
            if "rcsb_polymer_entity" in data:
                ec_lineage = data["rcsb_polymer_entity"].get("rcsb_ec_lineage", [])
                if ec_lineage:
                    ec_number = ec_lineage[0].get("id")
            
            # Get basic info
            info = {
                'pdb_id': pdb_id,
                'ec_number': ec_number,
                'name': data.get("rcsb_polymer_entity", {}).get("pdbx_description", ""),
                'organism': data.get("rcsb_entity_source_organism", [{}])[0].get("scientific_name", "") if data.get("rcsb_entity_source_organism") else ""
            }
            
            return info
    except:
        pass
    
    return None

def search_enzymes_by_keyword():
    """Search for enzymes using keyword search"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Search for common enzyme keywords
    enzyme_keywords = [
        "dehydrogenase",
        "kinase",
        "phosphatase",
        "protease",
        "lipase",
        "amylase",
        "oxidase",
        "reductase",
        "transferase",
        "hydrolase",
        "lyase",
        "isomerase",
        "ligase",
        "synthase",
        "peptidase"
    ]
    
    all_pdbs = set()
    
    for keyword in enzyme_keywords:
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "struct.title",
                            "operator": "contains_words",
                            "value": keyword
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "entity_poly.rcsb_entity_polymer_type",
                            "operator": "exact_match",
                            "value": "Protein"
                        }
                    }
                ]
            },
            "request_options": {
                "return_all_hits": False,
                "results_content_type": ["compact"],
                "rows": 50
            },
            "return_type": "entry"
        }
        
        try:
            response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
            if response.status_code == 200:
                data = response.json()
                pdbs = [result["identifier"] for result in data.get("result_set", [])]
                all_pdbs.update(pdbs)
                print(f"Found {len(pdbs)} structures for '{keyword}'")
        except:
            continue
        
        time.sleep(0.5)
    
    return list(all_pdbs)

def get_detailed_info(pdb_id):
    """Get detailed information about a PDB entry"""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            
            info = {
                'pdb_id': pdb_id,
                'title': data.get("struct", {}).get("title", ""),
                'resolution': None,
                'method': data.get("exptl", [{}])[0].get("method", "") if data.get("exptl") else "",
                'ec_number': None,
                'organism': ""
            }
            
            # Get resolution
            if "rcsb_entry_info" in data:
                res_list = data["rcsb_entry_info"].get("resolution_combined", [])
                if res_list:
                    info['resolution'] = res_list[0]
            
            # Get polymer entity info
            polymer_info = get_enzyme_info_simple(pdb_id)
            if polymer_info:
                info.update(polymer_info)
            
            return info
    except:
        pass
    
    return None

def main():
    """Main function to find enzyme structures"""
    
    print("Finding Enzyme Structures from PDB")
    print("=" * 50)
    
    # Load existing dataset to avoid duplicates
    existing_pdbs = set()
    if os.path.exists("enzyme_ph_dataset.csv"):
        df = pd.read_csv("enzyme_ph_dataset.csv")
        existing_pdbs = set(df['pdb_id'].values)
        print(f"Loaded {len(existing_pdbs)} existing PDB IDs")
    
    # Search by enzyme keywords
    print("\nSearching for enzyme structures by keywords...")
    enzyme_pdbs = search_enzymes_by_keyword()
    print(f"Found {len(enzyme_pdbs)} total enzyme structures")
    
    # Filter out existing ones
    new_pdbs = [pdb for pdb in enzyme_pdbs if pdb not in existing_pdbs]
    print(f"{len(new_pdbs)} are new structures")
    
    # Get detailed information
    enzyme_data = []
    print("\nGetting detailed information...")
    
    for i, pdb_id in enumerate(new_pdbs[:200]):  # Limit to 200 for speed
        if i % 20 == 0:
            print(f"Processing {i}/{min(len(new_pdbs), 200)}...")
        
        info = get_detailed_info(pdb_id)
        if info and info.get('ec_number'):
            # Categorize based on EC number
            ec_parts = info['ec_number'].split('.')
            if ec_parts:
                major_class = ec_parts[0]
                category_map = {
                    "1": "oxidoreductases",
                    "2": "transferases",
                    "3": "hydrolases",
                    "4": "lyases",
                    "5": "isomerases",
                    "6": "ligases"
                }
                info['category'] = category_map.get(major_class, "unclassified")
            else:
                info['category'] = "unclassified"
            
            # Clean enzyme name
            title_words = info['title'].split()
            enzyme_name = title_words[0].lower() if title_words else "enzyme"
            info['enzyme_name'] = enzyme_name.replace(",", "").replace(".", "")
            
            # Add placeholder pH
            info['ph_optimum'] = 7.0
            
            enzyme_data.append(info)
        
        time.sleep(0.1)
    
    if enzyme_data:
        # Create DataFrame
        new_df = pd.DataFrame(enzyme_data)
        
        # Save dataset
        columns = ['category', 'enzyme_name', 'ec_number', 'pdb_id', 'ph_optimum', 'title', 'organism', 'resolution', 'method']
        save_df = new_df[[col for col in columns if col in new_df.columns]]
        save_df.to_csv("expanded_enzyme_dataset.csv", index=False)
        
        # Create organized directory structure
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
        
        # Create summary
        summary = {
            'total_new_structures': len(save_df),
            'downloaded': downloaded,
            'unique_ec_numbers': save_df['ec_number'].nunique(),
            'categories': save_df['category'].value_counts().to_dict(),
            'organisms': save_df['organism'].value_counts().head(10).to_dict() if 'organism' in save_df else {}
        }
        
        with open("expansion_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"\n{'='*50}")
        print(f"Expansion complete!")
        print(f"New structures with EC numbers: {len(save_df)}")
        print(f"Downloaded: {downloaded}")
        print(f"Categories: {', '.join(save_df['category'].unique())}")
        print(f"\nFiles created:")
        print(f"  - expanded_enzyme_dataset.csv")
        print(f"  - expansion_summary.json")
        print(f"  - {base_dir}/ (organized PDB files)")
        print(f"\nNext: Run BRENDA matcher on expanded_enzyme_dataset.csv")
        print(f"{'='*50}")
    else:
        print("No new enzyme structures with EC numbers found")

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

if __name__ == "__main__":
    main()