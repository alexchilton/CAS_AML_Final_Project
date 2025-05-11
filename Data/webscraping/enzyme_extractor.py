#!/usr/bin/env python3
"""
Simplified script to find enzyme structures with complete EC numbers
"""

import requests
import pandas as pd
import time
import os
from Bio.PDB import PDBList
import json

def search_pdb_simple(max_per_type=50):
    """Simple search for enzyme structures with pagination"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    enzyme_names = [
        "protease", "kinase", "phosphatase", "dehydrogenase", 
        "oxidase", "lipase", "amylase", "synthase", "transferase",
        "hydrolase", "lyase", "isomerase", "ligase"
    ]
    
    all_pdbs = []
    
    for enzyme in enzyme_names:
        pdb_ids_for_enzyme = []
        start = 0
        
        while len(pdb_ids_for_enzyme) < max_per_type:
            query = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": enzyme
                    }
                },
                "request_options": {
                    "paginate": {
                        "start": start,
                        "rows": 10  # API seems to limit to 10 per request
                    }
                },
                "return_type": "entry"
            }
            
            try:
                response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
                if response.status_code == 200:
                    data = response.json()
                    results = data.get("result_set", [])
                    
                    if not results:  # No more results
                        break
                    
                    for result in results:
                        if len(pdb_ids_for_enzyme) < max_per_type:
                            pdb_ids_for_enzyme.append(result["identifier"])
                    
                    start += 10  # Move to next page
                else:
                    break
            except Exception as e:
                print(f"Error searching for {enzyme}: {e}")
                break
            
            time.sleep(0.2)  # Be nice to the API
        
        all_pdbs.extend(pdb_ids_for_enzyme)
        print(f"Found {len(pdb_ids_for_enzyme)} structures for '{enzyme}'")
        time.sleep(0.5)
    
    return list(set(all_pdbs))  # Remove duplicates

# def search_pdb_simple():
#     """Simple search for enzyme structures"""
#     url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
#     # Search for common enzymes using full text search
#     enzyme_names = [
#         "protease", "kinase", "phosphatase", "dehydrogenase", 
#         "oxidase", "lipase", "amylase", "synthase", "transferase",
#         "hydrolase", "lyase", "isomerase", "ligase"
#     ]
    
#     all_pdbs = []
    
#     for enzyme in enzyme_names:
#         query = {
#             "query": {
#                 "type": "terminal",
#                 "service": "full_text",
#                 "parameters": {
#                     "value": enzyme
#                 }
#             },
#             "return_type": "entry"
#         }
        
#         try:
#             response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
#             if response.status_code == 200:
#                 data = response.json()
#                 results = data.get("result_set", [])[:50]  # Get up to 50 per enzyme type
#                 for result in results:
#                     all_pdbs.append(result["identifier"])
#                 print(f"Found {len(results)} structures for '{enzyme}'")
#         except Exception as e:
#             print(f"Error searching for {enzyme}: {e}")
        
#         time.sleep(0.5)
    
#     return list(set(all_pdbs))  # Remove duplicates

def get_polymer_entity_info(pdb_id):
    """Get polymer entity information including EC number"""
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            
            # Extract EC number
            ec_lineage = data.get("rcsb_polymer_entity", {}).get("rcsb_ec_lineage", [])
            
            # Find the most complete EC number
            best_ec = None
            max_parts = 0
            
            for ec_info in ec_lineage:
                ec_id = ec_info.get("id", "")
                if ec_id and "." in ec_id:
                    parts = ec_id.split(".")
                    valid_parts = [p for p in parts if p and p != "-"]
                    if len(valid_parts) > max_parts:
                        max_parts = len(valid_parts)
                        best_ec = ec_id
            
            return {
                "ec_number": best_ec,
                "description": data.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
            }
    except:
        pass
    
    return None

def get_structure_info(pdb_id):
    """Get structure information"""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            
            info = {
                "pdb_id": pdb_id,
                "title": data.get("struct", {}).get("title", ""),
                "resolution": None,
                "method": "",
                "organism": ""
            }
            
            # Get resolution
            if "rcsb_entry_info" in data:
                res_list = data["rcsb_entry_info"].get("resolution_combined", [])
                if res_list:
                    info["resolution"] = res_list[0]
            
            # Get method
            if "exptl" in data and data["exptl"]:
                info["method"] = data["exptl"][0].get("method", "")
            
            # Get organism
            if "rcsb_entity_source_organism" in data and data["rcsb_entity_source_organism"]:
                info["organism"] = data["rcsb_entity_source_organism"][0].get("scientific_name", "")
            
            return info
    except:
        pass
    
    return None

def categorize_enzyme(ec_number):
    """Categorize enzyme based on EC number"""
    if not ec_number:
        return "unclassified"
    
    ec_parts = ec_number.split(".")
    if not ec_parts:
        return "unclassified"
    
    # Main EC classes
    main_classes = {
        "1": "oxidoreductases",
        "2": "transferases",
        "3": "hydrolases",
        "4": "lyases",
        "5": "isomerases",
        "6": "ligases"
    }
    
    # More specific categories
    if len(ec_parts) >= 2:
        ec_prefix = f"{ec_parts[0]}.{ec_parts[1]}"
        specific_categories = {
            "3.4": "proteases",
            "3.1.3": "phosphatases",
            "3.1.1": "lipases",
            "3.2": "glycosidases",
            "1.1": "dehydrogenases",
            "2.7": "kinases"
        }
        
        if ec_prefix in specific_categories:
            return specific_categories[ec_prefix]
    
    return main_classes.get(ec_parts[0], "unclassified")

# def main():
#     print("Finding Enzyme Structures with Complete EC Numbers")
#     print("=" * 50)
    
#     # Step 1: Search for enzyme structures
#     print("\nSearching for enzyme structures...")
#     pdb_ids = search_pdb_simple()
#     print(f"\nFound {len(pdb_ids)} candidate structures")


def main(max_structures_per_type=50):
    """Main function
    
    Args:
        max_structures_per_type: Maximum number of structures to retrieve per enzyme type
    """
    print("Finding Enzyme Structures with Complete EC Numbers")
    print("=" * 50)
    print(f"Maximum structures per enzyme type: {max_structures_per_type}")
    
    # Step 1: Search for enzyme structures
    print("\nSearching for enzyme structures...")
    pdb_ids = search_pdb_simple(max_per_type=max_structures_per_type)  # <-- This is the key change
    print(f"\nFound {len(pdb_ids)} candidate structures")
    if not pdb_ids:
        print("No structures found. Check your internet connection.")
        return
    
    # Step 2: Get complete information for each structure
    enzyme_data = []
    print("\nGetting enzyme information...")
    
    for i, pdb_id in enumerate(pdb_ids):
        if i % 50 == 0:
            print(f"Processing {i}/{len(pdb_ids)}...")
        
        # Get polymer entity info (EC number)
        polymer_info = get_polymer_entity_info(pdb_id)
        
        if polymer_info and polymer_info["ec_number"]:
            ec_number = polymer_info["ec_number"]
            
            # Check if EC number is reasonably complete (at least 3 parts)
            ec_parts = [p for p in ec_number.split(".") if p and p != "-"]
            if len(ec_parts) >= 3:
                # Get structure info
                struct_info = get_structure_info(pdb_id)
                
                if struct_info:
                    # Combine information
                    enzyme_entry = {
                        "pdb_id": pdb_id,
                        "ec_number": ec_number,
                        "category": categorize_enzyme(ec_number),
                        "enzyme_name": polymer_info["description"].lower().split()[0] if polymer_info["description"] else "enzyme",
                        "ph_optimum": 7.0,  # Placeholder
                        "title": struct_info["title"],
                        "organism": struct_info["organism"],
                        "resolution": struct_info["resolution"],
                        "method": struct_info["method"]
                    }
                    
                    enzyme_data.append(enzyme_entry)
                    
                    # Show progress
                    if len(enzyme_data) <= 5:
                        print(f"  {pdb_id}: EC {ec_number}")
        
        time.sleep(0.1)  # Rate limiting
    
    # Step 3: Save results
    if enzyme_data:
        df = pd.DataFrame(enzyme_data)
        
        print(f"\nFound {len(df)} structures with complete EC numbers")
        print(f"EC number distribution:")
        ec_lengths = df['ec_number'].apply(lambda x: len([p for p in x.split(".") if p and p != "-"]))
        print(ec_lengths.value_counts().sort_index())
        
        # Save to CSV
        df.to_csv("enzyme_dataset_complete_ec.csv", index=False)
        
        # Create output directory
        output_dir = "enzyme_structures_complete_ec"
        os.makedirs(output_dir, exist_ok=True)
        
        # Download PDB files
        print("\nDownloading PDB files...")
        downloaded = 0
        
        for idx, row in df.iterrows():
            category_dir = os.path.join(output_dir, row['category'])
            os.makedirs(category_dir, exist_ok=True)
            
            if download_pdb(row['pdb_id'], category_dir):
                downloaded += 1
                if downloaded % 10 == 0:
                    print(f"Downloaded {downloaded}/{len(df)}")
            
            time.sleep(0.5)
        
        print(f"\nComplete! Found {len(df)} enzymes with complete EC numbers")
        print(f"Files saved to: enzyme_dataset_complete_ec.csv")
        print(f"PDB files saved to: {output_dir}/")
    else:
        print("\nNo structures with complete EC numbers found")

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
    main(max_structures_per_type=20)