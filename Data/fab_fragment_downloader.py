#!/usr/bin/env python3
"""
Script to search for and download Fab fragment structures from the RCSB PDB database,
with a focus on approved therapeutic antibodies.
"""

import os
import json
import time
import requests
import csv
from tqdm import tqdm  # For progress bar (install with pip install tqdm)

# List of FDA-approved therapeutic antibodies (abbreviated list)
# This helps identify Fab fragments from marketed drugs
APPROVED_ANTIBODIES = [
    "adalimumab", "bevacizumab", "certolizumab", "denosumab", "eculizumab",
    "golimumab", "infliximab", "natalizumab", "omalizumab", "palivizumab",
    "ranibizumab", "rituximab", "tocilizumab", "trastuzumab", "ustekinumab",
    "vedolizumab", "alemtuzumab", "belimumab", "daratumumab", "durvalumab",
    "idarucizumab", "ipilimumab", "nivolumab", "obinutuzumab", "ocrelizumab",
    "ofatumumab", "pembrolizumab", "ramucirumab", "secukinumab", "siltuximab"
]

def search_fab_fragments(approved_only=False):
    """
    Search for Fab fragment structures in the PDB using the RCSB REST API.
    Handles pagination correctly to get all results, not just the first 10.
    
    Args:
        approved_only (bool): If True, only search for Fab fragments from approved antibodies
    
    Returns:
        list: List of PDB IDs for Fab fragments
    """
    print("Searching for Fab fragment structures in PDB...")
    
    # Base query for Fab fragments
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "struct.title",
                                "operator": "contains_phrase",
                                "value": "Fab"
                            }
                        },
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "struct_keywords.pdbx_keywords",
                                "operator": "contains_phrase",
                                "value": "Fab"
                            }
                        },
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "entity_poly.rcsb_entity_polymer_type",
                                "operator": "exact_match",
                                "value": "Antibody"
                            }
                        }
                    ]
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater",
                        "value": 1
                    }
                }
            ]
        },
        "return_type": "entry"
    }
    
    # If approved_only is True, add additional filter for approved antibodies
    if approved_only:
        approved_queries = []
        for antibody in APPROVED_ANTIBODIES:
            approved_queries.append({
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "struct.title",
                    "operator": "contains_phrase",
                    "value": antibody
                }
            })
        
        # Add the approved antibodies filter to the query
        query["query"]["nodes"].append({
            "type": "group",
            "logical_operator": "or",
            "nodes": approved_queries
        })
    
    # Send POST request to RCSB Search API
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # First request to get total count
    response = requests.post(url, json=query)
    
    if response.status_code != 200:
        raise Exception(f"Search failed with status code {response.status_code}: {response.text}")
    
    results = response.json()
    total_count = results.get("total_count", 0)
    print(f"Found a total of {total_count} Fab fragment structures.")
    
    # Initialize with first batch of results
    pdb_ids = [hit["identifier"] for hit in results.get("result_set", [])]
    
    # Use a different approach - query all results
    if total_count > 10:
        # Perform a new query requesting all results at once
        all_query = query.copy()
        # The trick: Add "return_all_hits" parameter instead of using pager
        all_query["request_options"] = {"return_all_hits": True}
        
        print("Retrieving all structures in a single request...")
        all_response = requests.post(url, json=all_query)
        
        if all_response.status_code != 200:
            print(f"Warning: Could not retrieve all structures at once. Error: {all_response.status_code}")
            print("Falling back to sequential retrieval method...")
            
            # Alternative approach: sequential fetching with pagination
            current_count = len(pdb_ids)
            while current_count < total_count:
                # Fetch next batch with URL parameters instead of JSON body pagination
                next_url = f"{url}?start={current_count}"
                next_response = requests.post(next_url, json=query)
                
                if next_response.status_code != 200:
                    print(f"Warning: Could not retrieve more results. Error: {next_response.status_code}")
                    break
                    
                next_results = next_response.json()
                next_ids = [hit["identifier"] for hit in next_results.get("result_set", [])]
                
                if not next_ids:
                    break  # No more results to fetch
                    
                pdb_ids.extend(next_ids)
                current_count = len(pdb_ids)
                print(f"Retrieved {current_count}/{total_count} structures...")
        else:
            # If the "return_all_hits" approach succeeded
            all_results = all_response.json()
            pdb_ids = [hit["identifier"] for hit in all_results.get("result_set", [])]
    
    print(f"Successfully retrieved {len(pdb_ids)} Fab fragment structures.")
    return pdb_ids

def get_pdb_metadata(pdb_id):
    """
    Get metadata for a PDB entry.
    
    Args:
        pdb_id (str): PDB ID
        
    Returns:
        dict: Metadata for the PDB entry
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Failed to get metadata for {pdb_id}: HTTP {response.status_code}")
        return {}
    
    return response.json()

def download_pdb_files(pdb_ids, output_dir="fab_pdbs"):
    """
    Download PDB files for the given PDB IDs.
    
    Args:
        pdb_ids (list): List of PDB IDs to download
        output_dir (str): Directory to save downloaded files
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # Base URL for downloading PDB files
    base_url = "https://files.rcsb.org/download/"
    
    print(f"Downloading {len(pdb_ids)} PDB files...")
    
    # Download each PDB file
    for pdb_id in tqdm(pdb_ids):
        # Construct the file URL and output path
        file_url = f"{base_url}{pdb_id}.pdb"
        output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
        
        # Skip if file already exists
        if os.path.exists(output_path):
            continue
        
        # Download the file
        try:
            response = requests.get(file_url)
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
            else:
                print(f"Failed to download {pdb_id}: HTTP {response.status_code}")
            
            # Be nice to the server with a small delay
            time.sleep(0.1)
        except Exception as e:
            print(f"Error downloading {pdb_id}: {str(e)}")
    
    print(f"Downloaded PDB files are saved in the '{output_dir}' directory.")

def create_summary_file(pdb_ids, output_dir="fab_pdbs"):
    """
    Create a summary file with detailed information about each PDB entry.
    
    Args:
        pdb_ids (list): List of PDB IDs
        output_dir (str): Directory to save the summary file
    """
    summary_path = os.path.join(output_dir, "fab_summary.csv")
    
    print("Creating detailed summary file...")
    
    with open(summary_path, 'w', newline='') as csvfile:
        fieldnames = ['PDB_ID', 'Title', 'Resolution', 'Release_Date', 'Deposition_Date', 
                      'Experimental_Method', 'Structure_Author', 'Keywords']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for pdb_id in tqdm(pdb_ids):
            try:
                metadata = get_pdb_metadata(pdb_id)
                if not metadata:
                    continue
                
                entry = metadata.get('rcsb_entry_info', {})
                struct = metadata.get('struct', {})
                
                # Write the data to the CSV file
                writer.writerow({
                    'PDB_ID': pdb_id,
                    'Title': struct.get('title', ''),
                    'Resolution': entry.get('resolution_combined', [0])[0] if entry.get('resolution_combined') else '',
                    'Release_Date': entry.get('release_date', ''),
                    'Deposition_Date': entry.get('deposit_date', ''),
                    'Experimental_Method': entry.get('experimental_method', ''),
                    'Structure_Author': struct.get('pdbx_descriptor', ''),
                    'Keywords': struct.get('pdbx_keywords', '')
                })
                
                # Small delay to be nice to the server
                time.sleep(0.1)
                
            except Exception as e:
                print(f"Error processing metadata for {pdb_id}: {str(e)}")
    
    print(f"Summary file created: {summary_path}")

def identify_therapeutic_relevance(pdb_ids, output_dir="fab_pdbs"):
    """
    Create a file identifying which Fab fragments are related to therapeutic antibodies.
    
    Args:
        pdb_ids (list): List of PDB IDs
        output_dir (str): Directory to save the output file
    """
    therapeutic_path = os.path.join(output_dir, "therapeutic_fabs.txt")
    
    print("Identifying therapeutically relevant Fab fragments...")
    
    therapeutic_pdb_ids = []
    
    for pdb_id in tqdm(pdb_ids):
        try:
            metadata = get_pdb_metadata(pdb_id)
            if not metadata:
                continue
            
            struct = metadata.get('struct', {})
            title = struct.get('title', '').lower()
            
            # Check if the PDB entry is related to any approved antibodies
            is_therapeutic = any(antibody.lower() in title for antibody in APPROVED_ANTIBODIES)
            
            if is_therapeutic:
                therapeutic_pdb_ids.append((pdb_id, title))
                
            # Small delay to be nice to the server
            time.sleep(0.1)
                
        except Exception as e:
            print(f"Error checking therapeutic relevance for {pdb_id}: {str(e)}")
    
    # Write the results to a file
    with open(therapeutic_path, 'w') as f:
        f.write(f"Found {len(therapeutic_pdb_ids)} therapeutically relevant Fab fragments\n\n")
        for pdb_id, title in therapeutic_pdb_ids:
            f.write(f"{pdb_id}: {title}\n")
    
    print(f"Therapeutic relevance file created: {therapeutic_path}")
    return therapeutic_pdb_ids

def main():
    """Main function to search and download Fab fragment structures."""
    print("Fab Fragment PDB Structure Downloader")
    print("=====================================")
    
    try:
        # Ask user if they want only approved therapeutic antibodies
        print("\nOptions:")
        print("1. Download all Fab fragment structures")
        print("2. Download only Fab fragments from approved therapeutic antibodies")
        
        choice = input("\nEnter your choice (1 or 2): ").strip()
        approved_only = (choice == "2")
        
        # Search for Fab fragment structures
        pdb_ids = search_fab_fragments(approved_only)
        
        if not pdb_ids:
            print("No Fab fragment structures found. Check your search criteria.")
            return
        
        # Download the PDB files
        output_dir = "fab_therapeutic_pdbs" if approved_only else "fab_pdbs"
        download_pdb_files(pdb_ids, output_dir)
        
        # Create a summary file with detailed information
        create_summary_file(pdb_ids, output_dir)
        
        # If downloading all Fab fragments, identify those with therapeutic relevance
        if not approved_only:
            therapeutic_pdb_ids = identify_therapeutic_relevance(pdb_ids, output_dir)
            print(f"Found {len(therapeutic_pdb_ids)} therapeutically relevant Fab fragments.")
        
        print("\nProcess completed successfully!")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()