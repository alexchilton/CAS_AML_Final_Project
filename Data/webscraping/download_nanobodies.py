#!/usr/bin/env python3
"""
Script to search for and download all nanobody structures from the RCSB PDB database.
Implements proper pagination to retrieve ALL matching structures.
"""

import os
import json
import time
import requests
from tqdm import tqdm  # For progress bar (install with pip install tqdm)

def search_nanobodies():
    """
    Search for nanobody structures in the PDB using the RCSB REST API.
    Uses advanced query to find structures containing nanobodies.
    Implements pagination to get ALL results, not just the first 10.
    """
    print("Searching for nanobody structures in PDB...")
    
    # Define the search query for nanobodies
    query = {
        "query": {
            "type": "group",
            "logical_operator": "or",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                        "operator": "exact_match",
                        "value": "Camelidae"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct_keywords.pdbx_keywords",
                        "operator": "contains_phrase",
                        "value": "nanobody"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": "nanobody"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": "VHH"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": "single domain antibody"
                    }
                }
            ]
        },
        "return_type": "entry"
    }
    
    # Send POST request to RCSB Search API
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # First request to get total count
    response = requests.post(url, json=query)
    
    if response.status_code != 200:
        raise Exception(f"Search failed with status code {response.status_code}: {response.text}")
    
    results = response.json()
    total_count = results.get("total_count", 0)
    print(f"Found a total of {total_count} nanobody structures.")
    
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
    
    print(f"Successfully retrieved {len(pdb_ids)} nanobody structures.")
    return pdb_ids

def download_pdb_files(pdb_ids, output_dir="nanobody_pdbs"):
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

def create_summary_file(pdb_ids, output_dir="nanobody_pdbs"):
    """
    Create a summary file listing all downloaded PDB IDs.
    """
    summary_path = os.path.join(output_dir, "nanobody_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"Total nanobody structures found: {len(pdb_ids)}\n\n")
        f.write("PDB IDs:\n")
        for pdb_id in sorted(pdb_ids):
            f.write(f"{pdb_id}\n")
    
    print(f"Summary file created: {summary_path}")

def main():
    """Main function to search and download nanobody structures."""
    print("Nanobody PDB Structure Downloader")
    print("=================================")
    
    try:
        # Search for nanobody structures
        pdb_ids = search_nanobodies()
        
        if not pdb_ids:
            print("No nanobody structures found. Check your search criteria.")
            return
        
        # Download the PDB files
        download_pdb_files(pdb_ids)
        
        # Create a summary file
        create_summary_file(pdb_ids)
        
        print("\nProcess completed successfully!")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()