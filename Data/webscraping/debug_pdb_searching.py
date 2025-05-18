#!/usr/bin/env python3
"""
Debug script to test PDB search functionality
"""

import requests
import json

def test_simple_search():
    """Test basic PDB search"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Very simple query for lysozyme
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "struct.title",
                "operator": "contains_words",
                "value": "lysozyme"
            }
        },
        "request_options": {
            "return_all_hits": False,
            "results_content_type": ["compact"],
            "rows": 5
        },
        "return_type": "entry"
    }
    
    print("Testing simple search for 'lysozyme'...")
    print(f"Query: {json.dumps(query, indent=2)}")
    
    try:
        response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
        print(f"Response status: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            print(f"Response: {json.dumps(data, indent=2)[:500]}...")  # First 500 chars
            
            if "result_set" in data:
                results = data["result_set"]
                print(f"\nFound {len(results)} results")
                for result in results[:3]:
                    print(f"  PDB ID: {result.get('identifier', 'N/A')}")
            else:
                print("No result_set in response")
        else:
            print(f"Error response: {response.text[:500]}")
    except Exception as e:
        print(f"Exception: {e}")

def test_ec_search():
    """Test EC number search"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Search for a specific EC number
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity.rcsb_ec_lineage.id",
                "operator": "exact_match",
                "value": "3.2.1.17"  # Lysozyme EC number
            }
        },
        "return_type": "entry"
    }
    
    print("\nTesting EC number search...")
    
    try:
        response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
        print(f"Response status: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            if "result_set" in data:
                print(f"Found {len(data['result_set'])} results")
            else:
                print("No results in response")
                print(f"Response: {json.dumps(data, indent=2)[:500]}...")
    except Exception as e:
        print(f"Exception: {e}")

def test_full_text_search():
    """Test full text search"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": "kinase"
            }
        },
        "return_type": "entry"
    }
    
    print("\nTesting full text search for 'kinase'...")
    
    try:
        response = requests.post(url, json=query, headers={"Content-Type": "application/json"})
        print(f"Response status: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            if "result_set" in data:
                print(f"Found {len(data['result_set'])} results")
                for result in data['result_set'][:3]:
                    print(f"  PDB ID: {result.get('identifier', 'N/A')}")
    except Exception as e:
        print(f"Exception: {e}")

def test_get_entry():
    """Test getting a specific PDB entry"""
    pdb_id = "1DPX"  # Known lysozyme structure
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    print(f"\nTesting direct entry retrieval for {pdb_id}...")
    
    try:
        response = requests.get(url)
        print(f"Response status: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            print(f"Title: {data.get('struct', {}).get('title', 'N/A')}")
            
            # Check for EC number
            polymer_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
            polymer_response = requests.get(polymer_url)
            
            if polymer_response.status_code == 200:
                polymer_data = polymer_response.json()
                ec_lineage = polymer_data.get("rcsb_polymer_entity", {}).get("rcsb_ec_lineage", [])
                if ec_lineage:
                    print(f"EC number: {ec_lineage[0].get('id', 'N/A')}")
                else:
                    print("No EC number found")
    except Exception as e:
        print(f"Exception: {e}")

def main():
    """Run all tests"""
    print("PDB Search Debug Tests")
    print("=" * 50)
    
    test_simple_search()
    test_ec_search()
    test_full_text_search()
    test_get_entry()
    
    print("\nDebug complete!")

if __name__ == "__main__":
    main()