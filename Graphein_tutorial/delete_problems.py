import os
import shutil

def delete_problematic_pdbs(problematic_file_path, pdb_dir, backup_dir=None):
    """
    Delete problematic PDB files listed in the text file.

    Parameters:
    problematic_file_path (str): Path to the text file containing problematic file names
    pdb_dir (str): Directory containing the PDB files
    backup_dir (str, optional): If provided, files will be moved here instead of deleted
    """
    # Check if the problematic files list exists
    if not os.path.exists(problematic_file_path):
        print(f"Error: File {problematic_file_path} does not exist.")
        return

    # Create backup directory if specified and doesn't exist
    if backup_dir and not os.path.exists(backup_dir):
        os.makedirs(backup_dir)
        print(f"Created backup directory: {backup_dir}")

    # Read the problematic files
    problematic_files = []
    with open(problematic_file_path, 'r') as f:
        # Skip the header lines (first two lines)
        lines = f.readlines()
        if len(lines) > 2:
            for line in lines[2:]:  # Skip first two lines
                if line.strip():  # Skip empty lines
                    # Extract just the filename part (before the colon)
                    parts = line.split(':')
                    if parts:
                        filename = parts[0].strip()
                        problematic_files.append(filename)

    if not problematic_files:
        print("No problematic files found to delete.")
        return

    print(f"Found {len(problematic_files)} problematic files to process.")

    # Process each problematic file
    deleted_count = 0
    for filename in problematic_files:
        file_path = os.path.join(pdb_dir, filename)

        if os.path.exists(file_path):
            try:
                if backup_dir:
                    # Move to backup instead of deleting
                    backup_path = os.path.join(backup_dir, filename)
                    shutil.move(file_path, backup_path)
                    print(f"Moved: {filename} to {backup_dir}")
                else:
                    # Delete the file
                    os.remove(file_path)
                    print(f"Deleted: {filename}")

                deleted_count += 1
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")
        else:
            print(f"Warning: File not found - {filename}")

    print(f"\nSummary: {deleted_count} of {len(problematic_files)} files were processed.")

    # Create a log file of deleted files
    log_file = "deleted_pdb_files.txt"
    with open(log_file, 'w') as f:
        f.write("Processed problematic PDB files:\n")
        f.write("="*50 + "\n")
        for filename in problematic_files:
            f.write(f"{filename}\n")

    print(f"Log created: {log_file}")

if __name__ == "__main__":
    # Default paths - update these as needed
    problematic_file_path = "problematic_pdb_files.txt"
    pdb_dir = os.path.expanduser("~/Downloads/nanobody_extracted")  # Update this path

    # Option 1: Delete files permanently
    # delete_problematic_pdbs(problematic_file_path, pdb_dir)

    # Option 2: Move files to backup directory instead of deleting
    backup_dir = "backup_problematic_pdbs"
    delete_problematic_pdbs(problematic_file_path, pdb_dir, backup_dir)