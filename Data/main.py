from processing.processor import ProteinProcessor

def main():
    processor = ProteinProcessor(
        folder_path="new_protein_data",
        output_path="new_protein_data_networkx_new",
        max_file_size_mb=25
    )
    processor.run()

if __name__ == "__main__":
    main()
