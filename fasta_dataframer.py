import pandas as pd  # Import pandas for data manipulation
import argparse  # Import argparse for command-line argument parsing
from Bio import SeqIO  # Import SeqIO from Biopython to read FASTA files
import os  # Import os to handle file paths

# Function to load the metadata file (CSV)
def load_metadata(metadata_file):
    print(f"Loading metadata from file: {metadata_file}")
    return pd.read_csv(metadata_file, low_memory=False)  # Load CSV metadata file into pandas DataFrame

# Function to process FASTA file and merge with metadata
def process_fasta_and_metadata(fasta_file, metadata_file, keyword, output_file):
    seq_list = []  # List to hold extracted SeqID, genome_id, name, and sequence

    # Process the FASTA file and extract required information
    print(f"Processing FASTA file: {fasta_file}")
    for record in SeqIO.parse(fasta_file, "fasta"):  # Read each entry from the FASTA file
        seq_id = record.id  # Extract SeqID (the ID of the sequence)
        
        # Split the description to get genome_id and protein name
        parts = record.description.split(" ")
        genome_id = parts[0].split("|")[1]  # Extract the genome_id from the description
        
        # Clean the genome_id by removing extra parts like .mat_peptide (we want just the numbers)
        genome_id_clean = ".".join(genome_id.split(".")[:2])  # Keep only the part before the first dot
        
        # Extract the protein name (description) and sequence
        name = " ".join(parts[1:]).split('[')[0].strip()  # Get everything after the first space
        sequence = str(record.seq)  # Convert sequence to string format
        
        # Filter by keyword (e.g., 'nsp1') to only keep relevant sequences
        if keyword.lower() in name.lower():
            seq_list.append((genome_id_clean, name, sequence))  # Add to list if the keyword matches

    print(f"Total number of filtered sequences (based on keyword '{keyword}'): {len(seq_list)}")

    # Create a DataFrame from the list of sequences
    seq_df = pd.DataFrame(seq_list, columns=["Genome_ID", "Protein_Name", "Sequence"])

    print("DataFrame from FASTA entries (filtered by keyword):")
    print(seq_df.head())  # Print the first few rows of the DataFrame

    # Load the metadata CSV file into a DataFrame
    metadata_df = load_metadata(metadata_file)

    # Clean metadata genome IDs (ensure they are strings and remove extra spaces)
    metadata_df["Genome ID"] = metadata_df["Genome ID"].astype(str).str.strip()

    # Debugging: Print the unique genome IDs in both the FASTA and metadata DataFrames
    print("Unique Genome IDs in FASTA DataFrame:")
    print(seq_df["Genome_ID"].unique())
    print("Unique Genome IDs in metadata DataFrame:")
    print(metadata_df["Genome ID"].unique())

    # Merge the two DataFrames on the cleaned genome ID
    merged_df = pd.merge(seq_df, metadata_df, left_on="Genome_ID", right_on="Genome ID", how="left")

    # Print the merged DataFrame and check for any NaN values after merging
    print("Merged DataFrame (first 5 rows):")
    print(merged_df.head())
    print("NaN values in merged DataFrame columns:")
    print(merged_df.isna().sum())

    # Save the final merged DataFrame to a CSV file
    merged_df.to_csv(output_file, index=False)
    print(f"Process complete! The merged output has been saved to {output_file}")

# Main function to handle command-line arguments
def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process FASTA and metadata files, and merge them into a CSV file.")
    parser.add_argument("-f", "--fasta_file", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-m", "--metadata_file", required=True, help="Path to the input metadata CSV file.")
    parser.add_argument("-k", "--keyword", required=True, help="Keyword to filter FASTA entries.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output CSV file.")
    args = parser.parse_args()
    
    # Call the processing function with the parsed arguments
    process_fasta_and_metadata(args.fasta_file, args.metadata_file, args.keyword, args.output_file)

# Entry point of the script
if __name__ == "__main__":
    main()

# Example command to run the script:
# python fasta_to_metadata.py -f combined_proteins.fasta -m BVBRC_genome.csv -k nsp1 -o merged_output.csv
