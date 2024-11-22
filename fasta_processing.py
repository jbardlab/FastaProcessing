import os
import re  # Import the regex module to clean up filenames because the fasta headers may contain weird characters
from Bio import SeqIO
import argparse  # Import argparse to handle command-line arguments

# Function to create the output directory where the extracted protein files will be saved
# If no directory is provided, it creates one based on the FASTA file name
def create_output_directory(output_dir, fasta_file):
    # If output directory is not specified, create one based on the FASTA file name
    if output_dir:
        output_directory = output_dir
    else:
        fasta_filename = os.path.splitext(os.path.basename(fasta_file))[0]
        output_directory = f"{fasta_filename}_output"
    
    os.makedirs(output_directory, exist_ok=True)  # Create directory if it doesn't exist, this is the same as mkdir command from command prompt
    return output_directory

# Function to load the FASTA file and return an iterator over parsed sequences
def load_fasta(fasta_file):
    return SeqIO.parse(fasta_file, "fasta")  # Parse the FASTA file format

# Function to search for the proteins using the fasta headers
# It takes the parsed sequences and a list of target proteins (keywords) to search for
def extract_proteins_by_header(fasta_sequences, nsp_proteins):
    extracted_proteins = {protein: [] for protein in nsp_proteins}  # Create a dictionary for each keyword
    
    # Loop through all the sequences in the FASTA file
    for seq_record in fasta_sequences:
        # Convert headers to lowercase to make the search case-insensitive
        header = seq_record.description.lower()
        
        # Loop through the nsp proteins to search for matches in the headers
        for protein in nsp_proteins:
            # If the target protein keyword is found in the header, extract this sequence
            if protein in header:
                extracted_proteins[protein].append(seq_record)  # Add the matching protein to the list for that keyword
    return extracted_proteins  # Return the dictionary of extracted proteins categorized by keyword

# This function removes invalid characters from the sequence ID for file saving
def clean_filename(filename):
    # Removes weird characters from the protein ID to make the filename valid
    return re.sub(r'[^\w\-_\.]', '_', filename)

# Function to save the extracted protein sequences into separate FASTA files and folders
def save_extracted_proteins(extracted_proteins, output_directory, nsp_proteins):
    combined_proteins = []  # List to store all extracted proteins for the combined file

    # Loop through each keyword (nsp1, nsp2, etc.)
    for protein in nsp_proteins:
        # Create a folder for each protein keyword inside the output directory
        protein_folder = os.path.join(output_directory, protein)
        os.makedirs(protein_folder, exist_ok=True)

        # Save each sequence in the respective folder
        for seq_record in extracted_proteins[protein]:
            cleaned_id = clean_filename(seq_record.id)
            output_path = os.path.join(protein_folder, f"{protein}_{cleaned_id}.fasta")
            
            # Save the sequence in FASTA format
            with open(output_path, "w") as output_handle:
                SeqIO.write(seq_record, output_handle, "fasta")
            print(f"Saved {protein.upper()} to {output_path}")
            combined_proteins.append(seq_record)  # Add to the combined list for the combined file

    # Save all proteins into a single file
    combined_file_path = os.path.join(output_directory, "combined_proteins.fasta")
    with open(combined_file_path, "w") as combined_handle:
        SeqIO.write(combined_proteins, combined_handle, "fasta")
    print(f"Saved all proteins to {combined_file_path}")

# Main function to process the FASTA file
def process_fasta_file(fasta_file, keywords, output_dir):
    # Create the output directory
    output_directory = create_output_directory(output_dir, fasta_file)
    
    # Load the FASTA file using SeqIO.parse
    fasta_sequences = load_fasta(fasta_file)
    
    # Extract protein sequences based on their headers
    extracted_proteins = extract_proteins_by_header(fasta_sequences, keywords)
    
    # Save each extracted protein sequence in separate files and combined file
    save_extracted_proteins(extracted_proteins, output_directory, keywords)

# Set up command-line argument parsing using argparse
def main():
    parser = argparse.ArgumentParser(description="Extract protein sequences from a FASTA file based on keywords.")
    
    # Define command-line arguments
    parser.add_argument("-f", "--fasta", required=True, help="Path to the input FASTA file")
    parser.add_argument("-k", "--keywords", required=True, nargs="+", help="List of keywords to search for in headers (e.g., nsp1 nsp2 nsp3 nsp4)")
    parser.add_argument("-o", "--output", help="Output directory (if not provided, it will be created based on the FASTA file name)")
    
    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Call the main function to process the FASTA file with the provided keywords and output directory
    process_fasta_file(args.fasta, args.keywords, args.output)

# If this script is run directly (not imported), execute the main function
if __name__ == "__main__":
    main()
