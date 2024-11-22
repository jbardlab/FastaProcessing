from Bio import SeqIO  # Import SeqIO from Biopython
import pandas as pd  # Import pandas for data manipulation
import argparse  # Import argparse for command-line argument parsing
import os  # Import os for handling file paths

# Function to load the metadata file (CSV or Excel)
def load_metadata(metadata_file):
    _, file_extension = os.path.splitext(metadata_file)

    if file_extension.lower() == ".csv":
        return pd.read_csv(metadata_file, low_memory=False)
    elif file_extension.lower() == ".xlsx":
        return pd.read_excel(metadata_file)
    else:
        raise ValueError("Unsupported file type. Please provide a CSV or Excel file for the metadata.")

# Main processing function
def process_fasta_and_metadata(fasta_file, metadata_file, output_file):
    seq_list = []  # List to hold extracted SeqID, genome_id, and name

    # Iterate through the FASTA entries
    for record in SeqIO.parse(fasta_file, "fasta"):
        parts = record.description.split(" ")  # Split the description
        seq_id = record.id  # Extract SeqID
        genome_id = parts[0].split("|")[1]  # Extract genome_id from the description

        # Extract the name from the description
        name = " ".join(parts[1:]).split('[')[0].strip()  # Get everything after the first space

        # Append the extracted data to the list
        seq_list.append((seq_id, genome_id, name))

    # Create a DataFrame from the list
    seq_df = pd.DataFrame(seq_list, columns=["SeqID", "genome_id", "Name"])

    # Debugging statement: Print the DataFrame after creation
    print("DataFrame from FASTA entries:")
    print(seq_df.head())  # Show the first few rows of seq_df

    # Clean genome_id to match 'Genome ID' format in the metadata DataFrame
    seq_df["genome_id"] = seq_df["genome_id"].apply(lambda x: ".".join(x.split(".")[:2]))

    # Load the metadata
    metadata_df = load_metadata(metadata_file)

    # Debugging statement: Print the metadata DataFrame
    print("Metadata DataFrame:")
    print(metadata_df.head())  # Show the first few rows of metadata_df

    # Ensure 'Genome ID' is treated as string
    metadata_df["Genome ID"] = metadata_df["Genome ID"].astype(str).str.strip()

    # Merge DataFrames on genome_id
    merged_df = pd.merge(seq_df, metadata_df, left_on="genome_id", right_on="Genome ID", how="left")

    # Debugging statement: Print the merged DataFrame before dropping the column
    print("Merged DataFrame before dropping 'genome_id':")
    print(merged_df.head())  # Show the first few rows of merged_df

    # Drop the genome_id column
    merged_df.drop(columns=['genome_id'], inplace=True)

    # Save the merged DataFrame to a TSV file
    merged_df.to_csv(output_file, sep="\t", index=False)

    # Final debugging statement: Print the final merged DataFrame
    print("Final merged DataFrame:")
    print(merged_df.head())  # Show the first few rows of the final merged_df
    print(f"Process complete! The merged output has been saved to {output_file}")

# Main function for argument parsing
def main():
    parser = argparse.ArgumentParser(description="Process FASTA and metadata files, and merge them into a TSV file.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file.")
    parser.add_argument("metadata_file", help="Path to the input metadata CSV or Excel file.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    args = parser.parse_args()
    
    # Call the processing function
    process_fasta_and_metadata(args.fasta_file, args.metadata_file, args.output_file)

# Entry point
if __name__ == "__main__":
    main()

# Command to run the script:
# python fasta_to_tsv.py combined_proteins.fasta BVBRC_genome.csv merged_output.tsv
