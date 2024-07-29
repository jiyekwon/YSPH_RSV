import os
import gzip
from Bio import SeqIO
import glob
from tqdm import tqdm

# Function to check if a sequence starts with any of the primers
def match_primers(sequence, primers):
    # Check the first 50 nucleotides
    start_seq = sequence[:50]
    # Check the last 50 nucleotides
    end_seq = sequence[-50:]
    
    # Check if either the start or end sequence matches any primer
    for primer in primers:
        if start_seq.startswith(primer) or end_seq.startswith(primer):
            return primer
    return None

# Function to process a single FASTQ file
def process_fastq_file(fastq_path, primers, output_dir):
    sample_name = os.path.basename(os.path.dirname(fastq_path))
    output_file = os.path.join(output_dir, f"{sample_name}_matched_reads.txt")

    # Debug statement to check file paths
    print(f"Processing file: {fastq_path}")
    print(f"Output file: {output_file}")

    try:
        with gzip.open(fastq_path, "rt") as handle, open(output_file, "w") as output_handle:
            for record in SeqIO.parse(handle, "fastq"):
                match = match_primers(str(record.seq)[:50], primers)
                if match:
                    output_handle.write(f"> {record.id}\n")
                    output_handle.write(f"{record.seq}\n")
                    output_handle.write(f"+\n")
                    output_handle.write(f"{record.letter_annotations['phred_quality']}\n")
                    output_handle.write(f"Matched Primer: {match}\n\n")
    except Exception as e:
        print(f"Error processing file {fastq_path}: {e}")

# Function to process all FASTQ files in a single folder
def process_fastq_files_in_folder(folder_path, primers, output_dir):
    # Updated glob pattern to find fastq.gz files in the specific folder
    fastq_files = glob.glob(os.path.join(folder_path, "*.fastq.gz"))
    
    # Debug statement to check found files
    print(f"Found FASTQ files in folder: {folder_path}")
    print(f"Files: {fastq_files}")

    # Process each FASTQ file in the folder
    for fastq_file in tqdm(fastq_files, desc=f"Processing folder: {folder_path}"):
        process_fastq_file(fastq_file, primers, output_dir)

# Function to process all folders within the main directory
def process_all_folders_in_directory(root_dir, primers, output_dir):
    # Iterate over each subfolder in the root directory
    for subfolder in glob.glob(os.path.join(root_dir, "*")):
        if os.path.isdir(subfolder):
            print(f"Processing folder: {subfolder}")
            process_fastq_files_in_folder(subfolder, primers, output_dir)

# List of primers to match
primers = [
   "RSVA_22_Left": "AGTTCCAACAAAAGAACAACAGACT",
    "RSVA_20_RIGHT": "CAGGACCTTGGATACGGCAAT",
    "RSVA_26_RIGHT": "CAGTAAGGAGTTTGCTCATGGC",
    "RSVAB_16_LEFT": "GGGCAAATGCAAACATGTCCAA",
    "RSVA_15_LEFT": "TCTGGCCTTACTTTACACTAATACACA",
    "RSVB_15_LEFT": "TTGGCCCTATTTTACACTAATACATATGA",
    "RSVAB_17_LEFT": "CTTCTCCAATCTGTCCGGAACT",
    "RSVA_16_RIGHT": "TTGTGGATTGTGGGGTTGACTC",
    "RSVB_16_RIGHT": "TTGGGTGATATTGTGGCTGAGT",
    "RSVAB_15_RIGHT": "AGATGATTGAGAGTGTCCCAGGT",
    "RSVAB_17_RIGHT": "TGATGGTTGGCTTTCCTGTAGG",
    "RSVBA_44_LEFT": "GCTGGACATTGGATTCTGATTATTCA",
    "RSVA_18_RIGHT": "AGATTGTGATGGGTACTCGGATG",
    "RSVB_45_LEFT": "AGTTTGCATAGAATAAAAGGTTGTCACA",
    "RSVAB_44_RIGHT": "TAACAACCCAAGGGCAAACTGT",
    "RSVB_45_RIGHT": "CCACTTATACAAAATTTAGGTTTGTTTCCA",
    "RSVAB_48_LEFT": "GGTGAAAATTTGACCATTCCTGCTA",
    "RSVAB_48_RIGHT": "GGCATGATGAAATTTTTGGTTCTTGA",
    "hRSVA_7_LEFT_mod": "TTATGAATGCCTATGGTGCAGGG",
    "RSVA_7_RIGHT": "CTTCTCCATGGAATTCAGGAGCA"
]

# Define the root directory containing multiple subfolders
root_directory = "/home/emd224/ycga_work/Raw_data/RSV004"

# Define the output directory for matched reads
output_directory = "/gpfs/gibbs/project/grubaugh/emd224/PGCOE/Analysis/Primer_Efficency/Primers"
os.makedirs(output_directory, exist_ok=True)

# Process all folders in the root directory
process_all_folders_in_directory(root_directory, primers, output_directory)

print("Processing complete.")
