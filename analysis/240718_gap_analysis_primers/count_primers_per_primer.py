import os
from collections import defaultdict

# Function to read and count primer occurrences from a single output file
def count_primers_in_file(file_path, primer_counts):
    try:
        with open(file_path, "r") as file:
            print(f"Reading file: {file_path}")  # Debug statement
            for line in file:
                if line.startswith("Matched Primer:"):
                    primer = line.strip().split(": ", 1)[1]
                    primer_counts[primer] += 1
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")

# Function to process all output files in the output directory
def count_primers_in_directory(output_dir):
    primer_counts = defaultdict(int)
    
    # Iterate over all TXT files in the output directory
    for file_name in os.listdir(output_dir):
        if file_name.endswith("_matched_reads.txt"):
            file_path = os.path.join(output_dir, file_name)
            count_primers_in_file(file_path, primer_counts)
    
    return primer_counts

# Function to write the summary table to a CSV file
def write_summary_to_csv(primer_counts, output_file):
    print(f"Writing summary to {output_file}")  # Debug statement
    with open(output_file, "w") as file:
        file.write("Primer,Count\n")
        for primer, count in primer_counts.items():
            file.write(f"{primer},{count}\n")
            print(f"{primer}: {count}")  # Debug statement

# Define the output directory containing matched reads files
output_directory = "/gpfs/gibbs/project/grubaugh/emd224/PGCOE/Analysis/Primer_Efficency/Primers"

# Define the output CSV file for the summary
summary_file = "/gpfs/gibbs/project/grubaugh/emd224/PGCOE/Analysis/Primer_Efficency/primer_summary.csv"

# Count primers in the output directory
primer_counts = count_primers_in_directory(output_directory)

# Write the summary to a CSV file
write_summary_to_csv(primer_counts, summary_file)

print("Primer count summary complete.")
