import sys
import pandas as pd
from Bio import SeqIO
import argparse



def read_positions_table(positions_file):
    # Reads the entire positions table using pandas
    df = pd.read_csv(positions_file, sep=None, engine='python')
    return df


def extract_alleles(fasta_file, positions):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    data = []
    for rec in records:
        alleles = [rec.seq[pos-1] for pos in positions]  # 1-based to 0-based
        data.append([rec.id] + alleles)
    columns = ['id'] + [f'pos_{p}' for p in positions]
    df = pd.DataFrame(data, columns=columns)
    return df

def calculate_allele_frequencies(df, positions):
    freq_dict = {}
    for pos in positions:
        col = f'pos_{pos}'
        counts = df[col].value_counts(normalize=True)
        freq_dict[col] = counts
    freq_df = pd.DataFrame(freq_dict).fillna(0)
    return freq_df

if __name__ == "__main__":
    #get command line arguments
    parser = argparse.ArgumentParser(description="Extract alleles and calculate frequencies from fasta and positions file.")
    parser.add_argument("-p", "--positions", required=True, help="Path to positions.txt file")
    parser.add_argument("-f", "--fasta", required=True, help="Path to alignments.fasta file")
    parser.add_argument("-o", "--output", required=False, help="Path to output file", default="output.txt")
    parser.add_argument("-s", "--samples", required=False, help="Path to samples.txt file (list of sample IDs)")
    args = parser.parse_args()

    positions_file = args.positions or None
    fasta_file = args.fasta or None
    outprefix = args.output
    samples_file = args.samples

    if not positions_file or not fasta_file:
        print("Usage: python pull_escapes.py -p <positions.txt> -f <alignments.fasta>")
        sys.exit(2)

    # Read sample IDs if provided
    sample_ids = None
    if samples_file:
        with open(samples_file) as f:
            sample_ids = set(line.strip() for line in f if line.strip())

    #get table of known escape mutations, pull list of positions
    muttable = read_positions_table(positions_file)
    positions = muttable['Site'].tolist()
    
    # Extract alleles from the fasta file for each position
    allelesdf = extract_alleles(fasta_file, positions)
    if sample_ids is not None:
        allelesdf = allelesdf[allelesdf['id'].isin(sample_ids)]

    #calculate allele frequencies for each locus
    freq_df = calculate_allele_frequencies(allelesdf, positions)
    print("\nAllele frequencies:")
    print(freq_df.round(4))

    # Iterate through the allele table, find alt alleles, print consequence
    results = []
    for i in range(allelesdf.shape[0]):
        sample_row = allelesdf.iloc[i]
        sample_id = sample_row['id']
        for j in range(muttable.shape[0]):
            col = f'pos_{muttable.at[j, "Site"]}'
            if sample_row[col] == muttable.at[j, "Alt"]:
                results.append({
                    'sample': sample_id,
                    'position': muttable.at[j, "Site"],
                    'allele': muttable.at[j, "Alt"],
                    'mutation': muttable.at[j, "Mutation"],
                    'affected_antibody': muttable.at[j, "Affected Antibody"],
                    'effect': muttable.at[j, "Effect"]
                })
    escapesdf = pd.DataFrame(results)
    print("\nSamples with Alt alleles:")
    print(escapesdf)

    # Save the results to the output file
    escapesdf.to_csv(outprefix + "_calls.txt", sep='\t', index=False)
    allelesdf.to_csv(outprefix + "_alleles.txt", sep='\t', index=False)
    freq_df.round(4).to_csv(outprefix + "_allele_frequencies.txt", sep='\t', index=True)