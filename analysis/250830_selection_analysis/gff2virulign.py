#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_gff_line(line):
    """Parse one GFF3 line into dict of fields."""
    chrom, source, feature, start, end, score, strand, phase, attrs = line.strip().split("\t")
    attr_dict = {}
    for field in attrs.split(";"):
        if "=" in field:
            k, v = field.split("=")
            attr_dict[k] = v
    return {
        "feature": feature,
        "start": int(start),
        "end": int(end),
        "strand": strand,
        "phase": int(phase) if phase != "." else 0,
        "name": attr_dict.get("Name", attr_dict.get("ID", "CDS")),
        "gene": attr_dict.get("gene", attr_dict.get("ID", "CDS"))
    }

def main():
    parser = argparse.ArgumentParser(description="Convert FASTA + GFF3 to VIRULIGN XML")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file")
    parser.add_argument("--gff", required=True, help="GFF3 annotation file")
    parser.add_argument("--output", required=True, help="Output VIRULIGN XML file")
    args = parser.parse_args()

    # Load reference sequence
    ref = SeqIO.read(args.fasta, "fasta")
    ref_seq = str(ref.seq).upper()

    genes = []
    with open(args.gff) as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue
            entry = parse_gff_line(line)
            if entry["feature"] != "CDS":
                continue

            frame = entry["phase"] + 1  # GFF3 phase → VIRULIGN frame
            start, end = entry["start"], entry["end"]

            if entry["strand"] == "-":
                # Flip coordinates to forward strand
                start, end = len(ref_seq) - end + 1, len(ref_seq) - start + 1
                # NOTE: VIRULIGN requires + strand reference. 
                # Here we just flip coordinates, assuming the reference seq is plus-strand oriented.

            genes.append((entry['gene'], start, end, frame))

    # Write XML
    with open(args.output, "w") as out:
        out.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        out.write("<root>\n")
        out.write(f'  <reference label="{ref.id}" sequence="{ref_seq}">\n')
        for name, start, end, frame in genes:
            out.write(f'    <gene name="{name}" start="{start}" stop="{end}" frame="{frame}"/>\n')
        out.write("  </reference>\n")
        out.write("</root>\n")

if __name__ == "__main__":
    main()
