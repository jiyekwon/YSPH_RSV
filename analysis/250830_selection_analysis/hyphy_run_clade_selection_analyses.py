
from Bio.Seq import Seq
import subprocess

from Bio import SeqIO
from Bio import Phylo

import pandas as pd
import json
import argparse
import tempfile
import os

tmpdir = None

def remove_zero_length_branches(tree):
    """
    Removes clades with a branch_length of 0 from a Biopython Tree object.
    This effectively prunes zero-length branches and creates polytomies.
    """
    # Create a list to store clades to be processed (children of zero-length clades)
    clades_to_process = []

    # Iterate through all clades in the tree
    for clade in tree.find_clades():
        # Check if the clade has children and a branch_length of 0
        if clade.branch_length == 0 and clade.clades:
            clades_to_process.append(clade)
    print(f"Removing {len(clades_to_process)} zero-length branches")
    # Process clades with zero-length branches
    for zero_clade in clades_to_process:
        parent = tree.get_path(zero_clade)[-2]  # Get the parent of the zero-length clade
        
        if parent:
            # Remove the zero-length clade from its parent's children
            parent.clades.remove(zero_clade)
            # Add the children of the zero-length clade directly to the parent
            parent.clades.extend(zero_clade.clades)

    return tree

def prune_tree_to_seqs(wgstree, fastaf, pruned_tree_file):
    """
    Prune a tree to only include tips present in the fasta file.
    Args:
        wgstree (str): Path to the input tree file (Newick format).
        fastaf (str): Path to the fasta file containing sequence names to keep.
        pruned_tree_file (str): Path to write the pruned tree.
    """
    seq_names = set(SeqIO.to_dict(SeqIO.parse(fastaf, "fasta")).keys())
    
    # Read the tree
    tree = Phylo.read(wgstree, "newick")
    tips = [c.name for c in tree.get_terminals()]
    # Find tips to prune (not in fasta)
    tips_to_prune = [term for term in tips if term not in seq_names]
    
    # Prune tips
    if(len(tips_to_prune) == len(tips)):
        print(len(tips), len(seq_names), 0)
        exit(1)
    else:
        for tip in tips_to_prune:
            tree.prune(tip)    
        #remove zero-length branches
        nozerotree = remove_zero_length_branches(tree) 

        # Write pruned tree
        print(len(tips), len(seq_names), len(tree.get_terminals()))

        Phylo.write(nozerotree, pruned_tree_file, "newick")
    return(len(tree.get_terminals()))

def filter_seqs_to_tree(fastaf, treefile, output_fasta):
    """
    Filter sequences in a FASTA file to only those present as tips in a tree file.
    Args:
        fastaf (str): Path to input FASTA file.
        treefile (str): Path to tree file (Newick format).
        output_fasta (str): Path to output filtered FASTA file.
    """
    # Get tip names from tree
    tree = Phylo.read(treefile, "newick")
    tip_names = set([term.name for term in tree.get_terminals()])
    # Parse FASTA and filter
    records = [rec for rec in SeqIO.parse(fastaf, "fasta")]
    nrecords = len(records)
    goodrecords = [rec for rec in records if rec.id in tip_names]
    # Write filtered FASTA
    with open(output_fasta, "w") as out:
        SeqIO.write(goodrecords, out, "fasta")

    print(nrecords, len(tip_names), len(goodrecords))
    return(len(goodrecords))

# %%
def prune_tree_to_samples(wgstree, samples, pruned_tree_file):
    """
    Prune a tree to only include tips present in the samples file.
    Args:
        wgstree (str): Path to the input tree file (Newick format).
        samples (str): Path to the samples file (one sample name per line).
        pruned_tree_file (str): Path to write the pruned tree.
    """
    with open(samples) as f:
        sample_names = set(line.strip() for line in f if line.strip())
    
    # Read the tree
    tree = Phylo.read(wgstree, "newick")
    tips = [c.name for c in tree.get_terminals()]
    # Find tips to prune (not in samples)
    tips_to_prune = [term for term in tips if term not in sample_names]
    
    # Prune tips
    if(len(tips_to_prune) == len(tips)):
        print(len(tips), len(sample_names), 0)
        exit(1)
    else:
        for tip in tips_to_prune:
            tree.prune(tip)    
        #remove zero-length branches
        nozerotree = remove_zero_length_branches(tree) 

        # Write pruned tree
        print(len(tips), len(sample_names), len(tree.get_terminals()))
        Phylo.write(nozerotree, pruned_tree_file, "newick")
    return(len(tree.get_terminals()))

def filter_seqs_to_samples(fastaf, samples, output_fasta):
    """
    Filter sequences in a FASTA file to only those present in the samples file.
    Args:
        fastaf (str): Path to input FASTA file.
        samples (str): Path to the samples file (one sample name per line).
        output_fasta (str): Path to output filtered FASTA file.
    """
    with open(samples) as f:
        sample_names = set(line.strip() for line in f if line.strip())
    
    # Parse FASTA and filter
    records = [rec for rec in SeqIO.parse(fastaf, "fasta")]
    nrecords = len(records)
    goodrecords = [rec for rec in records if rec.id in sample_names]
    # Write filtered FASTA
    with open(output_fasta, "w") as out:
        SeqIO.write(goodrecords, out, "fasta")

    print(nrecords, len(sample_names), len(goodrecords))
    return(len(goodrecords))

# %%
def remove_duplicate_sequences(fasta_file, output_file):
    """
    Reads a FASTA file, removes duplicate sequences (keeping the one with the lowest alphanumeric ID),
    and writes the unique sequences to output_file.
    """

    seq_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_str = str(record.seq)
        if seq_str not in seq_dict or record.id < seq_dict[seq_str].id:
            seq_dict[seq_str] = record

    with open(output_file, "w") as out:
        SeqIO.write(seq_dict.values(), out, "fasta")

# %%
def all_descendants(clade):
    """
    Recursively get all descendant clades of a given clade.
    """
    descendants = []
    for child_clade in clade.clades:
        descendants.append(child_clade)
        descendants.extend(all_descendants(child_clade))
    return descendants


#annotate tree with clade assignments
def annotate_tree(treefile,outtree,clades,clade=None,strategy="MRCA"):
    """
    Annotate a tree with clade assignments.
    Args:
        treefile (str): Path to the input tree file (Newick format).
        outtree (str): Path to write the annotated tree.
        clades (dict): Dictionary mapping clade names to lists of tip names.
        clade (str, optional): Specific clade to annotate. If None, all clades are annotated.
        strategy (str): Strategy for annotation. Options are "MRCA", "tips", "path", "all".
    """
    tree = Phylo.read(treefile, "newick")
    if clade is not None: cladenames = [clade]
    else: cladenames = clades.keys()
    
    for C in cladenames:
        cladelist = clades[C]
        print(C, len(clades[C]))
        tip_names = set([term.name for term in tree.get_terminals()])
        goodtips = [rec for rec in cladelist if rec in tip_names]
        print(f"Annotating clade {C} with {len(goodtips)} tips in length {len(tip_names)} tree")
        
        tag = '{Foreground}'
        mrca = tree.common_ancestor(goodtips)
        # Set the color and width of the branch leading to this clade
        cladestr = '{'+f'{C}'+'}'
        if  strategy == "MRCA":  # Tag MRCA
            mrca.name = f'{mrca.name}{cladestr}{tag}'
        if  strategy == "tips":  # Tag tips
            for C in [c for c in tree.get_terminals() if c.name in goodtips]:
                cladestr = '{'+f'{C}'+'}'
                C.name = f'{C.name}{cladestr}{tag}'
        if strategy == "path":  # Tag path from tips to MRCA
            pathnodes = {}
            for tip in goodtips:
                #append only new nodes in path to list of paths
                for n in get_path_to_ancestor(tree, tip, mrca):
                    if n.name not in pathnodes:
                        pathnodes[n.name] = n
            pathnodes = [pathnodes[p] for p in pathnodes]
            for p in pathnodes:
                cladestr = '{'+f'{C}'+'}'
                for n in p:
                    n.name = f'{n.name}{cladestr}{tag}'
            print(f"Marking {len(goodtips)} of {len(tree.get_terminals())} ({len(pathnodes)} nodes, mrca {mrca.name})")
        if  strategy == "all":  # Tag MRCA and all descendants
            mrca.name = f'{mrca.name}{cladestr}{tag}'
            descs = all_descendants(mrca)
            for dclade in descs:
                cladestr = '{'+f'{C}'+'}'
                dclade.name = f'{dclade.name}{cladestr}{tag}'
        # Draw the tree (e.g., to ASCII or using Matplotlib)
        #Phylo.draw_ascii(tree)
            print(f"Marking {len(goodtips)} of {len(tree.get_terminals())} ({len(descs)} nodes, mrca {mrca.name})")
    Phylo.write(tree, outtree, "newick")
    return len(goodtips)


# Define a function to find the path from a tip to an ancestor
def get_path_to_ancestor(tree, tip, ancestor):
    """
    Get the path from a tip to a specified ancestor in the tree.
    Args:
        tree (Bio.Phylo.BaseTree.Tree): The phylogenetic tree.
        tip (str): The name of the tip node.
        ancestor (Bio.Phylo.BaseTree.Clade): The ancestor clade.
    Returns:
        list: List of clades from the tip to the ancestor.
    """
    path_to_root = tree.get_path(tip)
    path_to_mrca = []
    # Iterate through the path from root and stop at the ancestor
    for clade in reversed(path_to_root):
        path_to_mrca.append(clade)
        if clade == ancestor:
            break
    return path_to_mrca




# %%
# Run RELAX using HyPhy on the codon alignment and IQ-TREE treefile


def run_relax(codontrimf, treefile, relaxout, tmpdir=None):
    """Run RELAX using HyPhy on the codon alignment and IQ-TREE treefile."""
    if tmpdir is None:
        errfile = f'{os.path.basename(relaxout)}.err'
        outfile = f'{os.path.basename(relaxout)}.out'
    else:
        errfile = os.path.join(tmpdir.name, f'{os.path.basename(relaxout)}.err')
        outfile = os.path.join(tmpdir.name, f'{os.path.basename(relaxout)}.out')
    command = [
        "hyphy", "relax",
        "--alignment", codontrimf,
        "--tree", treefile,
        "--output", relaxout,
        "--models", "Minimal",
        "--test", "Foreground"
    ]
    with open(errfile, "w") as errfile, open(outfile, "w") as outfile:
        print(" ".join(command))
        subprocess.run(command, check=True,stderr=errfile, stdout=outfile)

def parse_relax_json(relax_json_file):
    """
    Parse RELAX JSON output and extract all values within 'test results'.
    Returns a dictionary of the extracted values.
    """
    with open(relax_json_file) as f:
        data = json.load(f)
    test_results = data.get('test results')
    print(test_results)
    LR = test_results.get('LRT')
    p = test_results.get('p-value')
    K = test_results.get('relaxation or intensification parameter')
    return LR, p, K



def run_busted(codontrimf, treefile, bustedout, tmpdir=None):
    """Run BUSTED using HyPhy on the codon alignment and IQ-TREE treefile."""
    if tmpdir is None:
        errfile = f'{os.path.basename(bustedout)}.err'
        outfile = f'{os.path.basename(bustedout)}.out'
    else:
        errfile = os.path.join(tmpdir.name, f'{os.path.basename(bustedout)}.err')
        outfile = os.path.join(tmpdir.name, f'{os.path.basename(bustedout)}.out')
    command = [
        "hyphy", "busted",
        "--alignment", codontrimf,
        "--tree", treefile,
        "--output", bustedout,
        "--branches", "Foreground"
    ]
    with open(errfile, "w") as errfile, open(outfile, "w") as outfile:
        print(" ".join(command))
        #subprocess.run(command, check=True,stderr=errfile, stdout=outfile)
        subprocess.run(command, check=True,stderr=errfile)



def main():
    parser = argparse.ArgumentParser(description="Run HyPhy RELAX analysis on a codon alignment with clade annotation.")
    parser.add_argument("-a", "--codon_alignment", required=True, help="Path to codon alignment FASTA file")
    parser.add_argument("-C", "--clade_table", required=True, help="Path to clade assignment table (TSV with seqName,clade,qc.overallStatus)")
    parser.add_argument("-c", "--clade", required=True, help="Clade ID to analyze")
    parser.add_argument("-p", "--tree", required=True, help="Path to Newick tree file")
    parser.add_argument("-g", "--gene", required=True, help="Gene ID")
    parser.add_argument("-M", "--majorclade", action="store_true", help="Use major clade")
    parser.add_argument("-t", "--target", required=True, help="Target name (reference / virus name)")
    parser.add_argument("-o", "--output", required=False, help="Output prefix", default="hyphy")
    args = parser.parse_args()

    # Assign variables for downstream code
    codonf = args.codon_alignment       #nucleotide codon alignment
    nextcladeF = args.clade_table
    wgstree = args.tree
    gene = args.gene
    clade = args.clade
    target = args.target
    majorclade = args.majorclade
    prefix = args.output

    # Create a temporary directory for intermediate files
    tmpdir = tempfile.TemporaryDirectory()
    print(f"Temporary directory created at: {tmpdir.name}")

    # Optionally, you can use tmpdir.name to construct paths for temp files
    codonfuniq = os.path.join(tmpdir.name, os.path.basename(codonf).replace(".fasta", "_uniq.fasta"))
    wgstrim = os.path.join(tmpdir.name, os.path.basename(wgstree).replace(".nwk", "_uniq.nwk"))
    wgsclade = os.path.join(tmpdir.name, os.path.basename(wgstrim).replace(".nwk", f"_{clade}.nwk"))



    # The rest of the analysis code would go here, using these variables
    df_clades = pd.read_csv(nextcladeF, usecols=['seqName', 'clade', 'qc.overallStatus'], sep="\t")

    # Add 'majorclade' column by removing everything after the second '.' in 'clade'
    if majorclade==True:
        df_clades['majorclade'] = df_clades['clade'].str.split('.').str[:4].str.join('.')
        clades = df_clades.groupby('majorclade')['seqName'].apply(list).to_dict()
    else:
        clades = df_clades.groupby('clade')['seqName'].apply(list).to_dict()

    #make unique annotated tree and codon alignment
    prune_tree_to_seqs(wgstree, codonf, wgstrim)
    remove_duplicate_sequences(codonf, codonfuniq)
    prune_tree_to_seqs(wgstrim, codonfuniq, wgstrim)
    filter_seqs_to_tree(codonfuniq, wgstrim, codonfuniq)
    ntips = annotate_tree(wgstrim,wgsclade,clades,clade=clade,strategy="path")
    

    outcsv = f"{prefix}_{target}_{gene}_{clade}_relax.csv"               #RELAX output csv

    if(ntips >= 5):
        #hyphy outputs
        outjson = f"{prefix}_{target}_{gene}_{clade}_relax.json"               #RELAX output json
        run_relax(codonfuniq, wgsclade, outjson, tmpdir=tmpdir)
        LR, p, K = parse_relax_json(outjson)
        outjson = f"{prefix}_{target}_{gene}_{clade}_busted.json"               #BUSTED output json
        run_busted(codonfuniq, wgsclade, outjson, tmpdir=tmpdir)

        results = pd.DataFrame([[target, gene, clade, ntips, LR, p, K]], 
                                columns=["target","gene","clade","n_tips","LR","p","K"]) 
    else: 
        print(f"Skipping {target} {gene} {clade} with {ntips} tips")
        LR, p, K = None, None, None
        results = pd.DataFrame([[target, gene, clade, ntips, LR, p, K]], 
                                columns=["target","gene","clade","n_tips","LR","p","K"]) 
    results.to_csv(outcsv, index=False)





if __name__ == "__main__":
    main()

