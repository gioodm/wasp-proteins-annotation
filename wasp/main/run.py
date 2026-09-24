#!/usr/bin/env python3
# @author Giorgia Del Missier

import argparse, os, sys, subprocess, glob, shutil
from pathlib import Path
import numpy as np
import networkx as nx
import pandas as pd
import json

from .fs_besthits import get_besthits, save_besthits
from .parse_fs import parse_m8, save_json
from .generate_network import run_network_generation, save_network
from .retrieve_annotation import fetch_annotations
from .SAFE_enrichment import run_safe_analysis

# Set random seed for reproducibility
np.random.seed(0)

def show_help():
    help_message = """
Usage: python3 run.py [-h] (-t taxid | -f proteins.fasta | --structures proteins.tar) [-e evalue_threshold] [-b bitscore_threshold] [-n max_neighbours] [-s step] [-i iterations] [-m max_seq_id]

WASP (Whole-proteome Annotation through Structural homology Pipeline) performs a "structural BLAST" using AlphaFold models to better annotate the target taxid proteome.
Parameters:

    -h, --help                  show this help message and exit
    -t, --taxid                 NCBI taxonomy identifier to be analysed
    -f, --fasta                 FASTA file to analyse using Foldseek ProstT5
    -p, --structures            tar archive containing .cif.gz or .pdb.gz structures
    -e, --evalue_thr            set the evalue threshold (default: 10e-10)
    -b, --bitscore_thr          set the bitscore threshold (default: 50)
    -n, --max_n                 set the max number of neighbours (default: 10)
    -s, --step                  set step to add to max neighbours (n) in additional iterations (default: 10)
    -i, --iters                 set number of iterations to perform (default: 3)
    -m, --max_seq_id            set maximum sequence identity (0.0 to 1.0) to exclude sequence homologs (default: 1.0)

Examples:
    python3 run.py -t 559292
    python3 run.py -t 559292 -e 1e-50 -b 200 -n 5 -i 5
    python3 run.py -t 559292 -s 5
    python3 run.py -f proteins.fasta
    python3 run.py --structures proteins.tar
    """
    print(help_message)
 
def check_command(command):
    try:
        subprocess.run([command], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        print(f"{command} is required but it's not installed. Aborting.")
        sys.exit(1)

def main():
    if "-h" in sys.argv or "--help" in sys.argv:
        show_help()
        sys.exit(0)

    parser = argparse.ArgumentParser(description="WASP Pipeline", add_help=False)

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-t", "--taxid")
    input_group.add_argument("-f", "--fasta", type=Path)
    input_group.add_argument("-p", "--structures", type=Path)
    parser.add_argument("-e", "--eval_thr", type=float, default=1e-10)
    parser.add_argument("-b", "--bits_thr", type=int, default=50)
    parser.add_argument("-n", "--max_n", type=int, default=10)
    parser.add_argument("-s", "--step", type=int, default=10)
    parser.add_argument("-i", "--iters", type=int, default=3)
    parser.add_argument("-m", "--max_seq_id", type=float, default=1.0)

    args = parser.parse_args()

    input_file = args.fasta or args.structures
    if input_file and not input_file.is_file():
        parser.error(f"Input file not found: {input_file}")

    input_id = args.taxid or input_file.stem
    is_fasta_input = args.fasta is not None
    is_structure_input = args.structures is not None
    # FASTA/ProstT5 input has no coordinate (_ca) database to concatenate or subset.
    ca_suffixes = [] if is_fasta_input else ["_ca"]

    system_tmp = os.environ.get("TMPDIR", "/tmp")
    fs_tmp = f"{system_tmp}/fs_tmp_{input_id}"
    os.makedirs(fs_tmp, exist_ok=True)
    os.chmod(fs_tmp, 0o755)

    # Check required commands
    check_command("foldseek")
    if not (is_fasta_input or is_structure_input):
        check_command("gsutil")

    print("Setting required variables:")
    if is_fasta_input:
        print(f"\nSelected FASTA file is: {args.fasta}")
    elif is_structure_input:
        print(f"\nSelected structure archive is: {args.structures}")
    else:
        print(f"\nSelected taxid is: {args.taxid}")
    print(f"Selected evalue threshold is: {args.eval_thr}")
    print(f"Selected bitscore threshold is: {args.bits_thr}")
    print(f"Selected max neighbours is: {args.max_n}")
    print(f"Selected step is: {args.step}")
    print(f"Selected number of iterations is: {args.iters}")
    print(f"Selected maximum sequence identity threshold is: {args.max_seq_id}")

    ####---- DOWNLOADING FILES AND DATABASES ----####

    # Define directory names
    db_dir = "foldseek_dbs"
    prot_dir = "proteomes"
    results_dir = "results"
    taxid_dir = f"{results_dir}/{input_id}"
    protein_input = str(input_file) if (is_fasta_input or is_structure_input) else f"{prot_dir}/{input_id}.tar"
    source_db = f"{db_dir}/{input_id}"
    combined_db = f"{db_dir}/afdb50sp{input_id}"
    prostt5_weights = f"{db_dir}/prostt5"

    # Create directories
    for directory in [db_dir, prot_dir, results_dir, taxid_dir, f"{taxid_dir}/SAFE"]:
        os.makedirs(directory, exist_ok=True)

    print("\nDownloading required AlphaFold models:")

    # Download and prepare databases
    if not os.path.exists(f"{db_dir}/afdb50sp"):
        if not os.path.exists(f"{db_dir}/afdb50"):
            subprocess.run(["foldseek", "databases", "Alphafold/UniProt50-minimal", f"{db_dir}/afdb50", "tmp", "--remove-tmp-files", "1"])
        if not os.path.exists(f"{db_dir}/swissprot"):
            subprocess.run(["foldseek", "databases", "Alphafold/Swiss-Prot", f"{db_dir}/swissprot", "tmp", "--remove-tmp-files", "1"])

        # Merge the databases
        for suffix in ["", "_h", "_ss", "_ca"]:
            subprocess.run(["foldseek", "concatdbs", f"{db_dir}/afdb50{suffix}", f"{db_dir}/swissprot{suffix}", f"{db_dir}/afdb50sp{suffix}"])
    else:
        print("Foldseek databases already downloaded")

    # Check if results already exist
    if all(os.path.exists(f"{taxid_dir}/{input_id}{suffix}.m8") for suffix in ["", "_bh", "_norm", "_norm_bh"]):
        print("Input data already prepared and Foldseek results already generated")
    else:
        if is_fasta_input:
            if not os.path.exists(prostt5_weights):
                subprocess.run(["foldseek", "databases", "ProstT5", prostt5_weights, fs_tmp], check=True)
            if not os.path.exists(source_db):
                subprocess.run(["foldseek", "createdb", protein_input, source_db, "--prostt5-model", prostt5_weights], check=True)
        elif not is_structure_input:
            # Download and prepare the AlphaFold proteome.
            if not os.path.exists(protein_input):
                os.makedirs(f"{prot_dir}/{input_id}", exist_ok=True)
                subprocess.run(["gsutil", "-m", "cp", f"gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-{input_id}-*_v4.tar", prot_dir], check=True)

                for filename in os.listdir(prot_dir):
                    if filename.startswith(f"proteome-tax_id-{input_id}") and filename.endswith("_v4.tar"):
                        subprocess.run(["tar", "-xf", f"{prot_dir}/{filename}", "-C", f"{prot_dir}/{input_id}"], check=True)
                        os.remove(f"{prot_dir}/{filename}")

                for filename in os.listdir(f"{prot_dir}/{input_id}"):
                    if filename.endswith(".json.gz"):
                        os.remove(f"{prot_dir}/{input_id}/{filename}")

                subprocess.run(["tar", "-cf", protein_input, "-C", prot_dir, input_id], check=True)
                shutil.rmtree(f"{prot_dir}/{input_id}")
            if not os.path.exists(source_db):
                subprocess.run(["foldseek", "createdb", protein_input, source_db], check=True)

        if not os.path.exists(combined_db):
            for suffix in ["", "_h", "_ss"] + ca_suffixes:
                subprocess.run(["foldseek", "concatdbs", f"{db_dir}/afdb50sp{suffix}", f"{source_db}{suffix}", f"{combined_db}{suffix}"], check=True)
        else:
            print("Input database already prepared")

    ####---- RECIPROCAL BEST STRUCTURE HITS SEARCH ----####

    print("\nPerforming Reciprocal Best Structural Hits search:")

    if not os.path.exists(f"{taxid_dir}/{input_id}.m8") or not os.path.exists(f"{taxid_dir}/{input_id}_bh.m8"):
        # Perform foldseek searches
        prostt5_options = ["--prostt5-model", prostt5_weights] if is_fasta_input else []
        subprocess.run(["foldseek", "easy-search", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits",
                        protein_input, combined_db, f"{taxid_dir}/{input_id}.m8", fs_tmp, "--threads", "64", *prostt5_options], check=True)

        subprocess.run(["foldseek", "easy-search", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits",
                        protein_input, source_db, f"{taxid_dir}/{input_id}_norm.m8", fs_tmp, "--threads", "64",
                        "--exhaustive-search", "1", "--min-seq-id", "0.9", *prostt5_options], check=True)

        best_hits = get_besthits(f"{taxid_dir}/{input_id}.m8", 1, args.eval_thr, args.bits_thr, args.max_seq_id)
        save_besthits(best_hits, f"{taxid_dir}/{input_id}_bh.txt")

        subprocess.run(["foldseek", "prefixid", f"{combined_db}_h", f"{combined_db}.lookup", "--tsv", "--threads", "1"], check=True)
        
        subset_tsv = f"{db_dir}/subset{input_id}.tsv"
        subset_db = f"{db_dir}/subdb{input_id}"
        awk_command = ["awk", "NR == FNR {f[$1] = $1; next} $2 in f {print $1}", f"{taxid_dir}/{input_id}_bh.txt", f"{combined_db}.lookup"]
        with open(subset_tsv, "w") as output_file:
            subprocess.run(awk_command, stdout=output_file, check=True)

        for suffix in ["", "_ss"] + ca_suffixes:
            subprocess.run(["foldseek", "createsubdb", subset_tsv, f"{combined_db}{suffix}", f"{subset_db}{suffix}"], check=True)
        
        os.remove(subset_tsv)

        search_bh_prefix = f"{taxid_dir}/{input_id}_bh"
        subprocess.run(["foldseek", "search", subset_db, combined_db, search_bh_prefix, fs_tmp,"-a", "1", "--threads", "64"], check=True)

        subprocess.run(["foldseek", "convertalis", subset_db, combined_db, search_bh_prefix, f"{search_bh_prefix}.m8",
            "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits"], check=True)

        search_norm_prefix = f"{taxid_dir}/{input_id}_norm_bh"
        subprocess.run(["foldseek", "search", subset_db, subset_db, search_norm_prefix, fs_tmp, "-a", "1", "--threads", "64"], check=True)

        subprocess.run(["foldseek", "convertalis", subset_db, subset_db, search_norm_prefix, f"{search_norm_prefix}.m8",
            "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits"], check=True)

        # Clean up temporary files
        for path in glob.glob(os.path.join(db_dir, f"subdb{input_id}*")):
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)

        for path in glob.glob(os.path.join(db_dir, f"afdb50sp{input_id}*")):
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)

        for path in glob.glob(os.path.join(db_dir, f"{input_id}*")):
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)

        for root, dirs, files in os.walk(taxid_dir):
            for fname in files:
                if not fname.endswith(".txt") and not fname.endswith(".m8"):
                    fpath = os.path.join(root, fname)
                    os.remove(fpath)
    else:
        print("Foldseek results already generated")

    # Count the input proteins for the search summary.
    if is_fasta_input:
        with args.fasta.open() as fasta_file:
            psize = sum(line.startswith(">") for line in fasta_file)
    else:
        psize = subprocess.run(["tar", "-tvf", protein_input], capture_output=True, text=True, check=True)
        psize = len([line for line in psize.stdout.splitlines() if line.endswith(".gz")])

    queries = parse_m8(f"{taxid_dir}/{input_id}.m8", f"{taxid_dir}/{input_id}_norm.m8", args.eval_thr, args.bits_thr, args.max_seq_id)
    reciprocal_queries = parse_m8(f"{taxid_dir}/{input_id}_bh.m8", f"{taxid_dir}/{input_id}_norm_bh.m8", args.eval_thr, args.bits_thr, args.max_seq_id)
    
    # Saving the results to JSON files
    save_json(queries, f"{taxid_dir}/{input_id}.json")
    save_json(reciprocal_queries, f"{taxid_dir}/{input_id}_bh.json")

    print(f"Found significant hits for {len(queries)} out of {psize} proteins in the target organism")

    open(f"{taxid_dir}/{input_id}_nan.txt", "w").close()

    for j in range(1, args.iters + 1):
        print(f"\nPerforming iteration {j}")

        ####---- NETWORK GENERATION ----####

        print("\nCreating RBSH network and identifying clusters of homologs:")

        neighbours = args.max_n + (args.step * (j - 1))
        print(f"Selected max number of neighbours for iteration {j} is: {neighbours}\n")

        # Run the network generation
        G, all_queries, diff, clusters_sorted = run_network_generation(f"{taxid_dir}/{input_id}.json", f"{taxid_dir}/{input_id}_bh.json",
                                           f"{taxid_dir}/{input_id}_nan.txt", neighbours)
        # Save the network and cluster details
        save_network(G, diff, clusters_sorted, f"{taxid_dir}/{input_id}_clusters_iter{j}.txt", f"{taxid_dir}/{input_id}_edgelist_iter{j}.txt", f"{taxid_dir}/{input_id}_nan.txt")

        print(f"{len(diff)} proteins in the target organism had no RBSH hits... trying again with increased number of neighbours")

        # Print network statistics
        print(f"Network statistics generated using {len(all_queries) - len(diff)} RBSH hits:")
        print(f"Number of nodes: {len(G.nodes())}")
        print(f"Number of edges: {len(G.edges())}")
        print(f"Number of generated clusters: {nx.number_connected_components(G)}")
        
        ####---- ANNOTATION ----####

        print("\nAnnotating network's nodes")

        fetch_annotations(f"{taxid_dir}/{input_id}_clusters_iter{j}.txt", f"{taxid_dir}/{input_id}_annotation_iter{j}.txt")
        
        ####---- SAFE ENRICHMENT AND STATISTICS COMPUTATION ----####

        print("Performing SAFE analysis and computing statistics on new annotation\n")

        run_safe_analysis(f"{taxid_dir}/{input_id}_annotation_iter{j}.txt", taxid_dir, f"{input_id}_edgelist_iter{j}.txt", j, f"{taxid_dir}/{input_id}_norm.m8",
                  f"{taxid_dir}/{input_id}_NEWannotation", f"{taxid_dir}/{input_id}_barcharts", f"{taxid_dir}/{input_id}_nan.txt")

        if os.path.getsize(f"{taxid_dir}/{input_id}_nan.txt") == 0:
            print(f"All nan2nan proteins have been annotated. Stopping at iteration {j}.")
            break
        else:
            with open(f"{taxid_dir}/{input_id}_nan.txt", 'r') as file:
                line_count = sum(1 for line in file)

            if j < args.iters:
                print(f"Found {line_count} nan2nan IDs in total in the target proteome... proceeding to next iteration.")
            else:
                print(f"Found {line_count} nan2nan IDs in total in the target proteome.")

    print("\nWASP pipeline completed successfully!")


if __name__ == "__main__":
    main()
