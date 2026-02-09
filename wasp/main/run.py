#!/usr/bin/env python3
# @author Giorgia Del Missier

import argparse, os, sys, subprocess, glob, shutil
import numpy as np

from .fs_besthits import get_besthits, save_besthits
from .parse_fs import parse_m8, save_json
from .generate_network import run_network_generation, save_network
from .retrieve_annotation import fetch_annotations
from .SAFE_enrichment import run_safe_analysis

# Set random seed for reproducibility
np.random.seed(0)

def show_help():
    help_message = """
Usage: python3 run.py [-h] -t taxid [-e evalue_threshold] [-b bitscore_threshold] [-n max_neighbours] [-s step] [-i iterations]

WASP (Whole-proteome Annotation through Structural homology Pipeline) performs a "structural BLAST" using AlphaFold models to better annotate the target taxid proteome.
Parameters:

    -h, --help                  show this help message and exit
    -t, --taxid                  NCBI taxonomy identifier to be analysed (required)
    -e, --evalue_thr             set the evalue threshold (default: 10e-10)
    -b, --bitscore_thr           set the bitscore threshold (default: 50)
    -n, --max_n                  set the max number of neighbours (default: 10)
    -s, --step                   set step to add to max neighbours (n) in additional iterations (default: 10)
    -i, --iters                  set number of iterations to perform (default: 3)

Examples:
    python3 run.py -t 559292
    python3 run.py -t 559292 -e 1e-50 -b 200 -n 5 -i 5
    python3 run.py -t 559292 -s 5
    """
    print(help_message)
 
def check_command(command):
    try:
        subprocess.run([command], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        print(f"{command} is required but it's not installed. Aborting.")
        sys.exit(1)

def main():

    parser = argparse.ArgumentParser(description="WASP Pipeline", add_help=False)
    parser.add_argument("-h", "--help", action="store_true")
    parser.add_argument("-t", "--taxid", required=True)
    parser.add_argument("-e", "--eval_thr", type=float, default=1e-10)
    parser.add_argument("-b", "--bits_thr", type=int, default=50)
    parser.add_argument("-n", "--max_n", type=int, default=10)
    parser.add_argument("-s", "--step", type=int, default=10)
    parser.add_argument("-i", "--iters", type=int, default=3)

    args = parser.parse_args()

    if args.help:
        show_help()
        sys.exit(0)

    # Set environment variable for foldseek
    os.environ["PATH"] = f"{os.getcwd()}/foldseek/bin/:{os.environ['PATH']}"

    # Check required commands
    check_command("foldseek")
    check_command("gsutil")

    print("Setting required variables:")
    print(f"\nSelected taxid is: {args.taxid}")
    print(f"Selected evalue threshold is: {args.eval_thr}")
    print(f"Selected bitscore threshold is: {args.bits_thr}")
    print(f"Selected max neighbours is: {args.max_n}")
    print(f"Selected step is: {args.step}")
    print(f"Selected number of iterations is: {args.iters}")

    ####---- DOWNLOADING FILES AND DATABASES ----####

    # Define directory names
    db_dir = "foldseek_dbs"
    prot_dir = "proteomes"
    results_dir = "results"
    taxid_dir = f"{results_dir}/{args.taxid}"

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
    if all(os.path.exists(f"{taxid_dir}/{args.taxid}{suffix}.m8") for suffix in ["", "_bh", "_norm", "_norm_bh"]):
        print("AlphaFold models of selected organism already downloaded and Foldseek results already generated")
    else:
        # Download and prepare proteome
        if not os.path.exists(f"{prot_dir}/{args.taxid}.tar"):
            os.makedirs(f"{prot_dir}/{args.taxid}", exist_ok=True)
            subprocess.run(["gsutil", "-m", "cp", f"gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-{args.taxid}-*_v4.tar", prot_dir])
            
            for f in os.listdir(prot_dir):
                if f.startswith(f"proteome-tax_id-{args.taxid}") and f.endswith("_v4.tar"):
                    subprocess.run(["tar", "-xf", f"{prot_dir}/{f}", "-C", f"{prot_dir}/{args.taxid}"])
                    os.remove(f"{prot_dir}/{f}")
            
            for f in os.listdir(f"{prot_dir}/{args.taxid}"):
                if f.endswith(".json.gz"):
                    os.remove(f"{prot_dir}/{args.taxid}/{f}")
            
            subprocess.run(["tar", "-cf", f"{prot_dir}/{args.taxid}.tar", "-C", prot_dir, args.taxid])
            subprocess.run(["rm", "-r", f"{prot_dir}/{args.taxid}"])
        if not os.path.exists(f"{db_dir}/afdb50sp{args.taxid}"):
            subprocess.run(["foldseek", "createdb", f"{prot_dir}/{args.taxid}.tar", f"{db_dir}/{args.taxid}"])

            for suffix in ["", "_h", "_ss", "_ca"]:
                subprocess.run(["foldseek", "concatdbs", f"{db_dir}/afdb50sp{suffix}", f"{db_dir}/{args.taxid}{suffix}", f"{db_dir}/afdb50sp{args.taxid}{suffix}"])
        else:
            print("AlphaFold models of selected organism already downloaded")

    ####---- RECIPROCAL BEST STRUCTURE HITS SEARCH ----####

    print("\nPerforming Reciprocal Best Structural Hits search:")

    if not os.path.exists(f"{taxid_dir}/{args.taxid}.m8") or not os.path.exists(f"{taxid_dir}/{args.taxid}_bh.m8"):
        # Perform foldseek searches
        subprocess.run(["foldseek", "easy-search", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits",
                        f"{prot_dir}/{args.taxid}.tar", f"{db_dir}/afdb50sp{args.taxid}", f"{taxid_dir}/{args.taxid}.m8", "tmp", "--threads", "64"])

        subprocess.run(["foldseek", "easy-search", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits",
                        f"{prot_dir}/{args.taxid}.tar", f"{db_dir}/{args.taxid}", f"{taxid_dir}/{args.taxid}_norm.m8", "tmp", "--threads", "64",
                        "--exhaustive-search", "1", "--min-seq-id", "0.9"])

        best_hits = get_besthits(f"{taxid_dir}/{args.taxid}.m8", 1, str(args.eval_thr), str(args.bits_thr))
        save_besthits(best_hits, f"{taxid_dir}/{args.taxid}_bh.txt")

        subprocess.run(["foldseek", "prefixid", f"{db_dir}/afdb50sp{args.taxid}_h", f"{db_dir}/afdb50sp{args.taxid}.lookup", "--tsv", "--threads", "1"])
        
        awk_command = ["awk", "NR == FNR {f[$1] = $1; next} $2 in f {print $1}", f"{taxid_dir}/{args.taxid}_bh.txt", f"{db_dir}/afdb50sp{args.taxid}.lookup"]
        with open(f"{db_dir}/subset{args.taxid}.tsv", "w") as output_file:
            subprocess.run(awk_command, stdout=output_file, check=True)

        for suffix in ["", "_ss", "_ca"]:
            subprocess.run(["foldseek", "createsubdb", f"{db_dir}/subset{args.taxid}.tsv", f"{db_dir}/afdb50sp{args.taxid}{suffix}", f"{db_dir}/subdb{args.taxid}{suffix}"])
        
        os.remove(f"{db_dir}/subset{args.taxid}.tsv")

        subprocess.run(["foldseek", "search", f"{db_dir}/subdb{args.taxid}", f"{db_dir}/afdb50sp{args.taxid}", f"{taxid_dir}/{args.taxid}_bh", "tmp", "-a", "1", "--threads", "64"])
        subprocess.run(["foldseek", "convertalis", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits",
                        f"{db_dir}/subdb{args.taxid}", f"{db_dir}/afdb50sp{args.taxid}", f"{taxid_dir}/{args.taxid}_bh", f"{taxid_dir}/{args.taxid}_bh.m8"])

        subprocess.run(["foldseek", "search", f"{db_dir}/subdb{args.taxid}", f"{db_dir}/subdb{args.taxid}", f"{taxid_dir}/{args.taxid}_norm_bh", "tmp",
                        "-a", "1", "--threads", "64", "--exhaustive-search", "1", "--min-seq-id", "0.9"])
        subprocess.run(["foldseek", "convertalis", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits",
                        f"{db_dir}/subdb{args.taxid}", f"{db_dir}/subdb{args.taxid}", f"{taxid_dir}/{args.taxid}_norm_bh", f"{taxid_dir}/{args.taxid}_norm_bh.m8"])

        # Clean up temporary files
        for path in glob.glob(os.path.join(db_dir, f"subdb{args.taxid}*")):
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)

        for path in glob.glob(os.path.join(db_dir, f"afdb50sp{args.taxid}*")):
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)

        for path in glob.glob(os.path.join(db_dir, f"{args.taxid}*")):
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

    # Count the number of proteins in the proteome
    psize = subprocess.run(["tar", "-tvf", f"{prot_dir}/{args.taxid}.tar"], capture_output=True, text=True)
    psize = len([line for line in psize.stdout.splitlines() if line.endswith(".gz")])

    queries = parse_m8(f"{taxid_dir}/{args.taxid}.m8", f"{taxid_dir}/{args.taxid}_norm.m8", args.eval_thr, args.bits_thr)
    reciprocal_queries = parse_m8(f"{taxid_dir}/{args.taxid}_bh.m8", f"{taxid_dir}/{args.taxid}_norm_bh.m8", args.eval_thr, args.bits_thr)
    
    # Saving the results to JSON files
    save_json(queries, f"{taxid_dir}/{args.taxid}.json")
    save_json(reciprocal_queries, f"{taxid_dir}/{args.taxid}_norm.json")

    print(f"Found significant hits for {len(queries)} out of {psize} proteins in the target organism")

    open(f"{taxid_dir}/{args.taxid}_nan.txt", "w").close()

    for j in range(1, args.iters + 1):
        print(f"\nPerforming iteration {j}")

        ####---- NETWORK GENERATION ----####

        print("\nCreating RBSH network and identifying clusters of homologs:")

        neighbours = args.max_n + (args.step * (j - 1))
        print(f"Selected max number of neighbours for iteration {j} is: {neighbours}\n")

        # Run the network generation
        G, all_queries, diff, clusters_sorted = run_network_generation(f"{taxid_dir}/{args.taxid}.json", f"{taxid_dir}/{args.taxid}_norm.json", 
                                                                       f"{taxid_dir}/{args.taxid}_nan.txt", neighbours)
        # Save the network and cluster details
        save_network(G, all_queries, diff, clusters_sorted, f"{taxid_dir}/{args.taxid}_clusters_iter{j}.txt", f"{taxid_dir}/{args.taxid}_edgelist_iter{j}.txt")

        print(f"{len(diff)} proteins in the target organism had no RBSH hits... trying again with increased number of neighbours")

        # Print network statistics
        print(f"Network statistics generated using {len(all_queries) - len(diff)} RBSH hits:")
        print(f"Number of nodes: {len(G.nodes())}")
        print(f"Number of edges: {len(G.edges())}")
        print(f"Number of generated clusters: {nx.number_connected_components(G)}")
        
        ####---- ANNOTATION ----####

        print("\nAnnotating network's nodes")

        fetch_annotations(f"{taxid_dir}/{args.taxid}_clusters_iter{j}.txt", f"{taxid_dir}/{args.taxid}_annotation_iter{j}.txt")

        ####---- SAFE ENRICHMENT AND STATISTICS COMPUTATION ----####

        print("Performing SAFE analysis and computing statistics on new annotation\n")

        run_safe_analysis(f"{taxid_dir}/{args.taxid}_annotation_iter{j}.txt", taxid_dir, f"{args.taxid}_edgelist_iter{j}.txt", j, f"{taxid_dir}/{args.taxid}_norm.m8", 
                          f"{taxid_dir}/{args.taxid}_NEWannotation", f"{taxid_dir}/{args.taxid}_barcharts", f"{taxid_dir}/{args.taxid}_nan.txt")

        if os.path.getsize(f"{taxid_dir}/{args.taxid}_nan.txt") == 0:
            print(f"All nan2nan proteins have been annotated. Stopping at iteration {j}.")
            break
        else:
            with open(f"{taxid_dir}/{args.taxid}_nan.txt", 'r') as file:
                line_count = sum(1 for line in file)

            if j < args.iter:
                print(f"Found {line_count} nan2nan IDs in total in the target proteome... proceeding to next iteration.")
            else:
                print(f"Found {line_count} nan2nan IDs in total in the target proteome.")

    print("\nWASP pipeline completed successfully!")


if __name__ == "__main__":
    main()
