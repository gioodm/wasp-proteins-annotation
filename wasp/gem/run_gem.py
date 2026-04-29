#!/usr/bin/env python3
# @author Giorgia Del Missier

import argparse, os, sys, subprocess, glob, shutil
import numpy as np

from .rxn2uniprot import map_orphan_rxns
from .gap_filling import run_gap_filling

# Set random seed for reproducibility
np.random.seed(0)

def show_help():
    help_message = """
Usage: python run_GEM.py [-h] -t taxid -g gaps_file [-e evalue_threshold] [-b bitscore_threshold] [-tm tmscore]

This script uses foldseek structural alignment results to perform gap-filling in the Genome-scale Metabolic model of interest.

    -h, --help                  show this help message and exit
    -t, --taxid                  specify the GEM taxid (required)
    -g, --gaps_file              .tsv file containing the orphan reactions identified in the GEM (required)
    -e, --evalue_thr             set the evalue threshold (default: 10e-05)
    -b, --bitscore_thr           set the bitscore threshold (default: 50)
    -tm, --tmscore_thr           set the tmscore threshold (default: 0.5)

Examples:
    python3 run_gem.py -t 559292 -g 559292_gaps.txt
    python3 run_gem.py -t 559292 -g 559292_gaps.txt -e 1e-50 -b 1000 -tm 0.5
    python3 run_gem.py -t 559292 -g 559292_gaps.txt -tm 0.8
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

    parser = argparse.ArgumentParser(description='GEM gap-filling script', add_help=False)

    parser.add_argument("-t", "--taxid", type=str, required=True)
    parser.add_argument("-g", "--gaps_file", type=str, required=True)
    parser.add_argument("-e", "--eval_thr", type=float, default=1e-05)
    parser.add_argument("-b", "--bits_thr", type=int, default=50)
    parser.add_argument("-tm", "--tmscore_thr", type=float, default=0.5)

    args = parser.parse_args()
    
    system_tmp = os.environ.get("TMPDIR", "/tmp")
    fs_tmp = f"{system_tmp}/fs_tmp_{args.taxid}"
    os.makedirs(fs_tmp, exist_ok=True)
    os.chmod(fs_tmp, 0o755)

    # Check required commands
    check_command("foldseek")
    check_command("gsutil")

    print("Setting required variables:")
    print(f"\nSelected taxid is: {args.taxid}")
    print(f"File containing orphan reactions: {args.gaps_file}")
    print(f"Selected evalue threshold is: {args.eval_thr}")
    print(f"Selected bitscore threshold is: {args.bits_thr}")
    print(f"Selected tmscore threshold is: {args.tmscore_thr}")

    ####---- DOWNLOADING FILES AND DATABASES ----####

    # Define directory names
    db_dir = "foldseek_dbs"
    prot_dir = "proteomes"
    results_dir = "results"
    taxid_dir = f"{results_dir}/{args.taxid}"

    # Create directories
    for directory in [db_dir, prot_dir, results_dir, taxid_dir]:
        os.makedirs(directory, exist_ok=True)

    print("\nDownloading required AlphaFold models:")

    # Download and prepare databases
    if not os.path.exists(f"{db_dir}/swissprot"):
        subprocess.run(["foldseek", "databases", "Alphafold/Swiss-Prot", f"{db_dir}/swissprot", "tmp", "--remove-tmp-files", "1"])
    else:
        print("Foldseek database already downloaded")
    
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
    if not os.path.exists(f"{db_dir}/{args.taxid}"):
        subprocess.run(["foldseek", "createdb", f"{prot_dir}/{args.taxid}.tar", f"{db_dir}/{args.taxid}"])
    else:
        print("AlphaFold models of selected organism already downloaded")

    ####---- UNIPROT MAPPING ----####

    print("\nMapping orphan reaction to UniProt IDs:\n")

    map_orphan_rxns(args.gaps_file, f"{taxid_dir}/{args.taxid}_rxn2up.txt", f"{taxid_dir}/{args.taxid}_upIDs.txt")

    ####---- FOLDSEEK STRUCTURAL ALIGNMENT ----####

    print("\nPerforming foldseek alignment search:\n")

    subprocess.run(["foldseek", "prefixid", f"{db_dir}/swissprot_h", f"{db_dir}/swissprot.lookup", "--tsv", "--threads", "1"], check=True)
    awk_command = ["awk", "NR == FNR {f[$1] = $1; next} $2 in f {print $1}"] + [f"{taxid_dir}/{args.taxid}_upIDs.txt",f"{db_dir}/swissprot.lookup"]
    with open(f"{db_dir}/subset{args.taxid}.tsv", "w") as output_file:
        subprocess.run(awk_command, stdout=output_file, check=True)

    subprocess.run(["foldseek", "createsubdb", f"{db_dir}/subset{args.taxid}.tsv", f"{db_dir}/swissprot", f"{db_dir}/subdb{args.taxid}"], check=True)
    subprocess.run(["foldseek", "createsubdb", f"{db_dir}/subset{args.taxid}.tsv", f"{db_dir}/swissprot_ss", f"{db_dir}/subdb{args.taxid}_ss"], check=True)
    subprocess.run(["foldseek", "createsubdb", f"{db_dir}/subset{args.taxid}.tsv", f"{db_dir}/swissprot_ca", f"{db_dir}/subdb{args.taxid}_ca"], check=True)

    os.remove(f"{db_dir}/subset{args.taxid}.tsv")

    subprocess.run(["foldseek", "search", f"{db_dir}/subdb{args.taxid}", f"{db_dir}/{args.taxid}", f"{taxid_dir}/{args.taxid}", "tmp", "-a", "1", "--threads", "64"], check=True)
    subprocess.run(["foldseek", "convertalis", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits", 
                    f"{db_dir}/subdb{args.taxid}", f"{db_dir}/{args.taxid}", f"{taxid_dir}/{args.taxid}", f"{taxid_dir}/{args.taxid}.m8"], check=True)

    subprocess.run(["foldseek", "search", f"{db_dir}/subdb{args.taxid}", f"{db_dir}/subdb{args.taxid}", f"{taxid_dir}/{args.taxid}_db_allvsall", "tmp", "-a", "1", "--threads", "64"], check=True)
    subprocess.run(["foldseek", "convertalis", "--format-output", "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,alntmscore,evalue,bits", 
                    f"{db_dir}/subdb{args.taxid}", f"{db_dir}/subdb{args.taxid}", f"{taxid_dir}/{args.taxid}_db_allvsall", f"{taxid_dir}/{args.taxid}_db_allvsall.m8"], check=True)

    for file in glob.glob(os.path.join(db_dir, f"subdb{args.taxid}*")):
        try:
            os.remove(file)
        except IsADirectoryError:
            shutil.rmtree(file)

    for file in glob.glob(os.path.join(db_dir, f"{args.taxid}*")):
        try:
            os.remove(file)
        except IsADirectoryError:
            shutil.rmtree(file)

    for root, dirs, files in os.walk(taxid_dir):
        for fname in files:
            if not (fname.endswith(".txt") or fname.endswith(".m8")):
                file_path = os.path.join(root, fname)
                os.remove(file_path)

    ####---- GAP FILLING ----####

    print("Identifying best hits for each orphan reaction and performing gap-filling")
    run_gap_filling(f"{taxid_dir}/{args.taxid}.m8", f"{taxid_dir}/{args.taxid}_db_allvsall.m8", f"{taxid_dir}/{args.taxid}_rxn2up.txt", 
                    f"{taxid_dir}/{args.taxid}_hits.txt", args.eval_thr, args.bits_thr, args.tmscore_thr)

    print("\nGEM gap-filling complete!")

if __name__ == "__main__":
    main()