#!/usr/bin/env python3
# @author Giorgia Del Missie

import json

def parse_m8(fin, fnormalisation, eval_thr, bits_thr, max_seq_id):
    """
    Parse Foldseek output files and normalize bitscores.

    Parameters:
    - fin: Path to the input file with Foldseek results
    - fnormalisation: Path to the input file with Foldseek results for normalization
    - eval_thr: E-value threshold
    - bits_thr: Bitscore threshold

    Returns:
    - Dictionary of queries and their significant hits
    """
    queries = dict()
    normalisation = dict()

    # Read normalization file
    with open(fnormalisation) as fnorm:
        for line in fnorm:
            line = line.strip().split()
            query, target, evalue, bitscore = line[0], line[1], float(line[-2]), int(line[-1])

            # Extract query and target IDs from headers
            query = query.split("-")[1] if "-" in query else query
            query = query.split(".gz")[0][:-4] if ".gz" in query else query
            target = target.split("-")[1] if "-" in target else target
            target = target.split(".gz")[0][:-4] if ".gz" in target else target
            query = query.split(".pdb")[0] if ".pdb" in query else query
            target = target.split(".pdb")[0] if ".pdb" in target else target

            if query == target:
                if query not in normalisation:
                    normalisation[query] = bitscore

    # Read input file
    with open(fin) as fs_in:
        for line in fs_in:
            line = line.strip().split()
            query, target, evalue, bitscore, fident = line[0], line[1], float(line[-2]), int(line[-1]), float(line[4])

            # Extract query and target IDs from headers
            query = query.split("-")[1] if "-" in query else query
            query = query.split(".gz")[0][:-4] if ".gz" in query else query
            target = target.split("-")[1] if "-" in target else target
            target = target.split(".gz")[0][:-4] if ".gz" in target else target
            query = query.split(".pdb")[0] if ".pdb" in query else query
            target = target.split(".pdb")[0] if ".pdb" in target else target

            if query != target:
                if bitscore > bits_thr and evalue < eval_thr and fident <= max_seq_id:
                    bitscore = round(bitscore / normalisation.get(query, 1), 3)
                    try:
                        if (target, evalue, bitscore) not in queries[query]:
                            queries[query].append((target, evalue, bitscore))
                    except KeyError:
                        queries[query] = [(target, evalue, bitscore)]

    return queries

def save_json(data, filename):
    """
    Save the given data to a JSON file.
    """
    with open(filename, "w") as f_out:
        json.dump(data, f_out)

