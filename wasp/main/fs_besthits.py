#!/usr/bin/env python3.10
# @author: Giorgia Del Missier

def get_besthits(input_file, n, evalue_threshold, bitscore_threshold, max_seq_id):
    """
    Extract the best hit for each query based on evalue and bitscore thresholds.
    
    Parameters:
    - input_file: Path to the input .m8 file
    - n: The nth hit to extract for each query
    - evalue_threshold: The maximum evalue for a hit to be considered
    - bitscore_threshold: The minimum bitscore for a hit to be considered
    
    Returns:
    - A set of the best hit targets
    """
    queries = dict()

    with open(input_file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split()
            query, target, evalue, bitscore, fident = line[0], line[1], float(line[-2]), int(line[-1]), float(line[4])
            
            # Strip potential suffix from query and target identifiers
            query_stripped = query.split("-")[1] if "-" in query else query
            target_stripped = target.split("-")[1] if "-" in target else target
            query_stripped = query.split(".pdb")[0] if ".pdb" in query else query
            target_stripped = target.split(".pdb")[0] if ".pdb" in target else target
            query_stripped = query.split(".cif")[0] if ".cif" in query else query
            target_stripped = target.split(".cif")[0] if ".cif" in target else target
            
            # Ensure query and target are not the same and apply thresholds
            if query_stripped != target_stripped:
                if evalue <= evalue_threshold and bitscore >= bitscore_threshold and fident <= max_seq_id: 
                    # Store targets for each query
                    if query not in queries:
                        queries[query] = []
                    if target not in queries[query]:
                        queries[query].append(target)
    
    # Extract the nth hit for each query if it exists
    hits = [targets[n-1] for query, targets in queries.items() if len(targets) >= n]

    return set(hits)

def save_besthits(best_hits, output_file):
    """
    Save the best hits to a specified output file.
    """
    with open(output_file, "w") as f_out:
        for hit in best_hits:
            f_out.write(hit + "\n")
