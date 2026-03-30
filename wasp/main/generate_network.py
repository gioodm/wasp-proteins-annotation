#!/usr/bin/env python3.10
# @author: Giorgia Del Missier

import json
import numpy as np
import networkx as nx

# Set random seed for reproducibility
np.random.seed(0)

def create_network(all_queries, reciprocal_queries, max_hits):
    """
    Create a network graph containing clusters of proteins with similar structure.

    Parameters:
    - all_queries: Dictionary of queries and their hits from the input file
    - reciprocal_queries: Dictionary of reciprocal hits from the input_bh file
    - max_hits: Maximum number of neighbours to include for each query

    Returns:
    - A networkx Graph object representing the protein clusters
    """
    G = nx.Graph()

    # Add edges for each query and its top n neighbours
    for query, hits in all_queries.items():
        if hits[0][0] in reciprocal_queries:
            reci_hits = {e[0]: e[2] for e in reciprocal_queries[hits[0][0]][:max_hits]}
            if query in reci_hits:
                G.add_edge(query, hits[0][0], weight=hits[0][2])

                i = 1
                while i < len(hits) and i < max_hits:
                    if hits[i][0] in reci_hits:
                        G.add_edge(query, hits[i][0], weight=hits[i][2])
                        G.add_edge(hits[0][0], hits[i][0], weight=reci_hits[hits[i][0]])
                    i += 1

    return G

def run_network_generation(input_file, input_bh_file, nan_file, max_hits):
    """
    Main function to run the network generation process.
    """
    # Load data from input files
    with open(input_file) as fall, open(input_bh_file) as freci:
        all_queries = json.load(fall)
        reciprocal_queries = json.load(freci)

    # Load IDs of uncharacterized proteins
    with open(nan_file) as fnan:
        nan2nan = [line.strip() for line in fnan]

    # Filter queries to include only uncharacterized proteins
    if nan2nan:
        all_queries = {key: value for key, value in all_queries.items() if key in nan2nan}

    # Create the network graph
    G = create_network(all_queries, reciprocal_queries, max_hits)

    # Find the difference between all queries and the graph nodes
    diff = set(all_queries.keys()).difference(set(G.nodes()))

    # Sort the clusters by size (number of nodes)
    clusters_sorted = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]

    return G, all_queries, diff, clusters_sorted

def save_network(G, diff, clusters_sorted, output, edgelist, out_nan):
    """
    Save the network clusters and edge list to output files.
    """
    # Write nan2nan IDs to the output file
    with open(out_nan, "w") as foutnan:
        for i in diff:
            foutnan.write(i + "\n")

    # Sort and save clusters and their metrics
    columns = ["#Cluster", "Number of Nodes", "Number of Edges", "Average Degree", "Average Clustering Coefficient", "Node List"]
    with open(output, "w") as output_file, open(edgelist, "w") as output_edgelist:
        output_file.write("\t".join(columns) + "\n")

        for idx, cluster in enumerate(clusters_sorted):
            S = G.subgraph(cluster)
            avg_degree = round(sum(dict(S.degree(weight='weight')).values()) / len(S), 3)
            avg_clustering_coef = round(nx.average_clustering(S, weight='weight'), 3)

            for edge in S.edges.data("weight"):
                output_edgelist.write("\t".join(str(item) for item in edge) + "\n")

            fields = [idx, len(S.nodes()), len(S.edges()), avg_degree, avg_clustering_coef, str(cluster).strip("{ }")]
            output_file.write("\t".join(str(item) for item in fields) + "\n")
