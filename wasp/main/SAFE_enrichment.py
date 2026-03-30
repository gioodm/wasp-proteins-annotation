#!/usr/bin/env python3.10
# @author Giorgia Del Missier

import os
import sys
import json
import multiprocessing
import pandas as pd
import numpy as np
import altair as alt
from altair_saver import save
from pathlib import Path

# Get the directory where safe.py is located
#safepy_dir = str(Path(__file__).parent.resolve())
#sys.path.append(safepy_dir)

from wasp.main.safepy import safe

# Set random seed for reproducibility
np.random.seed(0)

def process_identifier(args_list):
    """Process protein identifiers for SAFE analysis.
    
    Args:
        args_list (tuple): A tuple containing the identifier, the dataframe, and the working directory.
    """
    i, df, wd = args_list
    subdf = df[['#Cluster', 'UniProt ID', i]]
    flattened_ids = subdf[i].dropna().astype(str).str.split(';').explode()
    final_matrix = pd.DataFrame()
    
    if not flattened_ids.empty:
        unique_ids = flattened_ids.dropna().unique()
        final_matrix = pd.DataFrame(columns=unique_ids)

        for index, row in subdf.iterrows():
            protein = row['UniProt ID']
            val = row[i]
            
            # Safely handle the value, ensuring it's treated as a string and split into an iterable list
            if pd.isna(val):
                ids = []
            else:
                ids = str(val).split(';')
                
            final_matrix.loc[protein] = [1 if j in ids else 0 for j in unique_ids]

    output_file = f"{wd}/SAFE/{i}_matrix.txt"
    final_matrix.to_csv(output_file, sep='\t')

def perform_safe_analysis(i, wd, radius, sf, enrichment_thr=0.05):
    """Perform SAFE analysis on the given identifier.
    
    Args:
        i (str): Identifier for the analysis.
        wd (str): Working directory for output files.
        radius (float): Neighborhood radius for SAFE analysis.
        sf (safe.Safe): SAFE analysis instance.
        enrichment_thr (float): Threshold for enrichment significance.
        
    Returns:
        tuple: The identifier and its significant results.
    """
    sf.load_attributes(attribute_file=f'SAFE/{i}_matrix.txt')
    sf.define_neighborhoods(neighborhood_radius=float(radius))
    sf.compute_pvalues(background='network', num_permutations=1000, processes=64)
    sf.print_output_files(output_dir=f'{wd}/SAFE/{i}_')

    df_tmp = pd.read_csv(f"{wd}/SAFE/{i}_node_properties_annotation.txt", sep="\t")
    
    significant = {}
    for index, row in df_tmp.iterrows():
        key = row['key']
        significant[key] = [
            (col, round(row[col], 3)) 
            for col in df_tmp.columns[3:] 
            if row[col] > -np.log10(enrichment_thr) and not col.startswith("Unnamed")
        ]
    
    return i, significant

def prepare_df(d):
    """Prepare a DataFrame from the provided data dictionary.
    
    Args:
        d (dict): A dictionary where keys are identifiers and values are dataframes.
        
    Returns:
        pd.DataFrame: A stacked DataFrame suitable for plotting.
    """
    df = pd.DataFrame()

    for i in d:
        tmp_df = pd.DataFrame(d[i]).transpose()
        tmp_df = tmp_df.stack().reset_index()
        tmp_df.columns = ['id', 'group', 'count']
        tmp_df['category'] = i

        df = pd.concat([df, tmp_df])

    return df

def make_barcharts(txid, coll, outfig, iteration):

    df_taxid = prepare_df(txid)
    df_collateral = prepare_df(coll)

    chart_taxid = alt.Chart(df_taxid).mark_bar().encode(

        # tell Altair which field to group columns on
        x=alt.X('group:N', title=None, sort=None),  # Add sort argument

        # tell Altair which field to use as Y values and how to calculate
        y=alt.Y('sum(count):Q',
            axis=alt.Axis(
                grid=False,
                title=None)),

        # tell Altair which field to use to use as the set of columns to be represented in each group
        column=alt.Column('id:N', sort=alt.SortField(field='id', order='ascending'), title=None),  # Add sort argument

        # tell Altair which field to use for color segmentation 
        color=alt.Color('category:N', sort=None,
            scale=alt.Scale(
                # make it look pretty with an enjoyable color pallet
                range=['#001219','#005F73','#94D2BD','#EE9B00','#BB3E03','#9B2226'],
            )),

        order=alt.Order(
        # Sort the segments of the bars by this field
        'category',
        sort='ascending')
            )\
            .configure_view(
            # remove grid lines around column clusters
            strokeOpacity=0    
    ).properties(
    width=150)

    chart_taxid = chart_taxid.configure_legend(labelFontSize=8)
    chart_taxid.save(f"{outfig}_taxid_iter{iteration}.png", ppi=200)

    chart_collateral = alt.Chart(df_collateral).mark_bar().encode(

        # tell Altair which field to group columns on
        x=alt.X('group:N', title=None, sort=None),  # Add sort argument

        # tell Altair which field to use as Y values and how to calculate
        y=alt.Y('sum(count):Q',
            axis=alt.Axis(
                grid=False,
                title=None)),

        # tell Altair which field to use to use as the set of columns to be represented in each group
        column=alt.Column('id:N', sort=alt.SortField(field='id', order='ascending'), title=None),  # Add sort argument

        # tell Altair which field to use for color segmentation 
        color=alt.Color('category:N', sort=None,
            scale=alt.Scale(
                # make it look pretty with an enjoyable color pallet
                range=['#001219','#005F73','#94D2BD','#EE9B00','#BB3E03','#9B2226'],
            )),

        order=alt.Order(
        # Sort the segments of the bars by this field
        'category',
        sort='ascending')
            )\
            .configure_view(
            # remove grid lines around column clusters
            strokeOpacity=0    
    ).properties(
    width=150)

    chart_collateral = chart_collateral.configure_legend(labelFontSize=8)
    chart_collateral.save(f"{outfig}_collateral_iter{iteration}.png", ppi=200)

def analyze_safe_results(df_anno, df_safe, ids, allp, output, wd, iteration, nan):

    taxid, collateral = {}, {}

    for key in ["total", "previously annotated", "not annotated", "nan2nan", "new annotation by WASP", "nan2nan in dark clusters"]:
        taxid[key], collateral[key] = {}, {}
        for identifier in ids:
            taxid[key][f'{ids.index(identifier) + 1}.{identifier}'] = {'total': 0, 'annotation': 0, 'nan2nan / new': 0, 'nan2nan (dark)': 0}
            collateral[key][f'{ids.index(identifier) + 1}.{identifier}'] = {'total': 0, 'annotation': 0, 'nan2nan / new': 0, 'nan2nan (dark)': 0}

    list_nan2nan = list()
    cols = ["UniProt ID"] + ids + ["Organism"]
    taxid_outdf, collateral_outdf = pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

    # Convert to sets/lists outside the loop for speed
    allp_ids = set(allp['UniProt ID'].tolist())
    
    nclusters = df_anno["#Cluster"].unique()
    for c in nclusters:
        
        subdf_anno = df_anno[df_anno["#Cluster"] == c]
        subdf_safe = df_safe[df_safe["UniProt ID"].isin(subdf_anno["UniProt ID"])]

        dark_anno = subdf_anno[ids].apply(lambda col: col.isna().all()).to_dict()
        dark_safe = subdf_safe[ids].apply(lambda col: col.isna().all()).to_dict()
        
        tmp_taxid, tmp_collateral = pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

        for i in ids:

            taxid['total'][f'{ids.index(i) + 1}.{i}']['total'] += len(subdf_anno[subdf_anno['UniProt ID'].isin(allp_ids)])
            collateral['total'][f'{ids.index(i) + 1}.{i}']['total'] += len(subdf_anno[~subdf_anno['UniProt ID'].isin(allp_ids)])

            if dark_anno[i] and dark_safe[i]:
                proteins = subdf_anno[["UniProt ID", "Organism"]]
                proteins_taxid = proteins[proteins['UniProt ID'].isin(allp_ids)]
                list_nan2nan.extend(proteins_taxid["UniProt ID"].tolist())

                taxid['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(proteins_taxid)
                taxid['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += len(proteins_taxid)
                taxid['nan2nan in dark clusters'][f'{ids.index(i) + 1}.{i}']['nan2nan (dark)'] += len(proteins_taxid)

                collateral['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(proteins) - len(proteins_taxid)
                collateral['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += len(proteins) - len(proteins_taxid)
                collateral['nan2nan in dark clusters'][f'{ids.index(i) + 1}.{i}']['nan2nan (dark)'] += len(proteins) - len(proteins_taxid)

            else:
                # Restored Vectorized Logic
                pre = subdf_anno[['UniProt ID', i, "Organism"]].reset_index(drop=True)
                new = subdf_safe[['UniProt ID', i]].reset_index(drop=True)
                
                # Align rows safely
                new = new.sort_values(by='UniProt ID', key=lambda x: x.map(dict(zip(pre['UniProt ID'], pre.index)))).reset_index(drop=True)

                pre_taxid = pre[pre['UniProt ID'].isin(allp_ids)]

                # Compute new annotation by WASP
                changes = pre[i].isna() & new[i].notna()
                changed_rows = new[changes].copy()
                
                # Compute nan2nan
                changes_n2n = pre[i].isna() & new[i].isna()
                unchanged_rows_n2n = new[changes_n2n].copy()

                changed_rows = changed_rows.merge(pre[changes][['UniProt ID', 'Organism']], on='UniProt ID', how='left')
                changed_rows_taxid = changed_rows[changed_rows['UniProt ID'].isin(allp_ids)]
                changed_rows_collateral = changed_rows[~changed_rows['UniProt ID'].isin(allp_ids)]

                unchanged_rows_n2n = unchanged_rows_n2n.merge(pre[changes_n2n][['UniProt ID', 'Organism']], on='UniProt ID', how='left')
                unchanged_rows_taxid_n2n = unchanged_rows_n2n[unchanged_rows_n2n['UniProt ID'].isin(allp_ids)]
                unchanged_rows_collateral_n2n = unchanged_rows_n2n[~unchanged_rows_n2n['UniProt ID'].isin(allp_ids)]
                
                list_nan2nan.extend(unchanged_rows_taxid_n2n["UniProt ID"].tolist())

                # Update Dictionaries
                taxid['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += (len(changed_rows_taxid) + len(unchanged_rows_taxid_n2n))
                taxid['new annotation by WASP'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(changed_rows_taxid)
                taxid['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(unchanged_rows_taxid_n2n)
                taxid['previously annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += len(pre_taxid) - (len(changed_rows_taxid) + len(unchanged_rows_taxid_n2n))

                collateral['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += (len(changed_rows_collateral) + len(unchanged_rows_collateral_n2n))
                collateral['new annotation by WASP'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(changed_rows_collateral)
                collateral['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(unchanged_rows_collateral_n2n)
                collateral['previously annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += (len(pre) - len(pre_taxid)) - (len(changed_rows_collateral) + len(unchanged_rows_collateral_n2n))

                # Append changed rows for Excel output
                if not changed_rows_taxid.empty:
                    tmp_taxid = pd.concat([tmp_taxid, changed_rows_taxid], ignore_index=True, sort=False)
                if not changed_rows_collateral.empty:
                    tmp_collateral = pd.concat([tmp_collateral, changed_rows_collateral], ignore_index=True, sort=False)

        # Restored Groupby & Aggregation for clean Excel formatting
        tmp_taxid = tmp_taxid.fillna('')
        tmp_collateral = tmp_collateral.fillna('')

        if not tmp_taxid.empty:
            tmp_taxid = tmp_taxid.groupby('UniProt ID').agg({
                        'Pfam': ''.join, 'PANTHER': ''.join, 'CATH': ''.join,
                        'EC number': ''.join, 'Rhea ID': ''.join, 'GO terms': ''.join,
                        'Organism': 'first'
                        }).reset_index()

        if not tmp_collateral.empty:
            tmp_collateral = tmp_collateral.groupby('UniProt ID').agg({
                        'Pfam': ''.join, 'PANTHER': ''.join, 'CATH': ''.join,
                        'EC number': ''.join, 'Rhea ID': ''.join, 'GO terms': ''.join,
                        'Organism': 'first'
                        }).reset_index()

        taxid_outdf = pd.concat([taxid_outdf, tmp_taxid], axis=0, ignore_index=True)
        collateral_outdf = pd.concat([collateral_outdf, tmp_collateral], axis=0, ignore_index=True)

    output_name = os.path.basename(output)
    
    # Save the output files safely
    taxid_excel = f"{wd}/taxid_{output_name}.xlsx"
    collateral_excel = f"{wd}/collateral_{output_name}.xlsx"

    def save_to_excel(df, filename, sheet_name):
        try:
            with pd.ExcelWriter(filename, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
                df.to_excel(writer, sheet_name=sheet_name, index=False)
        except (FileNotFoundError, ValueError):
            with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                df.to_excel(writer, sheet_name=sheet_name, index=False)

    save_to_excel(taxid_outdf, taxid_excel, f'Iteration_{iteration}')
    save_to_excel(collateral_outdf, collateral_excel, f'Iteration_{iteration}')
    
    if len(list_nan2nan) > 0:
        unique_nan = set(list_nan2nan)
        with open(nan, "a") as f_out:
            for prot_id in unique_nan:
                f_out.write(f"{prot_id}\n")
        print(f"Iteration {iteration}: {len(unique_nan)} IDs appended to {nan}")

    return taxid, collateral

def run_safe_analysis(input_file, wd, edgelist, iteration, identifiers, output, outfig, nan_file):
    """Run the SAFE analysis pipeline with the specified arguments."""

    if not os.path.exists(f"{wd}/SAFE"):
        os.makedirs(f"{wd}/SAFE")

    df_anno = pd.read_csv(input_file, sep="\t", header=0)
        
    all_proteins = list()
    with open(identifiers) as fup:
        for line in fup:
            line = line.strip().split()
            query = line[0]
            query = (lambda query : query.split("-")[1] if "-" in query else query)(query)
            query = (lambda query : query.split(".gz")[0][:-4] if ".gz" in query else query)(query)
            query = (lambda query : query.split(".pdb")[0] if ".pdb" in query else query)(query)
            query = (lambda query : query.split(".cif")[0] if ".cif" in query else query)(query)
            all_proteins.append(query)

    all_proteins = list(set(all_proteins))  # Ensure unique entries
    all_proteins = pd.DataFrame({'UniProt ID': all_proteins})

    safe_file = safe.SAFE(path_to_safe_data=f'{wd}/')
    safe_file.load_network(network_file=edgelist)

    ids = ["Pfam", "PANTHER", "CATH", "EC number", "Rhea ID", "GO terms"]

    with multiprocessing.Pool(processes=64) as pool:
        pool.map(process_identifier, [(i, df_anno, wd) for i in ids])
        results_list = pool.starmap(perform_safe_analysis, [(i, wd, 0.5, safe_file) for i in ids])
    
    # Rebuild safe_df from the multiprocessing results (restored from your original script)
    safe_enriched = dict(results_list)
    safe_df = pd.DataFrame.from_dict(safe_enriched, orient='index').transpose()
    for i in ids:
        # Sort and format the SAFE scores safely
        safe_df[i] = safe_df[i].apply(lambda x: sorted(x, key=lambda item: item[1], reverse=True) if isinstance(x, list) else x)
        safe_df[i] = safe_df[i].apply(lambda x: ';'.join(f'({item[0]}, {item[1]})' for item in x) if isinstance(x, list) and len(x) > 0 else np.nan)
    
    safe_df = safe_df.reset_index().rename(columns={'index': 'UniProt ID'})

    # Pass safe_df as the second argument
    taxid, collateral = analyze_safe_results(df_anno, safe_df, ids, all_proteins, output, wd, iteration, nan_file)
    make_barcharts(taxid, collateral, outfig, iteration)
