#%%

import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.preprocessing import MultiLabelBinarizer
from tabulate import tabulate


def strip_last_digit(ec_string):
    """Converts 1.2.3.4 to 1.2.3 and 1.2.3.- to 1.2.3"""
    if not ec_string or pd.isna(ec_string): return ""
    labels = set(str(ec_string).split(';'))
    stripped = {'.'.join(l.split('.')[:-1]) for l in labels if len(l.split('.')) > 1}
    return ';'.join(filter(None, stripped))

def calculate_metrics(true_df, pred_df):
    """Aligns predictions to true labels and calculates weighted metrics."""
    # Align pred_df to the UniProt IDs present in the ground truth
    aligned_pred = pd.merge(true_df[['UniProt ID']], pred_df, on='UniProt ID', how='left').fillna('')
    
    mlb = MultiLabelBinarizer()
    parse = lambda x: set(str(x).split(';')) if (isinstance(x, str) and x != '') else set()
    
    # Fit MLB on all possible labels present in both sets
    all_labels = pd.concat([true_df['EC number'], aligned_pred['EC number']]).apply(parse)
    mlb.fit(all_labels)
    
    y_true = mlb.transform(true_df['EC number'].apply(parse))
    y_pred = mlb.transform(aligned_pred['EC number'].apply(parse))
    
    return [
        precision_score(y_true, y_pred, average='weighted', zero_division=0),
        recall_score(y_true, y_pred, average='weighted', zero_division=0),
        f1_score(y_true, y_pred, average='weighted', zero_division=0)
    ]


results_table = []

datasets = [
    {
        "name": "New-392", 
        "gt": "new.csv", 
        "clean": "new_maxsep.csv", 
        "wasp": "new_annotated_taxid.xlsx", 
        "deepec": "new_DeepEC_Result.txt",
        "deepec_3d": "new_3digit_EC_prediction.txt",
        "eggnog": "new_eggnog.xlsx",
        "bh_saved": "new_bh_transfer_annotated.tsv" 
    },
    {
        "name": "Price-149", 
        "gt": "price.csv", 
        "clean": "price_maxsep.csv", 
        "wasp": "price_annotated_taxid.xlsx", 
        "deepec": "price_DeepEC_Result.txt",
        "deepec_3d": "price_3digit_EC_prediction.txt",
        "eggnog": "price_eggnog.xlsx",
        "mapping": "price_mapped.tsv",
        "bh_saved": "price_bh_transfer_annotated.tsv"
    }
]

for ds in datasets:
    print(f"\n--- Processing {ds['name']} ---")
    
    # 0. Load Mapping if it exists
    mapping_df = None
    if "mapping" in ds:
        mapping_df = pd.read_csv(ds["mapping"], sep="\t").drop_duplicates(['Entry'])

    def map_ids(df, id_column='UniProt ID'):
        """Helper to convert Entry names to UniProt IDs using the mapping file."""
        if mapping_df is None or df.empty: return df
        # Rename the provided ID column to 'Entry' to match mapping file
        df = df.rename(columns={id_column: 'Entry'})
        # Merge and keep only the real UniProt ID
        df = df.merge(mapping_df, on='Entry', how='left')
        return df[['UniProt ID', 'EC number']].dropna(subset=['UniProt ID'])

    # 1. Load Ground Truth
    gt = pd.read_csv(ds['gt'], sep="\t")
    if 'Entry' in gt.columns and 'UniProt ID' not in gt.columns:
        gt = gt.rename(columns={'Entry': 'UniProt ID'})
    
    # Map GT if needed (for Price-149)
    if mapping_df is not None:
        gt = map_ids(gt, id_column='UniProt ID')

    gt = gt.groupby('UniProt ID').agg({'EC number': ';'.join}).reset_index()
    print(f"Ground Truth loaded: {len(gt)} samples.")

    # 2. Load CLEAN & map
    clean_dict = {}
    try:
        with open(ds['clean']) as f:
            for line in f:
                p = line.strip().split(",")
                clean_dict[p[0]] = ";".join([x.split("/")[0][3:] for x in p[1:]])
        clean_df = pd.DataFrame(list(clean_dict.items()), columns=['UniProt ID', 'EC number'])
        clean_df = map_ids(clean_df)
    except: clean_df = pd.DataFrame(columns=['UniProt ID', 'EC number'])

    # 3. Load WASP & map
    try:
        xl = pd.ExcelFile(ds['wasp'])
        wasp_df = pd.concat([xl.parse(s) for s in xl.sheet_names])[['UniProt ID', 'EC number']].fillna('')
        wasp_df['EC number'] = wasp_df['EC number'].apply(lambda x: ';'.join([t.strip('()').split(', ')[0] for t in str(x).split(';') if '-' not in t and t]))
        wasp_df = wasp_df.groupby('UniProt ID').agg({'EC number': ';'.join}).reset_index()
    except: wasp_df = pd.DataFrame(columns=['UniProt ID', 'EC number'])

    # 4. Load DeepEC (4-Digit & 3-Digit)
    try:
        # Load 4-digit
        d_df = pd.read_csv(ds['deepec'], sep="\t")
        id_col = 'Query ID' if 'Query ID' in d_df.columns else 'query'
        pred_col = 'Predicted EC number' if 'Predicted EC number' in d_df.columns else 'EC number'
        deepec_df = d_df.rename(columns={id_col: 'UniProt ID', pred_col: 'EC number'})
        deepec_df['EC number'] = deepec_df['EC number'].astype(str).str.replace('EC:', '').replace('nan', '')
        deepec_df = map_ids(deepec_df.groupby('UniProt ID').agg({'EC number': ';'.join}).reset_index())

        # Load 3-digit specific file
        d3_raw = pd.read_csv(ds['deepec_3d'], sep="\t")
        id_col3 = 'Query ID' if 'Query ID' in d3_raw.columns else 'query'
        pred_col3 = 'Predicted EC number' if 'Predicted EC number' in d3_raw.columns else 'EC number'
        deepec_3d_df = d3_raw.rename(columns={id_col3: 'UniProt ID', pred_col3: 'EC number'})
        # Filter out "not predicted" entries
        deepec_3d_df = deepec_3d_df[~deepec_3d_df['EC number'].str.contains('not predicted', na=False, case=False)]
        deepec_3d_df['EC number'] = deepec_3d_df['EC number'].astype(str).str.replace('EC:', '').replace('nan', '')
        deepec_3d_df = map_ids(deepec_3d_df.groupby('UniProt ID').agg({'EC number': ';'.join}).reset_index())
    except: 
        deepec_df = pd.DataFrame(columns=['UniProt ID', 'EC number'])
        deepec_3d_df = pd.DataFrame(columns=['UniProt ID', 'EC number'])

    # 5. Load EggNOG & map
    try:
        f_egg = pd.ExcelFile(ds['eggnog'])
        eggnog_df = f_egg.parse(f_egg.sheet_names[-1], header=2).iloc[:-3]
        eggnog_df = eggnog_df[['query', 'EC']].rename(columns={'query': 'UniProt ID', 'EC': 'EC number'})
        eggnog_df = eggnog_df[eggnog_df['EC number'] != '-'].groupby('UniProt ID').agg({'EC number': ';'.join}).reset_index()
        eggnog_df = map_ids(eggnog_df)
    except: eggnog_df = pd.DataFrame(columns=['UniProt ID', 'EC number'])

    # 6. Load BH Transfer
    try:
        bh_df = pd.read_csv(ds['bh_saved'], sep="\t")[['UniProt ID', 'EC number']].fillna('')
    except: bh_df = pd.DataFrame(columns=['UniProt ID', 'EC number'])

    # 7. Create Hybrid (CLEAN + WASP)
    hybrid_df = pd.merge(clean_df, wasp_df, on='UniProt ID', how='outer', suffixes=('_c', '_w')).fillna('')
    hybrid_df['EC number'] = hybrid_df.apply(lambda r: ';'.join(set(filter(None, str(r['EC number_c']).split(';') + str(r['EC number_w']).split(';')))), axis=1)

    # -- Evaluation --
    models = {
        "WASP": wasp_df, "CLEAN": clean_df, "DeepEC": deepec_df, 
        "EggNOG": eggnog_df, "Hybrid (C+W)": hybrid_df, "BH Transfer": bh_df
    }

    for level in ["4th Digit", "3rd Digit"]:
        target_gt = gt.copy()
        if level == "3rd Digit":
            target_gt['EC number'] = target_gt['EC number'].apply(strip_last_digit)

        for name, pred in models.items():
            current_pred = pred.copy()
            
            if level == "3rd Digit" and name == "DeepEC" and not deepec_3d_df.empty:
                current_pred = deepec_3d_df.copy()
            elif level == "3rd Digit":
                current_pred['EC number'] = current_pred['EC number'].apply(strip_last_digit)
            
            stats = calculate_metrics(target_gt, current_pred)
            results_table.append([ds['name'], level, name] + [f"{s:.3f}" for s in stats])


headers = ["Dataset", "Level", "Model", "Precision", "Recall", "F1-Score"]
print("\n" + "="*80)
print("FINAL EC PREDICTION BENCHMARK")
print("="*80)
print(tabulate(results_table, headers=headers, tablefmt="grid"))
# %%
