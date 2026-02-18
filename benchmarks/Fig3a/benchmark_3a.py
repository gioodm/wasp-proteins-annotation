#%%
# 
import pandas as pd
import numpy as np
import re
from sklearn.metrics import f1_score, precision_score, recall_score
from sklearn.preprocessing import MultiLabelBinarizer

def get_metrics(true_labels_path, excel_path, ids):
    true_labels = pd.read_csv(true_labels_path)
    fsafe = pd.ExcelFile(excel_path)
    
    # Filter sheets for UniProt IDs present in truth set
    dfs = [fsafe.parse(s) for s in fsafe.sheet_names]
    dfs = [df[df['UniProt ID'].isin(true_labels['UniProt ID'])] for df in dfs]

    metrics = {lab: {'precision': [0.0], 'recall': [0.0], 'f1': [0.0]} for lab in ids}
    
    for i in range(1, len(dfs) + 1):
        for j in ids:
            # Aggregate and format predictions
            merged = pd.concat(dfs[:i])[["UniProt ID", j]].fillna('')
            merged = merged.groupby('UniProt ID').agg({j: '; '.join}).reset_index()
            
            merged[j] = merged[j].apply(lambda x: ';'.join(list(set([t.strip('()').split(', ')[0] for t in x.split('; ')]))))
            if j == 'GO terms':
                merged[j] = merged[j].apply(lambda x: ';'.join(['GO:' + g for g in re.findall(r'GO:(\d+)', x)]))

            # Align with true_labels (filling missing IDs with empty strings)
            merged = pd.merge(true_labels[['UniProt ID']], merged, on='UniProt ID', how='left').fillna('')

            # Multi-label Binarization
            mlb = MultiLabelBinarizer()
            # Combine sets to ensure identical column space for both true and pred
            y_true_sets = true_labels[j].apply(lambda x: set(x.split(';')) if x else set())
            y_pred_sets = merged[j].apply(lambda x: set(x.split(';')) if x else set())
            
            mlb.fit(pd.concat([y_true_sets, y_pred_sets]))
            y_true_bin = mlb.transform(y_true_sets)
            y_pred_bin = mlb.transform(y_pred_sets)

            # Scoring
            metrics[j]['precision'].append(round(precision_score(y_true_bin, y_pred_bin, average='weighted', zero_division=0), 3))
            metrics[j]['recall'].append(round(recall_score(y_true_bin, y_pred_bin, average='weighted', zero_division=0), 3))
            metrics[j]['f1'].append(round(f1_score(y_true_bin, y_pred_bin, average='weighted', zero_division=0), 3))
            
    return metrics


ids = ["Pfam", "PANTHER", "CATH", "EC number", "Rhea ID", "GO terms"]
datasets = {
    "E. coli (Struct)": ("83333_truelabs.csv", "83333_annotated_taxid.xlsx"),
    "S. cere (Struct)": ("559292_truelabs.csv", "559292_annotated_taxid.xlsx"),
    "E. coli (Seq)": ("83333_truelabs.csv", "83333_annotated_taxid_mmseqs2.xlsx"),
    "S. cere (Seq)": ("559292_truelabs.csv", "559292_annotated_taxid_mmseqs2.xlsx")
}

for name, paths in datasets.items():
    print(f"\n{'='*30}\n{name}\n{'='*30}")
    results = get_metrics(paths[0], paths[1], ids)
    
    # Track F1 arrays to calculate average later
    all_f1_arrays = []

    for label in ids:
        f1_scores = results[label]['f1']
        all_f1_arrays.append(f1_scores)
        print(f"{label} F1: {f1_scores}")

    # Calculate and print the average F1 across all descriptors for this dataset
    avg_f1 = np.mean(all_f1_arrays, axis=0).round(3).tolist()
    print(f"\n>>> OVERALL AVG F1: {avg_f1}")
# %%
