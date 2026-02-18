#%%

import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.preprocessing import MultiLabelBinarizer

def extract_ec_digits(ec_number, parts):
    """Extracts first N parts of an EC number string."""
    if not ec_number: return ""
    processed_ecs = []
    for ec in ec_number.split(';'):
        split_ec = ec.split('.')
        processed_ecs.append('.'.join(split_ec[:parts]))
    return ';'.join(set(processed_ecs))

def add_ec_levels(df, col='EC_4dig'):
    """Generates 1st, 2nd, and 3rd digit columns for a dataframe."""
    for i in range(1, 4):
        df[f'EC_{i}dig'] = df[col].apply(lambda x: extract_ec_digits(x, i))
    return df

def calculate_metrics(y_true_series, y_pred_series):
    """Calculates multi-label metrics for a specific EC level."""
    mlb = MultiLabelBinarizer()
    # Convert strings to sets of labels
    true_labels = y_true_series.apply(lambda x: set(x.split(';')) if x else set())
    pred_labels = y_pred_series.apply(lambda x: set(x.split(';')) if x else set())
    
    # Fit binarizer on all possible labels in this batch to ensure consistent columns
    mlb.fit(pd.concat([true_labels, pred_labels]))
    
    y_true = mlb.transform(true_labels)
    y_pred = mlb.transform(pred_labels)
    
    return {
        'precision': round(precision_score(y_true, y_pred, average='weighted', zero_division=0), 3),
        'recall': round(recall_score(y_true, y_pred, average='weighted', zero_division=0), 3),
        'f1': round(f1_score(y_true, y_pred, average='weighted', zero_division=0), 3)
    }

def evaluate_wasp_iterations(true_df, excel_path):
    """Processes iterations from Excel sheets and returns metrics history."""
    fsafe = pd.ExcelFile(excel_path)
    sheet_dfs = [fsafe.parse(name) for name in fsafe.sheet_names]
    
    history = []
    # Initial state (Iteration 0)
    for level in range(1, 5):
        history.append({'Iteration': 0, 'EC_Level': level, 'precision': 0.0, 'recall': 0.0, 'f1': 0.0})

    # Prepare ground truth - ensure columns are named correctly before the loop
    # We rename 'EC_1dig' -> 'EC_1dig_true', etc.
    true_labels_renamed = true_df.rename(columns={
        f'EC_{i}dig': f'EC_{i}dig_true' for i in range(1, 5)
    })
    
    for i in range(1, len(sheet_dfs) + 1):
        merged = pd.concat(sheet_dfs[:i])[["UniProt ID", 'EC number']].fillna('')
        merged = merged.groupby('UniProt ID').agg({'EC number': ';'.join}).reset_index()
        
        merged['EC_4dig_pred'] = merged['EC number'].apply(
            lambda x: ';'.join(list(set([t.strip('()').split(', ')[0] for t in x.split(';')])))
        )
        
        # Merge predictions onto the true labels
        combined = true_labels_renamed.merge(
            merged[['UniProt ID', 'EC_4dig_pred']], 
            on='UniProt ID', 
            how='left'
        )
        combined['EC_4dig_pred'] = combined['EC_4dig_pred'].fillna('')
        
        # Manually generate predicted level columns with the '_pred' suffix
        for level in range(1, 4):
            combined[f'EC_{level}dig_pred'] = combined['EC_4dig_pred'].apply(
                lambda x: extract_ec_digits(x, level)
            )
        
        # Now the columns like 'EC_1dig_true' and 'EC_1dig_pred' exist!
        for level in range(1, 5):
            m = calculate_metrics(
                combined[f'EC_{level}dig_true'], 
                combined[f'EC_{level}dig_pred']
            )
            m.update({'Iteration': i, 'EC_Level': level})
            history.append(m)
            
    return pd.DataFrame(history)


# 1. Process New-392 Dataset
true_new = pd.read_csv("new.csv", sep="\t")[["Entry", "EC number"]]
true_new = true_new.rename(columns={'Entry': 'UniProt ID', 'EC number': 'EC_4dig'})
true_new = add_ec_levels(true_new)
df_new = evaluate_wasp_iterations(true_new, "new_annotated_taxid.xlsx")
df_new['Dataset'] = 'New-392'

# 2. Process Price-149 Dataset
true_price_raw = pd.read_csv("price.csv", sep="\t")[["Entry", "EC number"]]
mapping_price = pd.read_csv("price_mapped.tsv", sep="\t").drop_duplicates(['Entry'])
true_price = true_price_raw.merge(mapping_price, on='Entry', how='left')[['UniProt ID', 'EC number']]
true_price = true_price.rename(columns={'EC number': 'EC_4dig'})
true_price = add_ec_levels(true_price)
df_price = evaluate_wasp_iterations(true_price, "price_annotated_taxid.xlsx")
df_price['Dataset'] = 'Price-149'

# 3. Consolidate Source Data
source_data = pd.concat([df_new, df_price], ignore_index=True)

# %%
