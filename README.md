# WASP: Protein Functional Annotation using AlphaFold structures

Welcome to the official repository for the paper *WASP: A pipeline for functional annotation based on AlphaFold structural models*!

**WASP**, **W**hole-proteome **A**nnotation via **S**tructural-homology **P**rediction, is a python-based software designed for comprehensive organism annotation at the whole-proteome level based on structural homology.

WASP is a user-friendly command-line tool that only requires the NCBI taxonomy ID of the organism of interest as an input. Using the computational speed of Foldseek [[1](https://doi.org/10.1038/s41587-023-01773-0)], WASP generates a graphical representation of reciprocal hits between the organism protein query and the AlphaFold database [[2](https://doi.org/10.1038/s41586-021-03819-2), [3](https://doi.org/10.1093/nar/gkab1061)], enabling downstream robust functional enrichment and statistical testing. WASP annotates uncharacterised proteins using multiple functional descriptors, including GO terms, Pfam domains, PANTHER family classification and CATH superfamilies, Rhea IDs and EC numbers. Additionally, WASP provides a module to map native proteins to orphan reactions in genome-scale models based on structural homology.

<p align="center">
<img src="logo.png" alt="drawing" width="800"/>
</p>

<!---
If you find WASP helpful in your research, please cite us:

    @article{,
      author = {},
      title = {},
      journal = {},
      volume = {},
      number = {},
      pages = {},
      year = {},
      doi = {},
      URL = {}
    }  
-->

## Table of Contents

1. [Installation](#1-installation)
   - [External requirements](#11-external-requirements)
   - [Quickstart](#12-quickstart)
2. [Usage](#2-usage)
   - [Whole-proteome annotation](#21-whole-proteome-annotation)
   - [GEM gap-filling module](#22-gem-gap-filling-module)
3. [Additional utils](#3-additional-utils)
   - [Download custom proteins](#31-downloading-custom-protein-structures-from-afdb)
   - [Predict custom strucrures](#32-predicting-structures-with-colabfold)
4. [References](#4-references)

## 1. Installation

### 1.1 External requirements

- Foldseek == 8
- gsutil (for accessing AlphaFold DB)

The manuscript results were obtained using Python 3.10.14 and Foldseek 8-ef4e960.

### 1.2 Quickstart

### ✅ Recommended: Install from PyPI

1. **Create and activate a virtual environment**:

```bash
python3 -m venv wasp-env
source wasp-env/bin/activate
```

or 

```bash
conda create -n wasp-env python==3.10 -y
conda activate wasp-env
```

2. **Install WASP**:

```bash
pip install wasp-proteins-annotation
```

Then, install Foldseek in the `/bin` of the project root. Follow the installation instructions at [Foldseek GitHub](https://github.com/steineggerlab/foldseek).

Example for a Linux AVX2 build:

```bash
# Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

```

Install `gsutil` and initialise the gcloud CLI, following instructions for your machine at [Google Cloud Storage Documentation](https://cloud.google.com/storage/docs/gsutil_install).


### 🛠️ For Development: Install from Source

1. **Clone the repository**:

```bash
git clone https://github.com/gioodm/wasp-proteins-annotation.git
cd WASP
```

2. **Create and activate a virtual environment**:

```bash
python3 -m venv wasp-env
source wasp-env/bin/activate
```

or 

```bash
conda create -n WASP python==3.10 -y
conda activate WASP
```

2. **Install the package locally**:

```bash
pip install .
```

Install Foldseek and `gsutil`.


## 2. Usage

### 2.1 Whole-proteome annotation

#### Usage

First, identify the NCBI taxonomy ID of your organism of interest. You can find the ID at [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy). You can use AlphaFold DB to check how many structures are linked to that ID by searching the same ID on the AFDB search bar.

Once you have identified the taxid, run the pipeline as follows (e.g., for organism S. cerevisiae S288c, taxid: 559292):

```bash
wasp-run -t 559292
```

This requires `gsutil` installed.    
Note: On the first run, the AlphaFold DB clustered at 50% will need to be downloaded using Foldseek. This can take some time depending on your machine's performance and internet speed. Ensure you have sufficient storage space to host the database. After the initial setup, WASP annotation will take up to 3 hours for iteration for a proteome of approximately 6000 proteins.

Additional parameters can be customised, including:

- `-e` evalue_threshold: set the evalue threshold (default: 10e-10)
- `-b` bitscore threshold: set the bitscore threshold (default: 50)
- `-n` max_neighbours: set the max number of neighbours (default: 10)
- `-s` step: set step to add to max neighbours (n) in additional iterations (default: 10)
- `-i` iterations: set number of iterations to perform (default: 3)
- `-m` max sequence identity: set maximum sequence identity (0.0 to 1.0) to exclude sequence homologs (default: 1.0)

Usage examples:

```sh
wasp-run -e 1e-50 -b 200 -n 5 -i 5 -t 559292
wasp-run -s 5 -t 559292
```

To use a custom dataset (e.g., a newly sequenced genome or a set of proteins from different organisms), create a tarred folder containing the protein structures (`.cif.gz` or `.pdb.gz` format) and place it in a folder called `proteomes/` within the WASP folder - example folder in `example_files/price.tar`. Then run WASP with:

```bash
wasp-run -t your_custom.tar
```

An example of the final output can be found at `example_files/price_annotated.xlsx`

### 2.2 GEM gap-filling module

#### Pre-processing

Some pre-processing steps are required to obtain a standardized input file for the WASP pipeline. Potential modifications to the Python scripts might be needed depending on the GEM format.

- `find-orphans`: identifies orphan reactions in the GEM (accepted extensions: `.xml`, `.sbml`, `.json`, and `.mat`) using the Python3 `cobrapy` module. Annotation in different formats (accepted: BiGG, Rhea, EC number, KEGG, PubMed, MetaNetX) present in the model is retrieved - input example at: `example_files/gap_filling/iYLI649.xml`, output example at: `example_files/gap_filling/iYLI649_orphans.txt`.

- `run-rxn2code`: each reaction annotation is mapped to the corresponding Rhea reaction ID and/or EC number when available. The output file includes: 1. reaction id, 2. reaction extended name, 3. Rhea IDs, 4. EC numbers.
If MetaNetX codes are present, the `reac_xref.tsv` file (retrieved from [MetaNetX](https://www.metanetx.org/mnxdoc/mnxref.html)) is needed - input example at: `example_files/gap_filling/iYLI649_orphans.txt`, output example at: `example_files/gap_filling/iYLI649_gaps.txt`.


The final `gaps_file.txt` must contain all 4 columns; if the reaction extended name is not present in the model, an empty column should be present. If no Rhea/EC IDs are identified, empty columns should be present instead.

#### Usage

```bash
wasp-gem [-h] [-e evalue_threshold] [-b bitscore_threshold] [-t tmscore] taxid gaps_file.txt
```

An example of the final output can be found at `example_files/gap_filling/iYLI649_hits.txt`

## 3 Additional utils

### 3.1 Downloading custom protein structures from AFDB
The afdb-download utility allows you to fetch protein structures from the AlphaFold Database (AFDB) using either a FASTA file containing protein sequences or a list of UniProt IDs.

#### Usage:

```bash
afdb-download [-h] --input_file INPUT_FILE --output_dir OUTPUT_DIR [--num_cores NUM_CORES]
```

### 3.2 Predicting structures with ColabFold
A Python script (`predict_structures.py`) is provided to run ColabFold predictions on a list of UniProt IDs.

Prerequisites:
- [ColabFold](https://github.com/sokrypton/ColabFold) must be installed and accessible in your PATH.
- Ensure your system meets hardware requirements (GPU recommended for large predictions).

#### Usage:

```bash
python3 predict_structures.py -i <uniprot_ids.txt> -o <output_dir> [-c <cores>]  
```

#### Best practices & recommendations

For large-scale predictions, combine proteins into a single FASTA file to reduce overhead, but split into smaller batches (50-100 proteins) to avoid memory issues. Use parallel processing (`-c 8` for 8 cores) and submit as a batch job (e.g., SLURM) for heavy workloads.

Pre-download FASTA sequences to avoid failures. Test with one protein first to verify setup. If predictions fail, reduce `--num-recycle` (default: 3) or check GPU memory.

For custom workflows, modify the script (e.g., add `--amber` or `--templates`). See the ColabFold docs for advanced options.

## 4 References

[1] van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, et al. Fast and accurate protein structure
search with Foldseek. Nature Biotechnology. 2023. doi: https://doi.org/10.1038/s41587-023-01773-0

[2] Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, et al. Highly accurate protein structure prediction
with AlphaFold. Nature. 2021;596(7873):583-9. doi: https://doi.org/10.1038/s41586-021-03819-2

[3] Varadi M, Bertoni D, Magana P, Paramval U, Pidruchna I, Radhakrishnan M, et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Research. 2023;52(D1):D368–D375. doi: http://dx.doi.org/10.1093/nar/gkad1011
