import os
import argparse
from pathlib import Path
from warnings import warn

import anndata as ad
import pandas as pd
import HTSeq




parser = argparse.ArgumentParser()
parser.add_argument('quant_files', nargs='+', type=Path)
parser.add_argument('-g', '--genes', type=Path, required=True)
parser.add_argument('--tpm', action='store_true', default=False)
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()


column_contents = {}
quant_col_key = 'TPM' if args.tpm else 'NumReads'
for file_path in args.quant_files:
    sample = os.path.split(file_path)[-1]
    column_contents[sample] = pd.read_csv(os.path.join(file_path, 'quant.sf'), sep='\t', index_col=0)[quant_col_key]

transcript_id_to_gene_name = {}
for entry in HTSeq.GFF_Reader(args.genes):
    if entry.type == 'transcript':
        transcript_id = entry.attr['transcript_id']
        if 'gene_name' in entry.attr:
            gene_name = entry.attr['gene_name']
        elif 'gene_id' in entry.attr:
            gene_name = entry.attr['gene_id']
        elif 'Parent' in entry.attr:
            gene_name = entry.attr['Parent'].removeprefix('gene-')
        else:
            print(f'Warning: no gene name or ID found for transcript {transcript_id}')
            continue
        transcript_id_to_gene_name[transcript_id] = gene_name

df_quant = pd.DataFrame(column_contents)

new_index = []
for transcript_id in df_quant.index:
    if transcript_id in transcript_id_to_gene_name:
        new_name = transcript_id_to_gene_name[transcript_id]
    elif transcript_id.split('.')[0] in transcript_id_to_gene_name:
        new_name = transcript_id_to_gene_name[transcript_id.split('.')[0]]
    else:
        new_name = transcript_id
    new_index.append(new_name)

df_quant.index = new_index
df_quant = df_quant.groupby(level=0).sum().T

adata = ad.AnnData(df_quant)
adata.write_h5ad(args.output)