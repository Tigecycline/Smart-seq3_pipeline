import os
import argparse
from pathlib import Path

import pandas as pd




parser = argparse.ArgumentParser()
parser.add_argument('umite_files', nargs='+', type=Path)
parser.add_argument('--filename_prefix', type=str, required=True)
parser.add_argument('--samplename_suffix', type=str, required=True)
args = parser.parse_args()

for uf in args.umite_files:
    df_count = pd.read_csv(uf, sep='\t', index_col='samples')
    df_count.index = [os.path.basename(sample).removesuffix(args.samplename_suffix) for sample in df_count.index]
    df_count.to_csv(uf, sep='\t')
    os.rename(uf, f'{os.path.dirname(uf)}/{args.filename_prefix}_{os.path.basename(uf)}')
