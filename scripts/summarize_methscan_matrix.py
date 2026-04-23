import argparse
import os
import pickle

import pandas as pd
import anndata as ad




parser = argparse.ArgumentParser()
parser.add_argument('mean_shrunken_residuals_csv')
parser.add_argument('-o', '--output', required=True)
args = parser.parse_args()


adata = ad.AnnData(pd.read_csv(args.mean_shrunken_residuals_csv, index_col=0))
adata.write_h5ad(args.output)
