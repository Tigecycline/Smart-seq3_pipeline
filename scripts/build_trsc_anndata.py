import pickle

import pandas as pd
import anndata as ad




def build_trsc_adata_from_umicount(umicount_file, parsed_gtf_pkl):
    df = pd.read_csv(umicount_file, sep='\t', index_col=0)
    with open(parsed_gtf_pkl, 'rb') as f:
        gene_id_to_name = {gene_id: gene_names[0] for gene_id, gene_names in pickle.load(f)[2].items()}

    df.rename(gene_id_to_name, axis=1, inplace=True)    
    df = df.T.groupby(level=0).sum().T

    adata = ad.AnnData(df)
    
    return adata








if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('umicount_tsv')
    parser.add_argument('-g', '--gtf', required=False)
    parser.add_argument('-o', required=True)
    args = parser.parse_args()

    adata = build_trsc_adata_from_umicount(args.umicount_tsv, args.gtf)
    adata.write(args.o, compression='gzip')