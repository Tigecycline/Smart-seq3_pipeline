import pandas as pd
import anndata as ad

import HTSeq




def build_trsc_adata_from_umicount(umicount_file, annotation_file=None):
    df = pd.read_csv(umicount_file, sep='\t', index_col=0)

    if annotation_file is not None:
        gtf = HTSeq.GFF_Reader(annotation_file)
        gene_id_to_name = {feature.name: feature.attr['gene_name'] for feature in gtf if feature.type == 'gene' and 'gene_name' in feature.attr.keys()}
        df.index = [gene_id_to_name.get(gene_id, gene_id) for gene_id in df.index]
        df = df.loc[df.index.value_counts() == 1] # drop duplicated genes
    
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