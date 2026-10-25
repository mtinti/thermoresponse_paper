"""Build the RNA-seq count matrix used for differential expression.

Edits applied before edgeR:

- stacks the TREU927 v68 and Tb427 BES count tables;
- drops the GFP-BLA marker and the 14 full-length BES VSGs (their reads
  are counted on the truncated VSG models instead);
- replaces the row of AGO1 (Tb927.10.10850) with fixed values: the
  observed WT counts and 1 for every ago1 null sample, because the gene is
  deleted in the mutants;
- keeps genes that have a GC content and length from bedtools nuc.

Usage:
    python scripts/01_rnaseq_matrix.py DATA_DIR OUT_CSV

DATA_DIR must contain featurecounts_tb927_v68.tsv.gz,
featurecounts_bes.tsv.gz, GC_content_927.txt.gz and GC_content_BES.txt.gz.
"""

import os
import sys

import pandas as pd

BES_VSGS_REPLACED = [
    'gene-Tb427.BES40.22', 'gene-Tb427.BES129.14', 'gene-Tb427.BES15.12',
    'gene-Tb427.BES28.15', 'gene-Tb427.BES98.12', 'gene-Tb427.BES65.13',
    'gene-Tb427.BES64.2', 'gene-Tb427.BES134.6', 'gene-Tb427.BES122.11',
    'gene-Tb427.BES29.9', 'gene-Tb427.BES56.13', 'gene-Tb427.BES10.10',
    'gene-Tb427.BES126.15', 'gene-Tb427.BES51.12']

AGO1 = 'Tb927.10.10850'
AGO1_ROW = (['Tb927_10_v5.1', '2635529', '2638240', '-', 2712,
             3418, 3683, 3518, 3110, 3276, 3201, 3862, 3895, 3922] + [1] * 18)


def load_gc(path):
    gc = pd.read_csv(path, sep='\t')
    gc = gc.drop_duplicates('4_usercol').set_index('4_usercol')
    gc = gc[['13_seq_len', '6_pct_gc']]
    gc.columns = ['length', 'gccontent']
    return gc


def main(data_dir, out_csv):
    tables = [pd.read_csv(os.path.join(data_dir, name), sep='\t', index_col=0)
              for name in ['featurecounts_tb927_v68.tsv.gz',
                           'featurecounts_bes.tsv.gz']]
    df = pd.concat(tables)
    df = df.drop('gene_GFP-BLA_MARKER')
    df['Chr'] = [c.split(';')[0] for c in df['Chr']]
    df = df.drop(BES_VSGS_REPLACED)
    df.loc[AGO1] = AGO1_ROW

    counts = df.iloc[:, 5:].dropna()
    gc = pd.concat([load_gc(os.path.join(data_dir, name))
                    for name in ['GC_content_927.txt.gz', 'GC_content_BES.txt.gz']])
    counts = counts.join(gc, how='inner')
    counts.to_csv(out_csv)
    print(f'{out_csv}: {counts.shape[0]} genes x {counts.shape[1] - 2} samples')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
