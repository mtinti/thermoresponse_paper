"""Collect the per-sample featureCounts tables into two count tables.

upstream/rnaseq_align_count.sh writes one counts.txt per sample and per
reference:

    <run_dir>/<sample>/res1/<sample>/counts.txt   TREU927 v68 genes
    <run_dir>/<sample>/res2/<sample>/counts.txt   Tb427 BES genes

This script joins the count column of every sample, in the order of
upstream/samples_rnaseq.tsv, and keeps the featureCounts annotation columns
(Chr, Start, End, Strand, Length). It does not change any value.

Usage:
    python upstream/collect_featurecounts.py RUN_DIR OUT_DIR
"""

import os
import sys

import pandas as pd

REFERENCES = {'res1': 'featurecounts_tb927_v68.tsv.gz',
              'res2': 'featurecounts_bes.tsv.gz'}


SAMPLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'samples_rnaseq.tsv')


def sample_order():
    return pd.read_csv(SAMPLES, sep='\t')['sample'].tolist()


def collect(run_dir, res, samples):
    table = None
    for sample in samples:
        path = os.path.join(run_dir, sample, res, sample, 'counts.txt')
        counts = pd.read_csv(path, sep='\t', comment='#', index_col=0)
        # The last column header is the BAM path in the job's TMPDIR.
        counts = counts.rename(columns={counts.columns[-1]: sample})
        table = counts if table is None else table.join(counts[[sample]])
    return table


def main(run_dir, out_dir):
    samples = sample_order()
    os.makedirs(out_dir, exist_ok=True)
    for res, name in REFERENCES.items():
        table = collect(run_dir, res, samples)
        table.to_csv(os.path.join(out_dir, name), sep='\t')
        print(f'{name}: {table.shape[0]} genes x {len(samples)} samples')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
