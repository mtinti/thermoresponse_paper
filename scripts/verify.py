"""Check the input data and compare rebuilt results with the published ones.

Usage:
    python scripts/verify.py inputs
        Every file in data/ against data/checksums.md5.
    python scripts/verify.py sequences DIR
        Sequence CSVs from 05_extract_sequences.py against data/sequences
        (expected: identical).
    python scripts/verify.py rnaseq TSV
        Table from 02_rnaseq_de.R against data/rnaseq/WT40_vs_WT34.txt
        (expected: same genes, a few class changes at |logFC| = 0.4).
    python scripts/verify.py proteomics CSV
        Table from 04_proteomics_de.R against
        data/proteomics/all_cellular_comparisons.csv.gz (expected: identical
        logFC for proteins without imputed values).
    python scripts/verify.py figure DIR
        AUCs and out-of-fold predictions from 07_figure_roc.py against the
        published figure and expected/ (expected: identical).

Each command exits with status 1 if the expectation is not met.
"""

import gzip
import hashlib
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
EXPECTED = os.path.join(ROOT, 'expected')

# AUCs printed on the published figure (expected/Heat_Shock_Prediction_Up_VS_Down.svg).
PUBLISHED_AUC = {('Proteomics', 'ALL'): 0.72, ('Proteomics', '3UTR'): 0.67,
                 ('Proteomics', '5UTR'): 0.70, ('Transcriptomics', 'ALL'): 0.76,
                 ('Transcriptomics', '3UTR'): 0.78, ('Transcriptomics', '5UTR'): 0.56}


def content_md5(path):
    opener = gzip.open if path.endswith('.gz') else open
    digest = hashlib.md5()
    with opener(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def report(ok, text):
    print(f'{"OK  " if ok else "FAIL"} {text}')
    return 0 if ok else 1


def check_inputs():
    failures = 0
    with open(os.path.join(DATA, 'checksums.md5')) as handle:
        for line in handle:
            if line.strip() and not line.startswith('#'):
                expected, name = line.split()
                failures += report(content_md5(os.path.join(DATA, name)) == expected, name)
    return failures


def check_sequences(out_dir):
    failures = 0
    for region in ['5utr', 'cds', '3utr']:
        name = f'trypanosoma_sequences_{region}.csv'
        same = (content_md5(os.path.join(out_dir, name))
                == content_md5(os.path.join(DATA, 'sequences', name + '.gz')))
        failures += report(same, f'{name} identical to data/sequences')
    return failures


def rnaseq_classes(table):
    """1 for Up, -1 for Down and 0 otherwise, for every gene in the table."""
    usable = ~table.index.str.startswith('gene') & (table['logCPM'] > 1)
    up = usable & (table['logFC'] > 0.4) & (table['FDR'] < 0.01)
    down = usable & (table['logFC'] < -0.4) & (table['FDR'] < 0.01)
    return up.astype(int) - down.astype(int)


def check_rnaseq(path):
    published = pd.read_csv(os.path.join(DATA, 'rnaseq', 'WT40_vs_WT34.txt'), sep='\t', index_col=0)
    rebuilt = pd.read_csv(path, sep='\t', index_col=0)
    failures = report(set(published.index) == set(rebuilt.index),
                      f'same {len(published)} genes tested')
    rebuilt = rebuilt.reindex(published.index)
    print(f'     logFC correlation {published["logFC"].corr(rebuilt["logFC"]):.6f}')
    changed = rnaseq_classes(published).ne(rnaseq_classes(rebuilt))
    print(f'     genes whose Up/Down class changes: {int(changed.sum())}')
    if changed.any():
        print(published.loc[changed, ['logFC', 'FDR']]
              .join(rebuilt[['logFC', 'FDR']], rsuffix='_rebuilt').to_string())
    return failures


def check_proteomics(path):
    published = pd.read_csv(os.path.join(DATA, 'proteomics', 'all_cellular_comparisons.csv.gz'), index_col=0)
    rebuilt = pd.read_csv(path, index_col=0).reindex(published.index)
    imputed = ((published['all_missing_in_one_condition'] == 1)
               | (published['missing_in_one_or more sample'] == 1))
    diff = (published['C34_C40_logFC'] - rebuilt['C34_C40_logFC']).abs()
    failures = report(diff[~imputed].max() < 1e-9,
                      f'C34_C40 logFC identical for {int((~imputed).sum())} proteins without imputed values')
    print(f'     {int(imputed.sum())} proteins with imputed values: largest logFC difference {diff[imputed].max():.3f}')
    return failures


def check_figure(out_dir):
    failures = 0
    auc = pd.read_csv(os.path.join(out_dir, 'auc_up_vs_down.tsv'), sep='\t')
    for _, row in auc.iterrows():
        published = PUBLISHED_AUC[(row['dataset'], row['features'])]
        failures += report(round(row['auc'], 2) == published,
                           f'{row["dataset"]} {row["features"]} AUC {row["auc"]:.4f} (published {published:.2f})')
    for dataset in ['proteomics', 'transcriptomics']:
        name = f'ML_results_Up_vs_Down_{dataset}.tsv'
        published = pd.read_csv(os.path.join(EXPECTED, name), sep='\t', index_col=0)
        rebuilt = pd.read_csv(os.path.join(out_dir, name), sep='\t', index_col=0)
        same = (list(published.index) == list(rebuilt.index)
                and (published['preds'].astype(int).values == rebuilt['preds'].values).all())
        failures += report(same, f'{name}: {len(rebuilt)} out-of-fold predictions identical to expected/')
    return failures


if __name__ == '__main__':
    commands = {'inputs': check_inputs, 'sequences': check_sequences, 'rnaseq': check_rnaseq,
                'proteomics': check_proteomics, 'figure': check_figure}
    if len(sys.argv) < 2 or sys.argv[1] not in commands:
        sys.exit(__doc__)
    sys.exit(1 if commands[sys.argv[1]](*sys.argv[2:]) else 0)
