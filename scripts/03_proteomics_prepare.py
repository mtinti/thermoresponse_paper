"""Prepare the cell-extract proteomics matrix for limma.

Steps, applied to the DIA-NN MaxLFQ protein quantities:

1. keep protein groups identified by at least 2 peptides;
2. keep the 9 cell-extract samples (34, 37 and 40 C, 3 replicates each)
   and drop proteins with no value in any of them;
3. loading normalisation: scale each sample so that its median equals the
   mean of the sample medians;
4. impute missing values, separately for each temperature:
   - a protein missing in all 3 replicates gets random values around the
     detection limit (median of the 10 smallest values of each sample,
     truncated normal with SD = 20% of that limit);
   - other missing values are estimated with scikit-learn IterativeImputer
     from the other proteins, and raised to the detection limit if below it;
5. log2.

The functions are adapted from nbdev-proteomics (github.com/mtinti/nbdev-proteomics,
commit 2ebeae76, Apache-2.0). The published table was made with an
unseeded draw for the detection-limit values. This version fixes the seed
and the order of conditions, so its output is the same on every run; only
the 140 imputed proteins differ from the published table.

Usage:
    python scripts/03_proteomics_prepare.py DATA_DIR OUT_DIR
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import truncnorm
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer

SEED = 1
MIN_PEPTIDES = 2


def load_quantities(data_dir):
    quant = pd.read_csv(os.path.join(data_dir, 'protein_groups.txt.gz'), sep='\t')
    # R table() output: the header has one field fewer than the rows.
    peptides = pd.read_csv(os.path.join(data_dir, 'peptides_counts.txt'), sep='\t')
    peptides = peptides.set_index('Var1')
    quant = quant.loc[peptides.index[peptides['Freq'] >= MIN_PEPTIDES]]
    quant = quant.replace(0, np.nan)
    # DIA-NN column names are full Windows paths to the .raw.dia files.
    quant.columns = [c.split('\\')[-1].replace('.dia', '') for c in quant.columns]
    mapping = pd.read_csv(os.path.join(data_dir, 'column_mapping.csv'))
    names = dict(zip(mapping['col_name'],
                     mapping['condition'] + '.' + mapping['replica'].astype(str)))
    return quant.rename(columns=names)[list(names.values())]


def normalise_loading(df):
    medians = df.median(axis=0)
    return df * (medians.mean() / medians)


def detection_limit(df):
    smallest = []
    for column in df.columns:
        values = df[column]
        smallest.extend(values[values > 0].nsmallest(10).tolist())
    return np.median(smallest)


def impute(df, conditions, rng):
    imputed = []
    for condition in sorted(set(conditions)):
        block = df.loc[:, [c for c, k in zip(df.columns, conditions) if k == condition]]
        limit = detection_limit(block)
        all_missing = block.isnull().all(axis=1)

        missing = block[all_missing]
        sd = 0.2 * limit
        a = (0 - limit) / sd
        draws = truncnorm.rvs(a=a, b=np.inf, loc=limit, scale=sd,
                              size=missing.shape, random_state=rng)
        missing = pd.DataFrame(draws, index=missing.index, columns=missing.columns)

        partial = block[~all_missing]
        estimated = IterativeImputer(random_state=0, imputation_order='roman').fit_transform(partial)
        partial = pd.DataFrame(np.maximum(estimated, limit),
                               index=partial.index, columns=partial.columns)

        imputed.append(pd.concat([missing, partial]).sort_index())
    return pd.concat(imputed, axis=1).reindex(columns=df.columns)


def main(data_dir, out_dir):
    quant = load_quantities(data_dir)
    cells = quant.iloc[:, :9].dropna(how='all')
    conditions = [c.split('.')[0] for c in cells.columns]
    imputed = impute(normalise_loading(cells), conditions, np.random.default_rng(SEED))
    os.makedirs(out_dir, exist_ok=True)
    np.log2(imputed).to_csv(os.path.join(out_dir, 'proteomics_log2_imputed.csv'))
    quant.to_csv(os.path.join(out_dir, 'proteomics_quantities.csv'))
    print(f'{out_dir}: {imputed.shape[0]} proteins x {imputed.shape[1]} samples')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
