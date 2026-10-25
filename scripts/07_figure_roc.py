"""ROC curves for predicting heat-shock Up vs Down genes from HyenaDNA embeddings.

Draws the two ROC panels of the figure "Heat Shock Prediction Up vs Down".

Gene classes, 40 C vs 34 C, wild-type cells:
- proteomics (limma): average expression > 18, adjusted P < 0.01,
  Up logFC > 0.3, Down logFC < -0.55;
- transcriptomics (edgeR): logCPM > 1, FDR < 0.01, |logFC| > 0.4.
Both lists keep only genes in data/gene_universe/ml_dataset_gene_ids.txt.

Features are the mean-pooled HyenaDNA embeddings (256 values per region).
"ALL" means 5'UTR + 3'UTR; the CDS embedding is not used. Scores are out-of-fold probabilities from 3-fold stratified
cross-validation of a random forest.

The published AUCs are reproduced exactly in the Docker image (linux/amd64,
scikit-learn 1.6.1). Other platforms or versions grow slightly different
forests and move the AUCs by up to 0.02.

Usage:
    python scripts/07_figure_roc.py [--proteomics CSV]
        [--transcriptomics TSV] [--out DIR] [--repeats N]

--repeats N also reruns the cross-validation with N different seeds and
reports the mean and standard deviation of each AUC.
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
MODEL = 'hyenadna-medium-450k-seqlen-hf'
FEATURE_SETS = {'ALL': ('5utr', '3utr'), '3UTR': ('3utr',), '5UTR': ('5utr',)}
COLOURS = {'ALL': 'magenta', '3UTR': 'cyan', '5UTR': 'orange'}


def load_embeddings():
    frames = []
    for region in ['cds', '5utr', '3utr']:
        ids = pd.read_csv(os.path.join(
            DATA, 'sequences', f'trypanosoma_sequences_{region}.csv.gz'))['gene_id']
        values = np.load(os.path.join(
            DATA, 'embeddings',
            f'embeddings_927_{MODEL}_trypanosoma_sequences_{region}.csv.npy'))
        frames.append(pd.DataFrame(
            values, index=ids.values,
            columns=[f'component_{i}_{region}' for i in range(values.shape[1])]))
    return frames[0].join(frames[1]).join(frames[2])


def gene_universe():
    with open(os.path.join(DATA, 'gene_universe', 'ml_dataset_gene_ids.txt')) as handle:
        return set(handle.read().split())


def label(up, down):
    labelled = pd.concat([up, down])
    labelled['Target'] = [1] * len(up) + [0] * len(down)
    return labelled[['Target']]


def proteomics_classes(path, genes):
    prot = pd.read_csv(path)
    # A protein group can list several proteins; each one gets the group's values.
    prot['protein'] = prot['Unnamed: 0'].str.split(';')
    prot = prot.explode('protein').reset_index(drop=True)
    prot['protein'] = [p.split(':')[0] for p in prot['protein']]
    prot = prot[prot['protein'].isin(genes)].set_index('protein')
    # limma reports 34 C minus 40 C; flip it so that Up means higher at 40 C.
    prot['logFC'] = -prot['C34_C40_logFC']
    prot = prot.sort_values('logFC')
    prot = prot[prot['C34_C40_AveExpr'] > 18]
    significant = prot['C34_C40_adj.P.Val'] < 0.01
    return label(prot[(prot['logFC'] > 0.3) & significant],
                 prot[(prot['logFC'] < -0.55) & significant])


def transcriptomics_classes(path, genes):
    rna = pd.read_csv(path, sep='\t', index_col=0)
    rna = rna[~rna.index.str.startswith('gene')]
    rna = rna[rna.index.isin(genes)].sort_values('logFC')
    rna = rna[rna['logCPM'] > 1]
    significant = rna['FDR'] < 0.01
    return label(rna[(rna['logFC'] > 0.4) & significant],
                 rna[(rna['logFC'] < -0.4) & significant])


def out_of_fold_scores(X, y, seed=42):
    model = RandomForestClassifier(n_estimators=50, random_state=seed,
                                   min_samples_leaf=2, max_depth=10,
                                   class_weight='balanced')
    folds = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    return cross_val_predict(model, X, y, cv=folds, method='predict_proba')[:, 1]


def evaluate(embeddings, classes, repeats):
    data = embeddings.join(classes, how='left').dropna(subset=['Target'])
    y = data['Target']
    results = {}
    for name, regions in FEATURE_SETS.items():
        columns = [c for c in data.columns
                   if c.startswith('component_') and c.rsplit('_', 1)[1] in regions]
        scores = out_of_fold_scores(data[columns], y)
        fpr, tpr, _ = roc_curve(y, scores)
        repeated = [roc_auc_score(y, out_of_fold_scores(data[columns], y, seed))
                    for seed in range(repeats)]
        results[name] = {'fpr': fpr, 'tpr': tpr, 'auc': roc_auc_score(y, scores),
                         'scores': scores, 'repeated': repeated}
    return y, results


def plot(panels, out_dir):
    fig, axes = plt.subplots(ncols=2, figsize=(10, 5))
    for ax, (title, (_, results)) in zip(axes, panels.items()):
        for name, r in results.items():
            ax.plot(r['fpr'], r['tpr'], c=COLOURS[name],
                    label=f'{name} (AUC = {r["auc"]:.2f})')
        ax.plot([0, 1], [0, 1], linestyle='--', c='gray', label='Random')
        ax.set_title(title)
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.legend()
    fig.suptitle('Heat Shock Prediction Up vs Down', fontsize=18)
    fig.tight_layout()
    for ext in ['png', 'svg']:
        fig.savefig(os.path.join(out_dir, f'Heat_Shock_Prediction_Up_VS_Down_ROC.{ext}'))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--proteomics', default=os.path.join(
        DATA, 'proteomics', 'all_cellular_comparisons.csv.gz'))
    parser.add_argument('--transcriptomics', default=os.path.join(
        DATA, 'rnaseq', 'WT40_vs_WT34.txt'))
    parser.add_argument('--out', default=os.path.join(ROOT, 'results', 'figures'))
    parser.add_argument('--repeats', type=int, default=0)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    plt.style.use('ggplot')
    embeddings = load_embeddings()
    genes = gene_universe()
    panels = {
        'Proteomics': evaluate(embeddings, proteomics_classes(args.proteomics, genes), args.repeats),
        'Transcriptomics': evaluate(embeddings, transcriptomics_classes(args.transcriptomics, genes), args.repeats),
    }
    plot(panels, args.out)

    rows = []
    for title, (y, results) in panels.items():
        predictions = pd.DataFrame({'Target': y.astype(int),
                                    'preds': (results['ALL']['scores'] > 0.5).astype(int)})
        predictions.index.name = 'gene_id'
        predictions.to_csv(os.path.join(
            args.out, f'ML_results_Up_vs_Down_{title.lower()}.tsv'), sep='\t')
        for name, r in results.items():
            row = {'dataset': title, 'features': name, 'n_up': int(y.sum()),
                   'n_down': int((y == 0).sum()), 'auc': round(r['auc'], 4)}
            if r['repeated']:
                row['auc_mean'] = round(float(np.mean(r['repeated'])), 4)
                row['auc_sd'] = round(float(np.std(r['repeated'], ddof=1)), 4)
            rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(args.out, 'auc_up_vs_down.tsv'), sep='\t', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
