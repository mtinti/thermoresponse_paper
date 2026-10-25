"""Compute mean-pooled HyenaDNA embeddings for one sequence file.

- Model LongSafari/hyenadna-medium-450k-seqlen-hf, revision 42dedd4d374e.
- Sequences are upper-cased and padded to 31,585 tokens, longer than every
  sequence here (the longest CDS is 24,645 nt), so nothing is truncated.
- The embedding of a sequence is the mean of the model's last hidden
  states over non-padding tokens: 256 values per sequence.

Batch size and GPU type change the result only by float rounding (at most
3e-6 per value between two GPU runs with batch sizes 10 and 20), which does
not change the figure. A GPU is needed in practice: about 10 minutes per
region on an NVIDIA T4.

Usage:
    python scripts/06_compute_embeddings.py SEQUENCES_CSV OUT_NPY [--batch-size N]
"""

import argparse

import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer

MODEL = 'LongSafari/hyenadna-medium-450k-seqlen-hf'
REVISION = '42dedd4d374eac0fb8168549e546a3472fbd27ae'
MAX_LENGTH = 31585


def embed(sequences, batch_size, device):
    tokenizer = AutoTokenizer.from_pretrained(MODEL, revision=REVISION, trust_remote_code=True)
    model = AutoModel.from_pretrained(MODEL, revision=REVISION, trust_remote_code=True)
    model = model.to(device).eval()
    out = []
    with torch.no_grad():
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start:start + batch_size]
            tokens = tokenizer(batch, padding='max_length', truncation=True,
                               max_length=MAX_LENGTH, return_tensors='pt')
            ids = tokens['input_ids'].to(device)
            hidden = model(ids)[0]
            mask = (ids != tokenizer.pad_token_id).unsqueeze(-1).to(hidden.dtype)
            pooled = (mask * hidden).sum(dim=1) / mask.sum(dim=1).clamp(min=1e-9)
            out.append(pooled.cpu())
    return torch.cat(out).numpy()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('sequences_csv')
    parser.add_argument('out_npy')
    parser.add_argument('--batch-size', type=int, default=20)
    args = parser.parse_args()
    sequences = pd.read_csv(args.sequences_csv)['sequence'].str.upper().tolist()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    embeddings = embed(sequences, args.batch_size, device)
    np.save(args.out_npy, embeddings)
    print(f'{args.out_npy}: {embeddings.shape}')


if __name__ == '__main__':
    main()
