"""Extract 5'UTR, CDS and 3'UTR sequences for the HyenaDNA embeddings.

Reads the curated TREU927 annotation (merged_fixed_6.gtf) and the
TriTrypDB-68 TREU927 genome, and writes one CSV per region with the
columns gene_id and sequence. Rules:

- keep genes with exactly one CDS and exactly one transcript;
- keep genes whose transcript extends more than 20 nt beyond the CDS at
  both ends, so both UTRs exist;
- keep genes whose CDS starts with ATG.

Row order follows the GTF, which is why the files start on chromosome 8.

Usage:
    python scripts/05_extract_sequences.py GTF GENOME_FASTA OUT_DIR
"""

import os
import re
import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

GTF_COLUMNS = ['seqid', 'source', 'type', 'start', 'end',
               'score', 'strand', 'phase', 'attributes']
GENE_ID = re.compile(r'gene_id.\"(.*?)\"')
MIN_UTR = 20


def load_single_isoform_genes(gtf_path):
    gtf = pd.read_csv(gtf_path, sep='\t', header=None, names=GTF_COLUMNS)
    gtf['gene_id'] = [GENE_ID.findall(a)[0] for a in gtf['attributes']]
    # One CDS row and one transcript row per gene, filtered in this order.
    for feature in ['CDS', 'transcript']:
        counts = gtf.loc[gtf['type'] == feature, 'gene_id'].value_counts()
        gtf = gtf[gtf['gene_id'].isin(counts[counts == 1].index)]
    return gtf


def feature_table(gtf, feature, suffix):
    table = gtf.loc[gtf['type'] == feature,
                    ['seqid', 'start', 'end', 'strand', 'gene_id']]
    table.columns = [c + suffix if c != 'gene_id' else c for c in table.columns]
    return table.set_index('gene_id')


def utr_coordinates(row):
    if row['strand_cds'] == '+':
        return (row['start_transcript'], row['start_cds'] - 1,
                row['end_cds'] + 1, row['end_transcript'])
    return (row['end_cds'] + 1, row['end_transcript'],
            row['start_transcript'], row['start_cds'] - 1)


def main(gtf_path, genome_path, out_dir):
    gtf = load_single_isoform_genes(gtf_path)
    genes = feature_table(gtf, 'CDS', '_cds').join(
        feature_table(gtf, 'transcript', '_transcript'), how='left')
    genes = genes[(genes['start_cds'] - genes['start_transcript']).abs() > MIN_UTR]
    genes = genes[(genes['end_cds'] - genes['end_transcript']).abs() > MIN_UTR]

    genome = {r.id: str(r.seq) for r in SeqIO.parse(genome_path, 'fasta')}
    records = []
    for gene_id, row in genes.iterrows():
        chrom = genome[row['seqid_transcript']]
        five_start, five_end, three_start, three_end = utr_coordinates(row)
        five = chrom[five_start - 1:five_end]
        three = chrom[three_start - 1:three_end]
        cds = chrom[row['start_cds'] - 1:row['end_cds']]
        if row['strand_cds'] == '-':
            five, three, cds = (str(Seq(s).reverse_complement())
                                for s in (five, three, cds))
        records.append((gene_id, five, cds, three))

    seqs = pd.DataFrame(records, columns=['gene_id', '5utr', 'cds', '3utr'])
    seqs = seqs[seqs['cds'].str.startswith('ATG')]
    os.makedirs(out_dir, exist_ok=True)
    for region in ['5utr', 'cds', '3utr']:
        out = os.path.join(out_dir, f'trypanosoma_sequences_{region}.csv')
        seqs[['gene_id', region]].rename(columns={region: 'sequence'}).to_csv(
            out, index=False)
        print(f'{out}: {len(seqs)} genes')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
