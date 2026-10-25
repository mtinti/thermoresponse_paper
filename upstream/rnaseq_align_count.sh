#!/usr/bin/env bash
# Align and count one paired-end RNA-seq sample.
#
# These are the commands that produced the featureCounts tables in
# data/rnaseq/, with site-specific paths and QC-only steps (FastQC,
# qualimap, bamCoverage, FastQ Screen) removed. Each sample was run twice:
#
#   reference   MAPQ   output table
#   tb927_v68   0      data/rnaseq/featurecounts_tb927_v68.tsv.gz
#   BES         2      data/rnaseq/featurecounts_bes.tsv.gz
#
# tb927_v68 is the TriTrypDB release 68 TREU927 genome extended with the
# kinetoplast maxicircle (M94286.1), Tb427 VSG genes and the selectable
# marker genes of the cell lines; the contig of each gene is in the Chr
# column of the count table. BES is the set of Tb427 bloodstream-form
# expression sites. Each reference needs a bowtie2 index (<prefix>.*.bt2)
# and a GTF (<prefix>.gtf).
#
# Tool versions: fastp 0.20.0, bowtie2 (very-sensitive-local), samtools,
# Picard 2.22.3, featureCounts (Subread) 1.6.4.
#
# Usage:
#   upstream/rnaseq_align_count.sh SAMPLE FASTQ_1 FASTQ_2 REFERENCE_PREFIX MAPQ OUT_DIR
#
# The per-sample counts.txt files are then combined with
# upstream/collect_featurecounts.py.

set -euo pipefail

if [ $# -ne 6 ]; then
    sed -n '2,24p' "$0"
    exit 1
fi

sample=$1; fastq_1=$2; fastq_2=$3; reference=$4; mapq=$5; out=$6
threads=${THREADS:-8}
mkdir -p "$out"

# Trim adapters and remove low-quality reads (fastp defaults).
fastp -i "$fastq_1" -I "$fastq_2" \
      -o "$out/${sample}_1.trimmed.fq.gz" -O "$out/${sample}_2.trimmed.fq.gz" \
      -h "$out/${sample}.fastp.html" -j "$out/${sample}.fastp.json"

# Align, keep reads at or above the MAPQ threshold, sort.
bowtie2 --very-sensitive-local -p "$threads" -x "$reference" \
        -1 "$out/${sample}_1.trimmed.fq.gz" -2 "$out/${sample}_2.trimmed.fq.gz" \
        2> "$out/${sample}.bowtie2.log" \
    | samtools view -bSu -q "$mapq" \
    | samtools sort -@ "$threads" -o "$out/${sample}.sorted.bam"

# Mark duplicates; they are kept and counted.
picard MarkDuplicates I="$out/${sample}.sorted.bam" O="$out/${sample}.marked.bam" \
       M="$out/${sample}.dupmetrics.txt" REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate
samtools index "$out/${sample}.marked.bam"

# Count read pairs on CDS features, per gene.
# -p -B: pairs with both ends mapped; -C: no chimeric pairs;
# -M: count multi-mapping reads; -O: count reads overlapping several genes.
featureCounts -p -B -C -M -O -T "$threads" -t CDS -g gene_id \
              -a "$reference.gtf" -o "$out/counts.txt" "$out/${sample}.marked.bam"
