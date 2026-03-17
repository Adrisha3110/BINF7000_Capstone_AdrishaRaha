#!/usr/bin/env python3

import os
import sys
import argparse
from collections import defaultdict, OrderedDict

# ------------------------------------------------------------
# FASTA helpers
# ------------------------------------------------------------

def read_fasta_entries(path):
    entries = []
    header = None
    seq_chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_chunks)))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        entries.append((header, "".join(seq_chunks)))

    return entries


def parse_block_header(header):
    if "-" not in header:
        return None, None, None, None, None

    species, rest = header.split("-", 1)

    if "/" not in rest:
        return species, None, None, None, None

    coord_part, pos_part = rest.split("/", 1)

    if "(" in coord_part:
        chrom = coord_part.split("(")[0]
        strand = coord_part.split("(")[1].split(")")[0]
    else:
        chrom = coord_part
        strand = "."

    if "-" not in pos_part:
        return species, chrom, strand, None, None

    s, e = pos_part.split("-", 1)

    try:
        start = int(s)
        end = int(e)
    except:
        start = None
        end = None

    return species, chrom, strand, start, end


# ------------------------------------------------------------
# GTF loading
# ------------------------------------------------------------

def parse_attributes(attr_str):
    attrs = {}
    for field in attr_str.split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            continue
        key, val = field.split(" ", 1)
        attrs[key] = val.strip().strip('"')
    return attrs


def load_ccds_gtf(gtf_path, chrom):

    transcripts = {}
    cds_intervals_by_chrom = defaultdict(list)

    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue

            g_chrom, source, feature, start, end, score, strand, frame, attrs = parts

            if g_chrom != chrom:
                continue

            if feature != "CDS":
                continue

            start = int(start)
            end = int(end)

            attr = parse_attributes(attrs)

            gene = attr.get("gene_name") or attr.get("gene_id")
            transcript = attr.get("transcript_id")

            if gene is None or transcript is None:
                continue

            key = (gene, transcript, g_chrom, strand)

            if key not in transcripts:
                transcripts[key] = {
                    "gene": gene,
                    "transcript": transcript,
                    "chrom": g_chrom,
                    "strand": strand,
                    "cds": []
                }

            transcripts[key]["cds"].append((start, end))
            cds_intervals_by_chrom[g_chrom].append((start, end, gene, transcript, strand))

    for key in transcripts:
        transcripts[key]["cds"].sort()

    for c in cds_intervals_by_chrom:
        cds_intervals_by_chrom[c].sort()

    return transcripts, cds_intervals_by_chrom


def find_transcripts(chrom, start, cds_intervals_by_chrom):

    hits = []

    if chrom not in cds_intervals_by_chrom:
        return hits

    for cds_start, cds_end, gene, tid, strand in cds_intervals_by_chrom[chrom]:

        if cds_start > start:
            break

        if cds_end < start:
            continue

        hits.append((gene, tid, strand))

    return hits


# ------------------------------------------------------------
# GAP REMOVAL (hg38)
# ------------------------------------------------------------

def remove_hg38_gaps(full_seqs, species_order, hg_species):

    hg_seq = full_seqs[hg_species]

    gap_positions = [i for i, b in enumerate(hg_seq) if b == "-"]

    if not gap_positions:
        return full_seqs

    cleaned = {}

    for sp in species_order:

        seq = list(full_seqs[sp])

        for pos in reversed(gap_positions):
            seq.pop(pos)

        cleaned[sp] = "".join(seq)

    return cleaned


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--gtf", required=True)
    parser.add_argument("--blocks-dir", required=True)
    parser.add_argument("--out-fasta-dir", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--hg-species", default="hg38")

    args = parser.parse_args()

    os.makedirs(args.out_fasta_dir, exist_ok=True)

    transcripts, cds_intervals = load_ccds_gtf(args.gtf, args.chrom)

    transcript_seqs = {}

    species_order = None

    files = sorted(f for f in os.listdir(args.blocks_dir) if f.endswith(".fa"))

    for fname in files:

        entries = read_fasta_entries(os.path.join(args.blocks_dir, fname))

        if species_order is None:
            species_order = []
            seen = set()

            for h, _ in entries:
                sp = h.split("-", 1)[0]
                if sp not in seen:
                    seen.add(sp)
                    species_order.append(sp)

        block = {}
        hg_info = None

        for h, seq in entries:

            sp, chrom, strand, start, end = parse_block_header(h)

            block[sp] = seq

            if sp == args.hg_species:
                hg_info = (chrom, start, end)

        if hg_info is None:
            continue

        chrom, start, end = hg_info

        if chrom != args.chrom:
            continue

        hits = find_transcripts(chrom, start, cds_intervals)

        for gene, tid, strand in hits:

            key = (gene, tid, chrom, strand)

            if key not in transcript_seqs:
                transcript_seqs[key] = []

            transcript_seqs[key].append((start, block))


    for key, segments in transcript_seqs.items():

        gene, tid, chrom, strand = key

        segments.sort(key=lambda x: x[0], reverse=(strand == "-"))

        full = {sp: [] for sp in species_order}

        for _, block in segments:
            for sp in species_order:
                seq = block.get(sp)
                if seq is None:
                    seq = "-" * len(block[args.hg_species])
                full[sp].append(seq)

        for sp in full:
            full[sp] = "".join(full[sp])

        full = remove_hg38_gaps(full, species_order, args.hg_species)

        name = f"{gene}_{tid}.fa"
        out = os.path.join(args.out_fasta_dir, name)

        with open(out, "w") as f:

            for sp in species_order:
                f.write(f">{sp}\n")

                seq = full[sp]

                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")


if __name__ == "__main__":
    main()
