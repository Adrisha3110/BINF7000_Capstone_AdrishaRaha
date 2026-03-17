#!/bin/bash

#SBATCH --job-name=chr_pad_blocks_to_genes_remove_hg38_gaps
#SBATCH --partition=short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=raha.a@northeastern.edu
#SBATCH --output=/courses/BINF7700.202630/students/raha.a/CAPSTONE/logs/%x_%j.log
#SBATCH --error=/courses/BINF7700.202630/students/raha.a/CAPSTONE/logs/%x_%j.err

set -euo pipefail

# -----------------------------------------
# Usage
# -----------------------------------------
# sbatch chr_pad_blocks_to_genes_remove_hg38_gaps.sh chrY
# sbatch chr_pad_blocks_to_genes_remove_hg38_gaps.sh chr21
# -----------------------------------------

CHR=$1

if [ -z "$CHR" ]; then
    echo "Usage: sbatch chr_pad_blocks_to_genes_remove_hg38_gaps.sh chrY"
    exit 1
fi

BASE=/courses/BINF7700.202630/students/raha.a/CAPSTONE

GTF=$BASE/Annotation/CCDS_hg38_one_id_per_gene.gtf
BLOCKS_DIR=$BASE/Results/${CHR}_chunks_padded
OUT_FASTA_DIR=$BASE/Results/gene_level_fastas/${CHR}

mkdir -p "$OUT_FASTA_DIR"

echo "[INFO] Chromosome      : $CHR"
echo "[INFO] Blocks dir      : $BLOCKS_DIR"
echo "[INFO] Output dir      : $OUT_FASTA_DIR"
echo "[INFO] GTF             : $GTF"

python $BASE/Data/Scripts/padded_blocks_to_genes_remove_hg38_gaps.py \
    --gtf "$GTF" \
    --blocks-dir "$BLOCKS_DIR" \
    --out-fasta-dir "$OUT_FASTA_DIR" \
    --chrom "$CHR" \
    --hg-species hg38

echo "[INFO] Completed chromosome: $CHR"
