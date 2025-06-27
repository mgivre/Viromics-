#!/bin/bash

# === PARAMÈTRES ===
INPUT="vOTUs.faa"
DB="diamond_db"
OUTPUT="diamond_all_vs_all.tsv"
EVALUE_THRESHOLD="0.05"
MATRIX="similarity_matrix.tsv"

# === 1. Construction de la base DIAMOND ===
diamond makedb --in "$INPUT" -d "$DB"

# === 2. Alignement all-vs-all ===
diamond blastp -d "$DB" -q "$INPUT" -o "$OUTPUT" \
  -f 6 qseqid sseqid evalue bitscore \
  --evalue "$EVALUE_THRESHOLD" \
  --more-sensitive \
  --threads 8

# === 3. Filtrage & construction d’une pseudo-matrice ===
awk -v threshold="$EVALUE_THRESHOLD" '$3 <= threshold { print $1 "\t" $2 "\t" $4 }' "$OUTPUT" \
| rev | sed 's/\t[[:digit:]]\+_/\t/' | rev \
| sed 's/_[[:digit:]]\+\t/\t/' \
| sort \
| ./hashsums \
| ./tree_bray > "$MATRIX"
