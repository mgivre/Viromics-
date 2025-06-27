!/bin/bash

# === PARAMÈTRES ===
INPUT="viral_analysis_results/prodigal_proteins.faa"
DB="mmseqs_db"
RESULT="mmseqs_result"
TMP="tmp_mmseqs"
OUTPUT="mmseqs_all_vs_all.m8"
EVALUE_THRESHOLD="0.05"
MATRIX="similarity_matrix_mm.tsv"

# === 1. Création base MMseqs ===
mmseqs createdb "$INPUT" "$DB"

# === 2. Alignement all-vs-all très sensible ===
mmseqs search "$DB" "$DB" "$RESULT" "$TMP" \
  --threads 8 \
  -e "$EVALUE_THRESHOLD" \
  -c 0.5 \
  --start-sens 7.5 \
  --min-seq-id 0.0 \
  --max-seqs 100000

# === 3. Export en format tabulaire style BLAST (format 6) ===
mmseqs convertalis "$DB" "$DB" "$RESULT" "$OUTPUT" \
  --format-output "query,target,evalue,bits"

# === 4. Filtrage & génération de la matrice de similarité ===
awk -v threshold="$EVALUE_THRESHOLD" '$3 <= threshold { print $1 "\t" $2 "\t" $4 }' "$OUTPUT" \
| rev | sed 's/\t[[:digit:]]\+_/\t/' | rev \
| sed 's/_[[:digit:]]\+\t/\t/' \
| sort \
| ./hashsums \
| ./tree_bray > "$MATRIX"