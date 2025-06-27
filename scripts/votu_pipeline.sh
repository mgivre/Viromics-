#!/bin/bash

# === PARAMÈTRES ===
ALIGNMENTS="votu_analysis/vOTUs_alignment.tsv"
CLEANED="filtered_clean.tsv"
MATRIX="vOTUs_braycurtis.mat"
TREE="vOTUs.nwk"
BITSCORE_THRESHOLD=50

# === ÉTAPE 1 : filtrage bitscore et extraction colonnes ===
echo "[INFO] Filtrage des alignements avec bitscore > $BITSCORE_THRESHOLD"
awk '{if ($11 <= 0.05) print $1 "\t" $2 "\t" $12}' votu_analysis/vOTUs_alignment.tsv > filtered.tsv

# === ÉTAPE 2 : nettoyage des identifiants (_xxx) ===
echo "[INFO] Nettoyage des identifiants (_xxx)"
cat filtered.tsv \
  | rev | sed 's/\t[[:digit:]]\+_/\t/' | rev \
  | sed 's/_[[:digit:]]\+\t/\t/' \
  > "$CLEANED"

# === ÉTAPE 3 : génération de la matrice Bray-Curtis ===
echo "[INFO] Calcul de la matrice Bray-Curtis"
python3 generate_bray_matrix.py "$CLEANED" "$MATRIX"

# === ÉTAPE 4 : génération de l’arbre avec RapidNJ ===
echo "[INFO] Construction de l’arbre phylogénétique avec RapidNJ"
rapidnj -i pd "$MATRIX" -o t > "$TREE"

echo "[INFO] Arbre généré : $TREE"
