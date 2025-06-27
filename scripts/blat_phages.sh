#!/bin/bash

# Fichier : cluster_genomes_blat_parallel.sh
# Chemin pour les scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Nombre de processus parallèles (modifiable)
NPROC=30

# 1. Préparer les dossiers
mkdir -p blat_output tmp_chunks

# 2. Fichier d'entrée (multi-fasta)
INPUT=14Apr2025_genomes.fa
BASENAME=$(basename "$INPUT" .fa)

# 3. Découpage du fichier en NPROC chunks
#echo "Découpage de $INPUT en $NPROC morceaux..."
#mkdir -p tmp_chunks tmp_chunks_grouped

# Extraire les séquences individuellement
#csplit -z -f tmp_chunks/chunk_ -b "%03d.fa" "$INPUT" '/^>/' '{*}' >/dev/null

# Lister tous les fichiers extraits
#FILES=(tmp_chunks/chunk_*.fa)
#TOTAL=${#FILES[@]}

# Répartition équilibrée dans NPROC fichiers multi-fasta
#for ((i=0; i<NPROC; i++)); do
    #OUTFILE="tmp_chunks_grouped/chunk_$(printf "%02d" $i).fa"
    #> "$OUTFILE"  # Vide d'abord
    #for ((j=i; j<TOTAL; j+=NPROC)); do
        #cat "${FILES[$j]}" >> "$OUTFILE"
    #done
#done

# 4. Lancer BLAT en parallèle (chaque chunk contre tout)
#echo "Lancement de BLAT en parallèle (max $NPROC jobs)..."
#parallel -j "$NPROC" blat "$INPUT" {} blat_output/{/.}.blat -out=blast8 ::: tmp_chunks_grouped/chunk_*.fa


#echo "BLAT terminé pour tous les morceaux."

# 5. Concaténer les résultats
#cat blat_output/chunk_*.blat > blat_output/${BASENAME}.blat

# 6. Calcul des longueurs des séquences - VERSION CORRIGÉE
echo "Calcul des longueurs des séquences..."

# Première étape : extraire les longueurs des séquences du fichier FASTA
cat "$INPUT" | "$SCRIPT_DIR/f2s" | "$SCRIPT_DIR/seqlengths" > blat_output/${BASENAME}.seq_lengths_full.tmp

# DIAGNOSTIC : Comparer les noms de séquences
echo "Diagnostic des noms de séquences..."
echo "Noms dans le fichier FASTA (premières 3) :"
head -3 blat_output/${BASENAME}.seq_lengths_full.tmp | cut -f1

echo "Noms dans le fichier BLAT (premières 3 uniques) :"
awk '$1 == $2 {print $1}' blat_output/${BASENAME}.blat | head -3

# Extraire seulement le premier mot des noms de séquences (avant le premier espace)
echo "Normalisation des noms de séquences..."
awk '{split($1, a, " "); print a[1] "\t" $2}' blat_output/${BASENAME}.seq_lengths_full.tmp > blat_output/${BASENAME}.seq_lengths.tmp

echo "Noms normalisés (premières 3) :"
head -3 blat_output/${BASENAME}.seq_lengths.tmp

# Debug : vérifier qu'on a bien des alignements auto dans le fichier BLAT
echo "Vérification des self-alignments dans le fichier BLAT..."
SELF_COUNT=$(awk '$1 == $2' blat_output/${BASENAME}.blat | wc -l)
echo "Nombre de self-alignments trouvés : $SELF_COUNT"

if [ "$SELF_COUNT" -eq 0 ]; then
    echo "ERREUR : Aucun self-alignment trouvé dans le fichier BLAT !"
    head -5 blat_output/${BASENAME}.blat
    exit 1
fi

# Deuxième étape : extraire les alignements auto (query == subject) avec leurs scores
echo "Extraction des self-alignments..."
awk '$1 == $2 {print $1 "\t" $12}' blat_output/${BASENAME}.blat > blat_output/${BASENAME}.self_alignments.tmp

echo "Premières lignes des self-alignments :"
head -3 blat_output/${BASENAME}.self_alignments.tmp

# Vérifier qu'on a des scores non-nuls
NONZERO_COUNT=$(awk '$2 > 0' blat_output/${BASENAME}.self_alignments.tmp | wc -l)
echo "Nombre de self-alignments avec score > 0 : $NONZERO_COUNT"

# Troisième étape : traiter les self alignments avec hashsums
cat blat_output/${BASENAME}.self_alignments.tmp | "$SCRIPT_DIR/hashsums" > blat_output/${BASENAME}.self_scores_raw.tmp

# Enlever la première ligne vide produite par hashsums
sed '1d' blat_output/${BASENAME}.self_scores_raw.tmp > blat_output/${BASENAME}.self_scores.tmp

echo "Premières lignes après hashsums :"
head -3 blat_output/${BASENAME}.self_scores.tmp

# Quatrième étape : joindre les longueurs avec les scores
echo "Jointure des données..."
cat blat_output/${BASENAME}.seq_lengths.tmp | "$SCRIPT_DIR/joincol" blat_output/${BASENAME}.self_scores.tmp > blat_output/${BASENAME}.lengths

# Vérifier que le fichier lengths a été créé correctement
if [[ ! -s blat_output/${BASENAME}.lengths ]]; then
    echo "Erreur : le fichier lengths est vide ou n'a pas été créé"
    exit 1
fi

echo "Fichier lengths créé avec $(wc -l < blat_output/${BASENAME}.lengths) lignes"
echo "Premières lignes du fichier lengths :"
head -3 blat_output/${BASENAME}.lengths

# Vérifier qu'on a des scores non-nuls dans le fichier final
FINAL_NONZERO=$(awk '$3 > 0' blat_output/${BASENAME}.lengths | wc -l)
echo "Séquences avec score > 0 dans le fichier final : $FINAL_NONZERO"

# Nettoyage des fichiers temporaires
rm -f blat_output/${BASENAME}.seq_lengths*.tmp blat_output/${BASENAME}.self_alignments.tmp blat_output/${BASENAME}.self_scores*.tmp

# 7. Identifier les chimères
echo "Identification des chimères..."

# Vérifier d'abord qu'on a des scores non-nuls
NONZERO_SCORES=$(awk '$3 > 0' blat_output/${BASENAME}.lengths | wc -l)
echo "Séquences avec des scores non-nuls : $NONZERO_SCORES"

if [ "$NONZERO_SCORES" -eq 0 ]; then
    echo "ATTENTION : Aucune séquence n'a de score BLAT > 0"
    echo "Création d'une liste de chimères vide..."
    touch blat_output/${BASENAME}.chimeras.list
else
    # Calcul normal des chimères (ratio score/longueur > 2.15)
    awk '$3 > 0 && $3/$2 > 2.15 {print $1}' blat_output/${BASENAME}.lengths > blat_output/${BASENAME}.chimeras.list
fi

CHIMERA_COUNT=$(wc -l < blat_output/${BASENAME}.chimeras.list)
echo "Nombre de chimères détectées : $CHIMERA_COUNT"

# 8. Filtrage et clustering
echo "Filtrage et clustering..."
cut -f1,2,12 blat_output/${BASENAME}.blat | "$SCRIPT_DIR/hashsums" | tail -n +2 \
    | "$SCRIPT_DIR/joincol" blat_output/${BASENAME}.chimeras.list \
    | awk '$NF == 0 {print $1 "\t" $2 "\t" $3}' \
    | "$SCRIPT_DIR/joincol" blat_output/${BASENAME}.lengths \
    | "$SCRIPT_DIR/joincol" blat_output/${BASENAME}.lengths 2 \
    | sort -k4,4nr -k1,1 \
    | awk '$3/$NF >= 0.90 {print $1 "\t" $2}' \
    | perl -lane 'unless (exists($clusters{$F[1]})) {$clusters{$F[1]} = $F[0]; print "$F[1]\t$F[0]"}' \
    > blat_output/OTUs.tsv

# 9. Extraction des séquences représentatives (OTUs)
echo "Extraction des séquences représentatives..."
cat "$INPUT" | "$SCRIPT_DIR/f2s" | "$SCRIPT_DIR/joincol" <(cut -f2 blat_output/OTUs.tsv) | awk '$NF == 1 {print $1 "\t" $2}' | "$SCRIPT_DIR/s2f" > blat_output/OTUs.fna

echo "Clustering terminé. Résultats disponibles dans le dossier blat_output/"

# Optionnel : nettoyage
rm -rf tmp_chunks tmp_chunks_grouped

echo "Statistiques finales :"
echo "- Nombre de séquences d'entrée : $(grep -c '^>' "$INPUT")"
echo "- Nombre de chimères détectées : $(wc -l < blat_output/${BASENAME}.chimeras.list)"
echo "- Nombre d'OTUs finaux : $(wc -l < blat_output/OTUs.tsv)"