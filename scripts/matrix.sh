#!/bin/bash

# Script pour convertir des alignements MMseqs2 en arbre phylogénétique
# Reproduit la pipeline: filtrage -> nettoyage -> hashsums -> tree_bray -> rapidnj
# Version corrigée avec diagnostics détaillés

export TMPDIR="$PWD/tmp_sort"
mkdir -p "$TMPDIR"

set -e  # Arrêter en cas d'erreur

# Paramètres par défaut
INPUT_FILE="full_results.tab"
EVALUE_THRESHOLD=0.05
OUTPUT_PREFIX="fasta36_new"
TEMP_DIR=$(mktemp -d)
DEBUG=false

# Fonction d'aide
usage() {
    echo "Usage: $0 -i <fichier_mmseqs2> [-e <seuil_evalue>] [-o <prefixe_sortie>] [-d] [-h]"
    echo ""
    echo "Options:"
    echo "  -i  Fichier d'alignements MMseqs2 (obligatoire)"
    echo "  -e  Seuil d'e-value (défaut: 0.05)"
    echo "  -o  Préfixe pour les fichiers de sortie (défaut: diamond_vOTUs)"
    echo "  -d  Mode debug (affiche les fichiers intermédiaires)"
    echo "  -h  Afficher cette aide"
    echo ""
    echo "Exemple: $0 -i alignments.tsv -e 0.01 -o mes_sequences"
    exit 1
}

# Nettoyage à la sortie
cleanup() {
    if [[ "$DEBUG" == "false" ]]; then
        rm -rf "$TEMP_DIR"
    else
        echo "Mode debug: fichiers temporaires conservés dans $TEMP_DIR"
    fi
}
trap cleanup EXIT

# Parsing des arguments
while getopts "i:e:o:dh" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        e) EVALUE_THRESHOLD="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        d) DEBUG=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Vérifications
if [[ -z "$INPUT_FILE" ]]; then
    echo "Erreur: Fichier d'entrée requis (-i)"
    usage
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Erreur: Fichier $INPUT_FILE introuvable"
    exit 1
fi

# Fonction pour trouver un outil
find_tool() {
    local tool="$1"
    local tool_path=""
    
    # 1. Répertoire courant
    if [[ -x "./$tool" ]]; then
        tool_path="./$tool"
    # 2. PATH système
    elif command -v "$tool" &> /dev/null; then
        tool_path="$tool"
    else
        return 1
    fi
    
    echo "$tool_path"
    return 0
}


echo "=== Pipeline MMseqs2 vers arbre phylogénétique ==="
echo "Fichier d'entrée: $INPUT_FILE"
echo "Seuil e-value: $EVALUE_THRESHOLD"
echo "Préfixe de sortie: $OUTPUT_PREFIX"
echo "Répertoire temporaire: $TEMP_DIR"
echo ""

# Analyser le format du fichier d'entrée
echo "Analyse du fichier d'entrée..."
echo "Premières lignes du fichier:"
head -3 "$INPUT_FILE"
echo ""

NUM_COLS=$(head -1 "$INPUT_FILE" | awk '{print NF}')
TOTAL_LINES=$(wc -l < "$INPUT_FILE")
echo "Nombre de colonnes détectées: $NUM_COLS"
echo "Nombre total de lignes: $TOTAL_LINES"
echo ""

# Vérifier le format (MMseqs2 standard = 12 colonnes)
if [[ $NUM_COLS -lt 11 ]]; then
    echo "Attention: Le fichier semble avoir moins de 11 colonnes"
    echo "Format MMseqs2 attendu: query target identity alnlen mismatches gaps qstart qend tstart tend evalue bitscore"
    echo "Continuons quand même..."
fi

# Étape 1: Filtrage et extraction des données pertinentes
echo "Étape 1: Filtrage des alignements (e-value <= $EVALUE_THRESHOLD)..."
FILTERED_FILE="$TEMP_DIR/filtered.tsv"

# Déterminer la colonne d'e-value et de bitscore selon le nombre de colonnes
if [[ $NUM_COLS -ge 12 ]]; then
    # Format standard MMseqs2: e-value=col11, bitscore=col12
    EVALUE_COL=11
    BITSCORE_COL=12
elif [[ $NUM_COLS -ge 11 ]]; then
    # Format sans bitscore: e-value=col11, utiliser identity comme score
    EVALUE_COL=11
    BITSCORE_COL=3
else
    echo "Erreur: Format de fichier non reconnu (moins de 11 colonnes)"
    exit 1
fi

echo "Utilisation: e-value=colonne $EVALUE_COL, score=colonne $BITSCORE_COL"

# Filtrage adaptatif
awk -v threshold="$EVALUE_THRESHOLD" -v eval_col="$EVALUE_COL" -v score_col="$BITSCORE_COL" '
    NF >= eval_col && $eval_col <= threshold && $eval_col != "" && $score_col != "" {
        # Éviter les auto-alignements
        if ($1 != $2) {
            print $1 "\t" $2 "\t" $score_col
        }
    }
' "$INPUT_FILE" > "$FILTERED_FILE"

if [[ ! -s "$FILTERED_FILE" ]]; then
    echo "Erreur: Aucun alignement ne passe le filtre e-value <= $EVALUE_THRESHOLD"
    echo "Vérifiez:"
    echo "1. Le seuil d'e-value (essayez une valeur plus élevée, ex: -e 1.0)"
    echo "2. Le format de votre fichier d'alignement"
    echo "3. La colonne d'e-value utilisée"
    
    # Diagnostics supplémentaires
    echo ""
    echo "Diagnostic: valeurs d'e-value dans le fichier:"
    awk -v eval_col="$EVALUE_COL" 'NF >= eval_col {print $eval_col}' "$INPUT_FILE" | sort -n | head -10
    exit 1
fi

FILTERED_COUNT=$(wc -l < "$FILTERED_FILE")
echo "  $FILTERED_COUNT alignements retenus"

if [[ $DEBUG == "true" ]]; then
    echo "Exemple d'alignements filtrés:"
    head -5 "$FILTERED_FILE"
    echo ""
fi

# Étape 2: Nettoyage des noms de séquences
echo "Étape 2: Nettoyage des noms de séquences..."
CLEANED_FILE="$TEMP_DIR/cleaned.tsv"

# Nettoyage plus robuste des suffixes numériques
sed -E 's/_[0-9]+(\t|\s)/\1/g' "$FILTERED_FILE" > "$CLEANED_FILE"

if [[ $DEBUG == "true" ]]; then
    echo "Exemple après nettoyage:"
    head -5 "$CLEANED_FILE"
    echo ""
fi

# Étape 3: Tri des données
echo "Étape 3: Tri des données..."
SORTED_FILE="$TEMP_DIR/sorted.tsv"
sort -k1,1 -k2,2 "$CLEANED_FILE" > "$SORTED_FILE"

# Étape 4: Calcul des sommes (hashsums)
echo "Étape 4: Calcul des sommes par paire..."
HASHSUMS_FILE="$TEMP_DIR/hashsums.tsv"

# Vérifier que hashsums fonctionne
if ! "$HASHSUMS_CMD" < "$SORTED_FILE" > "$HASHSUMS_FILE" 2>/dev/null; then
    echo "Erreur: Échec du calcul des hashsums"
    echo "Vérifiez le format des données d'entrée pour hashsums"
    
    if [[ $DEBUG == "true" ]]; then
        echo "Données envoyées à hashsums:"
        head -5 "$SORTED_FILE"
    fi
    exit 1
fi

if [[ ! -s "$HASHSUMS_FILE" ]]; then
    echo "Erreur: hashsums n'a produit aucun résultat"
    exit 1
fi

HASHSUMS_COUNT=$(wc -l < "$HASHSUMS_FILE")
echo "  $HASHSUMS_COUNT paires uniques"

if [[ $DEBUG == "true" ]]; then
    echo "Exemple de hashsums:"
    head -5 "$HASHSUMS_FILE"
    echo ""
fi

echo "Contenu de $HASHSUMS_FILE :"
head -5 "$HASHSUMS_FILE"

# Étape 5: Calcul de la matrice de distances (tree_bray)
echo "Étape 5: Calcul de la matrice de distances de Bray-Curtis..."
MATRIX_FILE="${OUTPUT_PREFIX}.mat"

# Vérifier que tree_bray fonctionne
if ! "$TREE_BRAY_CMD" < "$HASHSUMS_FILE" > "$MATRIX_FILE" 2>/dev/null; then
    echo "Erreur: Échec du calcul de la matrice de distances"
    echo "Vérifiez le format des données d'entrée pour tree_bray"
    
    if [[ $DEBUG == "true" ]]; then
        echo "Données envoyées à tree_bray:"
        head -5 "$HASHSUMS_FILE"
    fi
    exit 1
fi

if [[ ! -s "$MATRIX_FILE" ]]; then
    echo "Erreur: tree_bray n'a produit aucun résultat"
    exit 1
fi

# Analyser la matrice produite
echo "Analyse de la matrice générée:"
MATRIX_LINES=$(wc -l < "$MATRIX_FILE")
echo "  Nombre de lignes dans la matrice: $MATRIX_LINES"

# Vérifier le format de la matrice
FIRST_LINE=$(head -n1 "$MATRIX_FILE")
echo "  Première ligne: $FIRST_LINE"

# Tenter d'extraire le nombre de séquences
if [[ "$FIRST_LINE" =~ ^[[:space:]]*([0-9]+) ]]; then
    NUM_SEQUENCES="${BASH_REMATCH[1]}"
    echo "  Nombre de séquences: $NUM_SEQUENCES"
else
    echo "  Attention: Format de matrice non standard"
    NUM_SEQUENCES="inconnu"
fi

if [[ $DEBUG == "true" ]]; then
    echo "Premières lignes de la matrice:"
    head -5 "$MATRIX_FILE"
    echo ""
fi

# Étape 6: Construction de l'arbre phylogénétique
echo "Étape 6: Construction de l'arbre phylogénétique..."
TREE_FILE="${OUTPUT_PREFIX}.nwk"

# Tester différents formats de matrice pour rapidnj
echo "Test de rapidnj avec format phylip distance (-i pd)..."
if "$RAPIDNJ_CMD" -i pd "$MATRIX_FILE" > "$TREE_FILE" 2>/dev/null; then
    echo "  Succès avec format phylip distance"
elif "$RAPIDNJ_CMD" -i dm "$MATRIX_FILE" > "$TREE_FILE" 2>/dev/null; then
    echo "  Succès avec format distance matrix"
elif "$RAPIDNJ_CMD" "$MATRIX_FILE" > "$TREE_FILE" 2>/dev/null; then
    echo "  Succès avec format par défaut"
else
    echo "Erreur: Échec de la construction de l'arbre avec tous les formats testés"
    echo ""
    echo "Diagnostic de la matrice:"
    echo "Taille du fichier: $(ls -lh "$MATRIX_FILE" | awk '{print $5}')"
    echo "Premières lignes:"
    head -3 "$MATRIX_FILE"
    echo ""
    echo "Essayez de vérifier manuellement le format de votre matrice"
    exit 1
fi

# Vérifier que l'arbre a été généré
if [[ ! -s "$TREE_FILE" ]]; then
    echo "Erreur: L'arbre phylogénétique est vide"
    exit 1
fi

echo "  Arbre phylogénétique généré: $TREE_FILE"

# Résumé final
echo ""
echo "=== Pipeline terminée avec succès ==="
echo "Fichiers générés:"
echo "  - Matrice de distances: $MATRIX_FILE ($(wc -l < "$MATRIX_FILE") lignes)"
echo "  - Arbre phylogénétique: $TREE_FILE ($(wc -c < "$TREE_FILE") caractères)"
echo "  - Nombre de séquences: $NUM_SEQUENCES"
echo ""

# Statistiques de l'arbre
if command -v grep &> /dev/null; then
    LEAF_COUNT=$(grep -o "(" "$TREE_FILE" | wc -l)
    echo "Statistiques de l'arbre:"
    echo "  - Nœuds internes estimés: $LEAF_COUNT"
fi

if [[ $DEBUG == "true" ]]; then
    echo ""
    echo "Fichiers temporaires conservés dans: $TEMP_DIR"
    echo "Contenu:"
    ls -la "$TEMP_DIR"
fi
