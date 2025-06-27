#!/bin/bash

# Script pour convertir des alignements MMseqs2 en arbre phylogénétique
# Reproduit la pipeline: filtrage -> nettoyage -> hashsums -> tree_bray -> rapidnj
# Version pour utilisation native (sans Singularity)

set -e  # Arrêter en cas d'erreur

# Paramètres par défaut
INPUT_FILE="viral_analysis_results/mmseqs_alignments/mmseqs_all_vs_all.tab"
EVALUE_THRESHOLD=0.05
OUTPUT_PREFIX="mm_sensitive_vOTUs"

TEMP_DIR="/srv/scratch/givreex/tmp_mmseqs"
mkdir -p "$TEMP_DIR"

# Fonction d'aide
usage() {
    echo "Usage: $0 -i <fichier_mmseqs2> [-e <seuil_evalue>] [-o <prefixe_sortie>] [-h]"
    echo ""
    echo "Options:"
    echo "  -i  Fichier d'alignements MMseqs2 (obligatoire)"
    echo "  -e  Seuil d'e-value (défaut: 0.05)"
    echo "  -o  Préfixe pour les fichiers de sortie (défaut: vOTUs)"
    echo "  -h  Afficher cette aide"
    echo ""
    echo "Exemple: $0 -i alignments.tsv -e 0.01 -o mes_sequences"
    exit 1
}

# Nettoyage à la sortie
cleanup() {
    rm -rf "$TEMP_DIR"
}
trap cleanup EXIT

# Parsing des arguments
while getopts "i:e:o:h" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        e) EVALUE_THRESHOLD="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
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

# Vérifier la disponibilité des outils
# Chercher les outils dans l'ordre de priorité :
# 1. Répertoire courant
# 2. PATH système
# 3. Répertoires personnalisés (à ajuster selon votre installation)

MISSING_TOOLS=()

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
    # 3. Répertoires personnalisés (décommentez et ajustez si nécessaire)
    # elif [[ -x "/usr/local/bin/$tool" ]]; then
    #     tool_path="/usr/local/bin/$tool"
    # elif [[ -x "$HOME/bin/$tool" ]]; then
    #     tool_path="$HOME/bin/$tool"
    else
        MISSING_TOOLS+=("$tool")
        return 1
    fi
    
    echo "$tool_path"
    return 0
}

# Vérifier et définir les chemins des outils
echo "Vérification des outils requis..."

HASHSUMS_CMD=$(find_tool "hashsums") || true
TREE_BRAY_CMD=$(find_tool "tree_bray") || true
RAPIDNJ_CMD=$(find_tool "rapidnj") || true

# Vérifier si des outils sont manquants
if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    echo "Erreur: Les outils suivants sont manquants ou non exécutables:"
    for tool in "${MISSING_TOOLS[@]}"; do
        echo "  - $tool"
    done
    echo ""
    echo "Solutions possibles:"
    echo "1. Installez les outils manquants dans votre PATH"
    echo "2. Placez les exécutables dans le répertoire courant"
    echo "3. Modifiez le script pour pointer vers les bons chemins"
    echo ""
    echo "Pour rapidnj, vous pouvez l'installer avec:"
    echo "  - Ubuntu/Debian: sudo apt-get install rapidnj"
    echo "  - Conda: conda install -c bioconda rapidnj"
    echo "  - Depuis les sources: https://github.com/johnlees/rapidnj"
    exit 1
fi

echo "Outils trouvés:"
echo "  - hashsums: $HASHSUMS_CMD"
echo "  - tree_bray: $TREE_BRAY_CMD"
echo "  - rapidnj: $RAPIDNJ_CMD"
echo ""

echo "=== Pipeline MMseqs2 vers arbre phylogénétique ==="
echo "Fichier d'entrée: $INPUT_FILE"
echo "Seuil e-value: $EVALUE_THRESHOLD"
echo "Préfixe de sortie: $OUTPUT_PREFIX"
echo "Répertoire temporaire: $TEMP_DIR"
echo ""

# Étape 1: Filtrage et extraction des données pertinentes
echo "Étape 1: Filtrage des alignements (e-value <= $EVALUE_THRESHOLD)..."
FILTERED_FILE="$TEMP_DIR/filtered.tsv"

# Adapter selon le format MMseqs2 standard (12 colonnes)
# Colonnes: query target identity alnlen mismatches gaps qstart qend tstart tend evalue bitscore
# On filtre sur l'e-value (colonne 11) et on extrait query, target, bitscore (colonnes 1, 2, 12)
awk -v threshold="$EVALUE_THRESHOLD" '
    NF >= 12 {
        e = $11 + 0  # conversion explicite vers numérique
        if (e <= threshold) {
            print $1 "\t" $2 "\t" $12
        }
    }
' "$INPUT_FILE" > "$FILTERED_FILE"

if [[ ! -s "$FILTERED_FILE" ]]; then
    echo "Erreur: Aucun alignement ne passe le filtre e-value <= $EVALUE_THRESHOLD"
    exit 1
fi

echo "  $(wc -l < "$FILTERED_FILE") alignements retenus"

# Étape 2: Nettoyage des noms de séquences
echo "Étape 2: Nettoyage des noms de séquences..."
CLEANED_FILE="$TEMP_DIR/cleaned.tsv"

# Équivalent de: rev | sed 's/\t[[:digit:]]\+_/\t/' | rev | sed 's/_[[:digit:]]\+\t/\t/'
# Cela enlève les suffixes numériques des noms de séquences
sed 's/_[0-9]\+\t/\t/g' "$FILTERED_FILE" | sed 's/\t_[0-9]\+/\t/g' > "$CLEANED_FILE"

# Étape 3: Tri des données
echo "Étape 3: Tri des données..."
SORTED_FILE="$TEMP_DIR/sorted.tsv"
sort "$CLEANED_FILE" > "$SORTED_FILE"

# Étape 4: Calcul des sommes (hashsums)
echo "Étape 4: Calcul des sommes par paire..."
HASHSUMS_FILE="$TEMP_DIR/hashsums.tsv"
"$HASHSUMS_CMD" < "$SORTED_FILE" > "$HASHSUMS_FILE"

if [[ ! -s "$HASHSUMS_FILE" ]]; then
    echo "Erreur: Échec du calcul des hashsums"
    exit 1
fi

echo "  $(wc -l < "$HASHSUMS_FILE") paires uniques"

# Étape 5: Calcul de la matrice de distances (tree_bray)
echo "Étape 5: Calcul de la matrice de distances de Bray-Curtis..."
MATRIX_FILE="${OUTPUT_PREFIX}.mat"
"$TREE_BRAY_CMD" < "$HASHSUMS_FILE" > "$MATRIX_FILE"

if [[ ! -s "$MATRIX_FILE" ]]; then
    echo "Erreur: Échec du calcul de la matrice de distances"
    exit 1
fi

# Extraire le nombre de séquences de la matrice
NUM_SEQUENCES=$(head -n1 "$MATRIX_FILE" | cut -f2)
echo "  Matrice ${NUM_SEQUENCES}x${NUM_SEQUENCES} générée: $MATRIX_FILE"

# Étape 6: Construction de l'arbre phylogénétique
echo "Étape 6: Construction de l'arbre phylogénétique..."
TREE_FILE="${OUTPUT_PREFIX}.nwk"

if "$RAPIDNJ_CMD" -i pd "$MATRIX_FILE" > "$TREE_FILE"; then
    echo "  Arbre phylogénétique généré: $TREE_FILE"
else
    echo "Erreur: Échec de la construction de l'arbre avec rapidnj"
    exit 1
fi

# Résumé final
echo ""
echo "=== Pipeline terminée avec succès ==="
echo "Fichiers générés:"
echo "  - Matrice de distances: $MATRIX_FILE"
echo "  - Arbre phylogénétique: $TREE_FILE"
echo "  - Nombre de séquences: $NUM_SEQUENCES"
echo ""

# Optionnel: afficher quelques statistiques
if command -v file &> /dev/null; then
    echo "Informations sur les fichiers:"
    echo "  $(file "$MATRIX_FILE")"
    echo "  $(file "$TREE_FILE")"
fi

echo "Pipeline MMseqs2 -> arbre phylogénétique terminée!"