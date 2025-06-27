#!/bin/bash

# Script pour continuer l'analyse virale en remplaçant FASTA36 par MMseqs2
# Ce script utilise les fichiers déjà générés par le pipeline original
# Utilisation: ./viral_analysis_mmseqs2_continue.sh
# Utilisation singularity : singularity exec --bind $(pwd):/mnt virome_shah.sif /mnt/viral_analysis_mmseqs2_continue.sh

# Répertoire de travail courant (où se trouvent les scripts hashsums et tree_bray)
CURRENT_DIR=$(pwd)
WORKDIR="viral_analysis_results"
PRODIGAL_OUT="$WORKDIR/prodigal_proteins.faa"
MMSEQS_TMPDIR="$WORKDIR/mmseqs_tmp"
MMSEQS_DB="$WORKDIR/mmseqs_db"
MMSEQS_RESULT="$WORKDIR/mmseqs_result"
MMSEQS_OUT="$WORKDIR/mmseqs_alignments"
LOG_FILE="$WORKDIR/pipeline.log"
CHECKPOINT_FILE="$WORKDIR/checkpoint.txt"
ALIGNMENT_FILE="$MMSEQS_OUT/mmseqs_all_vs_all.tab"

# Vérifier si les scripts requis sont disponibles
check_scripts() {
    # Vérifier les scripts locaux
    if [ ! -x "$CURRENT_DIR/hashsums" ]; then
        echo "ERREUR: Le script 'hashsums' n'est pas trouvé ou n'est pas exécutable dans $CURRENT_DIR"
        echo "Assurez-vous que le script existe et est exécutable (chmod +x hashsums)"
        exit 1
    fi
    
    if [ ! -x "$CURRENT_DIR/tree_bray" ]; then
        echo "ERREUR: Le script 'tree_bray' n'est pas trouvé ou n'est pas exécutable dans $CURRENT_DIR"
        echo "Assurez-vous que le script existe et est exécutable (chmod +x tree_bray)"
        exit 1
    fi
    
    # Vérifier les outils système
    for cmd in mmseqs rapidnj; do
        if ! command -v $cmd >/dev/null 2>&1; then
            echo "ERREUR: La commande '$cmd' n'est pas trouvée dans le PATH"
            echo "Veuillez installer ce programme ou vous assurer qu'il est disponible dans le conteneur Singularity"
            exit 1
        fi
    done
}

# Vérifier que le répertoire de travail et les fichiers requis existent
check_prerequisites() {
    if [ ! -d "$WORKDIR" ]; then
        echo "ERREUR: Le répertoire $WORKDIR n'existe pas."
        echo "Veuillez d'abord exécuter le pipeline original jusqu'à l'étape Prodigal."
        exit 1
    fi
    
    if [ ! -s "$PRODIGAL_OUT" ]; then
        echo "ERREUR: Le fichier $PRODIGAL_OUT n'existe pas ou est vide."
        echo "Veuillez d'abord exécuter le pipeline original jusqu'à l'étape Prodigal."
        exit 1
    fi
    
    PROTEIN_COUNT=$(grep -c "^>" "$PRODIGAL_OUT")
    echo "Fichier Prodigal détecté: $PROTEIN_COUNT protéines"
}

# Création des répertoires de travail pour MMseqs2
mkdir -p "$MMSEQS_OUT"
mkdir -p "$MMSEQS_TMPDIR"

# Vérifier la disponibilité des scripts et outils nécessaires
check_scripts

# Vérifier les prérequis (fichiers existants du pipeline original)
check_prerequisites

# Fonction pour enregistrer les étapes complétées
mark_completed() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$LOG_FILE"
    echo "$1" >> "$CHECKPOINT_FILE"
}

# Fonction pour vérifier si une étape est déjà complétée
is_completed() {
    if [ -f "$CHECKPOINT_FILE" ]; then
        grep -q "^$1$" "$CHECKPOINT_FILE" && return 0
    fi
    return 1
}

echo "=== Continuation du pipeline viral avec MMseqs2 ==="
echo "Répertoire de travail: $WORKDIR"
echo "Fichier de protéines: $PRODIGAL_OUT ($PROTEIN_COUNT protéines)"

# Étape 3: Création de la base de données MMseqs2 et alignement all-vs-all
if ! is_completed "mmseqs_completed"; then
    echo "=== Alignement all-vs-all avec MMseqs2 ==="
    
    # Création de la base de données MMseqs2
    echo "Création de la base de données MMseqs2..."
    mmseqs createdb "$PRODIGAL_OUT" "$MMSEQS_DB"
    
    # Recherche all-vs-all avec MMseqs2
    echo "Lancement de la recherche all-vs-all..."
    mmseqs search "$MMSEQS_DB" "$MMSEQS_DB" "$MMSEQS_RESULT" "$MMSEQS_TMPDIR" \
        --threads 0 \
        --sensitivity 7.5 \
        -e 1e-5 \
        --max-seqs 1000 \
        --min-seq-id 0.3 \
        --cov-mode 1 \
        -c 0.8
    
    # Conversion en format tabulaire
    echo "Conversion des résultats en format tabulaire..."
    mmseqs convertalis "$MMSEQS_DB" "$MMSEQS_DB" "$MMSEQS_RESULT" "$ALIGNMENT_FILE" \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    
    # Vérifier si MMseqs2 a généré des alignements
    if [ ! -s "$ALIGNMENT_FILE" ]; then
        echo "MMseqs2 n'a pas généré d'alignements. Vérifiez vos paramètres."
        exit 1
    fi
    
    ALIGNMENT_COUNT=$(wc -l < "$ALIGNMENT_FILE")
    echo "Nombre d'alignements générés: $ALIGNMENT_COUNT"
    
    mark_completed "mmseqs_completed"
else
    echo "=== Alignement MMseqs2 déjà effectué, étape ignorée ==="
    ALIGNMENT_COUNT=$(wc -l < "$ALIGNMENT_FILE")
    echo "Nombre d'alignements (précédemment): $ALIGNMENT_COUNT"
fi

# Étape 4: Traitement des alignements et création de la matrice de distance
MATRIX_FILE="$WORKDIR/vOTUs.mat"
if ! is_completed "matrix_completed"; then
    echo "=== Traitement des alignements pour créer la matrice de distance ==="
    
    # Vérification de l'existence du fichier d'alignement
    if [ ! -f "$ALIGNMENT_FILE" ]; then
        echo "ERREUR : le fichier d'alignement $ALIGNMENT_FILE est introuvable."
        exit 1
    fi

    # Vérifier le chemin complet du fichier d'alignement
    ALIGNMENT_FULL_PATH=$(readlink -f "$ALIGNMENT_FILE")
    echo "Fichier d'alignement: $ALIGNMENT_FULL_PATH"
    
    # Création d'un lien symbolique avec chemin absolu
    (cd "$WORKDIR" && ln -sf "$ALIGNMENT_FULL_PATH" vOTUs.mmseqs)
    
    # Vérifier que le lien symbolique a été créé
    if [ ! -f "$WORKDIR/vOTUs.mmseqs" ]; then
        echo "ERREUR: Impossible de créer le lien symbolique vOTUs.mmseqs"
        exit 1
    fi
    
    echo "Filtrage des alignements (E-value <= 0.05)..."
    # Compter les alignements avant filtrage
    TOTAL_ALIGN=$(wc -l < "$ALIGNMENT_FILE")
    echo "Alignements totaux: $TOTAL_ALIGN"
    
    # Filtrage des alignements et création de la matrice
    # Adaptation du filtrage pour le format MMseqs2 (colonne 11 = E-value, colonne 12 = bit score)
    (cd "$WORKDIR" && \
    awk 'NF >= 12 && $11 <= 0.05 { print $1 "\t" $2 "\t" $12 }' vOTUs.mmseqs | \
    rev | sed 's/\t[[:digit:]]\+_/\t/' | rev | \
    sed 's/_[[:digit:]]\+\t/\t/' | sort | \
    "$CURRENT_DIR/hashsums" | "$CURRENT_DIR/tree_bray" > vOTUs.mat)
    
    # Compter les alignements après filtrage
    FILTERED_ALIGN=$(awk 'NF >= 12 && $11 <= 0.05' "$WORKDIR/vOTUs.mmseqs" | wc -l)
    echo "Alignements filtrés (E-value <= 0.05): $FILTERED_ALIGN"
    
    # Vérification de la création de la matrice
    if [ ! -s "$MATRIX_FILE" ]; then
        echo "Erreur lors de la création de la matrice de distance."
        echo "Vérification des étapes de filtrage..."
        
        # Diagnostic des problèmes possibles
        echo "Contenu des 5 premières lignes du fichier d'alignement:"
        head -5 "$WORKDIR/vOTUs.mmseqs"
        
        echo "Nombre de colonnes dans la première ligne:"
        head -1 "$WORKDIR/vOTUs.mmseqs" | awk '{print NF}'
        
        echo "Alignements avec E-value <= 0.05:"
        awk 'NF >= 12 && $11 <= 0.05 {count++} END {print count+0}' "$WORKDIR/vOTUs.mmseqs"
        
        exit 1
    fi
    
    MATRIX_LINES=$(wc -l < "$MATRIX_FILE")
    echo "Matrice de distance créée: $MATRIX_FILE ($MATRIX_LINES lignes)"
    mark_completed "matrix_completed"
else
    echo "=== Matrice de distance déjà créée, étape ignorée ==="
fi

# Étape 5: Génération de l'arbre phylogénétique avec rapidnj
TREE_FILE="$WORKDIR/vOTUs.nwk"
if ! is_completed "tree_completed"; then
    echo "=== Génération de l'arbre phylogénétique avec rapidnj ==="
    
    # Exécution de rapidnj
    (cd "$WORKDIR" && rapidnj -i pd vOTUs.mat > vOTUs.nwk)
    
    # Vérification de la création de l'arbre
    if [ ! -s "$TREE_FILE" ]; then
        echo "Erreur lors de la génération de l'arbre phylogénétique."
        exit 1
    fi
    
    echo "Arbre phylogénétique créé: $TREE_FILE"
    mark_completed "tree_completed"
else
    echo "=== Arbre phylogénétique déjà généré, étape ignorée ==="
fi

# Nettoyage des fichiers temporaires MMseqs2
if [ -d "$MMSEQS_TMPDIR" ]; then
    echo "=== Nettoyage des fichiers temporaires MMseqs2 ==="
    rm -rf "$MMSEQS_TMPDIR"
fi

# Comparaison avec les anciens résultats FASTA36 (si disponibles)
FASTA36_FILE="$WORKDIR/fasta36_alignments/fasta36_all_vs_all.tab"
if [ -f "$FASTA36_FILE" ]; then
    echo "=== Comparaison avec les résultats FASTA36 ==="
    FASTA36_COUNT=$(wc -l < "$FASTA36_FILE")
    MMSEQS_COUNT=$(wc -l < "$ALIGNMENT_FILE")
    echo "Alignements FASTA36 : $FASTA36_COUNT"
    echo "Alignements MMseqs2 : $MMSEQS_COUNT"
    echo "Ratio MMseqs2/FASTA36 : $(awk -v m=$MMSEQS_COUNT -v f=$FASTA36_COUNT 'BEGIN {printf "%.2f", m/f}')"
fi

# Statistiques sur les alignements
echo "=== Statistiques des alignements MMseqs2 ==="
if [ -s "$ALIGNMENT_FILE" ]; then
    ALIGNMENT_COUNT=$(wc -l < "$ALIGNMENT_FILE")
    echo "Nombre total d'alignements: $ALIGNMENT_COUNT"

    # Filtrer les auto-alignements
    NON_SELF_COUNT=$(awk '$1 != $2 { count++ } END { print count+0 }' "$ALIGNMENT_FILE")
    echo "Nombre d'alignements (hors auto-alignements): $NON_SELF_COUNT"
    
    # Statistiques sur les E-values
    echo "Distribution des E-values:"
    awk '$11 < 1e-10 { count_very_good++ } 
         $11 >= 1e-10 && $11 < 1e-5 { count_good++ } 
         $11 >= 1e-5 && $11 < 0.01 { count_moderate++ } 
         $11 >= 0.01 { count_weak++ } 
         END { 
             print "  E-value < 1e-10 (très significatif): " count_very_good+0
             print "  1e-10 ≤ E-value < 1e-5 (significatif): " count_good+0  
             print "  1e-5 ≤ E-value < 0.01 (modéré): " count_moderate+0
             print "  E-value ≥ 0.01 (faible): " count_weak+0
         }' "$ALIGNMENT_FILE"
    
    # Statistiques sur l'identité
    echo "Distribution de l'identité de séquence:"
    awk '$3 >= 90 { count_very_high++ }
         $3 >= 70 && $3 < 90 { count_high++ }
         $3 >= 50 && $3 < 70 { count_medium++ }
         $3 < 50 { count_low++ }
         END {
             print "  Identité ≥ 90% (très élevée): " count_very_high+0
             print "  70% ≤ Identité < 90% (élevée): " count_high+0
             print "  50% ≤ Identité < 70% (moyenne): " count_medium+0
             print "  Identité < 50% (faible): " count_low+0
         }' "$ALIGNMENT_FILE"
    
    # Résumé final
    echo -e "\n=== Résumé final ==="
    echo "Protéines analysées: $PROTEIN_COUNT"
    echo "Alignements MMseqs2: $ALIGNMENT_COUNT (dont $NON_SELF_COUNT hors auto-alignements)"
    
    # Vérification de l'existence de la matrice et de l'arbre
    if [ -f "$WORKDIR/vOTUs.mat" ]; then
        MATRIX_LINES=$(wc -l < "$WORKDIR/vOTUs.mat")
        echo "Matrice de distance générée: $MATRIX_LINES lignes"
    fi
    
    if [ -f "$WORKDIR/vOTUs.nwk" ]; then
        TREE_SIZE=$(stat -c %s "$WORKDIR/vOTUs.nwk")
        echo "Arbre phylogénétique généré: $TREE_SIZE octets"
        echo "Vous pouvez visualiser l'arbre avec un outil comme iTOL, FigTree ou en ligne sur https://itol.embl.de/"
    fi
    
    echo -e "\n=== Commandes utiles pour analyser les résultats ==="
    echo "Meilleurs alignements (E-value < 1e-10):"
    echo "awk '\$11 < 1e-10' $ALIGNMENT_FILE | sort -k11,11g | head -20"
    
    echo -e "\nAlignements par identité (>90%):"
    echo "awk '\$3 > 90' $ALIGNMENT_FILE | sort -k3,3nr | head -20"
    
    echo -e "\nAlignements les plus longs:"
    echo "awk '\$4 > 100' $ALIGNMENT_FILE | sort -k4,4nr | head -20"
    
else
    echo "Aucun fichier d'alignement trouvé ou fichier vide."
fi

echo -e "\n=== Pipeline MMseqs2 terminé avec succès ==="