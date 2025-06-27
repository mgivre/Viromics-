#!/bin/bash

# Script pour analyser les séquences virales avec Prodigal et FASTA36
# Utilisation: ./viral_analysis.sh input.tsv blat_output/OTUs.fna
# Utilisation singularity : singularity exec --bind $(pwd):/mnt virome_shah.sif /mnt/pipeline3.sh benchmark.tsv /mnt/blat_output/OTUs.fna
# Le script vérifie à chaque étape si le traitement a déjà été effectué

# Vérification des arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <fichier_tsv> <fichier_fasta_supplémentaire>"
    exit 1
fi


# Répertoire de travail courant (où se trouvent les scripts hashsums et tree_bray)
CURRENT_DIR=$(pwd)
INPUT_FILE="$1"
EXTRA_FASTA="$2"
WORKDIR="viral_analysis_results"
VIRAL_SEQS="$WORKDIR/viral_sequences.fasta"
PRODIGAL_OUT="$WORKDIR/prodigal_proteins.faa"
FASTA36_OUT="$WORKDIR/fasta36_alignments"
LOG_FILE="$WORKDIR/pipeline.log"
CHECKPOINT_FILE="$WORKDIR/checkpoint.txt"
ALIGNMENT_FILE="$FASTA36_OUT/fasta36_all_vs_all.tab"



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
    for cmd in prodigal fasta36 rapidnj; do
        if ! command -v $cmd >/dev/null 2>&1; then
            echo "ERREUR: La commande '$cmd' n'est pas trouvée dans le PATH"
            echo "Veuillez installer ce programme ou vous assurer qu'il est disponible dans le conteneur Singularity"
            exit 1
        fi
    done
}

# Création du répertoire de travail
mkdir -p "$WORKDIR"
mkdir -p "$FASTA36_OUT"

# Vérifier la disponibilité des scripts et outils nécessaires
check_scripts

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

# Comptage des séquences dans le fichier TSV
TOTAL_SEQUENCES=$(wc -l < "$INPUT_FILE")
TOTAL_SEQUENCES=$((TOTAL_SEQUENCES - 1))  # Soustraire l'en-tête
echo "Nombre total de séquences dans le fichier TSV: $TOTAL_SEQUENCES"

# Comptage des séquences virales sans les extraire
VIRAL_COUNT=$(awk -F'\t' 'NR > 1 && $3 == 1 {count++} END {print count}' "$INPUT_FILE")
echo "Nombre de séquences virales détectées (colonne 3 == 1): $VIRAL_COUNT"
echo "Pourcentage viral: $(awk -v v=$VIRAL_COUNT -v t=$TOTAL_SEQUENCES 'BEGIN {printf "%.2f%%", (v/t)*100}')"

# Étape 1: Extraction des séquences virales + ajout du fichier FNA
if ! is_completed "extraction_completed"; then
    echo "=== Extraction des séquences virales et ajout des séquences OTU ==="

    # Extraction des séquences virales depuis le TSV
    awk -F'\t' 'NR > 1 && $3 == 1 { print ">"$1"\n"$NF }' "$INPUT_FILE" > "$VIRAL_SEQS"

    # Compter les séquences extraites du TSV
    TSV_SEQ_COUNT=$(grep -c "^>" "$VIRAL_SEQS")
    echo "Séquences virales extraites du TSV : $TSV_SEQ_COUNT"

    # Ajouter les séquences depuis le fichier FNA
    if [ -s "$EXTRA_FASTA" ]; then
        cat "$EXTRA_FASTA" >> "$VIRAL_SEQS"
        FNA_SEQ_COUNT=$(grep -c "^>" "$EXTRA_FASTA")
        echo "Séquences dans $EXTRA_FASTA : $FNA_SEQ_COUNT"
    else
        echo "ERREUR : Le fichier $EXTRA_FASTA est introuvable ou vide"
        exit 1
    fi

    # Compter le total combiné
    TOTAL_SEQ_COUNT=$(grep -c "^>" "$VIRAL_SEQS")
    echo "Nombre total de séquences (TSV + FNA) : $TOTAL_SEQ_COUNT"

    if [ "$TOTAL_SEQ_COUNT" -eq 0 ]; then
        echo "Aucune séquence détectée. Vérifiez vos fichiers d'entrée."
        exit 1
    fi

    mark_completed "extraction_completed"
else
    echo "=== Séquences virales déjà extraites, étape ignorée ==="
    TOTAL_SEQ_COUNT=$(grep -c "^>" "$VIRAL_SEQS")
    echo "Nombre total de séquences virales (précédemment) : $TOTAL_SEQ_COUNT"
fi

# Étape 2: Prédiction des protéines avec Prodigal
if ! is_completed "prodigal_completed"; then
    echo "=== Prédiction des protéines avec Prodigal ==="
    # Exécution de Prodigal pour prédire les protéines
    prodigal -i "$VIRAL_SEQS" -a "$PRODIGAL_OUT" -p meta 

    # Vérifier si Prodigal a généré des protéines
    if [ ! -s "$PRODIGAL_OUT" ]; then
        echo "Prodigal n'a pas pu prédire de protéines. Vérifiez vos séquences."
        exit 1
    fi

    PROTEIN_COUNT=$(grep -c "^>" "$PRODIGAL_OUT")
    echo "Nombre de protéines prédites: $PROTEIN_COUNT"
    
    mark_completed "prodigal_completed"
else
    echo "=== Prédiction des protéines déjà effectuée, étape ignorée ==="
    PROTEIN_COUNT=$(grep -c "^>" "$PRODIGAL_OUT")
    echo "Nombre de protéines prédites (précédemment): $PROTEIN_COUNT"
fi


#!/bin/bash

run_fasta36_parallel() {
    echo "=== Début de l'alignement FASTA36 ==="
    
    # Vérification des variables d'environnement
    if [[ -z "$PRODIGAL_OUT" || -z "$FASTA36_OUT" || -z "$WORKDIR" ]]; then
        echo "❌ Erreur: Variables d'environnement manquantes (PRODIGAL_OUT, FASTA36_OUT, WORKDIR)"
        return 1
    fi
    
    # Vérification des fichiers d'entrée
    if [[ ! -f "$PRODIGAL_OUT" ]]; then
        echo "❌ Erreur: Fichier d'entrée $PRODIGAL_OUT non trouvé"
        return 1
    fi
    
    # Créer la base de données si elle n'existe pas
    if [[ ! -f "$WORKDIR/db.faa" ]]; then
        echo "→ Copie du fichier de protéines vers db.faa"
        cp "$PRODIGAL_OUT" "$WORKDIR/db.faa"
    fi
    
    # Configuration
    MAX_JOBS=4  # Réduit pour éviter la surcharge
    
    # Préparation des répertoires
    echo "→ Préparation des répertoires"
    mkdir -p "$FASTA36_OUT"/{chunks,output,logs}
    
    # Nettoyage
    rm -f "$FASTA36_OUT/chunks"/*.faa "$FASTA36_OUT/output"/*.tab "$FASTA36_OUT/logs"/*.log 2>/dev/null
    
    # Découpage simple du fichier
    echo "→ Découpage du fichier en $MAX_JOBS parties"
    total_seqs=$(grep -c "^>" "$PRODIGAL_OUT")
    echo "  Total de séquences: $total_seqs"
    
    if [[ $total_seqs -eq 0 ]]; then
        echo "❌ Aucune séquence trouvée dans $PRODIGAL_OUT"
        return 1
    fi
    
    seqs_per_chunk=$(( (total_seqs + MAX_JOBS - 1) / MAX_JOBS ))
    echo "  Séquences par chunk: $seqs_per_chunk"
    
    # Découpage avec awk (plus fiable que seqkit)
    awk -v seqs_per_chunk="$seqs_per_chunk" -v output_dir="$FASTA36_OUT/chunks" '
    BEGIN { file_num = 1; seq_count = 0; output_file = "" }
    /^>/ { 
        if (seq_count >= seqs_per_chunk && seq_count > 0) {
            close(output_file)
            file_num++
            seq_count = 0
        }
        if (seq_count == 0) {
            output_file = output_dir "/chunk_" sprintf("%02d", file_num) ".faa"
        }
        seq_count++
    }
    { print > output_file }
    END { if (output_file) close(output_file) }
    ' "$PRODIGAL_OUT"
    
    # Vérification du découpage
    chunk_files=("$FASTA36_OUT/chunks"/*.faa)
    if [[ ! -f "${chunk_files[0]}" ]]; then
        echo "❌ Erreur: Aucun chunk créé"
        return 1
    fi
    
    chunk_count=${#chunk_files[@]}
    echo "  ✔️ $chunk_count chunks créés"
    
    # Test de fasta36
    if ! command -v fasta36 >/dev/null 2>&1; then
        echo "❌ fasta36 non trouvé dans le PATH"
        return 1
    fi
    
    # Lancement des alignements
    echo "→ Lancement des alignements ($(date))"
    
    # Fonction pour traiter un chunk
    process_chunk() {
        local chunk_file="$1"
        local base=$(basename "$chunk_file" .faa)
        local output_file="$FASTA36_OUT/output/${base}.tab"
        local log_file="$FASTA36_OUT/logs/${base}.log"
        
        echo "  Traitement de $base..." > "$log_file"
        
        # Commande fasta36 simplifiée
        if fasta36 -m 8 -E 1e-5 "$chunk_file" "$WORKDIR/db.faa" > "$output_file" 2>> "$log_file"; then
            if [[ -s "$output_file" ]]; then
                local lines=$(wc -l < "$output_file")
                echo "  ✔️ $base terminé: $lines alignements"
                return 0
            else
                echo "  ⚠️ $base: aucun alignement trouvé"
                return 2
            fi
        else
            local exit_code=$?
            echo "  ❌ $base: erreur (code $exit_code)"
            return $exit_code
        fi
    }
    
    # Export pour l'utilisation dans les sous-processus
    export -f process_chunk
    export FASTA36_OUT WORKDIR
    
    # Traitement parallèle simplifié
    for chunk_file in "${chunk_files[@]}"; do
        process_chunk "$chunk_file" &
        
        # Limiter le nombre de jobs parallèles
        if (( $(jobs -r | wc -l) >= MAX_JOBS )); then
            wait -n  # Attendre qu'un job se termine
        fi
    done
    
    # Attendre tous les jobs restants
    wait
    
    # Fusion des résultats
    echo "→ Fusion des résultats"
    final_output="$FASTA36_OUT/fasta36_all_vs_all.tab"
    > "$final_output"
    
    success_count=0
    for output_file in "$FASTA36_OUT/output"/*.tab; do
        if [[ -f "$output_file" && -s "$output_file" ]]; then
            cat "$output_file" >> "$final_output"
            ((success_count++))
        fi
    done
    
    # Statistiques finales
    final_lines=$(wc -l < "$final_output" 2>/dev/null || echo "0")
    echo "=== Résumé ==="
    echo "  Chunks traités avec succès: $success_count/$chunk_count"
    echo "  Alignements totaux: $final_lines"
    echo "  Fichier final: $final_output"
    echo "  Fin: $(date)"
    
    if [[ $final_lines -gt 0 ]]; then
        echo "✔️ Alignement terminé avec succès"
        return 0
    else
        echo "❌ Aucun alignement trouvé"
        return 1
    fi
}

# Étape 3: Alignement pairwise avec FASTA36
if ! is_completed "fasta36_completed"; then
    echo "=== Alignement pairwise avec FASTA36 ==="
    
    # Vérifier que la base de données existe
    if [[ ! -f "$WORKDIR/db.faa" ]]; then
        echo "→ Création de la base de données"
        cp "$PRODIGAL_OUT" "$WORKDIR/db.faa"
    fi
    
    if run_fasta36_parallel; then
        mark_completed "fasta36_completed"
        echo "✔️ Étape FASTA36 terminée"
    else
        echo "❌ Échec de l'étape FASTA36"
        exit 1
    fi
else
    echo "=== Alignement déjà effectué, étape ignorée ==="
fi


# Étape 4: Traitement des alignements et création de la matrice de distance
# Fichiers de sortie spécifiques à FASTA36
MATRIX_FILE="$WORKDIR/vOTUs.fasta36.mat"
if ! is_completed "fasta36_matrix_completed"; then
    echo "=== Traitement des alignements pour créer la matrice de distance ==="
    
    # Création d'un lien symbolique pour simplifier le nom du fichier d'alignement
    if [ ! -f "$ALIGNMENT_FILE" ]; then
        echo "ERREUR : le fichier d'alignement $ALIGNMENT_FILE est introuvable."
        exit 1
    fi

    
    # Filtrage des alignements et création de la matrice
    (cd "$WORKDIR" && \
    awk 'NF >= 12 && $11 <= 0.05 { print $1 "\t" $2 "\t" $12 }' vOTUs.fasta36 | \
    rev | sed 's/\t[[:digit:]]\+_/\t/' | rev | \
    sed 's/_[[:digit:]]\+\t/\t/' | sort | \
    "$CURRENT_DIR/hashsums" | "$CURRENT_DIR/tree_bray" > vOTUs.fasta36.mat)
    
    # Vérification de la création de la matrice
    if [ ! -s "$MATRIX_FILE" ]; then
        echo "Erreur lors de la création de la matrice de distance."
        exit 1
    fi
    
    echo "Matrice de distance créée: $MATRIX_FILE"
    mark_completed "fasta36_matrix_completed"
else
    echo "=== Matrice de distance déjà créée, étape ignorée ==="
fi

# Étape 5: Génération de l'arbre phylogénétique avec rapidnj
TREE_FILE="$WORKDIR/vOTUs.fasta36.nwk"
if ! is_completed "fasta36_tree_completed"; then
    echo "=== Génération de l'arbre phylogénétique avec rapidnj ==="
    
    # Exécution de rapidnj
    (cd "$WORKDIR" && rapidnj -i pd vOTUs.fasta36.mat > vOTUs.fasta36.nwk)
    
    # Vérification de la création de l'arbre
    if [ ! -s "$TREE_FILE" ]; then
        echo "Erreur lors de la génération de l'arbre phylogénétique."
        exit 1
    fi
    
    echo "Arbre phylogénétique créé: $TREE_FILE"
    mark_completed "fasta36_tree_completed"
else
    echo "=== Arbre phylogénétique déjà généré, étape ignorée ==="
fi

# Statistiques sur les alignements
echo "=== Statistiques des alignements ==="
if [ -s "$ALIGNMENT_FILE" ]; then
    ALIGNMENT_COUNT=$(wc -l < "$ALIGNMENT_FILE")
    echo "Nombre total d'alignements: $ALIGNMENT_COUNT"

    # Filtrer les auto-alignements
    NON_SELF_COUNT=$(awk '$1 != $2 { count++ } END { print count }' "$ALIGNMENT_FILE")
    echo "Nombre d'alignements (hors auto-alignements): $NON_SELF_COUNT"
    
    # Résumé des statistiques
    echo -e "\n=== Résumé ==="
    echo "Séquences totales dans le TSV: $TOTAL_SEQUENCES"
    echo "Séquences virales: $VIRAL_COUNT ($TSV_SEQ_COUNT extraites)"
    echo "Protéines prédites: $PROTEIN_COUNT"
    echo "Alignements identifiés: $ALIGNMENT_COUNT (dont $NON_SELF_COUNT hors auto-alignements)"
    
    # Vérification de l'existence de la matrice et de l'arbre
    if [ -f "$WORKDIR/vOTUs.fasta36.mat" ]; then
        echo "Matrice de distance générée: $(wc -l < "$WORKDIR/vOTUs.fasta36.mat") lignes"
    fi
    
    if [ -f "$WORKDIR/vOTUs.fasta36.nwk" ]; then
        echo "Arbre phylogénétique généré: $(stat -c %s "$WORKDIR/vOTUs.fasta36.nwk") octets"
        echo "Vous pouvez visualiser l'arbre avec un outil comme iTOL, FigTree ou en ligne sur https://itol.embl.de/"
    fi
    
    echo -e "\nPour visualiser les meilleurs alignements (E-value < 1e-10):"
    echo "awk '\$11 < 1e-10' $ALIGNMENT_FILE | sort -k11,11g | head -20"
else
    echo "Aucun fichier d'alignement trouvé ou fichier vide."
fi