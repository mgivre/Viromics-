#!/bin/bash

WORK_DIR="./processed"
LOGS_DIR="./logs"
TRIMMED_SAMPLES="$LOGS_DIR/trimmed_samples.txt"

# Créer le dossier de logs s'il n'existe pas
mkdir -p "$LOGS_DIR"

# Initialiser/Vider le fichier de suivi si l'option est spécifiée
if [[ "$1" == "--reset" ]]; then
    echo "Réinitialisation du fichier de suivi des échantillons nettoyés..."
    > "$TRIMMED_SAMPLES"
fi

# Créer le fichier de suivi s'il n'existe pas
[ -f "$TRIMMED_SAMPLES" ] || touch "$TRIMMED_SAMPLES"

echo "Analyse du dossier $WORK_DIR..."
echo "Recherche des fichiers de séquences nettoyées..."

# Compteurs
count_pe=0
count_se=0

# Recherche des fichiers paired-end nettoyés
for r1_file in "$WORK_DIR"/*_clean_1.fastq*; do
    # Vérifier si le fichier existe (au cas où le glob ne trouve rien)
    if [ -f "$r1_file" ]; then
        # Extraire le nom de l'échantillon
        sample=$(basename "$r1_file" | sed 's/_clean_1\.fastq.*//')
        
        # Vérifier si le R2 correspondant existe aussi
        r2_file="${r1_file/_clean_1.fastq/_clean_2.fastq}"
        if [ -f "$r2_file" ]; then
            # Vérifier si l'entrée existe déjà dans le fichier de suivi
            if ! grep -q "^${sample}_PE$" "$TRIMMED_SAMPLES"; then
                echo "${sample}_PE" >> "$TRIMMED_SAMPLES"
                echo "Paired-end détecté et ajouté: $sample"
                ((count_pe++))
            else
                echo "Paired-end déjà enregistré: $sample"
            fi
        else
            echo "ATTENTION: Fichier $r1_file trouvé mais pas son correspondant R2"
        fi
    fi
done

# Recherche des fichiers single-end nettoyés
for s_file in "$WORK_DIR"/*_clean_singles.fastq*; do
    # Vérifier si le fichier existe (au cas où le glob ne trouve rien)
    if [ -f "$s_file" ]; then
        # Extraire le nom de l'échantillon
        sample=$(basename "$s_file" | sed 's/_clean_singles\.fastq.*//')
        
        # Vérifier si l'entrée existe déjà dans le fichier de suivi
        if ! grep -q "^${sample}_SE$" "$TRIMMED_SAMPLES"; then
            echo "${sample}_SE" >> "$TRIMMED_SAMPLES"
            echo "Single-end détecté et ajouté: $sample"
            ((count_se++))
        else
            echo "Single-end déjà enregistré: $sample"
        fi
    fi
done

# Afficher un résumé
echo ""
echo "===== Récapitulatif ====="
echo "Nouveaux échantillons paired-end ajoutés: $count_pe"
echo "Nouveaux échantillons single-end ajoutés: $count_se"
echo "Total des entrées dans $TRIMMED_SAMPLES: $(wc -l < "$TRIMMED_SAMPLES")"
echo "======================="

echo "Analyse terminée!"