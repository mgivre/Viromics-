#!/bin/bash

INPUT_DIR="./PRJEB46943"       # Tes données brutes
WORK_DIR="./processed"         # Données nettoyées
ASSEMBLY_DIR="./assemblies"    # Résultats d'assemblage
LOGS_DIR="./logs"              # Dossier pour les logs

THREADS=4
mkdir -p "$WORK_DIR" "$ASSEMBLY_DIR" "$LOGS_DIR"

# Fichiers de suivi séparés pour le trimming et l'assemblage
TRIMMED_SAMPLES="$LOGS_DIR/trimmed_samples.txt"
ASSEMBLED_SAMPLES="$LOGS_DIR/assembled_samples.txt"

# Création des fichiers de suivi s'ils n'existent pas
[ -f "$TRIMMED_SAMPLES" ] || touch "$TRIMMED_SAMPLES"
[ -f "$ASSEMBLED_SAMPLES" ] || touch "$ASSEMBLED_SAMPLES"
[ -f "$LOGS_DIR/assembly_stats.tsv" ] || echo -e "sample\ttotal_contigs\tcontigs_>1kb" > "$LOGS_DIR/assembly_stats.tsv"

# Liste des préfixes uniques des paires
shopt -s nullglob
for file in "$INPUT_DIR"/*_1.fastq.gz; do
    sample=$(basename "$file" | sed 's/_1\.fastq\.gz//')
    echo "Sample détecté : $sample"

    # Vérifier si l'échantillon a déjà été assemblé
    if grep -q "^$sample$" "$ASSEMBLED_SAMPLES"; then
        echo "L'échantillon $sample a déjà été assemblé. Passage au suivant."
        continue
    fi

    # Vérifier si l'assemblage final existe déjà et est complet
    if [ -f "$ASSEMBLY_DIR/${sample}.assembly/contigs.fasta" ] && [ -f "$ASSEMBLY_DIR/${sample}.assembly/scaffolds.fasta" ]; then
        echo "L'assemblage de $sample existe déjà et semble complet. Passage au suivant."
        # Marquer comme assemblé dans le fichier de suivi
        echo "$sample" >> "$ASSEMBLED_SAMPLES"
        continue
    fi

    R1="$INPUT_DIR/${sample}_1.fastq.gz"
    R2="$INPUT_DIR/${sample}_2.fastq.gz"

    # Chercher les fichiers single-end
    SINGLES=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "${sample}*.fastq.gz" ! -name "*_1.fastq.gz" ! -name "*_2.fastq.gz")

    # Nettoyage paired-end
    CLEAN_R1="$WORK_DIR/${sample}_clean_1.fastq"
    CLEAN_R2="$WORK_DIR/${sample}_clean_2.fastq"
    CLEAN_S="$WORK_DIR/${sample}_clean_singles.fastq"

    # Vérifier si l'échantillon a déjà été nettoyé pour paired-end
    TRIMMING_NEEDED=true
    if grep -q "^${sample}_PE$" "$TRIMMED_SAMPLES"; then
        echo "Le nettoyage paired-end pour $sample a déjà été effectué."
        if [ -f "$CLEAN_R1" ] && [ -f "$CLEAN_R2" ]; then
            echo "Fichiers nettoyés trouvés, on saute l'étape de nettoyage paired-end"
            TRIMMING_NEEDED=false
        else
            echo "ATTENTION: Fichiers nettoyés non trouvés malgré l'indication dans le fichier de suivi!"
        fi
    fi

    # Effectuer le nettoyage paired-end si nécessaire
    if [ "$TRIMMING_NEEDED" = true ]; then
        echo "Nettoyage paired-end"
        fastp \
            --in1 "$R1" \
            --in2 "$R2" \
            --out1 "$CLEAN_R1" \
            --out2 "$CLEAN_R2" \
            --qualified_quality_phred 13 \
            --length_required 32 \
            --cut_right \
            --cut_right_mean_quality 13 \
            --adapter_sequence CTGTCTCTTATACACATCT \
            --thread $THREADS \
            --json "$LOGS_DIR/${sample}_fastp.json" \
            --html "$LOGS_DIR/${sample}_fastp.html"
        
        # Vérifier si le nettoyage a réussi
        if [ -f "$CLEAN_R1" ] && [ -f "$CLEAN_R2" ]; then
            echo "Nettoyage paired-end réussi pour $sample"
            echo "${sample}_PE" >> "$TRIMMED_SAMPLES"
        else
            echo "ERREUR: Le nettoyage paired-end a échoué pour $sample"
            continue
        fi
    fi

    # Nettoyage single-end si présent
    if [[ -n "$SINGLES" ]]; then
        # Vérifier si l'échantillon a déjà été nettoyé pour single-end
        SE_TRIMMING_NEEDED=true
        if grep -q "^${sample}_SE$" "$TRIMMED_SAMPLES"; then
            echo "Le nettoyage single-end pour $sample a déjà été effectué."
            if [ -f "$CLEAN_S" ]; then
                echo "Fichier single-end nettoyé trouvé, on saute l'étape de nettoyage single-end"
                SE_TRIMMING_NEEDED=false
            else
                echo "ATTENTION: Fichier single-end nettoyé non trouvé malgré l'indication dans le fichier de suivi!"
            fi
        fi

        # Effectuer le nettoyage single-end si nécessaire
        if [ "$SE_TRIMMING_NEEDED" = true ]; then
            echo "Nettoyage single-end détecté : $SINGLES"
            cat $SINGLES | zcat | \
            fastp \
                --stdin \
                --stdout \
                --qualified_quality_phred 13 \
                --length_required 32 \
                --cut_right \
                --cut_right_mean_quality 13 \
                --adapter_sequence CTGTCTCTTATACACATCT \
                --thread $THREADS \
                --json "$LOGS_DIR/${sample}_singles_fastp.json" \
                --html "$LOGS_DIR/${sample}_singles_fastp.html" > "$CLEAN_S"
            
            # Vérifier si le nettoyage a réussi
            if [ -f "$CLEAN_S" ]; then
                echo "Nettoyage single-end réussi pour $sample"
                echo "${sample}_SE" >> "$TRIMMED_SAMPLES"
            else
                echo "ERREUR: Le nettoyage single-end a échoué pour $sample"
            fi
        fi
        USE_S1="--s1 $CLEAN_S"
    else
        USE_S1=""  # Aucun single-end à inclure
    fi

    # Créer un dossier temporaire unique pour l'assemblage si le précédent a échoué
    if [ -d "$ASSEMBLY_DIR/${sample}.assembly" ]; then
        echo "Un assemblage incomplet existe déjà pour $sample, on le renomme"
        mv "$ASSEMBLY_DIR/${sample}.assembly" "$ASSEMBLY_DIR/${sample}.assembly.$(date +%Y%m%d_%H%M%S).bak"
    fi

    # Assemblage SPAdes avec données hybrides
    echo "Assemblage SPAdes"
    spades.py \
        --meta \
        --only-assembler \
        -1 "$CLEAN_R1" \
        -2 "$CLEAN_R2" \
        $USE_S1 \
        -o "$ASSEMBLY_DIR/${sample}.assembly" \
        -t $THREADS

    # Vérifier si l'assemblage a réussi
    if [ -f "$ASSEMBLY_DIR/${sample}.assembly/contigs.fasta" ] && [ -f "$ASSEMBLY_DIR/${sample}.assembly/scaffolds.fasta" ]; then
        
        echo "Assemblage de $sample terminé avec succès"
        # Marquer comme assemblé dans le fichier de suivi
        echo "$sample" >> "$ASSEMBLED_SAMPLES"

          # Compter le nombre total de contigs assemblés (avant filtrage)
        N_TOTAL_CONTIGS=$(grep -c "^>" "$ASSEMBLY_DIR/${sample}.assembly/contigs.fasta")
        echo "Nombre total de contigs pour $sample : $N_TOTAL_CONTIGS"

        # Filtrage des contigs > 1kb
        FILTERED_CONTIGS="$ASSEMBLY_DIR/${sample}.assembly/contigs.filtered_1kb.fasta"
        echo "Filtrage des contigs > 1kb pour $sample"
        awk '/^>/ {if (seqlen >= 1000) print header "\n" seq; header=$0; seq=""; seqlen=0; next} 
             {seqlen += length($0); seq = seq $0} 
             END {if (seqlen >= 1000) print header "\n" seq}' \
            "$ASSEMBLY_DIR/${sample}.assembly/contigs.fasta" > "$FILTERED_CONTIGS"

        # Compter les contigs > 1kb
        if [ -s "$FILTERED_CONTIGS" ]; then
            N_CONTIGS_1KB=$(grep -c "^>" "$FILTERED_CONTIGS")
            echo "Contigs > 1kb filtrés avec succès pour $sample"
            echo "Nombre de contigs > 1kb pour $sample : $N_CONTIGS_1KB"
        else
            N_CONTIGS_1KB=0
            echo "ATTENTION: Aucun contig > 1kb trouvé pour $sample"
        fi

        # Écrire les statistiques dans un fichier log
        echo -e "$sample\t$N_TOTAL_CONTIGS\t$N_CONTIGS_1KB" >> "$LOGS_DIR/assembly_stats.tsv"

    else
        echo "ATTENTION: L'assemblage de $sample semble avoir échoué"
    fi

    echo "Fini : $sample"
done
shopt -u nullglob

echo "Tous les échantillons ont été traités."