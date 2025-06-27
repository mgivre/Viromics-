#!/bin/bash

# Fichier: cluster_genomes_vclust.sh

# 1. Créer les dossiers nécessaires
mkdir -p vclust_input vclust_output

# 2. Ne découper que si vclust_input est vide
if [ -z "$(ls -A vclust_input/*.fa 2>/dev/null)" ]; then
    echo "Découpage des génomes multi-fasta en fichiers individuels..."
    seqkit split -i -O vclust_input 14Apr2025_genomes.fa
else
    echo "Fichiers déjà présents dans vclust_input/, découpage sauté."
fi

# 3. Concaténer tous les .fa en un seul fichier pour vclust
cat vclust_input/*.fa > vclust_input/all_genomes.fa

# 4. Clustering avec vclust
echo "Clustering avec vclust à 95% d'identité..."

# Paramètres classiques (à ajuster selon ta config)
IDENTITY=0.95
THREADS=8

vclust -i vclust_input/all_genomes.fa \
       -o vclust_output \
       --min-identity $IDENTITY \
       --threads $THREADS

# 5. Extraire les séquences représentatives (vclust génère généralement un fichier OTU)

REPSEQ_FILE="vclust_output/representatives.fasta"

if [ -f "$REPSEQ_FILE" ]; then
    echo "Extraction des génomes représentatifs terminée. Fichier : $REPSEQ_FILE"
else
    echo "Erreur : fichier des représentants non trouvé dans vclust_output/"
    exit 1
fi

echo "Déréplication terminée. Résultats dans vclust_output/"
