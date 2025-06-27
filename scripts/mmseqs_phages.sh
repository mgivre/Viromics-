#!/bin/bash

# Fichier: cluster_genomes_mmseqs.sh

# 1. Créer les dossiers nécessaires
mkdir -p mmseqs_input mmseqs_tmp mmseqs_output

# 2. Ne découper que si mmseqs_input est vide
if [ -z "$(ls -A mmseqs_input/*.fa 2>/dev/null)" ]; then
    echo "Découpage des génomes multi-fasta en fichiers individuels..."
    seqkit split -i -O mmseqs_input 14Apr2025_genomes.fa
else
    echo "Fichiers déjà présents dans mmseqs_input/, découpage sauté."
fi

# 3. Concaténer tous les .fa en un seul fichier pour MMseqs2
cat mmseqs_input/*.fa > mmseqs_input/all_genomes.fa

# 4. Déréplication avec MMseqs2
echo "Déréplication avec MMseqs2 à 95% d'identité..."
mmseqs createdb mmseqs_input/all_genomes.fa mmseqs_output/genomesDB
#mmseqs cluster mmseqs_output/genomesDB mmseqs_output/clusters mmseqs_tmp --min-seq-id 0.95 -c 0.9 --threads 4
mmseqs cluster mmseqs_output/genomesDB mmseqs_output/clusters mmseqs_tmp \
  --min-seq-id 0.95 -c 0.8 --cov-mode 1 --threads 8
mmseqs createtsv mmseqs_output/genomesDB mmseqs_output/genomesDB mmseqs_output/clusters mmseqs_output/cluster_results.tsv

# 5. Extraire les séquences représentatives
echo "Extraction des génomes représentatifs..."
mmseqs result2repseq mmseqs_output/genomesDB mmseqs_output/clusters mmseqs_output/repseqDB
mmseqs convert2fasta mmseqs_output/repseqDB mmseqs_output/OTUs_mmseqs.fna

echo "Déréplication terminée. Résultats dans mmseqs_output/"
