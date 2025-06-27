#!/bin/bash

# Fichier: cluster_genomes_drep.sh

# 1. Créer les dossiers s'ils n'existent pas
mkdir -p drep_input
mkdir -p drep_output

# 2. Ne découper que si drep_input est vide
if [ -z "$(ls -A drep_input/*.fa 2>/dev/null)" ]; then
    echo "Découpage des génomes multi-fasta en fichiers individuels..."
    seqkit split -i -O drep_input 14Apr2025_genomes.fa
else
    echo "Fichiers déjà présents dans drep_input/, découpage sauté."
fi

# 3. Lancer (ou relancer) dRep
echo "Lancement de dRep..."
dRep dereplicate drep_output -g drep_input/*.fa -pa 0.95 -sa 0.95 -p 30
