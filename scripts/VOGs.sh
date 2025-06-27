#!/bin/bash

# Script de traitement des vOTUs
# Auteur: Script généré pour le traitement des fichiers vOTUs.faa et vOTUs.fasta36

set -e  # Arrêter le script en cas d'erreur

# Vérification de la présence des fichiers requis
echo "Vérification des fichiers d'entrée..."
for file in votu_analysis/vOTUs.faa full_results.tab; do
    if [[ ! -f "$file" ]]; then
        echo "Erreur: Le fichier $file n'existe pas"
        exit 1
    fi
done

# Vérification de la présence des scripts requis
echo "Vérification des scripts requis..."
for script in f2s seqlengths joincol mcl; do
    if ! command -v "$script" &> /dev/null; then
        echo "Erreur: Le script/commande $script n'est pas disponible dans le PATH"
        exit 1
    fi
done

echo "Tous les fichiers et scripts requis sont présents."
echo

# Étape 1: Génération des longueurs de séquences
echo "Étape 1: Génération du fichier de longueurs de séquences..."
cat otu_analysis/vOTUs.faa | f2s | seqlengths > vOTUs.faa.lengths

if [[ ! -s vOTUs.faa.lengths ]]; then
    echo "Erreur: Le fichier vOTUs.faa.lengths n'a pas été créé ou est vide"
    exit 1
fi

echo "Fichier vOTUs.faa.lengths créé avec succès ($(wc -l < vOTUs.faa.lengths) lignes)"

# Étape 2: Traitement principal avec pipeline
echo "Étape 2: Traitement principal avec pipeline de filtrage..."
cat full_results.tab | \
joincol vOTUs.faa.lengths | \
joincol vOTUs.faa.lengths 2 | \
awk '{
    if ($13 < $14) {
        s = ($3*$4)/($13*100) - 0.75;
        if (s < 0) {
            s = 0
        }
    } else {
        s = 0
    }
    print $1 "\t" $2 "\t" $11 "\t" $13/$14 "\t" ($8-$7)/(2*$13)+($10-$9)/(2*$14) "\t" ($7+$8-$13)/$13-($9+$10-$14)/$14 "\t" s
}' | \
awk '{if ($3 <= 0.05) print}' | \
awk '{if ($5 >= 0.4) print}' | \
awk '{
    if (sqrt(log($4)^2) - (-0.0181/($5-0.32)+0.23) + sqrt($6^2) <= 0.15 + $7) 
        print $1 "\t" $2
}' | \
mcl - -o - --abc | \
awk '{
    j++; 
    for (i = 1; i <= NF; i++) {
        print $i "\t" j
    }
}' > vOTUs.VOGs.tsv

# Vérification du fichier de sortie
if [[ -s vOTUs.VOGs.tsv ]]; then
    echo "Succès: Fichier vOTUs.VOGs.tsv créé avec $(wc -l < vOTUs.VOGs.tsv) lignes"
    echo "Aperçu des premières lignes:"
    head -5 vOTUs.VOGs.tsv
else
    echo "Attention: Le fichier vOTUs.VOGs.tsv est vide ou n'a pas été créé"
    echo "Cela peut être normal si aucune séquence ne passe les filtres"
fi

echo
echo "Traitement terminé!"
echo "Fichiers générés:"
echo "  - vOTUs.faa.lengths"
echo "  - vOTUs.VOGs.tsv"