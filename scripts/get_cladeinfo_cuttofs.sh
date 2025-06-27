#!/bin/bash

# Dossier contenant les fichiers familles
FAMILY_DIR="families"
# Fichier de l'arbre
TREE="fasta36_rooted_3_clean.nwk"
# Fichier résultat
OUTFILE="cladeinfo_fasta36_3_results.tsv"

# En-tête du fichier résultat
echo -e "Family\tCutoff" > "$OUTFILE"

# Boucle sur chaque fichier txt dans FAMILY_DIR
for cladefile in "$FAMILY_DIR"/*.txt; do
    family=$(basename "$cladefile" .txt)
    echo "Processing family: $family"

    # Exécution de treetool --cladeinfo
    # Ici on suppose que la sortie contient une ligne avec "cutoff" et un nombre, à adapter si besoin
    cutoff=$(treetool --cladeinfo="$cladefile" "$TREE") #| grep -i "cutoff" | head -n1 | grep -oE "[0-9.]+")
    
    if [ -z "$cutoff" ]; then
        echo -e "$family\tNO_CUTOFF_FOUND" >> "$OUTFILE"
        echo "Warning: cutoff not found for $family"
    else
        echo -e "$family\t$cutoff" >> "$OUTFILE"
    fi
done

echo "All done. Results saved in $OUTFILE"
