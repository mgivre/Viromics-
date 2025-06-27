#!/bin/bash

output_dir="output_analysis"
output_file="viral_sequence_ids.txt"

# Vider le fichier de sortie s'il existe
> "$output_file"

for genomad_dir in "$output_dir"/seq*_genomad; do
    if [ -d "$genomad_dir" ]; then
        seq_id=$(basename "$genomad_dir" | cut -d'_' -f1)
        virus_fna="$genomad_dir/${seq_id}_summary/${seq_id}_virus.fna"
        
        if [ -f "$virus_fna" ] && [ -s "$virus_fna" ]; then
            echo "$seq_id" >> "$output_file"
        fi
    fi
done

echo "[DONE] IDs des séquences virales extraits dans $output_file"
