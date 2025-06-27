import os
import re

# === Paramètres ===
annotation_file = "14Apr2025_itol_family_annotations.txt"
tree_file = "fasta36_cleaned.nwk"
MIN_TAXONS = 20  # seuil minimal de taxons pour écrire le fichier
output_dir = "families"  # dossier de sortie

# Création du dossier de sortie s'il n'existe pas
os.makedirs(output_dir, exist_ok=True)

# Lecture de l'arbre et extraction des labels (sans guillemets)
with open(tree_file, "r") as tree_in:
    tree_data = tree_in.read()
tree_data_clean = tree_data.replace('"', '').replace("'", '')

labels = re.findall(r'[\(,]([^:(),]+):', tree_data_clean)
tree_taxa = set(labels)
print(f"[✓] {len(tree_taxa)} taxons trouvés dans l’arbre.")

family_to_taxa = {}
inside_data = False

with open(annotation_file, "r") as f:
    for line in f:
        line = line.strip()
        if line == "DATA":
            inside_data = True
            continue
        if not inside_data or not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        taxon_id = parts[0]
        family = parts[2].strip()

        family_to_taxa.setdefault(family, []).append(taxon_id)

# Écriture des fichiers dans le dossier output_dir pour familles avec au moins MIN_TAXONS présents dans l’arbre
for family, taxon_list in family_to_taxa.items():
    present = [taxon for taxon in taxon_list if taxon in tree_taxa]
    if len(present) < MIN_TAXONS:
        print(f"[!] Famille '{family}' ignorée (seulement {len(present)} taxons dans l’arbre, < {MIN_TAXONS})")
        continue
    filename = os.path.join(output_dir, f"{family}.txt")
    with open(filename, "w") as out:
        out.write("\n".join(present) + "\n")
    print(f"[✓] {len(present)} taxons écrits dans {filename}")