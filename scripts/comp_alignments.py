import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3
from itertools import combinations

# === 1. Définir les fichiers d’alignement ===

alignment_files = {
    "diamond": "diamond_all_vs_all_12cols.tsv",
    "mmseqs2": "vOTUs_alignment.tsv"
}

# === 2. Charger les données d’alignement ===

def load_alignments(path):
    df = pd.read_csv(path, sep="\t", header=None)
    df.columns = [
        "query", "target", "identity", "length",
        "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"
    ]
    df["pair"] = df["query"] + "||" + df["target"]
    return df

alignments = {name: load_alignments(path) for name, path in alignment_files.items()}

# === 3. Statistiques de base ===

print("\n--- Statistiques de base ---")
for name, df in alignments.items():
    print(f"{name}: {len(df)} alignments, "
          f"mean identity = {df['identity'].mean():.2f}, "
          f"mean e-value = {df['evalue'].mean():.2e}")

# === 4. Overlap des paires alignées ===

def compute_overlap(df1, df2):
    set1, set2 = set(df1["pair"]), set(df2["pair"])
    inter = set1 & set2
    jaccard = len(inter) / len(set1 | set2)
    return len(inter), jaccard

print("\n--- Overlap entre outils (alignements partagés) ---")
names = list(alignments.keys())
for i, j in combinations(names, 2):
    n_common, jaccard = compute_overlap(alignments[i], alignments[j])
    print(f"{i} vs {j}: {n_common} alignments communs, Jaccard = {jaccard:.3f}")

# === 5. Venn diagram (2 ou 3 outils max) ===

plt.figure(figsize=(6,6))
sets = [set(df["pair"]) for df in alignments.values()]
labels = list(alignments.keys())

if len(sets) == 2:
    venn2(sets, set_labels=labels)
elif len(sets) == 3:
    venn3(sets, set_labels=labels)
plt.title("Overlap des alignements")
plt.savefig("venn_alignments.png")
plt.close()

# === 6. Visualisation distributions ===

for metric in ["identity", "evalue", "bitscore"]:
    plt.figure(figsize=(8, 5))
    for name, df in alignments.items():
        sns.kdeplot(df[metric], label=name, fill=True, common_norm=False)
    plt.xlabel(metric)
    plt.ylabel("Densité")
    plt.title(f"Distribution de {metric}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"distribution_{metric}.png")
    plt.close()
