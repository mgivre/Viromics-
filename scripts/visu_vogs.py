import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Charger le fichier
df = pd.read_csv("vOTUs.VOGs.tsv", sep="\t", names=["protein", "cluster"])

# Extraire le contig : tout avant le dernier "_"
df['contig'] = df['protein'].apply(lambda x: '_'.join(x.split('_')[:-1]))

# === 3. Sélectionner les top 50 contigs les plus riches ===
top_contigs = df['contig'].value_counts().head(150).index

# === 4. Sélectionner les top 50 clusters les plus fréquents ===
top_clusters = df['cluster'].value_counts().head(150).index

# === 5. Filtrer le DataFrame ===
filtered_df = df[df['contig'].isin(top_contigs) & df['cluster'].isin(top_clusters)]

# === 6. Grouper et créer la table contig × cluster ===
heatmap_data = filtered_df.groupby(['contig', 'cluster']).size().unstack(fill_value=0)

# === 7. Afficher la heatmap ===
plt.figure(figsize=(14, 10))
#sns.heatmap(heatmap_data, cmap='viridis', linewidths=0.5, linecolor='gray')
sns.heatmap(heatmap_data, cmap='viridis', cbar=True,
            xticklabels=False, yticklabels=False,
            linecolor=None, linewidths=0)
plt.title("Heatmap : Top 150 vOTUs × Top 150 VOGs")
plt.xlabel("VOGs")
plt.ylabel("vOTUs")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("heatmap_top150_vOTUs_VOGs.png", dpi=300)
plt.show()