import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

# Dictionnaire : nom de l’outil -> chemin du fichier .csv
tools = {
    "FASTA36": "fasta36_vOTUs.mat",
    "MMseqs2_sens7.5": "mm_sensitive_vOTUs.mat",
    "DIAMOND": "diamond_vOTUs.mat"
}

# Afficher une heatmap pour chaque outil
for name, filepath in tools.items():
    
    # 1. Lire la matrice depuis le fichier CSV
    df = pd.read_csv(filepath, index_col=0)
    
    print(f"{name} shape: {df.shape}")
    print(f"Colonnes: {df.columns[:5]}")
    print(f"Lignes: {df.index[:5]}")


    # 2. Vérification basique
    assert df.shape[0] == df.shape[1], f"Matrice non carrée pour {name}"
    
    # 3. Clustering hiérarchique basé sur les distances
    linkage_matrix = linkage(squareform(df.values), method="average")
    
    # 4. Générer la heatmap avec clustering
    g = sns.clustermap(df,
                       row_linkage=linkage_matrix,
                       col_linkage=linkage_matrix,
                       cmap="mako",  # ou 'viridis', 'coolwarm', etc.
                       figsize=(12, 10),
                       xticklabels=False,
                       yticklabels=False)

    g.fig.suptitle(f"Heatmap – {name}", y=1.02)
    plt.show()
