import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet
import matplotlib

# Empêche les problèmes d'affichage en environnement sans GUI
matplotlib.use('Agg')


def load_pairs(filepath, name):
    """Charge toutes les paires query||target du premier fichier"""
    if not os.path.exists(filepath):
        print(f"❌ Fichier non trouvé : {filepath}", file=sys.stderr)
        return set()
    
    print(f"📥 Lecture de {name} depuis {filepath}")
    pairs = set()
    chunk_size = 10**6

    for chunk in pd.read_csv(filepath, sep="\t", header=None, usecols=[0, 1], dtype=str, chunksize=chunk_size):
        pairs.update(chunk[0] + "||" + chunk[1])

    print(f"✅ {name} : {len(pairs):,} paires extraites")
    return pairs


def count_overlap(set1, filepath, name2):
    """Compare le fichier donné à set1 et compte sans stocker"""
    common = 0
    only_in_2 = 0
    total = 0

    print(f"📥 Streaming de {name2} et comptage...")

    for chunk in pd.read_csv(filepath, sep="\t", header=None, usecols=[0, 1], dtype=str, chunksize=10**6):
        for pair in chunk[0] + "||" + chunk[1]:
            total += 1
            if pair in set1:
                common += 1
            else:
                only_in_2 += 1

    return total, common, only_in_2


def plot_upset(common, only1, only2):
    """Génère un UpSet plot à partir des comptes"""
    from upsetplot import UpSet
    import pandas as pd
    import matplotlib.pyplot as plt

    print("📊 Génération du UpSet plot...")

    data = pd.Series(
        [common, only1, only2],
        index=pd.MultiIndex.from_tuples(
            [
                (True, True),
                (True, False),
                (False, True)
            ],
            names=["mmseqs2", "fasta36"]
        )
    )

    os.makedirs("plots", exist_ok=True)
    upset = UpSet(data, subset_size='count', show_counts=True)
    upset.plot()
    plt.suptitle("UpSet Plot - Recouvrement des paires alignées")
    plt.savefig("plots/upset_plot.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("✅ UpSet plot enregistré dans plots/upset_plot.png")


def main():
    name1, file1 = "mmseqs2(7.5)", "mmseqs_sensitive_result.m8"
    name2, file2 = "fasta36", "full_results.tab"

    set1 = load_pairs(file1, name1)
    total2, common, only2 = count_overlap(set1, file2, name2)
    only1 = len(set1) - common

    print("\n📊 Résumé :")
    print(f"🔹 {name1} : {len(set1):,} paires")
    print(f"🔹 {name2} : {total2:,} paires (streamé)")
    print(f"✅ En commun : {common:,}")
    print(f"➖ Uniques {name1} : {only1:,}")
    print(f"➖ Uniques {name2} : {only2:,}")

    plot_upset(common, only1, only2)


if __name__ == "__main__":
    main()
