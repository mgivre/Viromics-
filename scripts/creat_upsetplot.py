import pandas as pd
from matplotlib import pyplot as plt
from upsetplot import UpSet, from_memberships

# Définir les données
data = from_memberships(
    [
        ('mmseqs2',),
        ('fasta36',),
        ('mmseqs2', 'fasta36'),
    ],
    data=[
        87717687,
        133842537,
        566041221,
    ]
)

# Tracer le UpSet plot
upset = UpSet(data, show_counts=True, show_percentages=True)
upset.plot()
plt.suptitle("UpSet plot des alignements  mmseqs2(7.5) vs fasta36")
plt.savefig("upsetplot_mmseqs2_vs_fasta36.png", dpi=300, bbox_inches='tight')
plt.show()
