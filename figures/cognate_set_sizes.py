"""Plot cognate sets per meaning class in data/cpacific.nex (ascertainment column excluded).

Usage: python3 figures/cognate_set_sizes.py   (needs matplotlib)
"""
import collections
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = Path(__file__).parent
nex = (here.parent / "data" / "cpacific.nex").read_text()

# Character labels: each meaning class starts with a "<name>_0ascertainment" column.
labels = re.findall(r"^\s*\d+\s+(\S+?),?\s*$", nex.split("CHARSTATELABELS")[1].split(";")[0], re.M)
sites = collections.Counter()  # meaning class -> number of cognate sets
cls = None
for label in labels:
    if label.endswith("_0ascertainment"):
        cls = label
    else:
        sites[cls] += 1
assert sum(sites.values()) == len(labels) - len(set(l for l in labels if l.endswith('_0ascertainment')))

counts = collections.Counter(sites.values())
bands = [(1, 4, "1-4", "#0072B2"), (5, 9, "5-9", "#E69F00"), (10, 14, "10-14", "#009E73")]  # Okabe-Ito; = BEAUti sizes 1-5, 6-10, 11-15 (+1 ascertainment column)
assert all(any(lo <= k <= hi for lo, hi, *_ in bands) for k in counts), "size outside bands"

fig, ax = plt.subplots(figsize=(5.5, 3.8))
for lo, hi, name, colour in bands:
    ks = [k for k in sorted(counts) if lo <= k <= hi]
    ax.bar(ks, [counts[k] for k in ks], color=colour, width=0.7,
           label=f"{name}: {sum(counts[k] for k in ks)} classes")
ax.set(xlabel="Cognate sets per meaning class", ylabel="Meaning classes",
       title=f"{len(sites)} meaning classes", xticks=range(min(counts), max(counts) + 1))
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
ax.set_axisbelow(True)
ax.yaxis.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(here / "cognate_set_sizes.png", dpi=150)
