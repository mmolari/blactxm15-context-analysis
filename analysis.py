# %%

import pypangraph as pp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.cluster.hierarchy import linkage, dendrogram
from collections import defaultdict

# %%
# load pangraph
pan = pp.Pangraph.load_json("results/pangraph.json")
bdf = pan.to_blockstats_df()
is_core = bdf["core"].to_dict()
block_length = bdf["len"].to_dict()
strains = pan.strains()

# save block info
print("block info")
print(bdf.to_markdown())
bdf.to_csv("results/block_info.csv")

# %%

# load mash distance matrix
D = pd.read_csv("results/mash.csv")
D = D.pivot(index="si", columns="sj", values="mash_dist")

# perform hierarchical clustering

Z = linkage(D, method="ward")

# plot dendrogram
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
dendrogram(Z, labels=D.index, ax=ax)
str_order = [s.get_text() for s in ax.get_xticklabels()]
str_order_inv = str_order[::-1]
plt.show()

# plot mash distance matrix
fig, ax = plt.subplots(1, 1, figsize=(10, 9))
mapp = ax.matshow(D.loc[str_order_inv, str_order_inv], cmap="coolwarm")
plt.colorbar(mapp, ax=ax, label="mash distance", shrink=0.8)
plt.xticks(np.arange(len(str_order_inv)), str_order_inv, rotation=90)
plt.yticks(np.arange(len(str_order_inv)), str_order_inv)
plt.tight_layout()
plt.savefig("figs/mash_dist_mat.png", facecolor="white")
plt.show()

# %%


# draw pangraph
cmap = iter(mpl.cm.get_cmap("Pastel1")(np.arange(len(bdf))))
block_color = defaultdict(lambda: next(cmap))
singleton_blocks = bdf.index[bdf["count"] == 1].to_list()
for block in singleton_blocks:
    block_color[block] = "lightgray"


def plot_path(ax, y, Bs, Ss, B_len):
    x = 0
    for bid, s in zip(Bs, Ss):
        ax.barh(
            y,
            B_len[bid],
            left=x,
            color=block_color[bid],
            edgecolor="black",
        )
        color = "black" if s else "red"
        ax.text(
            x + B_len[bid] / 2,
            y,
            bid[:4],
            ha="center",
            va="center",
            fontsize=8,
            color=color,
        )
        x += B_len[bid]


fig, ax = plt.subplots(1, 1, figsize=(10, 10))

for n, iso in enumerate(str_order):
    blocks = pan.paths[iso].block_ids
    strands = pan.paths[iso].block_strands
    plot_path(ax, n, blocks, strands, block_length)

ax.set_yticks(np.arange(len(str_order)))
ax.set_yticklabels(str_order)
ax.set_xlabel("position (bp)")
plt.tight_layout()
plt.savefig("figs/pangraph.png", facecolor="white")
plt.show()

# %%
