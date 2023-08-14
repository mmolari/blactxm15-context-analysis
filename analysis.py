# %%

import pypangraph as pp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
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
D_sq = squareform(D)
Z = linkage(D_sq, method="complete")

# plot dendrogram
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
dendrogram(Z, labels=D.index, ax=ax)
str_order = [s.get_text() for s in ax.get_xticklabels()]
str_order_inv = str_order[::-1]
# plt.show()
plt.close(fig)

# plot mash distance matrix
fig, ax = plt.subplots(1, 1, figsize=(10, 9))
mapp = ax.matshow(D.loc[str_order_inv, str_order_inv], cmap="coolwarm")
plt.colorbar(mapp, ax=ax, label="mash distance", shrink=0.8)
plt.xticks(np.arange(len(str_order_inv)), str_order_inv, rotation=90)
plt.yticks(np.arange(len(str_order_inv)), str_order_inv)
plt.tight_layout()
plt.savefig("figs/mash_dist_mat.png", facecolor="white")
# plt.show()
plt.close(fig)

# %%


# pick colormap for blocks
N_nonsingleton_blocks = (bdf["count"] > 1).sum()
if N_nonsingleton_blocks <= 8:
    cmap = iter(plt.get_cmap("Pastel1")(np.arange(N_nonsingleton_blocks)))
elif N_nonsingleton_blocks <= 20:
    cmap = iter(plt.get_cmap("tab20")(np.arange(N_nonsingleton_blocks)))
else:
    cmap = plt.get_cmap("rainbow")(np.linspace(0, 1, N_nonsingleton_blocks))
    np.random.shuffle(cmap)
    cmap = iter(cmap)
block_color = defaultdict(lambda: next(cmap))

singleton_blocks = bdf.index[bdf["count"] == 1].to_list()
for block in singleton_blocks:
    block_color[block] = "white"

# align on the right -> evaluate maximum path length
L_max = 0
for path in pan.paths:
    Bs = path.block_ids
    Ls = [block_length[b] for b in Bs]
    L_max = max(L_max, sum(Ls))


def plot_path(ax, y, Bs, Ss, B_len):
    x = L_max
    for bid, s in zip(Bs[::-1], Ss[::-1]):
        l = B_len[bid]
        x -= l
        ax.barh(
            y,
            l,
            left=x,
            color=block_color[bid],
            edgecolor="black",
        )
        color = "black" if s else "red"
        ax.text(
            x + l / 2,
            y,
            bid[:4],
            ha="center",
            va="center",
            fontsize=8,
            color=color,
        )


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
# plt.show()
plt.close(fig)

# %%
