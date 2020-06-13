#!/usr/bin/env python3

import sys
import pyslim, tskit
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections as cs

usage = """
Makes a plot of *lineages* moving back for the given number of
time steps at a randomly chosen set of positions on the genome
from the chromosomes of the given number of diploid individuals.
Usage:
    {} (num indivs) (num positions) (script name) [KEY=VALUE [KEY=VALUE]]
where the KEY=VALUE pairs get passed to SLiM.
""".format(sys.argv[0])

if len(sys.argv) < 4:
    raise ValueError(usage)

try:
    num_indivs = int(sys.argv[1])
except ValueError:
    raise ValueError("First argument should be numeric (a number of individuals' lineages to plot).")
try:
    num_positions = int(sys.argv[2])
except ValueError:
    raise ValueError("Second argument should be numeric (a number of positions to plot lineages at).")
script = sys.argv[3]
kwargs = {}
for kv in sys.argv[4:]:
    k, v = kv.split("=")
    kwargs[k] = v

if 'seed' not in kwargs:
    kwargs['seed'] = np.random.randint(10000)

kwargs['seed'] = int(kwargs['seed'])
np.random.seed(kwargs['seed'])

def plot_lineages(ts, ax, children, positions, max_time):
    """
    A plot of the lineages ancestral to the given children
    at the given positions.
    """
    path_dict = sps.get_lineages(ts, children, positions, max_time)
    locs = ts.individual_locations
    xmax = np.ceil(max(locs[:,0]))
    ymax = np.ceil(max([np.max(u[:, 1]) for u in path_dict.values()]))
    ax.set_xlabel("geographic position")
    ax.set_ylabel("time ago")
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    colormap = lambda x: plt.get_cmap("cool")(x/max(ts.individual_ages))
    treecolors = {p : plt.get_cmap("viridis")(x)
                  for p, x in zip(positions, np.linspace(0, 1, len(positions)))}
    paths = []
    pathcolors = []
    for u, p in path_dict:
        paths.append(path_dict[(u, p)])
        pathcolors.append(treecolors[p])
    lc = cs.LineCollection(paths, linewidths=0.5, colors=pathcolors)
    ax.add_collection(lc)


treefile = sps.run_slim(script = script, **kwargs)
patchfile = treefile + ".landscape"
outbase = ".".join(treefile.split(".")[:-1])

ts = pyslim.load(treefile)
patchdata = np.loadtxt(patchfile)
patchtimes = patchdata[:, 0]
patchvals = patchdata[:, 1:]
patches = sps.patch_polygons(patchtimes, patchvals, ts.slim_generation)

today = ts.individuals_alive_at(0)
has_parents = ts.has_individual_parents()
max_time = np.max(ts.individual_times[has_parents])
if len(today) < num_indivs:
    raise ValueError(f"Not enough individuals: only {len(today)} alive today!")

fig, ax = plt.subplots(figsize=(9, 9))
plot_lineages(ts, ax,
              children=np.random.choice(today, num_indivs, replace=False),
              positions=np.random.randint(0, ts.sequence_length - 1, num_positions),
              max_time=max_time)
ax.add_collection(patches)
fig.savefig(outbase + ".lineages.png")
plt.close(fig)
