#!/usr/bin/env python3

import sys
import pyslim, tskit
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections as cs
from matplotlib.patches import Polygon

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

np.random.seed(kwargs['seed'])

def plot_lineage(ts, children, positions):
    """
    A plot of the lineages ancestral to the given children
    at the given positions.
    """
    locs = ts.individual_locations
    inds = ts.individuals_alive_at(0)
    # will record here tuples of the form (time, x)
    nodes = np.concatenate([ts.individual(i).nodes for i in children])
    node_times = ts.tables.nodes.time
    node_indivs = ts.tables.nodes.individual
    paths = []
    for p in positions:
        tree = ts.at(p)
        for u in nodes:
            out = [np.array([locs[node_indivs[u], 0], node_times[u]])]
            u = tree.parent(u)
            while u is not tskit.NULL:
                uind = node_indivs[u]
                if uind is tskit.NULL:
                    break
                out.append(np.array([locs[uind, 0], node_times[u]]))
                u = tree.parent(u)
            paths.append(np.row_stack(out))

    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(111)
    ax.set_xlabel("position")
    ax.set_ylabel("time ago")
    xmax = np.ceil(max(locs[:,0]))
    ymax = np.ceil(max(node_times))
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    colormap = lambda x: plt.get_cmap("cool")(x/max(ts.individual_ages))
    treecolors = [plt.get_cmap("viridis")(x) for x in np.linspace(0, 1, len(positions))]
    pathcolors = []
    for c in treecolors:
        pathcolors.extend([c] * num_positions)

    lc = cs.LineCollection(paths, linewidths=0.5, colors=pathcolors)
    ax.add_collection(lc)
    return fig

def patch_polygons(times, vals):
    active = []
    finished = []
    for t, x in zip(times, vals):
        lefts = np.where(np.diff(np.concatenate([[0], x])) > 0)[0]
        rights = np.where(np.diff(np.concatenate([x, [0]])) < 0)[0]
        new = np.array([True for _ in lefts])
        for k, u in enumerate(active):
            try:
                l = u[-1][1]
                j = np.where(np.abs(lefts - l) <= 1)[0][0]
                new[j] = False
                u.append([t, lefts[j], rights[j]])
            except IndexError:
                # becomes inactive!
                finished.append(u)
                del active[k]
        for j in np.where(new)[0]:
            active.append([[t, lefts[j], rights[j]]])
    finished.extend(active)
    patches = []
    for u in finished:
        v = np.array(u)
        v[:, 0] = ts.slim_generation - v[:, 0]
        patches.append(Polygon(np.row_stack([v[:, [1,0]], v[::-1, [2,0]]])))
    return cs.PatchCollection(patches, alpha=0.4)


treefile = sps.run_slim(script = script, **kwargs)
patchfile = treefile + ".landscape"
outbase = ".".join(treefile.split(".")[:-1])

ts = pyslim.load(treefile)
patchdata = np.loadtxt(patchfile)
patchtimes = patchdata[:, 0]
patchvals = patchdata[:, 1:]
patches = patch_polygons(patchtimes, patchvals)

today = np.where(ts.individual_times == 0)[0]
fig = plot_lineage(ts, 
                   np.random.choice(today, num_indivs), 
                   np.random.randint(0, ts.sequence_length - 1, num_positions))
ax = fig.axes[0]
ax.add_collection(patches)
fig.savefig(outbase + ".lineages.png")
plt.close(fig)
