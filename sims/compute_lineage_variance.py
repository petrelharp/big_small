#!/usr/bin/env python3

import sys
import numpy as np
import pyslim, tskit
import scipy
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

usage = """
Estimates the variance per unit time of the linages in the provided tree sequences.
Usage:
    {} [file [file]]
""".format(sys.argv[0])

if len(sys.argv) < 2:
    raise ValueError(usage)

tsfiles = sys.argv[1:]

for treefile in tsfiles:
    outbase = ".".join(treefile.split(".")[:-1])
    outfile = outbase + ".variances.txt"
    plotfile = outbase + ".variances.png"
    print(f"Writing variances in {treefile} to {outfile} and plotting to {plotfile}")
    ts = pyslim.load(treefile)

    today = ts.individuals_alive_at(0)
    has_parents = ts.has_individual_parents()
    has_parents_nodes = np.append(has_parents, False)[ts.tables.nodes.individual]

    # inR[i,j] will give the proportion of the genome
    # that node j inherits from individual i -- so,
    # all of the nonzero entries of inR[:,has_parents_nodes] should be 1.0,
    # and each column should have just one nonzero
    R = sps.node_relatedness_matrix(ts)
    R /= ts.sequence_length
    N = sps.individual_node_matrix(ts)
    inR = N @ R
    assert(np.allclose(np.sum(inR[:,has_parents_nodes], axis=0), 1.0))
    _, j, v = scipy.sparse.find(inR)
    assert(np.allclose(v[has_parents_nodes[j]], 1.0))
    ind_locs = ts.individual_locations
    node_locs = ts.individual_locations[ts.tables.nodes.individual]
    times = ts.individual_times
    # TODO: need to compute var/dt separately for each parent-child link
    dt = (times @ inR - ts.tables.nodes.time)
    dx2 = (ind_locs[:,0] @ inR - node_locs[:,0]) ** 2
    vardt = dx2 / dt
    assert(np.min(dt[has_parents_nodes]) > 0)

    with open(outfile, 'w') as f:
        print("\t".join(["mean", "stdev", "2.5%", "25%", "50%", "75%", "97.5%"]), file=f)
        print("\t".join(map(str,
            [np.mean(vardt[has_parents_nodes]), np.std(vardt[has_parents_nodes])]
            + list(np.quantile(vardt[has_parents_nodes], [.025, .25, .5, .75, .975])))), file=f)


    kde = scipy.stats.gaussian_kde(np.row_stack([dt[has_parents_nodes], np.sqrt(dx2[has_parents_nodes])]))
    X, Y = np.meshgrid(
            np.linspace(0.0, np.max(dt[has_parents_nodes]), 51),
            np.linspace(0.0, np.max(np.sqrt(dx2[has_parents_nodes])), 51))
    Z = kde([X.flatten(), Y.flatten()])
    Z.shape = X.shape
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.scatter(dt[has_parents_nodes], np.sqrt(dx2[has_parents_nodes]), s=5)
    ax.set_xlabel("dt")
    ax.set_ylabel("|dx|")
    ax.contour(X, Y, Z,
               colors='r',
               alpha=0.95,
               zorder=-1)

    fig.savefig(plotfile)
    plt.close(fig)

print("Done.\n")
