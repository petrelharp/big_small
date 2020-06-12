#!/usr/bin/env python3

import sys
import numpy as np
import pyslim, tskit
import spatial_slim as sps

usage = """
Estimates the variance per unit time of the linages in the provided tree sequences.
Usage:
    {} [file [file]]
""".format(sys.argv[0])

if len(sys.argv) < 2:
    raise ValueError(usage)

tsfiles = sys.argv[1:]

def quantile(x, q):
    if hasattr(np, 'quantile'):
        return list(np.quantile(x, q))
    else:
        x.sort()
        return [x[int(len(x) * qq)] for qq in q]

for treefile in tsfiles:
    outbase = ".".join(treefile.split(".")[:-1])
    outfile = outbase + ".variances.txt"
    print(f"Writing variances in {treefile} to {outfile}")
    ts = pyslim.load(treefile)

    today = ts.individuals_alive_at(0)
    has_parents = ts.has_individual_parents()

    R = sps.individual_relatedness_matrix(ts)
    R /= 2 * ts.sequence_length
    locs = ts.individual_locations
    times = ts.individual_times
    dt = (times @ R - times)
    dx = (locs[:,0] @ R - locs[:,0])
    var = dx **2 / dt
    assert(np.min(dt[has_parents]) > 0)

    with open(outfile, 'w') as f:
        print("\t".join(["mean", "stdev", "2.5%", "25%", "50%", "75%", "97.5%"]), file=f)
        print("\t".join(map(str,
            [np.mean(var[has_parents]), np.std(var[has_parents])]
            + quantile(var[has_parents], [.025, .25, .5, .75, .975]))), file=f)

print("Done.\n")
