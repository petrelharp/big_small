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

for treefile in tsfiles:
    print("Doing {}.\n".format(treefile))
    ts = pyslim.load(treefile)
    outbase = ".".join(treefile.split(".")[:-1])

    today = ts.individuals_alive_at(0)
    has_parents = ts.has_individual_parents()

    R = sps.individual_relatedness_matrix(ts)
    R /= 2 * ts.sequence_length
    locs = ts.individual_locations
    times = ts.individual_times
    dt = (times @ R - times)
    dx = (locs[:,0] @ R - locs[:,0])
    dy = (locs[:,1] @ R - locs[:,1])
    var = (dx **2 + dy ** 2) / dt
    assert(np.min(dt[has_parents]) > 0)

    if len(today) < num_indivs:
        raise ValueError(f"Not enough individuals: only {len(today)} alive today!")

    path_dict = sps.get_lineages(ts, 
                                 np.random.choice(today, num_indivs, replace=False),
                                 np.linspace(0, ts.sequence_length - 1, num_positions),
                                 max_time)
    path_keys = list(path_dict.keys())

    ntimes = 40
    times = np.linspace(0, max_time, ntimes)
    x = np.zeros((ntimes, len(path_dict)))
    for i, t in enumerate(times):
        for j, k in enumerate(path_keys):
            u = path_dict[k]
            x[i, j] = u[np.max(np.where(u[:, 1] <= t)[0]), 0]

    dx = x[1:, :] - np.row_stack([x[0, :]] * (ntimes - 1))
    variances = np.concatenate([[0], np.var(dx, axis=1)])

    np.savetxt(outbase + ".variances.txt", np.column_stack([times, variances]), header="time\tvariance")

print("Done.\n")
