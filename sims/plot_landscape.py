#!/usr/bin/env python3

import sys
import numpy as np
import spatial_slim as sps

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.collections as cs

usage = """
Makes a plot of a 1D "landscape" changing in time from the provided files,
which should have no header and be whitespace-separated binary matrices, with
the first column giving the time step.
Usage:
    {} [file.landscape [file.landscape]]

Also plots statistics of 
""".format(sys.argv[0])

if len(sys.argv) < 2:
    raise ValueError(usage)

for patchfile in sys.argv[1:]:
    outfile = patchfile + ".png"
    statsfile = patchfile + ".lengths.png"
    patchdata = np.loadtxt(patchfile)
    times = patchdata[:, 0]
    vals = patchdata[:, 1:]

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlabel("geographic position")
    ax.set_ylabel("time ago")
    ax.set_title(patchfile)
    ax.set_xlim(-1, vals.shape[1] + 1)
    ax.set_ylim(0, np.max(times) - np.min(times))
    patches = sps.patch_polygons(times, vals)
    ax.add_collection(patches)
    fig.savefig(outfile)
    plt.close(fig)

    lengths = sps.patch_lengths(times, vals, breaks=np.linspace(min(times), max(times)+1, 11))
    has_lengths = np.sum(lengths, axis=0)
    rates = lengths[:, 1:max(np.where(has_lengths > 0)[0])]
    for k in range(rates.shape[1]):
        with np.errstate(invalid='ignore'):
            rates[:, k] /= np.sum(rates[:, k:], axis=1)
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlabel("patch length")
    ax.set_ylabel("hazard rate")
    ax.set_title(patchfile)
    ax.plot(rates.T)
    fig.savefig(statsfile)
    plt.close(fig)


