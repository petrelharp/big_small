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

""".format(sys.argv[0])

if len(sys.argv) < 2:
    raise ValueError(usage)

for patchfile in sys.argv[1:]:
    outfile = patchfile + ".png"
    patchdata = np.loadtxt(patchfile)
    patchtimes = patchdata[:, 0]
    patchvals = patchdata[:, 1:]

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlabel("geographic position")
    ax.set_ylabel("time ago")
    ax.set_title(patchfile)
    ax.set_xlim(-1, patchvals.shape[1] + 1)
    ax.set_ylim(0, np.max(patchtimes) - np.min(patchtimes))
    patches = sps.patch_polygons(patchtimes, patchvals)
    ax.add_collection(patches)
    fig.savefig(outfile)
    plt.close(fig)

