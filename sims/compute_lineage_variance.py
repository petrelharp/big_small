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

do_local = False

def local_var(ts, W):
    """
    Returns the vectors of dx and dt across all parent-offspring pairs.
    """
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
    dt = (times @ inR - ts.tables.nodes.time)
    dx = (ind_locs[:,0] @ inR - node_locs[:,0])
    # for periodic boundaries
    dx[dx > W / 2] -= W
    dx[dx < - W / 2] += W
    assert(np.min(dt[has_parents_nodes]) > 0)
    return dx[has_parents_nodes], dt[has_parents_nodes]


def global_var(ts, W, num_targets, max_n):
    """
    Picks num_targets random individuals alive today
    and returns (max_n, num_targets) arrays for each of
    dx, dt, and var, giving the mean displacement, mean dt,
    and mean (dx**2/dt) across level-n ancestors.
    Also computes the mean parent-child (dx**2/dt) long those lineages.
    """
    today = ts.individuals_alive_at(0)
    has_parents = ts.has_individual_parents()
    num_targets = min(num_targets, ts.num_individuals)
    targets = np.random.choice(np.where(has_parents)[0], num_targets, replace=False)
    print(f" choosing {num_targets} individuals out of {len(today)} alive (and {ts.num_individuals} total).")

    # R[i,j] will give the proportion of the genome
    # that individual j inherits from individual i -- so,
    # most of the nonzero entries of inR[:,has_parents_nodes] should be 0.5,
    # with occasional 1.0s, and columns should sum to 1.0
    R = sps.individual_relatedness_matrix(ts)
    # R.tocsr()
    R /= 2 * ts.sequence_length
    assert(np.allclose(np.sum(R[:,has_parents], axis=0), 1.0))
    Ru = scipy.sparse.csc_matrix(
            (np.repeat(1.0, num_targets),
             (targets,
              np.arange(num_targets))), shape=(ts.num_individuals, num_targets))
    x = ts.individual_locations[:, 0].reshape((1, ts.num_individuals))
    x0 = x[0, targets]
    t = ts.individual_times.reshape((1,ts.num_individuals))
    t0 = t[0, targets]
    # to get the parent-child mean dx**2/dt, we need to compute for each k
    #  sum_ij R[i,j] (xi - xj)^2 / (ti - tj) Ru[j,k]
    # so we precompute pcR[j] = sum_i R[ij] (xi - xj)**2 / (ti - tj)
    i, j, v = scipy.sparse.find(R)
    v *= (x[0,i] - x[0,j])**2 / (t[0,i] - t[0,j])
    pcR = np.sum(scipy.sparse.csr_matrix((v, (i, j)), shape=R.shape), axis=0)
    # set up output
    dx = np.repeat(np.nan, max_n * num_targets).reshape((max_n, num_targets))
    dt = np.repeat(np.nan, max_n * num_targets).reshape((max_n, num_targets))
    var = np.repeat(np.nan, max_n * num_targets).reshape((max_n, num_targets))
    pc_var = np.repeat(np.nan, max_n * num_targets).reshape((max_n, num_targets))
    for j in range(max_n):
        Ru = R.dot(Ru)
        dx[j, :] = x @ Ru - x0
        dt[j, :] = t @ Ru - t0
        pc_var[j, :] = pcR @ Ru
        totals = np.array(np.sum(Ru[has_parents, :], axis=0)).reshape((num_targets,))
        badones = np.array(~np.isclose(totals, 1)).reshape((num_targets,))
        dx[j, badones] = np.nan
        dt[j, badones] = np.nan
        pc_var[j, badones] = np.nan
        for k in range(num_targets):
            if np.isclose(totals[k], 1):
                with np.errstate(invalid='ignore', divide='ignore'):
                    z = (x - x0[k])** 2 / (t - t0[k])
                a = Ru[:, k]
                znan = ~np.isfinite(z)
                assert(np.sum(znan @ a) == 0)
                z[znan] = 0
                var[j, k] = z @ a
        if np.all(badones):
            print("All remaining probabilities less than 1:", totals)
            print(f" stopping at generation {j}.")
            break
    return dx, dt, var, pc_var

for treefile in tsfiles:
    outbase = ".".join(treefile.split(".")[:-1])
    outfile = outbase + ".variances.txt"
    local_outfile = outbase + ".local_variances.txt"
    local_plotfile = outbase + ".local_variances.png"
    print(f"Writing variances in {treefile} to {outfile} and {local_outfile} and plotting to {local_plotfile}")
    ts = pyslim.load(treefile)
    # get the width of the region
    patchfile = treefile + ".landscape"
    patchrows = np.loadtxt(patchfile, max_rows=2)
    W = patchrows.shape[1]

    # global statistics
    dx, dt, var, pc_var = global_var(ts, W, num_targets=50, max_n=1000)
    with open(outfile, 'w') as f:
        print("\t".join(["n", "mean_dt", "sd_dt", "mean_var", "sd_var", "2.5%", "25%", "50%", "75%", "97.5%", "mean_pc_var", "sd_pc_var"]), file=f)
        for n in range(dt.shape[0]):
            print("\t".join(map(str,
                [n, np.nanmean(dt[n,:]), np.nanstd(dt[n,:]), np.nanmean(var[n,:]), np.nanstd(var[n,:])]
                + list(np.nanquantile(var[n,:], [.025, .25, .5, .75, .975]))
                + [np.nanmean(pc_var[n,:]), np.nanstd(pc_var[n,:])])), file=f)

    if do_local:
        # local statistics
        ldx, ldt = local_var(ts, W)
        vardt = ldx ** 2 / ldt
        with open(local_outfile, 'w') as f:
            print("\t".join(["mean", "stdev", "2.5%", "25%", "50%", "75%", "97.5%"]), file=f)
            print("\t".join(map(str,
                [np.mean(vardt), np.std(vardt)]
                + list(np.quantile(vardt, [.025, .25, .5, .75, .975])))), file=f)

        kde = scipy.stats.gaussian_kde(np.row_stack([ldt, np.abs(ldx)]))
        X, Y = np.meshgrid(
                np.linspace(0.0, np.max(ldt), 51),
                np.linspace(0.0, np.max(np.abs(ldx)), 51))
        Z = kde([X.flatten(), Y.flatten()])
        Z.shape = X.shape
        fig, ax = plt.subplots(figsize=(9, 9))
        ax.scatter(ldt, np.abs(ldx), s=5)
        ax.set_xlabel("dt")
        ax.set_ylabel("|dx|")
        ax.contour(X, Y, Z,
                   colors='r',
                   alpha=0.95)

        fig.savefig(local_plotfile)
        plt.close(fig)
    else:
        print("Skipping local stats.")

print("Done.\n")
