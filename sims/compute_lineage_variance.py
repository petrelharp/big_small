#!/usr/bin/env python3

import sys
import warnings
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
num_targets = 50
max_n = 1000
# num_targets = 5
# max_n = 10
# How often to record output
init_steps = 40
num_steps = 25
# whether to do lengthy error checking
debug = False

def dW(a, b, W):
    # shorter displacement in the W-circle
    # ..(b-W)..|.a......b..|.(a+W)...
    # ..(a-W)..|.b......a..|.(b+W)...
    return (a - b + W / 2) % W - W / 2

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
    dx = dW(ind_locs[:,0] @ inR, node_locs[:,0], W)
    assert(np.min(dt[has_parents_nodes]) > 0)
    return dx[has_parents_nodes], dt[has_parents_nodes]


def global_var(ts, W, num_targets, record_steps):
    """
    Picks num_targets random individuals alive today
    and returns (len(record_steps), num_targets) arrays for each of
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
    v *= dW(x[0,i], x[0,j], W)**2 / (t[0,i] - t[0,j])
    pcR = np.sum(scipy.sparse.csr_matrix((v, (i, j)), shape=R.shape), axis=0)
    # set up output
    nsteps = len(record_steps)
    dx = np.repeat(np.nan, nsteps * num_targets).reshape((nsteps, num_targets))
    dt = np.repeat(np.nan, nsteps * num_targets).reshape((nsteps, num_targets))
    var = np.repeat(np.nan, nsteps * num_targets).reshape((nsteps, num_targets))
    pc_var = np.repeat(np.nan, nsteps * num_targets).reshape((nsteps, num_targets))
    for j in range(max(record_steps) + 1):
        print(j)
        if j in record_steps:
            row = list(record_steps).index(j)
            pc_var[row, :] = pcR @ Ru
        Ru = R.dot(Ru)
        if j in record_steps:
            row = list(record_steps).index(j)
            totals = np.array(np.sum(Ru[has_parents, :], axis=0)).reshape((num_targets,))
            badones = np.array(~np.isclose(totals, 1)).reshape((num_targets,))
            dt[row, :] = t @ Ru - t0
            pc_var[row, badones] = np.nan
            dt[row, badones] = np.nan
            for k in range(num_targets):
                if np.isclose(totals[k], 1):
                    xdiff = dW(x, x0[k], W)
                    with np.errstate(invalid='ignore', divide='ignore'):
                        z = xdiff**2 / (t - t0[k])
                    a = Ru[:, k]
                    dx[row, k] = xdiff @ a
                    znan = ~np.isfinite(z)
                    assert(np.sum(znan @ a) == 0)
                    z[znan] = 0
                    var[row, k] = z @ a
            if np.all(badones):
                print("All remaining probabilities less than 1:", totals)
                print(f" stopping at generation {j}.")
                break
    if debug:
        test_dx, test_dt, test_var, test_pc_var = naive_global_var(ts, W, targets, record_steps)
        if not (np.allclose(dx, test_dx, equal_nan=True)
                and np.allclose(dt, test_dt, equal_nan=True)
                and np.allclose(var, test_var, equal_nan=True)
                and np.allclose(pc_var, test_pc_var, equal_nan=True)):
            print("Something wrong happened.")
    return dx, dt, var, pc_var


def naive_global_var(ts, W, targets, record_steps):
    """
    Naive implementation of global_var( ) that should always agree (but be slower);
    used for testing.
    """
    num_targets = len(targets)
    has_parents = ts.has_individual_parents()
    indiv_nodes = np.isin(ts.tables.nodes.individual, np.where(has_parents)[0])
    nsteps = len(record_steps)
    dx = np.repeat(np.nan, nsteps * num_targets).reshape((nsteps, num_targets))
    dt = np.repeat(np.nan, nsteps * num_targets).reshape((nsteps, num_targets))
    var = np.repeat(np.nan, nsteps * num_targets).reshape((nsteps, num_targets))
    pc_var = np.repeat(np.nan, nsteps * num_targets).reshape((nsteps, num_targets))
    edges = ts.tables.edges
    for j, target in enumerate([ts.individual(t) for t in targets]):
        inds = [target]
        done = False
        for n in range(max(record_steps) + 1):
            this_dx = []
            this_dt = []
            this_var = []
            this_pc_var = []
            next_inds = []
            for ind in inds:
                ind_edges = np.logical_and(
                                np.isin(edges.child, ind.nodes),
                                indiv_nodes[edges.parent])
                total = np.sum(edges.right[ind_edges] - edges.left[ind_edges])
                if not np.isclose(total, 2 * ts.sequence_length):
                    print('done:', j, n, total/(2*ts.sequence_length))
                    done = True
                    break
                parent_ids = list(set([ts.node(u).individual for u in edges.parent[ind_edges]]))
                if len(parent_ids) == 1:
                    parent_ids = [parent_ids[0], parent_ids[0]]
                parents = [ts.individual(i) for i in parent_ids]
                pc_dx = np.array([dW(p.location[0], ind.location[0], W) for p in parents])
                pc_dt = np.array([p.time - ind.time for p in parents])
                adx = np.array([dW(p.location[0], target.location[0], W) for p in parents])
                adt = np.array([p.time - target.time for p in parents])
                this_dx.append(np.mean(adx))
                this_dt.append(np.mean(adt))
                this_var.append(np.mean(adx**2 / adt))
                this_pc_var.append(np.mean(pc_dx**2 / pc_dt))
                next_inds += parents
            if done:
                break
            inds = next_inds
            if n in record_steps:
                row = list(record_steps).index(n)
                dx[row, j] = np.mean(this_dx)
                dt[row, j] = np.mean(this_dt)
                var[row, j] = np.mean(this_var)
                pc_var[row, j] = np.mean(this_pc_var)
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

    # which time steps to record at
    record_steps = [n for n in range(max_n) if n < init_steps or n % num_steps == 0]
    
    # global statistics
    dx, dt, var, pc_var = global_var(ts, W, num_targets=num_targets, record_steps=record_steps)
    with open(outfile, 'w') as f:
        print("\t".join(["n", "mean_dt", "sd_dt", "mean_var", "sd_var", "2.5%", "25%", "50%", "75%", "97.5%", "mean_pc_var", "sd_pc_var"]), file=f)
        for j, n in enumerate(record_steps):
            with warnings.catch_warnings():
                # np.nanmean throws warnings when whole rows are nan
                warnings.simplefilter("ignore", category=RuntimeWarning)
                print("\t".join(map(str,
                    [n, np.nanmean(dt[j,:]), np.nanstd(dt[j,:]), np.nanmean(var[j,:]), np.nanstd(var[j,:])]
                    + list(np.nanquantile(var[j,:], [.025, .25, .5, .75, .975]))
                    + [np.nanmean(pc_var[j,:]), np.nanstd(pc_var[j,:])])), file=f)

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
