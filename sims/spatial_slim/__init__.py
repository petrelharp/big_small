import os, subprocess
import scipy.stats
import numpy as np
import matplotlib
import tskit

# pop_width = 8.0, 
# numgens = 100, 
# sigma = 1.0, 
# K = 5.0,

def run_slim(script, seed = 23, 
             **kwargs):
    scriptbase = "_".join(script.split(".")[:-1])
    if not os.path.isdir(scriptbase):
        os.mkdir(scriptbase)
    kwstrings = [str(u) + "_" + str(kwargs[u]) for u in kwargs]
    base = os.path.join(scriptbase, 
            "_".join(["run"] + kwstrings + ["seed_" + str(seed)]))
    treefile = base + ".trees"
    if os.path.isfile(treefile):
        print(treefile, "already exists.")
    else:
        logfile = base + ".log"
        slim_command = ["slim", "-s {}".format(seed)]
        slim_command += ["-d {}={}".format(k, v) for k, v in kwargs.items()]
        slim_command += ["-d \"OUTPATH='{}'\"".format(treefile), script]
        print(" ".join(slim_command))
        with open(logfile, "w") as log:
            subprocess.call(" ".join(slim_command), shell=True, stdout=log)
        if not os.path.isfile(treefile):
            raise ValueError("SLiM failed to produce output file {}".format(treefile))
    return(treefile)

# four_colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"]
four_colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]  # Dark2
four_markers = ["v", "<", "^", ">"]


def plot_density(ts, time, ax, scatter=True, alpha=0.8):
    """
    Plot a 2D kernel density estimate of the population density.
    """

    locs = ts.individual_locations[:,:2]
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    alive = ts.individuals_alive_at(time)
    tlocs = locs.T
    kde = scipy.stats.gaussian_kde(tlocs)
    X, Y = np.meshgrid(
            np.linspace(0.0, xmax, 51),
            np.linspace(0.0, ymax, 51))
    Z = kde([X.flatten(), Y.flatten()])
    Z.shape = X.shape
    if scatter:
        ax.scatter(locs[alive, 0], locs[alive, 1],
                   s=10,
                   alpha=0.5,
                   c='black',
                   marker="o",
                   edgecolors='none')
    ax.contour(X, Y, Z,
               colors='c',
               alpha=alpha,
               zorder=-1)

def get_lineages(ts, children, positions, max_time):
    """
    A plot of the lineages ancestral to the given children
    at the given positions.
    """
    locs = ts.individual_locations
    # will record here tuples of the form (time, x)
    nodes = np.concatenate([ts.individual(i).nodes for i in children])
    node_times = ts.tables.nodes.time
    # careful: some are tskit.NULL
    node_indivs = ts.tables.nodes.individual
    paths = {}
    for p in positions:
        tree = ts.at(p)
        for u in nodes:
            out = [np.array([locs[node_indivs[u], 0], node_times[u]])]
            u = tree.parent(u)
            while u is not tskit.NULL:
                uind = node_indivs[u]
                if uind is tskit.NULL or ts.node(u).time > max_time:
                    break
                out.append(np.array([locs[uind, 0], node_times[u]]))
                u = tree.parent(u)
            paths[(u, p)] = np.row_stack(out)
    return paths
