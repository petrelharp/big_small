import numpy as np
import matplotlib.collections as cs
import matplotlib.patches

def get_patches(times, vals, periodic=True):
    """
    Produces an iterator over the slices that returns the patches
    at that time.
    """
    n = vals.shape[1]
    for t, x in zip(times, vals):
        lefts = np.where(np.diff(np.concatenate([[0], x])) > 0)[0]
        rights = np.where(np.diff(np.concatenate([x, [0]])) < 0)[0]
        if periodic and lefts[0] == 0 and rights[-1] == n-1:
            lefts = lefts[1:]
            rights = np.append(rights[1:-1], rights[0])
        assert(len(lefts) == len(rights))
        yield (t, lefts, rights)


def patch_lengths(times, vals, breaks=None):
    """
    Returns the array of patch lengths.
    """
    n = vals.shape[1]
    if breaks is None:
        breaks = [min(times), max(times) + 1]
    counts = np.zeros((len(breaks) - 1, n))
    lengths = []
    k = 0
    for t, lefts, rights in get_patches(times, vals):
        if t >= breaks[k+1]:
            counts[k] += np.bincount(lengths, minlength=n)
            lengths = []
            k += 1
        lengths.extend(1 + (rights - lefts) % n)
    counts[k] += np.bincount(lengths, minlength=n)
    return counts


def match_patches(patches, t, lefts, rights, tolerance):
    """
    Given a list of lists of patches ("patches"), append each of the new patches
    defined by t, lefts, rights to the closest-matching ones, possibly creating
    new lists of patches.
    """
    new_matches = [[] for _ in lefts]
    patch_matches = [[] for _ in patches]
    for j, (l, r) in enumerate(zip(lefts, rights)):
        for k, (a, u) in enumerate(patches):
            if a:
                ll = u[-1][1]
                lr = u[-1][2]
                dx = np.abs(l - ll) + np.abs(r - lr)
                if dx <= tolerance:
                    new_matches[j].append((dx, ll - lr, k))
                    patch_matches[k].append((dx, l - r, j))
        new_matches[j].sort()
    matches = np.repeat(-1, len(new_matches))
    for j, m in enumerate(new_matches):
        if len(m) > 0:
            u = m.pop(0)
            matches[j] = u[-1]
    num_matches = np.bincount(matches + 1, minlength=len(patches)+1)[1:]
    while np.any(num_matches > 1):
        for k in np.where(num_matches > 1)[0]:
            xm = [x for x in patch_matches[k] if matches[x[2]] == k]
            xm.sort()
            for _, _, j in xm[1:]:
                if len(new_matches[j]) > 0:
                    _, _, k = new_matches[j].pop(0)
                    matches[j] = k
                else:
                    matches[j] = -1
        num_matches = np.bincount(matches + 1, minlength=len(patches)+1)[1:]
    for j, k in enumerate(matches):
        if k == -1:
            patches.append([True, [(t, lefts[j], rights[j])]])
        else:
            patches[k][1].append((t, lefts[j], rights[j]))
    for k in np.where(num_matches == 0)[0]:
        patches[k][0] = False


def patch_polygons(times, vals, max_time=None, tolerance=5):
    if max_time is None:
        max_time = np.max(times)
    patches = []
    for t, lefts, rights in get_patches(times, vals, periodic=False):
        match_patches(patches, t, lefts, rights, tolerance)
    polygons = []
    print(f"Found {len(patches)} patches.")
    for _, u in patches:
        v = np.array(u)
        v[:, 0] = max_time - v[:, 0]
        v[:, 1] -= 0.5
        v[:, 2] += 0.5
        if v.shape[0] == 1:
            t, x, y = v[0]
            v = np.row_stack([
                [t-0.5, (x+y)/2, (x+y)/2],
                v,
                [t+0.5, (x+y)/2, (x+y)/2]])
        polygons.append(matplotlib.patches.Polygon(np.row_stack([v[:, (1,0)], v[::-1, (2,0)]])))
    return cs.PatchCollection(polygons, alpha=0.4)
