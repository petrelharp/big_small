import numpy as np
import matplotlib.collections as cs
import matplotlib.patches

def patch_polygons(times, vals, max_time=None, tolerance=10):
    if max_time is None:
        max_time = np.max(times)
    active = []
    finished = []
    for t, x in zip(times, vals):
        lefts = np.where(np.diff(np.concatenate([[0], x])) > 0)[0]
        rights = np.where(np.diff(np.concatenate([x, [0]])) < 0)[0]
        new = np.array([True for _ in lefts])
        for k, u in enumerate(active):
            try:
                l = u[-1][1] # last left coord
                r = u[-1][2] # last left coord
                j = np.where(
                        np.logical_and(
                            np.abs(lefts - l) <= tolerance,
                            np.abs(rights - r) <= tolerance)
                        )[0][0]
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
    print(f"Found {len(finished)} patches.")
    for u in finished:
        v = np.array(u)
        v[:, 0] = max_time - v[:, 0]
        v[:, 1] -= 0.5
        v[:, 2] += 0.5
        patches.append(matplotlib.patches.Polygon(np.row_stack([v[:, (1,0)], v[::-1, (2,0)]])))
    return cs.PatchCollection(patches, alpha=0.4)
