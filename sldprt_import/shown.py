# SPDX-License-Identifier: GPL-3.0-or-later
"""Which bodies of the exact geometry the original software shows: those found in the display mesh of the file.

The display mesh is stored per body, each face with the Parasolid id of its face (reader.bodies). A body of the exact
geometry is shown when a display body has its face ids; bodies with the same ids (mirror or pattern copies) are told
apart by their position. Nothing is guessed: when the display mesh does not fit the exact geometry, None is returned
and the caller keeps every body visible (or uses its own test)."""


def _gap(a, b):
    """Distance between two boxes (low, high): sum of the distances of their corners."""
    return sum(abs(a[0][c] - b[0][c]) + abs(a[1][c] - b[1][c]) for c in range(3))


def _size(b):
    return sum((b[1][c] - b[0][c]) ** 2 for c in range(3)) ** 0.5


def shown(display, bodies):
    """display: [(face ids, box)] from reader.bodies; bodies: {key: (face ids, box or None)} of the exact geometry.

    Returns the set of keys shown, or None when the display mesh cannot be trusted for this part."""
    if not display or not bodies:
        return None
    every = set().union(*(ids for ids, _b in bodies.values()))
    seen = [i for ids, _b in display for i in ids]
    if not seen or sum(1 for i in seen if i in every) < 0.98 * len(seen):
        return None   # ids that are not those of this geometry (other configuration, older file…)
    pairs = []
    for g, (ids, box) in enumerate(display):
        if not ids:
            continue   # a body of one face: its id is not stored with it
        for key, (eids, ebox) in bodies.items():
            if ids <= eids:
                pairs.append((_gap(box, ebox) if ebox else float('inf'), g, key))
    pairs.sort(key=lambda p: p[0])
    used_g, out = set(), set()
    for d, g, key in pairs:
        if g in used_g or key in out:
            continue
        if d == float('inf'):   # no position to choose between bodies with the same faces: all of them kept
            out.update(k for _d, g2, k in pairs if g2 == g)
        else:
            out.add(key)
        used_g.add(g)
    # display bodies of a single face (its id is not stored): the exact body at the same place, or every body of that
    # size that could be it
    for g, (ids, box) in enumerate(display):
        if ids or g in used_g:
            continue
        near = sorted((_gap(box, eb) if eb else float('inf'), k) for k, (e, eb) in bodies.items() if k not in out and len(e) <= 2)
        if near and near[0][0] < 0.05 * _size(box) + 0.001 and (len(near) == 1 or near[1][0] > 4 * near[0][0] + 0.001):
            out.add(near[0][1])
        else:
            out.update(k for _d, k in near)
    return out
