# SPDX-License-Identifier: GPL-3.0-or-later
"""Finer mesh from the stored one: each triangle split in four, new points bent onto the surface with the stored normals
(cubic Hermite curve on each edge). Along an edge between two CAD faces with different normals, the curve follows the
intersection of both (tangent = n1 × n2), computed once for both faces: no crack, sharp edges stay sharp, flat faces stay flat."""
import math

Q = 1e6   # clé de position : au micromètre


def _key(p):
    return (round(p[0] * Q), round(p[1] * Q), round(p[2] * Q))


def _sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def _add(a, b): return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
def _mul(a, s): return (a[0] * s, a[1] * s, a[2] * s)
def _dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def _cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _unit(a):
    n = math.sqrt(_dot(a, a))
    return _mul(a, 1 / n) if n > 1e-12 else None


def _tangent(e, normals):
    """Unit tangent at an end of edge e: along the crease of two different normals, else e flattened onto the surface."""
    if len(normals) >= 2:
        t = _unit(_cross(normals[0], normals[1]))
        if t is not None and abs(_dot(normals[0], normals[1])) < 0.9999:
            return t if _dot(t, e) >= 0 else _mul(t, -1)
    n = normals[0]
    return _unit(_sub(e, _mul(n, _dot(e, n))))


def _midpoint(p0, p1, n0s, n1s):
    e = _sub(p1, p0)
    c = math.sqrt(_dot(e, e))
    mid = _mul(_add(p0, p1), 0.5)
    if c == 0:
        return mid
    t0, t1 = _tangent(e, n0s), _tangent(e, n1s)
    if t0 is None or t1 is None:
        return mid
    return _add(mid, _mul(_sub(t0, t1), c / 8))


def _level(F):
    edge_faces = {}   # arête (clés de position) -> [(face, normale au 1er bout, normale au 2e bout)]
    for fi, (vs, ns, ts) in enumerate(F):
        if ns is None:
            continue
        for tri in ts:
            for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                ka, kb = _key(vs[a]), _key(vs[b])
                if ka > kb:
                    ka, kb, a, b = kb, ka, b, a
                lst = edge_faces.setdefault((ka, kb), [])
                if not any(x[0] == fi for x in lst):
                    lst.append((fi, ns[a], ns[b], vs[a], vs[b]))
    mids = {}
    for k, lst in edge_faces.items():
        _fi, n0, n1, p0, p1 = lst[0]
        n0s = [n0] + [x[1] for x in lst[1:] if _dot(x[1], n0) < 0.9999][:1]
        n1s = [n1] + [x[2] for x in lst[1:] if _dot(x[2], n1) < 0.9999][:1]
        mids[k] = _midpoint(p0, p1, n0s, n1s)
    out = []
    for vs, ns, ts in F:
        vs, ns = list(vs), list(ns) if ns else None
        local = {}
        def mid(a, b):
            k = (a, b) if a < b else (b, a)
            if k not in local:
                ka, kb = _key(vs[a]), _key(vs[b])
                p = mids.get((ka, kb) if ka <= kb else (kb, ka)) if ns else None
                vs.append(p if p is not None else _mul(_add(vs[a], vs[b]), 0.5))
                if ns:
                    ns.append(_unit(_add(ns[a], ns[b])) or ns[a])
                local[k] = len(vs) - 1
            return local[k]
        nt = []
        for a, b, c in ts:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            nt += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        out.append((vs, ns, nt))
    return out


def refine(F, levels):
    for _ in range(levels):
        F = _level(F)
    return F
