# SPDX-License-Identifier: GPL-3.0-or-later
"""Meshing of the B-rep faces (Blender's constrained Delaunay triangulation in each face's (u, v)).

Each edge is sampled once (shared by its two faces: no crack); each face is triangulated in its (u, v) plane with its loops as
constraints and a grid of inner points where the surface bends, then carried back onto the surface."""
import math

from mathutils import Vector
from mathutils.geometry import delaunay_2d_cdt

from . import brep as B


class FaceMesh:
    __slots__ = ("verts", "tris", "normals", "face", "body", "fallback")

    def __init__(self, face, body):
        self.verts, self.tris, self.normals, self.face, self.body, self.fallback = [], [], [], face, body, False


def _sense(c):
    return -1 if c == b'-' else 1


class Mesher:
    def __init__(self, model, tol=0.0002, max_edge=None):
        self.M, self.N, self.tol = model, model.N, tol
        self.max_edge = max_edge
        self.edges = {}
        self.stats = {"faces": 0, "fallback": 0, "errors": []}

    # ------------------------------------------------------------ arêtes

    def _vertex_pt(self, vx):
        n = self.N.get(vx)
        return list(self.N[n.f[5]].f[5]) if n is not None else None

    def _edge_curve(self, edge, fin=None):
        c = self.M.curve(edge.f[6])
        if c is None and fin is not None:
            c = self.M.curve(self.N[fin].f[7])
        return c

    def _curve_sense(self, edge):
        n = self.N.get(edge.f[6])
        return _sense(n.f[6]) if n is not None and n.f and len(n.f) > 6 and isinstance(n.f[6], bytes) else 1

    def edge_points(self, ei):
        """3D points of edge ei from its start vertex to its end vertex (edge direction = its positive fin)."""
        if ei in self.edges:
            return self.edges[ei]
        e = self.N[ei]
        fin = e.f[3]
        fins = [fin, self.N[fin].f[5]] if fin in self.N else []
        pos = next((f for f in fins if f in self.N and self.N[f].f[9] == b'+'), fin)
        # une fin va du sommet de la fin précédente (départ) à son propre sommet (arrivée)
        bwd = self.N[pos].f[3] if pos in self.N else None
        a = self._vertex_pt(self.N[bwd].f[4]) if bwd in self.N else None
        b = self._vertex_pt(self.N[pos].f[4]) if pos in self.N else None
        cands = [self.M.curve(e.f[6])] + [self.M.curve(self.N[f].f[7]) for f in fins if f in self.N]
        pts = None
        for curve in cands:
            if curve is None:
                continue
            if isinstance(curve, B.Trimmed):   # courbe taillée : échantillonnée entre ses paramètres, si elle est cohérente
                if a is None:
                    continue
                if curve.ok():
                    n = self._count(curve)
                    pts = curve.points(a, n)
                    pts[0], pts[-1] = a, b
                    break
                base = curve.basis   # paramètres du fichier inutilisables (tracé d'intersection…) : projection des bouts
                if isinstance(base, (B.Polyline, B.Line, B.Circle, B.NurbsCurve)):
                    pts = self._sample(base, a, b, self._curve_sense(e))
                    if self._consistent(pts, base):
                        break
                    pts = None
                continue
            pts = self._sample(curve, a, b, self._curve_sense(e))
            break
        if pts is None:
            pts = [a, b] if a is not None and b is not None else []
        self.edges[ei] = (pts, pos)
        return self.edges[ei]

    def _sample(self, c, a, b, csense):
        tol = self.tol
        closed = a is None or b is None or (a == b)
        if c is None:
            return [a, b] if a is not None and b is not None else []
        if isinstance(c, B.Line) or (not closed and isinstance(c, B.Line)):
            return [a, b]
        per = getattr(c, 'periodic', None)
        if closed and a is None:   # courbe fermée sans sommet : cercle complet
            if per:
                n = max(8, int(math.ceil(per / c.step(tol))))
                ts = [k * per / n * csense for k in range(n)]
                return [c.eval(t) for t in ts] + [None]   # None : boucle fermée
            return []
        ta = c.project(a)
        tb = c.project(b) if not closed else ta
        if per:
            if csense > 0:
                while tb <= ta + (1e-12 if not closed else 0.5 * per):
                    tb += per
                if closed:
                    tb = ta + per
            else:
                while tb >= ta - (1e-12 if not closed else 0.5 * per):
                    tb -= per
                if closed:
                    tb = ta - per
        step = c.step(tol)
        if isinstance(c, B.Polyline):
            lo, hi = sorted((ta, tb))
            inner = [k for k in range(int(math.floor(lo)) + 1, int(math.ceil(hi)))]
            if ta > tb:
                inner.reverse()
            mids = [c.pts[k] for k in inner]
        else:
            n = max(1, int(math.ceil(abs(tb - ta) / step))) if step else 1
            n = min(n, 2000)
            mids = [c.eval(ta + (tb - ta) * k / n) for k in range(1, n)]
        return [a] + mids + [b]

    def _consistent(self, pts, curve):
        """Sampled points stay near the curve and the polyline is not much longer than the straight line… unless closed."""
        if len(pts) < 3 or pts[-1] is None:
            return True
        chord = B.length(B.sub(pts[-1], pts[0]))
        L = sum(B.length(B.sub(q, p)) for p, q in zip(pts, pts[1:]))
        return chord < 1e-12 or L < 20 * chord + 1e-6

    def _count(self, tc):
        """Number of segments for a trimmed curve: chord deviation under the tolerance (measured on 16 samples)."""
        pts = tc.points(tc.p1, 16)
        dev = 0.0
        for k in range(1, 16, 2):
            dev = max(dev, B.length(B.sub(pts[k], B.lerp(pts[k - 1], pts[k + 1], 0.5))))
        if dev < 1e-12:
            return 1
        seg = 8 * math.sqrt(dev / self.tol)   # flèche ∝ (longueur)² : on raffine les 8 segments mesurés
        return max(1, min(int(math.ceil(seg)), 2000))

    def loop_points(self, fins):
        """3D boundary of a loop (closed, first point not repeated)."""
        out = []
        for fn in fins:
            fin = self.N[fn]
            ei = fin.f[6]
            if ei not in self.N:
                continue
            pts, pos = self.edge_points(ei)
            if not pts:
                continue
            closed_ring = pts[-1] is None
            seg = pts[:-1] if closed_ring else pts
            if fn != pos:
                seg = list(reversed(seg))
                if closed_ring:
                    seg = [seg[-1]] + seg[:-1]
            out += seg if not out else (seg[1:] if out and out[-1] == seg[0] else seg)
        if len(out) > 1 and out[0] == out[-1]:
            out.pop()
        return out

    # ------------------------------------------------------------ faces

    def mesh_face(self, face):
        try:
            return self._mesh_face(face)
        except Exception as e:   # noqa: BLE001 — une face qui échoue ne doit pas arrêter l'import
            self.stats["errors"].append((face.index, "face", f"{type(e).__name__}: {e}"))
            self.stats["fallback"] += 1
            fm = FaceMesh(face.index, self.M.body_of(face))
            fm.fallback = True
            return fm

    def _mesh_face(self, face):
        fm = FaceMesh(face.index, self.M.body_of(face))
        surf = self.M.surface(face.f[7])
        loops = [self.loop_points(fins) for fins in self.M.loops(face)]
        loops = [l for l in loops if len(l) >= 2]
        self.stats["faces"] += 1
        if not loops:
            return fm
        fsense = _sense(face.f[8]) * self._surf_sense(face.f[7])
        try:
            if surf is None:
                raise B.GeomError("surface not supported")
            self._cdt(fm, surf, loops, fsense)
        except Exception as e:   # noqa: BLE001 — repli : face plane au mieux
            self.stats["errors"].append((face.index, self.N[face.f[7]].name, f"{type(e).__name__}: {e}"))
            self._fallback(fm, loops, fsense)
        return fm

    def _surf_sense(self, si):
        n = self.N.get(si)
        return _sense(n.f[6]) if n is not None and n.f and isinstance(n.f[6], bytes) else 1

    def _cdt(self, fm, surf, loops, fsense):
        pu, pv = getattr(surf, 'per_u', None), getattr(surf, 'per_v', None)
        uvl = []
        for l in loops:
            uv = [list(surf.project(p)) for p in l]
            for k in range(1, len(uv)):   # déroulé continu des directions périodiques
                if pu:
                    uv[k][0] = B.wrap(uv[k][0], uv[k - 1][0], pu)
                if pv:
                    uv[k][1] = B.wrap(uv[k][1], uv[k - 1][1], pv)
            uvl.append(uv)
        pts3, pts2, polys = [], [], []

        def wraps(uv, i, per):
            if not per:
                return 0
            d = uv[0][i] - B.wrap(uv[-1][i], uv[0][i], per) + (uv[-1][i] - uv[0][i])
            return int(round((uv[-1][i] - uv[0][i] + (B.wrap(uv[0][i], uv[-1][i], per) - uv[-1][i])) / per)) if per else 0

        # boucles qui font le tour (bande de cylindre, calotte) : coupées sur une couture
        wrap_dir = None
        winding = []
        for uv in uvl:
            w = [0, 0]
            for i, per in ((0, pu), (1, pv)):
                if per:
                    end = B.wrap(uv[0][i], uv[-1][i], per)   # le premier point, ramené à côté du dernier
                    w[i] = int(round((end - uv[0][i]) / per))
            winding.append(w)
            if w[0] and wrap_dir is None:
                wrap_dir = 0
            elif w[1] and wrap_dir is None:
                wrap_dir = 1
        if wrap_dir is None:
            # recentrer les trous sur la boucle extérieure (même période)
            ref = max(range(len(uvl)), key=lambda k: _area(uvl[k]))
            cu = sum(p[0] for p in uvl[ref]) / len(uvl[ref])
            cv = sum(p[1] for p in uvl[ref]) / len(uvl[ref])
            for k, uv in enumerate(uvl):
                if k == ref:
                    continue
                mu = sum(p[0] for p in uv) / len(uv)
                mv = sum(p[1] for p in uv) / len(uv)
                du = (B.wrap(mu, cu, pu) - mu) if pu else 0.0
                dv = (B.wrap(mv, cv, pv) - mv) if pv else 0.0
                for p in uv:
                    p[0] += du
                    p[1] += dv
            for l, uv in zip(loops, uvl):
                base = len(pts2)
                pts3 += l
                pts2 += uv
                polys.append(list(range(base, base + len(uv))))
        else:
            i, per = wrap_dir, (pu if wrap_dir == 0 else pv)
            j = 1 - i
            wl = [k for k, w in enumerate(winding) if w[i]]
            others = [k for k, w in enumerate(winding) if not w[i]]
            s = uvl[wl[0]][0][i]
            chains = []
            for k in wl:
                c3, c2 = _cut(loops[k], uvl[k], i, s, per, surf)
                chains.append((c3, c2))
            if len(chains) == 2:
                (a3, a2), (b3, b2) = chains
                ring3 = a3 + list(reversed(b3))
                ring2 = a2 + list(reversed(b2))
            elif len(chains) == 1:   # calotte : fermée par le pôle (ligne v = ±pi r)
                a3, a2 = chains[0]
                pole = self._pole(surf, a2, j, fsense)
                end2 = [a2[-1][0], a2[-1][1]]
                beg2 = [a2[0][0], a2[0][1]]
                end2[j] = pole
                beg2[j] = pole
                p_end, p_beg = surf.eval(*end2), surf.eval(*beg2)
                ring3 = a3 + [p_end, p_beg]
                ring2 = a2 + [end2, beg2]
            else:
                raise B.GeomError(f"{len(chains)} wrapping loops")
            pts3 += ring3
            pts2 += ring2
            polys.append(list(range(len(ring2))))
            lo = s
            for k in others:
                uv = uvl[k]
                m = sum(p[i] for p in uv) / len(uv)
                d = (B.wrap(m, lo + per / 2, per) - m)
                base = len(pts2)
                pts3 += loops[k]
                pts2 += [[p[0] + (d if i == 0 else 0), p[1] + (d if i == 1 else 0)] for p in uv]
                polys.append(list(range(base, base + len(uv))))
        # points intérieurs
        us = [p[0] for p in pts2]
        vs = [p[1] for p in pts2]
        box = (min(us), max(us), min(vs), max(vs))
        su, sv = surf.steps(self.tol, box)
        nb = len(pts2)
        if su or sv:
            su = su or (box[1] - box[0])
            sv = sv or (box[3] - box[2])
            su = max(su, (box[1] - box[0]) / 300, 1e-7)
            sv = max(sv, (box[3] - box[2]) / 300, 1e-7)
            if su < box[1] - box[0] or sv < box[3] - box[2]:
                near = {}
                cell = max(min(su, sv) * 0.6, 1e-9)
                for p in pts2:
                    near[(int(p[0] // cell), int(p[1] // cell))] = True
                nu = int((box[1] - box[0]) / su)
                nv = int((box[3] - box[2]) / sv)
                for a in range(1, nu + 1):
                    for b in range(1, nv + 1):
                        u = box[0] + (box[1] - box[0]) * a / (nu + 1)
                        v = box[2] + (box[3] - box[2]) * b / (nv + 1)
                        ci, cj = int(u // cell), int(v // cell)
                        if any((ci + x, cj + y) in near for x in (-1, 0, 1) for y in (-1, 0, 1)):
                            continue
                        pts2.append([u, v])
                        pts3.append(None)
        edges = []
        for poly in polys:
            edges += [(poly[k], poly[(k + 1) % len(poly)]) for k in range(len(poly))]
        span = max(box[1] - box[0], box[3] - box[2], 1e-9)
        out = delaunay_2d_cdt([Vector(p) for p in pts2], edges, polys, 1, span * 1e-9, True)
        overts, ofaces, orig = out[0], out[2], out[3]
        inside = _Inside([[pts2[k] for k in poly] for poly in polys])   # trous et concavités retirés ici (le mode 3 varie selon la version)
        ofaces = [f for f in ofaces if inside(sum(overts[k].x for k in f) / len(f), sum(overts[k].y for k in f) / len(f))]
        used = {k for f in ofaces for k in f}   # seuls les sommets des triangles gardés (pas les points hors de la face)
        remap = {}
        for k, ov in enumerate(overts):
            if k not in used:
                continue
            src = [o for o in orig[k] if o < len(pts3)]
            p3 = pts3[src[0]] if src and pts3[src[0]] is not None else surf.eval(ov.x, ov.y)
            remap[k] = len(fm.verts)
            fm.verts.append(list(p3))
            fm.normals.append(B.mul(surf.normal(ov.x, ov.y), fsense))
        for f in ofaces:
            for t in range(1, len(f) - 1):
                a, b, c = remap[f[0]], remap[f[t]], remap[f[t + 1]]
                pa, pb, pc = fm.verts[a], fm.verts[b], fm.verts[c]
                nrm = B.cross(B.sub(pb, pa), B.sub(pc, pa))
                ref = B.add(B.add(fm.normals[a], fm.normals[b]), fm.normals[c])
                fm.tris.append((a, b, c) if B.dot(nrm, ref) >= 0 else (a, c, b))

    def _pole(self, surf, chain2, j, fsense):
        if isinstance(surf, B.Sphere):
            half = math.pi / 2 * surf.r
            m = sum(p[j] for p in chain2) / len(chain2)
            return half if m >= 0 else -half
        raise B.GeomError("single wrapping loop on a surface without pole")

    def _fallback(self, fm, loops, fsense):
        """Face as flat as its boundary: triangulated in its mean plane (Newell normal)."""
        fm.fallback = True
        self.stats["fallback"] += 1
        allp = [p for l in loops for p in l]
        n = [0.0, 0.0, 0.0]
        for l in loops:
            for k in range(len(l)):
                a, b = l[k], l[(k + 1) % len(l)]
                n[0] += (a[1] - b[1]) * (a[2] + b[2])
                n[1] += (a[2] - b[2]) * (a[0] + b[0])
                n[2] += (a[0] - b[0]) * (a[1] + b[1])
        n = B.unit(n)
        x = B.unit(B.cross(n, [1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0]))
        y = B.cross(n, x)
        pts2, polys = [], []
        for l in loops:
            base = len(pts2)
            pts2 += [Vector((B.dot(p, x), B.dot(p, y))) for p in l]
            polys.append(list(range(base, base + len(l))))
        edges = [(poly[k], poly[(k + 1) % len(poly)]) for poly in polys for k in range(len(poly))]
        try:
            out = delaunay_2d_cdt(pts2, edges, polys, 1, 1e-9, True)
        except Exception:   # noqa: BLE001
            return
        inside = _Inside([[(pts2[k].x, pts2[k].y) for k in poly] for poly in polys])
        out = (out[0], out[1], [f for f in out[2] if inside(sum(out[0][k].x for k in f) / len(f), sum(out[0][k].y for k in f) / len(f))], out[3])
        remap = {}
        used = {k for f in out[2] for k in f}
        for k, ov in enumerate(out[0]):
            if k not in used:
                continue
            src = [o for o in out[3][k] if o < len(allp)]
            remap[k] = len(fm.verts)
            fm.verts.append(list(allp[src[0]]) if src else B.add(B.mul(x, ov.x), B.mul(y, ov.y)))
            fm.normals.append(n)
        for f in out[2]:
            for t in range(1, len(f) - 1):
                fm.tris.append((remap[f[0]], remap[f[t]], remap[f[t + 1]]))

    def mesh_all(self):
        return [self.mesh_face(f) for f in self.M.faces()]


class _Inside:
    """Even-odd point-in-polygons test (all loops together), edges binned by v for speed."""

    def __init__(self, polys, bins=64):
        segs = []
        for poly in polys:
            for k in range(len(poly)):
                a, b = poly[k], poly[(k + 1) % len(poly)]
                if a[1] != b[1]:
                    segs.append((a[0], a[1], b[0], b[1]))
        vs = [v for s in segs for v in (s[1], s[3])] or [0.0, 1.0]
        self.v0, self.v1 = min(vs), max(vs)
        self.n = bins
        self.h = (self.v1 - self.v0) / bins or 1.0
        self.b = [[] for _ in range(bins)]
        for s in segs:
            lo, hi = sorted((s[1], s[3]))
            for i in range(max(0, int((lo - self.v0) / self.h)), min(bins - 1, int((hi - self.v0) / self.h)) + 1):
                self.b[i].append(s)

    def __call__(self, u, v):
        if v < self.v0 or v > self.v1:
            return False
        inside = False
        for x1, y1, x2, y2 in self.b[min(self.n - 1, int((v - self.v0) / self.h))]:
            if (y1 > v) != (y2 > v) and u < x1 + (v - y1) * (x2 - x1) / (y2 - y1):
                inside = not inside
        return inside


def _area(uv):
    a = 0.0
    for k in range(len(uv)):
        x1, y1 = uv[k]
        x2, y2 = uv[(k + 1) % len(uv)]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2


def _cut(pts3, uv, i, s, per, surf):
    """Wrapping loop cut on the seam u = s (or v): chain going from s to s + per (3D and uv)."""
    n = len(uv)
    seq2 = [list(p) for p in uv] + [None]
    end = list(uv[0])
    end[i] = B.wrap(uv[0][i], uv[-1][i], per)
    seq2[-1] = end
    seq3 = list(pts3) + [pts3[0]]
    if seq2[-1][i] < seq2[0][i]:   # toujours dans le sens croissant
        seq2.reverse()
        seq3.reverse()
    # ramener la séquence pour que s soit dans [début, début + per)
    off = per * math.floor((seq2[0][i] - s) / per)
    for p in seq2:
        p[i] -= off
    # trouver le passage par s + per (ou s) et faire commencer la chaîne là
    target = s + per if seq2[0][i] > s else s
    for k in range(len(seq2) - 1):
        a, b = seq2[k], seq2[k + 1]
        if (a[i] - target) * (b[i] - target) <= 0 and a[i] != b[i]:
            f = (target - a[i]) / (b[i] - a[i])
            c2 = [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f]
            c2[i] = target
            c3 = B.lerp(seq3[k], seq3[k + 1], f)
            if f <= 1e-12:
                c3, c2 = seq3[k], [a[0], a[1]]
                c2[i] = target
                chain2 = seq2[k + 1:-1] + seq2[:k + 1]
                chain3 = seq3[k + 1:-1] + seq3[:k + 1]
            else:
                chain2 = seq2[k + 1:-1] + seq2[:k + 1]
                chain3 = seq3[k + 1:-1] + seq3[:k + 1]
            # chaîne : du passage (en target) jusqu'au passage (en target + per), coordonnées déroulées
            out2 = [c2]
            out3 = [c3]
            last = target
            for p2, p3 in zip(chain2, chain3):
                q = list(p2)
                q[i] = B.wrap(q[i], last, per)
                if q[i] < target - 1e-12:
                    q[i] += per
                out2.append(q)
                out3.append(p3)
                last = q[i]
            e2 = list(c2)
            e2[i] = target + per
            out2.append(e2)
            out3.append(c3)
            base = out2[0][i]
            for p in out2:
                p[i] -= base - s
            return out3, out2
    raise B.GeomError("seam not crossed")
