# SPDX-License-Identifier: GPL-3.0-or-later
"""Meshing of the B-rep faces (Blender's constrained Delaunay triangulation in each face's (u, v)).

Each edge is sampled once (shared by its two faces: no crack); each face is triangulated in its (u, v) plane with its loops as
constraints and a grid of inner points where the surface bends, then carried back onto the surface."""
import math

from mathutils import Vector
from mathutils.geometry import delaunay_2d_cdt

from . import brep as B


class FaceMesh:
    __slots__ = ("verts", "tris", "normals", "face", "body", "fallback", "used")

    def __init__(self, face, body):
        self.verts, self.tris, self.normals, self.face, self.body, self.fallback = [], [], [], face, body, False
        self.used = {}   # arête -> version de ses points utilisée par ce maillage


def _sense(c):
    return -1 if c == b'-' else 1


DEBUG = False
GRID_ASPECT = 4       # une maille de la grille de départ est au plus 4 fois plus longue que large
REFINE_ROUNDS = 6     # passes d'affinage d'une face (chaque passe coupe en deux les arêtes trop loin de la surface)
REFINE_MAX = 60000    # points par face au plus
EDGE_DEPTH = 7        # une arête de bord est coupée en deux au plus 7 fois (128 segments par segment de départ)


class Mesher:
    def __init__(self, model, tol=0.0002, max_edge=None):
        self.M, self.N, self.tol = model, model.N, tol
        self.max_edge = max_edge
        self.edges = {}
        self.changed = {}    # arête -> nombre de fois où des points lui ont été ajoutés après coup
        self._face_edges = ()
        self._adj = None
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
                if a is None:   # arête sans sommet : le cercle entier de la courbe de base
                    if getattr(curve.basis, 'periodic', None):
                        pts = self._sample(curve.basis, None, None, self._curve_sense(e))
                        curve = curve.basis
                        if pts:
                            break
                        pts = None
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
        exact = isinstance(curve, (B.Line, B.Circle, B.Ellipse)) and pts is not None and len(pts) > 0   # déjà au pas de la tolérance
        if not exact and len(pts) >= 2:
            pts = self._refine_edge(ei, pts)
        self.edges[ei] = (pts, pos)
        return self.edges[ei]

    def _adjacent(self, ei):
        """Surfaces of the faces that share edge ei."""
        if self._adj is None:
            self._adj = {}
            for face in self.M.faces():
                surf = self.M.surface(face.f[7])
                if surf is None:
                    continue
                for fins in self.M.loops(face):
                    for fn in fins:
                        e = self.N[fn].f[6]
                        lst = self._adj.setdefault(e, [])
                        if surf not in lst:
                            lst.append(surf)
        return self._adj.get(ei, [])

    def _on_surfaces(self, p, surfs, hints):
        """p brought onto the surfaces (their intersection when there are two) by alternate projections; hints: for each
        surface the (u, v) of a point nearby, or None. Returns (point, the new (u, v) of each surface)."""
        q, uv = p, list(hints)
        for _ in range(2 if len(surfs) > 1 else 1):
            for k, s in enumerate(surfs):
                uv[k] = s.project(q, uv[k]) if getattr(s, 'hint', False) else s.project(q)
                q = s.eval(*uv[k])
        return q, uv

    def _refine_edge(self, ei, pts):
        """Boundary edge given by a few points (chart of an intersection, spline): points added where the straight segments
        leave the surfaces of the two faces by more than the tolerance. The same points serve both faces: no gap."""
        surfs = [s for s in self._adjacent(ei) if getattr(s, 'refine', True)]   # (congés : trop coûteux, l'autre face suffit)
        if all(isinstance(s, B.Plane) for s in surfs):   # between flat faces: straight
            return pts
        ring = pts[-1] is None
        base = pts[:-1] if ring else pts
        tol = self.tol

        def split(a, b, ha, hb, depth):
            L2 = sum((a[c] - b[c]) ** 2 for c in range(3))
            if depth >= EDGE_DEPTH or L2 < 16 * tol * tol:
                return []
            m = [(a[c] + b[c]) / 2 for c in range(3)]
            hm = [None if x is None or y is None else ((x[0] + y[0]) / 2, (x[1] + y[1]) / 2) for x, y in zip(ha, hb)]
            try:
                q, hq = self._on_surfaces(m, surfs, hm)
            except Exception:   # noqa: BLE001 — projection impossible : segment laissé droit
                return []
            d2 = sum((q[c] - m[c]) ** 2 for c in range(3))
            if d2 <= tol * tol or d2 > 0.0625 * L2:   # dans la tolérance, ou projection partie ailleurs
                return []
            t = sum((q[c] - a[c]) * (b[c] - a[c]) for c in range(3)) / L2
            if not 0.25 < t < 0.75:   # point qui glisse le long de l'arête (deux faces tangentes) : il sortirait de l'ordre
                return []
            return split(a, q, ha, hq, depth + 1) + [list(q)] + split(q, b, hq, hb, depth + 1)

        def where(pt):   # (u, v) of a point of the edge on each spline surface (the others need no hint)
            return [s.project(pt) if getattr(s, 'hint', False) and pt is not None else None for s in surfs]
        out = []
        n = len(base)
        try:
            hs = [where(pt) for pt in base]
        except Exception:   # noqa: BLE001
            return pts
        for k in range(n if ring else n - 1):
            a, b = base[k], base[(k + 1) % n]
            out.append(a)
            if a is not None and b is not None:
                out += split(a, b, hs[k], hs[(k + 1) % n], 0)
        if not ring:
            out.append(base[-1])
        return out + ([None] if ring else [])

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
            if isinstance(c, B.Polyline) and len(c.pts) > 3 and B.length(B.sub(c.pts[0], c.pts[-1])) < 1e-7:   # tracé refermé
                ring = [list(q) for q in c.pts[:-1]]
                return (ring if csense > 0 else ring[::-1]) + [None]
            if hasattr(c, 't0') and B.length(B.sub(c.eval(c.t0), c.eval(c.t1))) < 1e-7:   # spline refermée
                n = max(16, int(math.ceil((c.t1 - c.t0) / (c.step(tol) or (c.t1 - c.t0) / 16))))
                ring = [c.eval(c.t0 + (c.t1 - c.t0) * k / n) for k in range(n)]
                return (ring if csense > 0 else ring[::-1]) + [None]
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

    def mesh_face(self, face, other_side=False):
        """Mesh of one face; other_side: on a ring of a surface periodic both ways (torus), take the other band between its two
        boundary circles (the two are possible and only the stored mesh tells which is right)."""
        try:
            return self._mesh_face(face, other_side)
        except Exception as e:   # noqa: BLE001 — une face qui échoue ne doit pas arrêter l'import
            self.stats["errors"].append((face.index, "face", f"{type(e).__name__}: {e}"))
            self.stats["fallback"] += 1
            fm = FaceMesh(face.index, self.M.body_of(face))
            fm.fallback = True
            return fm

    def _mesh_face(self, face, other_side=False):
        self.stats["faces"] += 1
        surf = self.M.surface(face.f[7])
        self._face_edges = [self.N[fn].f[6] for fins in self.M.loops(face) for fn in fins if self.N[fn].f[6] in self.N]
        for attempt in range(6):
            fm = FaceMesh(face.index, self.M.body_of(face))
            loops = [self.loop_points(fins) for fins in self.M.loops(face)]
            loops = [l for l in loops if len(l) >= 2]
            fm.used = {e: self.changed.get(e, 0) for e in self._face_edges}
            if not loops:
                return fm
            fsense = _sense(face.f[8]) * self._surf_sense(face.f[7])
            try:
                if surf is None:
                    raise B.GeomError("surface not supported")
                self._cdt(fm, surf, loops, fsense, other_side, split=attempt < 5)
            except _Split:   # bords recoupés : la face est refaite avec les nouveaux points
                continue
            except Exception as e:   # noqa: BLE001 — repli : face plane au mieux
                self.stats["errors"].append((face.index, self.N[face.f[7]].name, f"{type(e).__name__}: {e}"))
                self._fallback(fm, loops, fsense)
            return fm
        return fm

    def stale(self, fm):
        """True when one of the edges of this face mesh got more points afterwards (for a neighbour): to mesh again."""
        return any(self.changed.get(e, 0) != v for e, v in fm.used.items())

    def _split_chords(self, chords):
        """Boundary chords that cross another one in the face's (u, v) (two arcs meeting at a sharp corner, a chord cutting
        across the next edge): each is cut at its middle brought back onto the curve, for both faces of the edge.
        chords: pairs of 3D points. Returns True when an edge got a point."""
        K = lambda p: (round(p[0] * 1e8), round(p[1] * 1e8), round(p[2] * 1e8))
        done = False
        for a, b in chords:
            ka, kb = K(a), K(b)
            for ei in self._face_edges:
                if ei not in self.edges:
                    continue
                pts, pos = self.edges[ei]
                ring = bool(pts) and pts[-1] is None
                base = pts[:-1] if ring else pts
                ks = [K(p) for p in base]
                n = len(ks)
                hit = next((k for k in range(n if ring else n - 1) if {ks[k], ks[(k + 1) % n]} == {ka, kb}), None)
                if hit is None:
                    continue
                p, q = base[hit], base[(hit + 1) % n]
                L = B.length(B.sub(p, q))
                if L < 0.5 * self.tol or self.changed.get(ei, 0) >= 8:   # déjà court ou déjà recoupée souvent : le croisement
                    break                                                 # vient d'ailleurs (projection en (u, v) dégénérée)
                mid = [(p[c] + q[c]) / 2 for c in range(3)]
                surfs = self._adjacent(ei)
                try:
                    m, _uv = self._on_surfaces(mid, surfs, [None] * len(surfs))
                    m = list(m)
                except Exception:   # noqa: BLE001
                    break
                if B.length(B.sub(m, mid)) > 0.5 * L:   # projection partie ailleurs
                    break
                base = base[:hit + 1] + [m] + base[hit + 1:]
                self.edges[ei] = (base + ([None] if ring else []), pos)
                self.changed[ei] = self.changed.get(ei, 0) + 1
                done = True
                break
        return done

    def _surf_sense(self, si):
        n = self.N.get(si)
        return _sense(n.f[6]) if n is not None and n.f and isinstance(n.f[6], bytes) else 1

    def _cdt(self, fm, surf, loops, fsense, other_side=False, split=False):
        if isinstance(surf, B.Sphere) and len(loops) == 1 and len(loops[0]) >= 3:
            surf = self._cap(surf, loops[0], fsense, other_side)
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
            for k in wl:   # coupée sur un de ses propres points (le plus près de la couture) : aucun point ajouté sur l'arête,
                sk = min((B.wrap(q[i], s, per) for q in uvl[k]), key=lambda c: abs(c - s))   # qui reste celle de la face voisine
                c3, c2 = _cut(loops[k], uvl[k], i, sk, per, surf)
                chains.append((c3, c2))
            if len(chains) == 2:
                (a3, a2), (b3, b2) = chains
                pj = pu if j == 0 else pv
                if other_side and pj:   # band between the two circles going round the other way
                    sh = -pj if sum(p[j] for p in b2) / len(b2) > sum(p[j] for p in a2) / len(a2) else pj
                    b2 = [[p[0] + (sh if j == 0 else 0), p[1] + (sh if j == 1 else 0)] for p in b2]
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
            if not isinstance(surf, (B.Cylinder, B.Cone)):   # (droites le long de l'axe : une rangée suffit)
                su, sv = min(su, GRID_ASPECT * sv), min(sv, GRID_ASPECT * su)   # mailles pas trop allongées : sinon des diagonales qui coupent la courbure
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
        inside = _Inside([[pts2[k] for k in poly] for poly in polys])   # trous et concavités retirés ici (le mode 3 varie selon la version)
        known = {}   # point 3D de chaque point (u, v) : celui du bord quand c'est un point du bord, sinon la surface

        def at(k, ov, orig):
            key = (ov.x, ov.y)
            if key not in known:
                src = [o for o in orig[k] if o < len(pts3)]
                known[key] = list(pts3[src[0]]) if src and pts3[src[0]] is not None else list(surf.eval(ov.x, ov.y))
            return known[key]

        tol2 = self.tol * self.tol
        good = set()
        # (u, v) recentrés quand la face est loin de l'origine devant sa taille (voir plus bas) ; sinon tels quels : avant Blender 5.2
        # la triangulation donne d'autres résultats sur des points déplacés
        ox, oy = (box[0] + box[1]) / 2, (box[2] + box[3]) / 2
        if max(abs(ox), abs(oy)) < 100 * max(box[1] - box[0], box[3] - box[2], 1e-9):
            ox = oy = 0.0
        if split:   # bord qui se recoupe en (u, v) (cordes de deux arcs qui se touchent dans un angle aigu) : recoupé avant
            bad = _crossings(pts2, polys, min(0.05 * self.tol, 1e-3 * span))
            chords = [(pts3[x], pts3[y]) for x, y in bad if pts3[x] is not None and pts3[y] is not None]
            if chords and self._split_chords(chords):
                raise _Split()
        for rnd in range(REFINE_ROUNDS + 1):
            # (u, v) ramenés près de 0 : les Vector sont en simple précision, et loin de l'origine (plan à 20 m) des points
            # voisins d'un bord presque droit se confondaient — triangles plats écartés, points du bord perdus, fentes
            out = delaunay_2d_cdt([Vector((p[0] - ox, p[1] - oy)) for p in pts2], edges, polys, 1, span * 1e-9, True)
            ofaces, orig = out[2], out[3]
            if rnd == 0 and split:
                # points créés par le maillage sur un bord (deux cordes du bord qui se croisent) : la face voisine ne les a
                # pas, ce serait une fente ; on recoupe ces cordes et on recommence
                crossed = [k for k, o in enumerate(orig) if not any(i < len(pts2) for i in o)]
                if crossed:
                    bad = set(crossed)
                    chords = []
                    for (i, j), oe in zip(out[1], out[4]):
                        if (i in bad or j in bad) and oe:
                            for ie in oe:
                                if ie < len(edges):
                                    x, y = edges[ie]
                                    if pts3[x] is not None and pts3[y] is not None:
                                        chords.append((pts3[x], pts3[y]))
                    if chords and self._split_chords(chords):
                        raise _Split()
            overts = [_UV(pts2[o[0]] if o and o[0] < len(pts2) else (v.x + ox, v.y + oy)) for v, o in zip(out[0], orig)]
            ofaces = [f for f in ofaces if inside(sum(overts[k].x for k in f) / len(f), sum(overts[k].y for k in f) / len(f))]
            if rnd == REFINE_ROUNDS or len(pts2) > REFINE_MAX or not getattr(surf, 'refine', True):
                break
            # arêtes intérieures dont le milieu s'écarte de la surface de plus que la tolérance : un point de plus au milieu
            count = {}
            for f in ofaces:
                for t in range(len(f)):
                    e = (f[t], f[(t + 1) % len(f)])
                    e = e if e[0] < e[1] else (e[1], e[0])
                    count[e] = count.get(e, 0) + 1
            more, worst = [], []
            for (i, j), n in count.items():
                if n < 2:   # bord de la face : ses points viennent de l'arête, partagée avec la face voisine
                    continue
                a, b = overts[i], overts[j]
                ek = (a.x, a.y, b.x, b.y) if (a.x, a.y) < (b.x, b.y) else (b.x, b.y, a.x, a.y)
                if ek in good:   # déjà dans la tolérance à une passe précédente
                    continue
                pa, pb = at(i, a, orig), at(j, b, orig)
                um, vm = (a.x + b.x) / 2, (a.y + b.y) / 2
                q = surf.eval(um, vm)
                d2 = sum((q[c] - (pa[c] + pb[c]) / 2) ** 2 for c in range(3))
                if d2 <= tol2:
                    good.add(ek)
                elif d2 < 0.25 * sum((pa[c] - pb[c]) ** 2 for c in range(3)):   # (écart absurde : projection ratée, on n'insiste pas)
                    more.append([um, vm])
                    if DEBUG:
                        worst.append(d2 ** 0.5)
                        if d2 > 2.5e-5 and len(worst) < 400 and rnd >= 4:
                            print('   DEBUG arête', (round(a.x, 5), round(a.y, 4)), (round(b.x, 5), round(b.y, 4)), 'L mm', round(sum((pa[c] - pb[c]) ** 2 for c in range(3)) ** 0.5 * 1000, 1), 'pa', [round(x * 1000) for x in pa], 'surf(a)', [round(x * 1000) for x in surf.eval(a.x, a.y)], 'orig', orig[i][:3], len(pts3))
                    known[(um, vm)] = list(q)
            if DEBUG:
                print('  DEBUG passe', rnd, 'points', len(pts2), 'ajouts', len(more), 'pire écart mm', round(max([0.0] + worst) * 1000, 4))
            if not more:
                break
            pts2 += more
            pts3 += [None] * len(more)
        used = {k for f in ofaces for k in f}   # seuls les sommets des triangles gardés (pas les points hors de la face)
        remap = {}
        for k, ov in enumerate(overts):
            if k not in used:
                continue
            remap[k] = len(fm.verts)
            fm.verts.append(at(k, ov, orig))
            fm.normals.append(B.mul(surf.normal(ov.x, ov.y), fsense))
        for f in ofaces:
            for t in range(1, len(f) - 1):
                a, b, c = remap[f[0]], remap[f[t]], remap[f[t + 1]]
                pa, pb, pc = fm.verts[a], fm.verts[b], fm.verts[c]
                nrm = B.cross(B.sub(pb, pa), B.sub(pc, pa))
                ref = B.add(B.add(fm.normals[a], fm.normals[b]), fm.normals[c])
                fm.tris.append((a, b, c) if B.dot(nrm, ref) >= 0 else (a, c, b))

    def _cap(self, surf, loop, fsense, other_side):
        """Sphere face inside one loop: the same sphere with its axis through the middle of the face, so that the loop goes
        round the axis and the face is the cap at the top. The side is the one on the left of the loop seen from outside the
        face (the loop of a face runs anticlockwise around it)."""
        w = [0.0, 0.0, 0.0]
        n = len(loop)
        for k in range(n):
            w = B.add(w, B.cross(B.sub(loop[k], surf.c), B.sub(loop[(k + 1) % n], surf.c)))
        if B.length(w) < 1e-18:
            return surf
        side = fsense * surf.sign * (-1.0 if other_side else 1.0)
        a = B.mul(B.unit(w), side)
        x = B.unit(B.cross(a, [1.0, 0.0, 0.0] if abs(a[0]) < 0.9 else [0.0, 1.0, 0.0]))
        cap = B.Sphere(surf.c, surf.r, a, x)
        cap.sign, cap.top = surf.sign, True
        return cap

    def _pole(self, surf, chain2, j, fsense):
        if isinstance(surf, B.Sphere) and getattr(surf, 'top', False):
            return math.pi / 2 * surf.r
        if isinstance(surf, B.Sphere):
            half = math.pi / 2 * surf.r
            m = sum(p[j] for p in chain2) / len(chain2)
            return half if m >= 0 else -half
        if isinstance(surf, B.Cone) and j == 1 and abs(surf.k) > 1e-12:   # pointe du cône : rayon nul
            return -surf.r / surf.k
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
        o = allp[0]   # près de 0 (simple précision des Vector)
        for l in loops:
            base = len(pts2)
            pts2 += [Vector((B.dot(B.sub(p, o), x), B.dot(B.sub(p, o), y))) for p in l]
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
            fm.verts.append(list(allp[src[0]]) if src else B.add(o, B.add(B.mul(x, ov.x), B.mul(y, ov.y))))
            fm.normals.append(n)
        for f in out[2]:
            for t in range(1, len(f) - 1):
                fm.tris.append((remap[f[0]], remap[f[t]], remap[f[t + 1]]))

    def mesh_all(self):
        return [self.mesh_face(f) for f in self.M.faces()]


class _Split(Exception):
    """Raised when boundary chords were cut: the face is meshed again."""


class _UV:
    """A CDT output point back in the face's (u, v), in double precision."""
    __slots__ = ("x", "y")

    def __init__(self, p):
        self.x, self.y = p[0], p[1]


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


def _crossings(pts, polys, eps):
    """Boundary segments (pairs of point indices) that cross or touch a segment not next to them (closer than eps)."""
    segs = [(poly[k], poly[(k + 1) % len(poly)]) for poly in polys for k in range(len(poly))]
    if len(segs) < 4:
        return []
    L = sum(math.hypot(pts[b][0] - pts[a][0], pts[b][1] - pts[a][1]) for a, b in segs) / len(segs)
    cell = max(L, eps * 4, 1e-12)
    grid = {}
    for n, (a, b) in enumerate(segs):
        (x0, x1), (y0, y1) = sorted((pts[a][0], pts[b][0])), sorted((pts[a][1], pts[b][1]))
        i0, i1, j0, j1 = int((x0 - eps) // cell), int((x1 + eps) // cell), int((y0 - eps) // cell), int((y1 + eps) // cell)
        if (i1 - i0 + 1) * (j1 - j0 + 1) > 400:   # segment très long devant les autres : cases grossières
            i0, i1, j0, j1 = i0, min(i1, i0 + 19), j0, min(j1, j0 + 19)
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                grid.setdefault((i, j), []).append(n)

    def dist(p, a, b):
        dx, dy = b[0] - a[0], b[1] - a[1]
        d2 = dx * dx + dy * dy
        t = 0.0 if d2 == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / d2))
        return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)

    def orient(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    out, seen = set(), set()
    for lst in grid.values():
        for x in range(len(lst)):
            for y in range(x + 1, len(lst)):
                m, n = lst[x], lst[y]
                if (m, n) in seen:
                    continue
                seen.add((m, n))
                (a, b), (c, d) = segs[m], segs[n]
                if len({a, b, c, d}) < 4:   # segments qui se suivent
                    continue
                pa, pb, pc, pd = pts[a], pts[b], pts[c], pts[d]
                # deux segments qui partagent un bout (point de couture ou de pôle en double) ne se touchent pas : seul un
                # vrai croisement compte alors
                shared = min(math.hypot(p[0] - q[0], p[1] - q[1]) for p in (pa, pb) for q in (pc, pd)) < eps
                if (orient(pa, pb, pc) * orient(pa, pb, pd) < 0 and orient(pc, pd, pa) * orient(pc, pd, pb) < 0) or \
                        (not shared and min(dist(pc, pa, pb), dist(pd, pa, pb), dist(pa, pc, pd), dist(pb, pc, pd)) < eps):
                    out.add(segs[m])
                    out.add(segs[n])
    return sorted(out)


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
            if f >= 1 - 1e-9:   # la couture passe par le point suivant : coupée sur ce point même (un point recalculé
                k = (k + 1) % n  # différerait de l'arrondi près : la face voisine ne l'aurait pas)
                a = seq2[k]
                f = 0.0
            c2 = [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f]
            c2[i] = target
            c3 = B.lerp(seq3[k], seq3[k + 1], f)
            if f <= 1e-9:
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
