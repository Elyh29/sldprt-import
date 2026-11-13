# SPDX-License-Identifier: GPL-3.0-or-later
"""B-rep model from the Parasolid nodes: faces with their surface and loops, edges with their curve; geometry evaluation.

Pure Python (lists of 3 floats, metres). Each surface gives eval(u, v), project(p) -> (u, v) and normal(u, v); (u, v) are
roughly in metres so that a 2D triangulation of them is well shaped. Each curve gives eval(t) and project(p) -> t."""
import bisect
import math
import struct

UNSET = -3.0e13


# ---------------------------------------------------------------- vecteurs

def add(a, b): return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
def sub(a, b): return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
def mul(a, s): return [a[0] * s, a[1] * s, a[2] * s]
def dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
def length(a): return math.sqrt(dot(a, a))


def unit(a):
    n = length(a)
    return mul(a, 1.0 / n) if n > 1e-300 else [0.0, 0.0, 1.0]


def lerp(a, b, s): return [a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s]


def wrap(a, ref, period):
    """a moved by whole periods to be the closest to ref."""
    return a + period * round((ref - a) / period)


class GeomError(Exception):
    pass


# ---------------------------------------------------------------- NURBS

class Nurbs:
    """B-spline (rational or not) of one or two parameters: knots expanded from (distinct knots, multiplicities)."""

    @staticmethod
    def knots(vals, mults):
        out = []
        for v, m in zip(vals, mults):
            out += [v] * m
        return out

    @staticmethod
    def span(U, p, n, t):
        """Index i with U[i] <= t < U[i+1], p <= i <= n - 1 (n = number of control points)."""
        if t >= U[n]:
            return n - 1
        if t <= U[p]:
            return p
        return min(max(bisect.bisect_right(U, t) - 1, p), n - 1)

    @staticmethod
    def basis(U, p, i, t):
        N = [1.0] + [0.0] * p
        left, right = [0.0] * (p + 1), [0.0] * (p + 1)
        for j in range(1, p + 1):
            left[j] = t - U[i + 1 - j]
            right[j] = U[i + j] - t
            saved = 0.0
            for r in range(j):
                den = right[r + 1] + left[j - r]
                tmp = N[r] / den if den != 0 else 0.0
                N[r] = saved + right[r + 1] * tmp
                saved = left[j - r] * tmp
            N[j] = saved
        return N


class NurbsCurve:
    def __init__(self, deg, U, cps, rational, dim):
        self.p, self.U, self.cps, self.rat, self.dim = deg, U, cps, rational, dim
        self.n = len(cps)
        self.t0, self.t1 = U[deg], U[self.n]
        self._samples = None

    def eval(self, t):
        t = min(max(t, self.t0), self.t1)
        i = Nurbs.span(self.U, self.p, self.n, t)
        N = Nurbs.basis(self.U, self.p, i, t)
        acc = [0.0] * len(self.cps[0])
        for k in range(self.p + 1):
            c = self.cps[i - self.p + k]
            for d in range(len(acc)):
                acc[d] += N[k] * c[d]
        if self.rat:
            w = acc[-1] or 1.0
            acc = [a / w for a in acc[:-1]]
        return acc

    def project(self, pt):
        if self._samples is None:
            m = max(32, self.n * 8)
            ts = [self.t0 + (self.t1 - self.t0) * k / m for k in range(m + 1)]
            self._samples = [(t, self.eval(t)) for t in ts]
        best = min(self._samples, key=lambda s: sum((s[1][d] - pt[d]) ** 2 for d in range(len(pt))))
        t = best[0]
        h = (self.t1 - self.t0) / (len(self._samples) - 1)
        for _ in range(40):   # recherche dorée autour du meilleur échantillon
            a, b = max(self.t0, t - h), min(self.t1, t + h)
            for _k in range(20):
                c1, c2 = a + (b - a) * 0.382, a + (b - a) * 0.618
                d1 = sum((x - y) ** 2 for x, y in zip(self.eval(c1), pt))
                d2 = sum((x - y) ** 2 for x, y in zip(self.eval(c2), pt))
                if d1 < d2:
                    b = c2
                else:
                    a = c1
            t = (a + b) / 2
            break
        return t


class NurbsSurface:
    def __init__(self, du, dv, U, V, grid, rational):
        self.pu, self.pv, self.U, self.V, self.grid, self.rat = du, dv, U, V, grid, rational
        self.nu, self.nv = len(grid), len(grid[0])
        self.u0, self.u1, self.v0, self.v1 = U[du], U[self.nu], V[dv], V[self.nv]
        self._samp = None
        # échelle : longueur moyenne d'une ligne iso, pour un (u, v) à peu près en mètres
        self.su = self._len(True) / max(self.u1 - self.u0, 1e-12)
        self.sv = self._len(False) / max(self.v1 - self.v0, 1e-12)

    def _len(self, along_u):
        tot = 0.0
        for f in (0.25, 0.5, 0.75):
            prev = None
            for k in range(9):
                s = k / 8
                p = self.raw(self.u0 + (self.u1 - self.u0) * s, self.v0 + (self.v1 - self.v0) * f) if along_u else \
                    self.raw(self.u0 + (self.u1 - self.u0) * f, self.v0 + (self.v1 - self.v0) * s)
                if prev is not None:
                    tot += length(sub(p, prev))
                prev = p
        return max(tot / 3, 1e-9)

    def raw(self, u, v):   # hors du domaine : le dernier polynôme est prolongé (Parasolid étend les surfaces)
        iu = Nurbs.span(self.U, self.pu, self.nu, u)
        iv = Nurbs.span(self.V, self.pv, self.nv, v)
        Nu = Nurbs.basis(self.U, self.pu, iu, u)
        Nv = Nurbs.basis(self.V, self.pv, iv, v)
        dim = len(self.grid[0][0])
        acc = [0.0] * dim
        for a in range(self.pu + 1):
            row = self.grid[iu - self.pu + a]
            for b in range(self.pv + 1):
                w = Nu[a] * Nv[b]
                c = row[iv - self.pv + b]
                for d in range(dim):
                    acc[d] += w * c[d]
        if self.rat:
            w = acc[3] or 1.0
            return [acc[0] / w, acc[1] / w, acc[2] / w]
        return acc[:3]

    def eval(self, u, v):
        return self.raw(u / self.su, v / self.sv)

    def normal(self, u, v):
        h = 1e-6
        p = self.eval(u, v)
        return unit(cross(sub(self.eval(u + h, v), p), sub(self.eval(u, v + h), p)))

    def steps(self, tol, box):
        return _curv_step(self, tol, box)

    def project(self, p):
        if self._samp is None:   # grille de départ : environ 3 points par travée de nœuds
            mu = min(max(3 * (self.nu - self.pu), 8), 90)
            mv = min(max(3 * (self.nv - self.pv), 8), 90)
            self._samp = [(self.u0 + (self.u1 - self.u0) * i / mu, self.v0 + (self.v1 - self.v0) * j / mv) for i in range(mu + 1) for j in range(mv + 1)]
            self._samp = [(u, v, self.raw(u, v)) for u, v in self._samp]
        eu, ev = (self.u1 - self.u0) * 0.25, (self.v1 - self.v0) * 0.25
        return _multi_newton(self.raw, p, self._samp, (self.u0 - eu, self.u1 + eu, self.v0 - ev, self.v1 + ev), (self.su, self.sv))


def _multi_newton(f, p, samp, box, scale, starts=3):
    d2 = lambda q: (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 + (q[2] - p[2]) ** 2
    best, bd = None, 1e300
    for u, v, _q in sorted(samp, key=lambda e: d2(e[2]))[:starts]:
        uu, vv = _newton(f, p, u, v, box)
        d = d2(f(uu, vv))
        if d < bd:
            best, bd = (uu, vv), d
    return best[0] * scale[0], best[1] * scale[1]


def _newton(f, p, u, v, box, scale=(1.0, 1.0), it=30):
    """(u, v) minimising |f(u, v) - p| (Gauss-Newton, finite differences), returned multiplied by scale."""
    u0, u1, v0, v1 = box
    hu, hv = (u1 - u0) * 1e-6, (v1 - v0) * 1e-6
    for _ in range(it):
        q = f(u, v)
        r = sub(p, q)
        fu = mul(sub(f(u + hu, v), q), 1 / hu)
        fv = mul(sub(f(u, v + hv), q), 1 / hv)
        a, b, c = dot(fu, fu), dot(fu, fv), dot(fv, fv)
        det = a * c - b * b
        if abs(det) < 1e-30:
            break
        ru, rv = dot(fu, r), dot(fv, r)
        du, dv = (c * ru - b * rv) / det, (a * rv - b * ru) / det
        u, v = min(max(u + du, u0), u1), min(max(v + dv, v0), v1)
        if abs(du) < (u1 - u0) * 1e-10 and abs(dv) < (v1 - v0) * 1e-10:
            break
    return u * scale[0], v * scale[1]


# ---------------------------------------------------------------- courbes

class Line:
    periodic = None

    def __init__(self, o, d):
        self.o, self.d = list(o), unit(d)

    def eval(self, t): return add(self.o, mul(self.d, t))
    def project(self, p): return dot(sub(p, self.o), self.d)
    def step(self, tol): return None   # droite : aucun point intermédiaire


class Circle:
    periodic = 2 * math.pi

    def __init__(self, c, n, x, r):
        self.c, self.n, self.x, self.r = list(c), unit(n), unit(x), r
        self.y = cross(self.n, self.x)

    def eval(self, t):
        return add(self.c, add(mul(self.x, self.r * math.cos(t)), mul(self.y, self.r * math.sin(t))))

    def project(self, p):
        d = sub(p, self.c)
        return math.atan2(dot(d, self.y), dot(d, self.x)) % (2 * math.pi)

    def step(self, tol):
        return _angle_step(self.r, tol)


class Ellipse(Circle):
    def __init__(self, c, n, x, a, b):
        super().__init__(c, n, x, a)
        self.b = b

    def eval(self, t):
        return add(self.c, add(mul(self.x, self.r * math.cos(t)), mul(self.y, self.b * math.sin(t))))

    def project(self, p):
        d = sub(p, self.c)
        return math.atan2(dot(d, self.y) / self.b, dot(d, self.x) / self.r) % (2 * math.pi)

    def step(self, tol):
        return _angle_step(max(self.r, self.b) ** 2 / max(min(self.r, self.b), 1e-9), tol)


class Polyline:
    """Intersection curve given by its chart points."""
    periodic = None
    closed = False

    def __init__(self, pts):
        self.pts = pts

    def eval(self, t):
        if self.closed:
            t %= (len(self.pts) - 1)
        i = min(max(int(math.floor(t)), 0), len(self.pts) - 2)
        return lerp(self.pts[i], self.pts[i + 1], t - i)

    def project(self, p):
        best, bt = 1e300, 0.0
        for i in range(len(self.pts) - 1):
            a, b = self.pts[i], self.pts[i + 1]
            ab = sub(b, a)
            L = dot(ab, ab)
            s = min(max(dot(sub(p, a), ab) / L, 0.0), 1.0) if L > 0 else 0.0
            q = lerp(a, b, s)
            d = dot(sub(p, q), sub(p, q))
            if d < best:
                best, bt = d, i + s
        return bt

    def step(self, tol):
        return 1.0   # un point par point du tracé


class CurveOnSurface:
    def __init__(self, surf, bcurve):
        self.s, self.b = surf, bcurve
        self.periodic = None
        self.t0, self.t1 = bcurve.t0, bcurve.t1
        self._samples = None

    def eval(self, t):
        u, v = self.b.eval(t)[:2]
        return self.s.raw(u, v) if hasattr(self.s, 'raw') else self.s.eval(u, v)

    def project(self, p):
        if self._samples is None:
            m = 200
            self._samples = [(self.t0 + (self.t1 - self.t0) * k / m) for k in range(m + 1)]
            self._samples = [(t, self.eval(t)) for t in self._samples]
        t = min(self._samples, key=lambda s: dot(sub(s[1], p), sub(s[1], p)))[0]
        h = (self.t1 - self.t0) / 200
        a, b = max(self.t0, t - h), min(self.t1, t + h)
        for _ in range(30):
            c1, c2 = a + (b - a) * 0.382, a + (b - a) * 0.618
            if length(sub(self.eval(c1), p)) < length(sub(self.eval(c2), p)):
                b = c2
            else:
                a = c1
        return (a + b) / 2

    def step(self, tol):
        return (self.t1 - self.t0) / 64


class Trimmed:
    """Part of a basis curve between two parameters, with its two end points (given by the file)."""
    periodic = None

    def __init__(self, basis, p1, p2, t1, t2):
        self.basis, self.p1, self.p2, self.t1, self.t2 = basis, list(p1), list(p2), t1, t2

    def ok(self, tol=1e-5):
        try:
            return length(sub(self.basis.eval(self.t1), self.p1)) < tol and length(sub(self.basis.eval(self.t2), self.p2)) < tol
        except Exception:   # noqa: BLE001
            return False

    def points(self, a, n):
        """n + 1 points from the end nearest to a to the other one."""
        t1, t2 = (self.t1, self.t2) if length(sub(a, self.p1)) <= length(sub(a, self.p2)) else (self.t2, self.t1)
        return [self.basis.eval(t1 + (t2 - t1) * k / n) for k in range(n + 1)]


class Curve3D(NurbsCurve):
    def step(self, tol):
        return (self.t1 - self.t0) / max(8, self.n * 2)


def _angle_step(r, tol):
    if r <= tol:
        return math.pi / 2
    return min(2 * math.acos(max(-1.0, 1 - tol / r)), math.pi / 4)


# ---------------------------------------------------------------- surfaces

class Plane:
    per_u = per_v = None

    def __init__(self, o, n, x):
        self.o, self.n, self.x = list(o), unit(n), unit(x)
        self.y = cross(self.n, self.x)

    def eval(self, u, v): return add(self.o, add(mul(self.x, u), mul(self.y, v)))
    def project(self, p): d = sub(p, self.o); return dot(d, self.x), dot(d, self.y)
    def normal(self, u, v): return self.n
    def steps(self, tol, box): return None, None


class Cylinder:
    per_v = None

    def __init__(self, o, a, r, x):
        self.o, self.a, self.r, self.x = list(o), unit(a), r, unit(x)
        self.y = cross(self.a, self.x)
        self.per_u = 2 * math.pi * r

    def eval(self, u, v):
        t = u / self.r
        return add(self.o, add(mul(self.a, v), add(mul(self.x, self.r * math.cos(t)), mul(self.y, self.r * math.sin(t)))))

    def project(self, p):
        d = sub(p, self.o)
        return math.atan2(dot(d, self.y), dot(d, self.x)) * self.r, dot(d, self.a)

    def normal(self, u, v):
        t = u / self.r
        return add(mul(self.x, math.cos(t)), mul(self.y, math.sin(t)))

    def steps(self, tol, box): return _angle_step(self.r, tol) * self.r, None


class Cone:
    per_v = None

    def __init__(self, o, a, r, s, c, x):
        self.o, self.a, self.r, self.sn, self.cs, self.x = list(o), unit(a), r, s, c, unit(x)
        self.y = cross(self.a, self.x)
        self.k = s / c if c else 0.0   # variation du rayon par unité de hauteur
        self.per_u = 2 * math.pi * max(r, 1e-9)

    def radius(self, h): return self.r + h * self.k

    def eval(self, u, v):
        t = u / max(self.r, 1e-9)
        rr = self.radius(v)
        return add(self.o, add(mul(self.a, v), add(mul(self.x, rr * math.cos(t)), mul(self.y, rr * math.sin(t)))))

    def project(self, p):
        d = sub(p, self.o)
        return math.atan2(dot(d, self.y), dot(d, self.x)) * max(self.r, 1e-9), dot(d, self.a)

    def normal(self, u, v):
        t = u / max(self.r, 1e-9)
        radial = add(mul(self.x, math.cos(t)), mul(self.y, math.sin(t)))
        return unit(sub(mul(radial, self.cs), mul(self.a, self.sn)))

    def steps(self, tol, box):
        rmax = max(abs(self.radius(box[2])), abs(self.radius(box[3])), 1e-9)
        return _angle_step(rmax, tol) * max(self.r, 1e-9) * min(1.0, max(self.r, 1e-9) / rmax), None


class Sphere:
    per_v = None

    def __init__(self, c, r, a, x):
        self.c, self.r, self.a, self.x = list(c), r, unit(a), unit(x)
        self.y = cross(self.a, self.x)
        self.per_u = 2 * math.pi * r

    def eval(self, u, v):
        t, f = u / self.r, v / self.r
        return add(self.c, mul(add(mul(self.a, math.sin(f)), mul(add(mul(self.x, math.cos(t)), mul(self.y, math.sin(t))), math.cos(f))), self.r))

    def project(self, p):
        d = sub(p, self.c)
        h = dot(d, self.a)
        return math.atan2(dot(d, self.y), dot(d, self.x)) * self.r, math.asin(max(-1.0, min(1.0, h / max(length(d), 1e-300)))) * self.r

    def normal(self, u, v): return unit(sub(self.eval(u, v), self.c))
    def steps(self, tol, box): s = _angle_step(self.r, tol) * self.r; return s, s


class Torus:
    def __init__(self, c, a, R, r, x):
        self.c, self.a, self.R, self.r, self.x = list(c), unit(a), R, r, unit(x)
        self.y = cross(self.a, self.x)
        self.per_u = 2 * math.pi * max(R, 1e-9)
        self.per_v = 2 * math.pi * r

    def eval(self, u, v):
        t, f = u / max(self.R, 1e-9), v / self.r
        radial = add(mul(self.x, math.cos(t)), mul(self.y, math.sin(t)))
        return add(self.c, add(mul(radial, self.R + self.r * math.cos(f)), mul(self.a, self.r * math.sin(f))))

    def project(self, p):
        d = sub(p, self.c)
        t = math.atan2(dot(d, self.y), dot(d, self.x))
        radial = add(mul(self.x, math.cos(t)), mul(self.y, math.sin(t)))
        f = math.atan2(dot(d, self.a), dot(d, radial) - self.R)
        return t * max(self.R, 1e-9), f * self.r

    def normal(self, u, v):
        t, f = u / max(self.R, 1e-9), v / self.r
        radial = add(mul(self.x, math.cos(t)), mul(self.y, math.sin(t)))
        return add(mul(radial, math.cos(f)), mul(self.a, math.sin(f)))

    def steps(self, tol, box):
        return _angle_step(self.R + self.r, tol) * max(self.R, 1e-9), _angle_step(self.r, tol) * self.r


class Generic:
    """Surface given by raw(s, t) on a box: (u, v) = (s * su, t * sv); projection by Gauss-Newton."""
    per_u = per_v = None

    def __init__(self, raw, box, samples=20):
        self.raw, self.box = raw, box
        s0, s1, t0, t1 = box
        self._samp = [(s0 + (s1 - s0) * i / samples, t0 + (t1 - t0) * j / samples) for i in range(samples + 1) for j in range(samples + 1)]
        self._samp = [(s, t, raw(s, t)) for s, t in self._samp]
        def iso(along):
            tot = 0.0
            for f in (0.25, 0.5, 0.75):
                pts = [raw(s0 + (s1 - s0) * k / 8, t0 + (t1 - t0) * f) if along else raw(s0 + (s1 - s0) * f, t0 + (t1 - t0) * k / 8) for k in range(9)]
                tot += sum(length(sub(b, a)) for a, b in zip(pts, pts[1:]))
            return max(tot / 3, 1e-9)
        self.su = iso(True) / max(s1 - s0, 1e-12)
        self.sv = iso(False) / max(t1 - t0, 1e-12)

    def eval(self, u, v): return self.raw(u / self.su, v / self.sv)

    def project(self, p):
        return _multi_newton(self.raw, p, self._samp, self.box, (self.su, self.sv))

    def normal(self, u, v):
        h = 1e-5
        p = self.eval(u, v)
        return unit(cross(sub(self.eval(u + h, v), p), sub(self.eval(u, v + h), p)))

    def steps(self, tol, box):
        return _curv_step(self, tol, box)


class Swept:
    """Section curve moved along a vector: eval(a, b) = section(a) + b × d."""
    per_u = per_v = None

    def __init__(self, sec, d):
        self.sec, self.d = sec, list(d)
        self.dl = max(length(d), 1e-12)
        self.su = 1.0
        if hasattr(sec, 't0'):
            pts = [sec.eval(sec.t0 + (sec.t1 - sec.t0) * k / 16) for k in range(17)]
            self.su = max(sum(length(sub(b, a)) for a, b in zip(pts, pts[1:])), 1e-9) / max(sec.t1 - sec.t0, 1e-12)
        elif isinstance(sec, Circle):
            self.su = sec.r
        if isinstance(sec, Circle):
            self.per_u = 2 * math.pi * self.su
        elif hasattr(sec, 't0') and length(sub(sec.eval(sec.t0), sec.eval(sec.t1))) < 1e-9:
            self.per_u = (sec.t1 - sec.t0) * self.su   # section fermée : surface périodique

    def raw(self, a, b):
        if self.per_u and hasattr(self.sec, 't0'):
            a = self.sec.t0 + (a - self.sec.t0) % (self.sec.t1 - self.sec.t0)
        return add(self.sec.eval(a), mul(self.d, b))
    def eval(self, u, v): return self.raw(u / self.su, v / self.dl)

    def project(self, p):
        b = 0.0
        for _ in range(6):
            a = self.sec.project(sub(p, mul(self.d, b)))
            b = dot(sub(p, self.sec.eval(a)), self.d) / (self.dl * self.dl)
        return a * self.su, b * self.dl

    def normal(self, u, v):
        h = 1e-5
        p = self.eval(u, v)
        return unit(cross(sub(self.eval(u + h, v), p), sub(self.eval(u, v + h), p)))

    def steps(self, tol, box): return _curv_step(self, tol, box)


class Spun:
    """Profile curve turned about an axis: eval(a, θ)."""

    def __init__(self, prof, base, ax):
        self.prof, self.base, self.ax = prof, list(base), unit(ax)
        if hasattr(prof, 't0'):
            t0, t1 = prof.t0, prof.t1
        elif isinstance(prof, Circle):
            t0, t1 = 0.0, 2 * math.pi
        else:
            t0, t1 = -1.0, 1.0
        self.t0, self.t1 = t0, t1
        pts = [prof.eval(t0 + (t1 - t0) * k / 32) for k in range(33)]
        self.samples = [(t0 + (t1 - t0) * k / 32, self._hr(q)) for k, q in enumerate(pts)]
        self.su = max(sum(length(sub(b, a)) for a, b in zip(pts, pts[1:])), 1e-9) / max(t1 - t0, 1e-12)
        self.rho = max(max(r for _t, (_h, r, _v) in self.samples), 1e-9)
        self.per_u, self.per_v = None, 2 * math.pi * self.rho

    def _hr(self, q):
        d = sub(q, self.base)
        h = dot(d, self.ax)
        rv = sub(d, mul(self.ax, h))
        return h, length(rv), rv

    def raw(self, a, th):
        h, r, rv = self._hr(self.prof.eval(a))
        e1 = unit(rv) if r > 1e-12 else [0.0, 0.0, 0.0]
        e2 = cross(self.ax, e1)
        return add(self.base, add(mul(self.ax, h), add(mul(e1, r * math.cos(th)), mul(e2, r * math.sin(th)))))

    def eval(self, u, v): return self.raw(u / self.su, v / self.rho)

    def project(self, p):
        h, r, rv = self._hr(p)
        a = min(self.samples, key=lambda s: (s[1][0] - h) ** 2 + (s[1][1] - r) ** 2)[0]
        step = (self.t1 - self.t0) / 32
        lo, hi = max(self.t0, a - step), min(self.t1, a + step)
        for _ in range(30):
            c1, c2 = lo + (hi - lo) * 0.382, lo + (hi - lo) * 0.618
            f = lambda c: (lambda hh: (hh[0] - h) ** 2 + (hh[1] - r) ** 2)(self._hr(self.prof.eval(c)))
            if f(c1) < f(c2):
                hi = c2
            else:
                lo = c1
        a = (lo + hi) / 2
        _h, _r, qv = self._hr(self.prof.eval(a))
        e1 = unit(qv); e2 = cross(self.ax, e1)
        th = math.atan2(dot(rv, e2), dot(rv, e1))
        return a * self.su, th * self.rho

    def normal(self, u, v):
        h = 1e-5
        p = self.eval(u, v)
        return unit(cross(sub(self.eval(u + h, v), p), sub(self.eval(u, v + h), p)))

    def steps(self, tol, box): return _curv_step(self, tol, box)


class Blend:
    """Rolling-ball fillet: ball centre along the spine, arc between the feet on the two supporting surfaces."""
    per_u = per_v = None

    def __init__(self, spine, s1, s2, r):
        self.sp, self.s1, self.s2, self.r = spine, s1, s2, abs(r)
        pts = spine.pts
        self.su = max(sum(length(sub(b, a)) for a, b in zip(pts, pts[1:])) / max(len(pts) - 1, 1), 1e-9)
        if length(sub(pts[0], pts[-1])) < 1e-9:
            self.per_u = (len(pts) - 1) * self.su   # contour fermé
            spine.closed = True
        mid = len(pts) // 2
        _c, d1, d2, ang = self._frame(mid)
        self.sv = max(self.r * ang, 1e-9)

    def _foot(self, s, c):
        return s.eval(*s.project(c))

    def _frame(self, t):
        if self.sp.closed:
            t %= (len(self.sp.pts) - 1)
        c = self.sp.eval(t)
        i = min(max(int(math.floor(t)), 0), len(self.sp.pts) - 2)
        T = unit(sub(self.sp.pts[i + 1], self.sp.pts[i]))   # le centre reste dans le plan de la section
        for _ in range(8):   # centre recalé : à r des deux surfaces (le tracé du fichier n'a que quelques points)
            f1, f2 = self._foot(self.s1, c), self._foot(self.s2, c)
            n1, n2 = sub(c, f1), sub(c, f2)
            l1, l2 = length(n1), length(n2)
            if l1 < 1e-12 or l2 < 1e-12:
                break
            n1, n2 = mul(n1, 1 / l1), mul(n2, 1 / l2)
            e1, e2 = self.r - l1, self.r - l2
            n1, n2 = sub(n1, mul(T, dot(n1, T))), sub(n2, mul(T, dot(n2, T)))
            a11, a22, g = dot(n1, n1), dot(n2, n2), dot(n1, n2)
            den = a11 * a22 - g * g
            if den < 1e-12:
                break
            a, b = (e1 * a22 - g * e2) / den, (e2 * a11 - g * e1) / den
            c = add(c, add(mul(n1, a), mul(n2, b)))
            if abs(e1) < 1e-9 and abs(e2) < 1e-9:
                break
        d1 = unit(sub(self._foot(self.s1, c), c))
        d2 = unit(sub(self._foot(self.s2, c), c))
        ang = math.acos(max(-1.0, min(1.0, dot(d1, d2))))
        return c, d1, d2, ang

    def eval(self, u, v):
        c, d1, d2, ang = self._frame(u / self.su)
        f = v / self.sv
        if ang < 1e-9:
            return add(c, mul(d1, self.r))
        a = ang * f   # arc de d1 vers d2 (interpolation sphérique)
        w1, w2 = math.sin(ang - a) / math.sin(ang), math.sin(a) / math.sin(ang)
        return add(c, mul(add(mul(d1, w1), mul(d2, w2)), self.r))

    def project(self, p):
        t = self.sp.project(p)
        c, d1, d2, ang = self._frame(t)
        d = unit(sub(p, c))
        a = math.acos(max(-1.0, min(1.0, dot(d1, d))))
        if dot(cross(d1, d), cross(d1, d2)) < 0:
            a = -a
        return t * self.su, (a / ang if ang > 1e-9 else 0.0) * self.sv

    def normal(self, u, v):
        c = self._frame(u / self.su)[0]
        return unit(sub(self.eval(u, v), c))

    def steps(self, tol, box): return self.su, _angle_step(self.r, tol) * self.r


def _curv_step(surf, tol, box):
    """Grid steps (u, v) from the bending of the surface measured on a coarse grid."""
    u0, u1, v0, v1 = box
    ku = kv = 0.0
    for i in range(1, 6):
        for j in range(1, 6):
            u = u0 + (u1 - u0) * i / 6
            v = v0 + (v1 - v0) * j / 6
            hu, hv = (u1 - u0) / 12 or 1e-6, (v1 - v0) / 12 or 1e-6
            p = surf.eval(u, v)
            ku = max(ku, length(add(sub(surf.eval(u + hu, v), p), sub(surf.eval(u - hu, v), p))) / hu ** 2)
            kv = max(kv, length(add(sub(surf.eval(u, v + hv), p), sub(surf.eval(u, v - hv), p))) / hv ** 2)
    su = math.sqrt(8 * tol / ku) if ku > 1e-12 else None
    sv = math.sqrt(8 * tol / kv) if kv > 1e-12 else None
    return su, sv


# ---------------------------------------------------------------- lecture des nœuds

def _doubles(raw):
    return list(struct.unpack(f'>{len(raw) // 8}d', raw))


def _shorts(raw):
    return list(struct.unpack(f'>{len(raw) // 2}h', raw))


class Model:
    """Faces of all the bodies of a partition, ready to be meshed."""

    def __init__(self, nodes):
        self.N = nodes
        self._surf, self._curve = {}, {}

    # -- courbes
    def curve(self, i):
        if i in self._curve:
            return self._curve[i]
        try:
            c = self._make_curve(i)
        except (KeyError, IndexError, struct.error, ZeroDivisionError, ValueError):
            c = None   # lien cassé ou donnée inattendue : courbe ignorée (l'arête devient droite)
        self._curve[i] = c
        return c

    def _make_curve(self, i):
        n = self.N.get(i)
        c = None
        if n is None:
            c = None
        elif n.type == 30:
            c = Line(n.f[7], n.f[8])
        elif n.type == 31:
            c = Circle(n.f[7], n.f[8], n.f[9], n.f[10])
        elif n.type == 32:
            c = Ellipse(n.f[7], n.f[8], n.f[9], n.f[10], n.f[11])
        elif n.type == 133:
            b = self.curve(n.f[7])
            c = Trimmed(b, n.f[8], n.f[9], n.f[10], n.f[11]) if b is not None and n.f[10] > UNSET and n.f[11] > UNSET else b
        elif n.type == 38:
            ch = self.N[n.f[9]]
            pts = [list(struct.unpack_from('>3d', ch.raw, 24 * k)) for k in range(len(ch.raw) // 24)]   # raw = points du tracé
            c = Polyline(pts)
        elif n.type == 134:
            c = self._nurbs_curve(n, Curve3D)
        elif n.type == 137:
            s = self.surface(n.f[7], raw=True)
            b = self._nurbs_curve(self.N[n.f[8]], NurbsCurve)
            c = CurveOnSurface(s, b) if s is not None and b is not None else None
        return c

    def _nurbs_curve(self, bc, cls):
        if bc.type != 134:
            return None
        nc, data = self.N[bc.f[7]], bc.f[8]
        r = nc.raw
        deg, nv, dim, nk = struct.unpack_from('>hihi', r, 0)
        rational = r[15]   # degré 2, n points 4, dim 2, n nœuds 4, type 1, périodique 1, fermée 1, rationnelle 1, forme 1
        pv, pm, pk = nc.f[17], nc.f[18], nc.f[19]
        vals = _doubles(self.N[pv].raw)
        cps = [vals[k * dim:(k + 1) * dim] for k in range(nv)]
        U = Nurbs.knots(_doubles(self.N[pk].raw), _shorts(self.N[pm].raw))
        return cls(deg, U, cps, rational, dim)

    # -- surfaces
    def surface(self, i, raw=False):
        key = (i, raw)
        if key in self._surf:
            return self._surf[key]
        try:
            s = self._make_surface(i)
        except (KeyError, IndexError, struct.error, ZeroDivisionError, ValueError):
            s = None   # face maillée à plat en secours
        self._surf[key] = s
        return s

    def _make_surface(self, i):
        n = self.N.get(i)
        s = None
        t = n.type if n else None
        if t == 50:
            s = Plane(n.f[7], n.f[8], n.f[9])
        elif t == 51:
            s = Cylinder(n.f[7], n.f[8], n.f[9], n.f[10])
        elif t == 52:
            s = Cone(n.f[7], n.f[8], n.f[9], n.f[10], n.f[11], n.f[12])
        elif t == 53:
            s = Sphere(n.f[7], n.f[8], n.f[9], n.f[10])
        elif t == 54:
            s = Torus(n.f[7], n.f[8], n.f[9], n.f[10], n.f[11])
        elif t == 124:
            s = self._nurbs_surface(n)
        elif t == 60:   # décalée : surface de base + d × normale
            base = self.surface(n.f[9])
            d = n.f[10] * (-1 if n.f[6] == b'-' else 1)   # sens « - » : décalage du côté opposé à la normale
            if isinstance(base, Plane):
                s = Plane(add(base.o, mul(base.n, d)), base.n, base.x)
            elif isinstance(base, Cylinder):
                s = Cylinder(base.o, base.a, base.r + d, base.x)
            elif isinstance(base, Sphere):
                s = Sphere(base.c, base.r + d, base.a, base.x)
            elif isinstance(base, Torus):
                s = Torus(base.c, base.a, base.R, base.r + d, base.x)
            elif base is not None:
                fb = base.raw if hasattr(base, 'raw') else base.eval
                box = getattr(base, 'box', None)
                if box:
                    eu, ev = (box[1] - box[0]) * 0.25, (box[3] - box[2]) * 0.25
                    box = (box[0] - eu, box[1] + eu, box[2] - ev, box[3] + ev)
                def raw_(a, b, fb=fb, d=d):
                    p = fb(a, b); h = 1e-6
                    nn = unit(cross(sub(fb(a + h, b), p), sub(fb(a, b + h), p)))
                    return add(p, mul(nn, d))
                if box:
                    s = Generic(raw_, box)
                    if hasattr(base, 'su'):   # départ : projection sur la surface de base (plus fiable)
                        def proj(p, s=s, base=base, raw_=raw_):
                            a, b = base.project(p)
                            a, b = _newton(raw_, p, a / base.su, b / base.sv, s.box)
                            return a * s.su, b * s.sv
                        s.project = proj
        elif t == 67:   # balayée
            sec = self.curve(n.f[7])
            if sec is not None:
                s = Swept(sec, n.f[8])
        elif t == 68:   # révolution
            prof = self.curve(n.f[7])
            if prof is not None:
                s = Spun(prof, n.f[8], n.f[9])
        elif t == 56:   # congé à bille roulante
            spine = self.curve(n.f[10])
            s1, s2 = self.surface(n.f[8]), self.surface(n.f[9])
            if isinstance(spine, Polyline) and s1 is not None and s2 is not None and n.f[11] > UNSET:
                s = Blend(spine, s1, s2, n.f[11])
        if s is not None and not hasattr(s, 'raw'):
            s.raw = s.eval
        return s

    def _nurbs_surface(self, n):
        r = self.N[n.f[7]].raw
        du, dv, nu, nv = struct.unpack_from('>hhii', r, 2)
        rational = r[24]   # … n nœuds u/v (16-23), rationnelle 24, fermée u/v 25-26, forme 27, dim 28
        dim = struct.unpack_from('>h', r, 28)[0]
        pv, pmu, pmv, pku, pkv = self.N[n.f[7]].f[30:35]
        vals = _doubles(self.N[pv].raw)
        grid = [[vals[(a * nv + b) * dim:(a * nv + b + 1) * dim] for b in range(nv)] for a in range(nu)]
        U = Nurbs.knots(_doubles(self.N[pku].raw), _shorts(self.N[pmu].raw))
        V = Nurbs.knots(_doubles(self.N[pkv].raw), _shorts(self.N[pmv].raw))
        s = NurbsSurface(du, dv, U, V, grid, rational)
        s.box = (s.u0, s.u1, s.v0, s.v1)
        return s

    # -- topologie
    def faces(self):
        return [n for n in self.N.values() if n.type == 14]

    def loops(self, face):
        """[[fin, …]] of the face, fins in order."""
        out, lp, guard = [], face.f[5], 0
        while lp in self.N and guard < 100000:
            fins, f0, fn = [], self.N[lp].f[2], self.N[lp].f[2]
            while fn in self.N:
                fins.append(fn)
                fn = self.N[fn].f[2]
                if fn == f0 or len(fins) > 100000:
                    break
            out.append(fins)
            lp = self.N[lp].f[4]
            guard += 1
        return out

    def body_of(self, face):
        sh = self.N.get(face.f[6])
        return sh.f[2] if sh else None
