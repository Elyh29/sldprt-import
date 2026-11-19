# SPDX-License-Identifier: GPL-3.0-or-later
"""Reads the display mesh stored in a .sldprt file (2015 and later): no CAD kernel, standard library only."""
import collections
import struct
import zlib

MARK = b'\x14\x00\x06\x00\x08\x00'
STRIPS = struct.pack('<3i', 4, 8, 2)
VEC3 = struct.pack('<3i', 12, 100, 2)


class FormatError(Exception):
    pass


def _name(b):
    return bytes(((x << 4) | (x >> 4)) & 0xff for x in b).decode('latin1')


def sections(data):
    """{name: bytes} of the compressed sections of the file."""
    out, i = {}, data.find(MARK)
    while i >= 0:
        if i + 26 > len(data):   # marker in the last bytes of a cut file
            break
        _crc, _a, size, _full, nl = struct.unpack_from('<IIIII', data, i + 6)
        if 0 < nl < 300:
            p = i + 26 + nl
            try:
                out.setdefault(_name(data[i + 26:p]), zlib.decompressobj(-15).decompress(data[p:p + size + 16]))
            except zlib.error:
                pass
        i = data.find(MARK, i + 1)
    return out


def faces(dl):
    """[(vertices, normals, strip lengths)] of each face of the display list (metres, Y up)."""
    out, pos = [], 0
    while True:
        j = dl.find(STRIPS, pos)
        if j < 0:
            return out
        pos = j + 4
        if j < 8:
            continue
        ntri, nst = struct.unpack_from('<2i', dl, j - 8)
        k = struct.unpack_from('<i', dl, j + 12)[0]
        if k != nst or not 0 < k < 1000000 or j + 16 + 4 * k > len(dl):
            continue
        lens = struct.unpack_from(f'<{k}i', dl, j + 16)
        nv = sum(lens)
        if min(lens) < 3 or nv - 2 * k != ntri:
            continue
        q = j + 16 + 4 * k
        if dl[q:q + 12] != VEC3 or struct.unpack_from('<i', dl, q + 12)[0] != nv or q + 16 + 12 * nv > len(dl):
            continue
        v = struct.unpack_from(f'<{3 * nv}f', dl, q + 16)
        r = q + 16 + 12 * nv
        n = struct.unpack_from(f'<{3 * nv}f', dl, r + 16) if dl[r:r + 12] == VEC3 and r + 16 + 12 * nv <= len(dl) else None
        out.append((v, n, lens))
        pos = r


def _records(dl):
    """[(start, end, vertices)] of the face records of the display list, as faces() finds them."""
    out, pos = [], 0
    while True:
        j = dl.find(STRIPS, pos)
        if j < 0:
            return out
        pos = j + 4
        if j < 10:
            continue
        ntri, nst = struct.unpack_from('<2i', dl, j - 8)
        k = struct.unpack_from('<i', dl, j + 12)[0]
        if k != nst or not 0 < k < 1000000 or j + 16 + 4 * k > len(dl):
            continue
        lens = struct.unpack_from(f'<{k}i', dl, j + 16)
        nv = sum(lens)
        if min(lens) < 3 or nv - 2 * k != ntri:
            continue
        q = j + 16 + 4 * k
        if dl[q:q + 12] != VEC3 or struct.unpack_from('<i', dl, q + 12)[0] != nv or q + 16 + 12 * nv > len(dl):
            continue
        r = q + 16 + 12 * nv
        out.append((j - 10, r, struct.unpack_from(f'<{3 * nv}f', dl, q + 16)))
        pos = r


BODY_CLASS, FACE_CLASS = b'uoTempBodyTessData_c', b'uoTempFaceTessData_c'


def bodies(path):
    """Bodies shown by the display mesh: [(Parasolid ids of their faces, (low, high) corners in metres, Z up)].

    The display list holds one body record (its number of faces), then the records of its faces; the Parasolid id of a
    face is the last field of its record, just before the next face (none for the last face of a body).
    Returns [] when the layout is not the one recognised (no guess)."""
    with open(path, 'rb') as f:
        data = f.read()
    dl = sections(data).get('Contents/DisplayLists') if data[:8] != b'\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1' else None
    if not dl:
        return []
    recs = _records(dl)
    i0, i1 = dl.find(BODY_CLASS), dl.find(FACE_CLASS)
    if len(recs) < 2 or i0 < 0 or i1 != i0 + 28 or recs[0][0] != i1 + 18:
        return []
    ftag = collections.Counter(dl[r[0]:r[0] + 2] for r in recs[1:]).most_common(1)[0][0]
    btag = struct.pack('<H', struct.unpack('<H', ftag)[0] - 2)   # class of the body, its first object, then class of the face
    starts = [0] + [i for i in range(1, len(recs)) if dl[recs[i][0] - 4:recs[i][0] - 2] == btag]
    counts = [struct.unpack_from('<h', dl, i0 + 20)[0]] + [struct.unpack_from('<h', dl, recs[i][0] - 2)[0] for i in starts[1:]]
    ends = starts[1:] + [len(recs)]
    for k, (a, b, n) in enumerate(zip(starts, ends, counts)):
        if b - a != n and not (k == len(starts) - 1 and b - a < n):   # the last body may be cut with the file
            return []
    out = []
    for a, b in zip(starts, ends):
        ids, lo, hi = set(), [float('inf')] * 3, [float('-inf')] * 3
        for i in range(a, b):
            if i + 1 < b:   # the id of the last face of a body is not at a fixed place: left out
                ids.add(struct.unpack_from('<i', dl, recs[i + 1][0] - 4)[0])
            v = recs[i][2]
            xs, ys, zs = v[0::3], v[1::3], v[2::3]
            for c, (x0, x1) in enumerate(((min(xs), max(xs)), (-max(zs), -min(zs)), (min(ys), max(ys)))):   # Y up -> Z up
                lo[c], hi[c] = min(lo[c], x0), max(hi[c], x1)
        out.append((ids, (lo, hi)))
    return out


def triangles(v, n, lens):
    """Triangles (vertex indices) of the strips, facing the way the stored normals say; degenerate ones left out."""
    s = 0
    for L in lens:
        for t in range(L - 2):
            a, b, c = (s + t, s + t + 1, s + t + 2) if t % 2 == 0 else (s + t + 1, s + t, s + t + 2)
            pa, pb, pc = v[3 * a:3 * a + 3], v[3 * b:3 * b + 3], v[3 * c:3 * c + 3]
            if pa == pb or pb == pc or pa == pc:
                continue
            if n:
                u = [pb[i] - pa[i] for i in range(3)]
                w = [pc[i] - pa[i] for i in range(3)]
                cr = (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0])
                nm = [n[3 * a + i] + n[3 * b + i] + n[3 * c + i] for i in range(3)]
                if sum(cr[i] * nm[i] for i in range(3)) < 0:
                    b, c = c, b
            yield a, b, c
        s += L


def read_faces(path):
    """[(vertices, normals, triangles)] per CAD face: metres, Z up, local vertex indices."""
    with open(path, 'rb') as f:
        data = f.read()
    if data[:8] == b'\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1':
        raise FormatError("File saved before 2015 (older format): not supported yet")
    dl = sections(data).get('Contents/DisplayLists')
    if not dl:
        raise FormatError("No display mesh in this file (not a part saved with its graphics data?)")
    F = faces(dl)
    if not F:
        raise FormatError("Display mesh not recognised in this file")
    out = []
    for v, n, lens in F:
        vs = [(v[i], -v[i + 2], v[i + 1]) for i in range(0, len(v), 3)]   # Y up -> Z up
        ns = [(n[i], -n[i + 2], n[i + 1]) for i in range(0, len(n), 3)] if n else None
        out.append((vs, ns, list(triangles(v, n, lens))))
    return out


def read(path, levels=0):
    """(vertices, triangles, CAD face of each triangle), the mesh refined `levels` times (each one: 4 triangles per triangle)."""
    F = read_faces(path)
    if levels:
        from . import refine
        F = refine.refine(F, levels)
    verts, tris, owner = [], [], []
    for fi, (vs, _ns, ts) in enumerate(F):
        base = len(verts)
        verts += vs
        tris += [(base + a, base + b, base + c) for a, b, c in ts]
        owner += [fi] * len(ts)
    return verts, tris, owner
