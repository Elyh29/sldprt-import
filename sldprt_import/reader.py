# SPDX-License-Identifier: GPL-3.0-or-later
"""Reads the display mesh stored in a .sldprt file (2015 and later): no CAD kernel, standard library only."""
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
