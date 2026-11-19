# SPDX-License-Identifier: GPL-3.0-or-later
"""Reader of the Parasolid partition stored in a .sldprt (binary neutral transmit, v28): topology and geometry.

Decoded on 12/11/2026 from a real part: big-endian; each node = type (2 bytes), at the first node of a type a schema mark
(0xff = as the base schema, else n + change codes C/D/I/A + 'Z'), a length (4 bytes) for variable nodes, its index (2 bytes),
then its fields. Pointers are indexes (2 bytes), 1 = none. Doubles -3.14158e13 = unset."""
import re
import struct
import zlib

from .reader import FormatError, sections

# type -> (bytes after the index, bytes per element for variable nodes or 0)
SIZES = {
    68: (161, 0),   # SPUN_SURF
    101: (37, 0),   # PARTITION (racine)
    12: (75, 0),
    13: (20, 0),
    14: (35, 0),
    15: (12, 0),
    16: (28, 0),
    17: (19, 0),
    18: (24, 0),
    19: (17, 0),
    29: (36, 0),
    30: (63, 0),
    31: (95, 0),
    32: (103, 0),
    38: (25, 0),
    40: (52, 24),
    41: (2, 24),
    45: (0, 8),
    50: (87, 0),
    51: (95, 0),
    52: (111, 0),
    53: (95, 0),
    54: (103, 0),
    56: (62, 0),
    59: (19, 0),
    60: (35, 0),
    67: (49, 0),
    70: (28, 0),
    74: (10, 2),
    79: (0, 1),
    80: (32, 1),
    81: (16, 2),
    82: (0, 4),
    83: (0, 8),
    84: (0, 1),
    124: (19, 0),
    125: (85, 0),
    126: (40, 0),
    127: (0, 2),
    128: (0, 8),
    133: (81, 0),
    134: (19, 0),
    135: (3, 0),
    136: (23, 0),
    137: (29, 0),
    141: (8, 0)
}
# fields: d int (4), p pointer (2, or 4 above 32767), f double, c byte, n short, v vector (3 doubles)
LAYOUT = {
    12: "dppppppffpppcpcc" + "pppppppp" + "dppppd",   # BODY (75 octets avec des pointeurs courts)
    13: "dpppppppp",          # SHELL     id, attr, body, next, face, edge, vertex, region, front face
    14: "dpfpppppcppppp",     # FACE      id, attr, tol, next, prev, loop, shell, surface, sense, …, front shell
    15: "dpppp",              # LOOP      id, attr, fin, face, next
    16: "dpfppppppp",         # EDGE      id, attr, tol, fin, prev, next, curve, …, owner
    17: "pppppppppc",         # FIN       attr, loop, forward, backward, vertex (end), other, edge, curve, next at vertex, sense
    18: "dpppppfp",           # VERTEX    id, attr, fin, prev, next, point, tol, owner
    19: "dpppppcp",           # REGION
    29: "dppppv",             # POINT     id, attr, owner, next, prev, position
    30: "dpppppcvv",          # LINE      …, sense, point, direction
    31: "dpppppcvvvf",        # CIRCLE    …, centre, normal, x axis, radius
    32: "dpppppcvvvff",       # ELLIPSE   …, centre, normal, x axis, major, minor
    38: "dpppppcppppp",       # INTERSECTION …, surfaces, chart, limits
    40: "ffdffff",            # CHART (+ points)
    41: "cc",                 # LIMIT (+ points)
    45: "",                   # BSPLINE_VERTICES (doubles)
    50: "dpppppcvvv",         # PLANE     …, point, normal, x axis
    51: "dpppppcvvfv",        # CYLINDER  …, point, axis, radius, x axis
    52: "dpppppcvvfffv",      # CONE      …, point, axis, radius, sin, cos, x axis
    53: "dpppppcvfvv",        # SPHERE    …, centre, radius, axis, x axis
    54: "dpppppcvvffv",       # TORUS     …, centre, axis, major, minor, x axis
    56: "dpppppccpppffffpppp",  # BLENDED_EDGE (rolling ball)
    59: "dpppppcnp",          # BLEND_BOUND
    60: "dpppppcccpff",       # OFFSET_SURF
    67: "dpppppcpvf",         # SWEPT_SURF
    68: "dpppppcpvvvvffvf",   # SPUN_SURF …, profile, base, axis, start, end, start parm, end parm, x axis, tolerance
    70: "dccpppdddpp",        # LIST
    74: "ddp",                # POINTER_LIST_BLOCK (+ pointers)
    79: "",                   # string
    80: "ppd" + "c" * 24,     # ATTRIB_DEF (+ field types)
    81: "dpppppp",            # ATTRIBUTE (+ field pointers)
    82: "", 83: "", 84: "",   # int / real / char values
    85: "", 86: "", 87: "", 88: "", 89: "", 98: "",   # point / vector / direction / axis / tag / unicode values (attributes)
    100: "dppp" + "f" * 9 + "vfdfff",   # TRANSFORM  id, owner, next, prev, rotation (9), translation, scale, flag, 3 unset
    101: "pppppppcppddpdpp",  # PARTITION (racine)
    124: "dpppppcpp",         # B_SURFACE …, nurbs, data
    125: "f" * 8 + "c" * 13 + "pppp",   # SURFACE_DATA
    126: "c" * 30 + "ppppp",  # NURBS_SURF (en-tête de 30 octets, puis vertices, mult u, mult v, nœuds u, nœuds v)
    127: "", 128: "",         # KNOT_MULT (shorts) / KNOT_SET (doubles)
    133: "dpppppcpvvff",      # TRIMMED_CURVE …, basis, point 1, point 2, parm 1, parm 2
    134: "dpppppcpp",         # B_CURVE
    135: "cp",                # CURVE_DATA
    136: "c" * 17 + "ppp",    # NURBS_CURVE (en-tête de 17 octets, puis vertices, mult, nœuds)
    137: "dpppppcpppf",       # SP_CURVE
    141: "pppp",              # GEOMETRIC_OWNER
}
# nœuds de longueur variable : type de leurs éléments
VAR = {40: "v", 41: "v", 45: "f", 74: "p", 79: "c", 80: "c", 81: "p", 82: "d", 83: "f", 84: "c", 85: "v", 86: "v", 87: "v", 88: "a",
       89: "d", 98: "n", 127: "n", 128: "f"}
_ELT = {"v": 24, "f": 8, "c": 1, "d": 4, "n": 2, "a": 48}   # bytes of one element of a variable node
SIZES = LAYOUT   # compatibilité
_SZ = {"d": 4, "f": 8, "c": 1, "n": 2, "v": 24}
_FMT = {"d": ">i", "f": ">d", "c": ">c", "n": ">h", "v": ">3d"}
NAMES = {12: "body", 13: "shell", 14: "face", 15: "loop", 16: "edge", 17: "fin", 18: "vertex", 29: "point", 30: "line", 31: "circle",
         32: "ellipse", 50: "plane", 51: "cylinder", 52: "cone", 53: "sphere", 54: "torus", 56: "blend", 60: "offset", 67: "swept", 68: "spun",
         124: "bspline_surface", 134: "bspline_curve", 38: "intersection", 133: "trimmed_curve", 137: "sp_curve"}


def partition_stream(path, config=0):
    """Decompressed Parasolid stream of a configuration of the .sldprt (the deltas that follow are left out)."""
    with open(path, "rb") as f:
        secs = sections(f.read())
    raw = secs.get(f"Contents/Config-{config}-Partition")
    if not raw:
        raise FormatError("No Parasolid partition in this file")
    # en-tête : taille (4) + empreinte (16) ; puis des blocs [taille décompressée (4), taille compressée (4), zlib] jusqu'à un bloc vide
    q, parts = 20, []
    while q + 8 <= len(raw):
        usize, csize = struct.unpack_from("<II", raw, q); q += 8
        if usize == 0 or csize == 0:
            break
        try:
            parts.append(zlib.decompress(raw[q:q + csize]))
        except zlib.error as e:
            raise FormatError(f"Parasolid partition: damaged block ({e})")
        q += csize
    data = b"".join(parts)
    if not data.startswith(b"PS"):
        raise FormatError("Not a Parasolid binary stream")
    return data


def body_streams(path, config=0):
    """[(name, stream)] of the bodies a part keeps per feature (newer files, imported bodies): each one a compressed Parasolid
    transmit file, after the UTF-16 name of its feature. Empty when the file has none."""
    with open(path, "rb") as f:
        secs = sections(f.read())
    raw = secs.get(f"Config-{config}-FeatureBodies/LocalBodies") or b""
    out, i = [], 0
    while True:
        z = _ZLIB.search(raw, i)
        if not z:
            return out
        p = z.start()
        try:
            o = zlib.decompressobj()
            data = o.decompress(raw[p:])
        except zlib.error:
            i = p + 1
            continue
        if not (o.eof and data.startswith(b"PS")):
            i = p + 1
            continue
        names = _NAME.findall(raw[max(0, p - 120):p])
        i = len(raw) - len(o.unused_data)
        # gros corps : blocs [taille décompressée (4), taille compressée (4), zlib] à la suite, comme la partition
        if p >= 8 and struct.unpack_from("<I", raw, p - 8)[0] == len(data):
            parts = [data]
            while i + 8 <= len(raw):
                usize, csize = struct.unpack_from("<II", raw, i)
                if not usize or not csize or i + 8 + csize > len(raw):
                    break
                try:
                    more = zlib.decompress(raw[i + 8:i + 8 + csize])
                except zlib.error:
                    break
                if len(more) != usize or more.startswith(b"PS"):
                    break
                parts.append(more)
                i += 8 + csize
            data = b"".join(parts)
        out.append((names[-1].decode("utf-16-le") if names else f"Body {len(out) + 1}", data))


_ZLIB = re.compile(rb"\x78[\x01\x5e\x9c\xda]")
_NAME = re.compile(rb"(?:[\x20-\x7e]\x00){2,}")


# Schema marks of the node types that differ from the base schema (13006) in the streams the table LAYOUT was written for
# (modeller 28): C = field of the base schema kept, D = dropped, a name = field added. Any other stream gives its own marks:
# its fields are read with them and put back in the LAYOUT order, so the rest of the reader does not depend on the version.
V28 = {
    12: "C C C C C C C C C C C C C D owner C C C C C C C C C boundary_mesh index_map_offset index_map node_id_index_map "
        "schema_embedding_map child lowest_node_id",
    19: "C C C C C C C owner",
    41: "C term_use C",
    70: "C list_type notransmit C C C D C C D finger_index finger_block C",
    74: "C index_map_offset C C",
    101: "C C C C C C C C C attdef_list C C mesh index_map_offset index_map schema_embedding_map",
}
_KINDS = {"p": "p", "d": "d", "f": "f", "c": "c", "u": "c", "l": "c", "n": "n", "w": "n", "v": "v", "h": "v", "i": "ff", "b": "ffffff"}
_DEFAULT = {"p": 1, "d": 0, "f": 0.0, "c": b"\x00", "n": 0, "v": (0.0, 0.0, 0.0)}


def _mark(x, q):
    """Schema mark at the first node of a type: (tokens, next offset); tokens = None (as the base schema) or a list of
    'C', 'D' and (name, kinds, count) for an added field."""
    if x[q] == 0xFF:
        return None, q + 1
    q += 1
    out = []
    while True:
        c = chr(x[q]); q += 1
        if c == "Z":
            return out, q
        if c in "IA":
            n = x[q]; name = x[q + 1:q + 1 + n].decode("latin1"); q += 1 + n
            cls, count = struct.unpack_from(">HH", x, q); q += 4
            kind = "p"
            if cls == 0:
                n = x[q]; kind = x[q + 1:q + 1 + n].decode("latin1"); q += 1 + n
            out.append((name, kind, count))
        elif c in "CD":
            out.append(c)
        else:
            raise FormatError("Parasolid schema mark not understood")


def _declare(x, q):
    """Full description of a node type at its first node: ((kinds of its fixed fields, kind of its variable part or ""), next
    offset). Layout: number of fields, name, title, then per field its name, class, count and (class 0) its kind; one byte
    tells whether the last field has a variable length."""
    nf = x[q]; q += 1
    for _ in range(2):   # name and title
        q += 1 + x[q]
    kinds = []
    for _ in range(nf):
        q += 1 + x[q]
        cls, count = struct.unpack_from(">HH", x, q); q += 4
        kind = "p"
        if cls == 0:
            n = x[q]; kind = x[q + 1:q + 1 + n].decode("latin1"); q += 1 + n
        k = "".join(_KINDS.get(c, "?") for c in kind)
        if "?" in k or len(k) != 1:
            raise FormatError(f"Parasolid field type « {kind} » not known")
        kinds.append((k, count))
    var = x[q]; q += 1
    if var and kinds:
        last, _count = kinds.pop()
        return ("".join(k * max(c, 1) for k, c in kinds), last), q
    return ("".join(k * max(c, 1) for k, c in kinds), ""), q


def _plan(t, tokens):
    """How to read the fields of node type t in a stream whose mark is tokens: [(kinds to read, place in LAYOUT or None)];
    None when the stream uses the LAYOUT order as it is."""
    lay = LAYOUT[t]
    slots = len(lay) + (1 if t in VAR else 0)   # the variable part of a node counts as its last field
    ref = V28[t].split() if t in V28 else []
    place, b, j = {}, 0, 0
    for tok in ref:
        if tok == "C":
            place[("b", b)] = j; b += 1; j += 1
        elif tok == "D":
            b += 1
        else:
            place[("n", tok)] = j; j += 1
    while j < slots:   # a mark stops after its last change: the base fields that follow are kept
        place[("b", b)] = j; b += 1; j += 1
    nbase = b
    if tokens is None:
        tokens = []
    names = [tk if isinstance(tk, str) else tk[0] for tk in tokens]
    while names and names[-1] == "C":
        names.pop()
    short = list(ref)
    while short and short[-1] == "C":
        short.pop()
    if names == short:
        return None
    plan, b = [], 0
    for tok in list(tokens) + ["C"] * nbase:   # (the trailing base fields; the loop stops when they are all read)
        if tok == "C":
            if b >= nbase:
                break
            pos = place.get(("b", b)); b += 1
            if pos is None:
                raise FormatError(f"Parasolid node type {t}: field of the base schema not known")
            plan.append((lay[pos] if pos < len(lay) else "", pos))
        elif tok == "D":
            b += 1
        else:
            name, kind, count = tok
            kinds = "".join(_KINDS.get(k, "?") for k in kind) * max(count, 1)
            if "?" in kinds:
                raise FormatError(f"Parasolid node type {t}: field type « {kind} » not known")
            pos = place.get(("n", name))
            plan.append((kinds, pos if pos is not None and len(kinds) == 1 else None))
    return plan


class Node:
    __slots__ = ("type", "index", "f", "raw")

    def __init__(self, t, idx, f, raw):
        self.type, self.index, self.f, self.raw = t, idx, f, raw

    @property
    def name(self):
        return NAMES.get(self.type, str(self.type))


def _ptr(x, q):
    """Pointer (index): 2 bytes up to 32767; above, a negative short h then an unsigned short l: -h × 32767 + l."""
    h = struct.unpack_from(">h", x, q)[0]
    if h >= 0:
        return h, q + 2
    lo = struct.unpack_from(">H", x, q + 2)[0]
    return -h * 32767 + lo, q + 4


def parse(x):
    """{index: Node} of the whole stream."""
    m = x.find(b"SCH_")
    if m < 0:
        raise FormatError("Parasolid schema not found")
    p = m + x[m - 1] + 6   # nom du schéma, puis 6 octets avant le 1er nœud (la partition)
    plans, nodes, declared = {}, {}, {}
    while p < len(x) - 4:
        t = struct.unpack_from(">H", x, p)[0]; q = p + 2
        if t == 1:
            break   # fin du flux
        if t in declared or (t not in LAYOUT and x[q] != 0xFF):
            if t not in declared:   # type the base schema does not have: the stream describes it in full
                declared[t], q = _declare(x, q)
            kinds, var = declared[t]
            n = 0
            if var:
                n = struct.unpack_from(">I", x, q)[0]; q += 4
            idx, q = _ptr(x, q)
            for ch in kinds:
                if ch == "p":
                    _v, q = _ptr(x, q)
                else:
                    q += _SZ[ch]
            for _k in range(n if var == "p" else 0):
                _v, q = _ptr(x, q)
            if var and var != "p":
                q += n * _SZ[var]
            if q > len(x):
                raise FormatError("Parasolid stream cut short")
            p = q   # (nothing of these nodes is used yet)
            continue
        lay = LAYOUT.get(t)
        if lay is None:
            raise FormatError(f"Parasolid node type {t} not supported yet (at byte {p})")
        if t not in plans:
            tokens, q = _mark(x, q)
            plans[t] = _plan(t, tokens)
        plan = plans[t]
        n = None
        el = VAR.get(t)
        if el:
            n = struct.unpack_from(">I", x, q)[0]; q += 4
        idx, q = _ptr(x, q)
        start = q
        if plan is None:
            f = []
            for ch in lay:
                if ch == "p":
                    v, q = _ptr(x, q)
                    f.append(v)
                else:
                    v = struct.unpack_from(_FMT[ch], x, q); q += _SZ[ch]
                    f.append(v if ch == "v" else v[0])
        else:   # another version of the schema: fields read in its order, stored in the LAYOUT order
            f = [_DEFAULT[ch] for ch in lay]
            for kinds, pos in plan:
                for ch in kinds:
                    if ch == "p":
                        v, q = _ptr(x, q)
                    else:
                        v = struct.unpack_from(_FMT[ch], x, q); q += _SZ[ch]
                        v = v if ch == "v" else v[0]
                    if pos is not None and pos < len(lay):
                        f[pos] = v
        rawq = q
        if el == "p":
            items = []
            for _k in range(n):
                v, q = _ptr(x, q)
                items.append(v)
            f.append(items)
        elif el:
            q += n * _ELT[el]
        nodes[idx] = Node(t, idx, f, x[rawq:q] if el and el != "p" else x[start:q])
        p = q
    # contrôle : les faces, boucles et fins doivent pointer vers des nœuds lus (sinon le flux a été mal lu)
    link = {14: 5, 15: 2, 17: 6}   # face -> boucle, boucle -> fin, fin -> arête
    bad = sum(1 for n in nodes.values() if n.type in link and n.f[link[n.type]] > 1 and n.f[link[n.type]] not in nodes)
    if bad:
        raise FormatError(f"Parasolid stream read incompletely ({bad} broken links)")
    return nodes
