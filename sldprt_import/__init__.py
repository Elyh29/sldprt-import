# SPDX-License-Identifier: GPL-3.0-or-later
"""SLDPRT Import: File > Import > SLDPRT (.sldprt), or drop the file in the 3D View."""
import math
import os
import time

import bmesh
import bpy
from bpy_extras.io_utils import ImportHelper

from . import reader


LEVELS = (("COARSE", "Coarse", "Half the triangles of the file (lighter, rounds get flatter)"),
          ("FILE", "As Saved", "The mesh stored in the file"),
          ("FINE", "Fine", "4 times more triangles: rounds and holes bent onto the surface with the stored normals"),
          ("VERY_FINE", "Very Fine", "16 times more triangles (heavy on large parts)"))
REFINE = {"COARSE": 0, "FILE": 0, "FINE": 1, "VERY_FINE": 2}


def _decimate(ob, ratio):
    mod = ob.modifiers.new("Coarse", 'DECIMATE')
    mod.ratio = ratio
    mod.use_collapse_triangulate = True
    deps = bpy.context.evaluated_depsgraph_get()
    me = bpy.data.meshes.new_from_object(ob.evaluated_get(deps), depsgraph=deps)
    ob.modifiers.remove(mod)
    old, ob.data = ob.data, me
    me.name = old.name
    bpy.data.meshes.remove(old)


def load(context, path, merge=True, smooth_angle=math.radians(30), face_ids=False, detail="FILE"):
    verts, tris, owner = reader.read(path, REFINE[detail])
    name = os.path.splitext(os.path.basename(path))[0]
    me = bpy.data.meshes.new(name)
    me.from_pydata(verts, [], tris)
    if face_ids:
        attr = me.attributes.new("cad_face", 'INT', 'FACE')
        attr.data.foreach_set("value", owner)
    me.validate(clean_customdata=False)
    if merge:
        bm = bmesh.new()
        bm.from_mesh(me)
        bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-6)
        bm.to_mesh(me)
        bm.free()
    me.shade_smooth()
    me.set_sharp_from_angle(angle=smooth_angle)
    ob = bpy.data.objects.new(name, me)
    context.collection.objects.link(ob)
    ob.location = context.scene.cursor.location
    if detail == "COARSE":
        _decimate(ob, 0.5)
        ob.data.shade_smooth()
        ob.data.set_sharp_from_angle(angle=smooth_angle)
    for o in context.selected_objects:
        o.select_set(False)
    ob.select_set(True)
    context.view_layer.objects.active = ob
    return ob


def _two_sided(mesher, face):
    """A face on a surface periodic both ways (torus, ring): its band between two circles can be taken either way round."""
    surf = mesher.M.surface(face.f[7])
    if surf is None:
        return False
    n = len(mesher.M.loops(face))
    if type(surf).__name__ == "Sphere":   # a cap inside one loop: the small side or the large one
        return n == 1
    return bool(getattr(surf, "per_u", None)) and bool(getattr(surf, "per_v", None)) and n == 2


def _gap(bvh, fm):
    """Largest distance from the points of a face mesh to the display mesh (same axes as the display mesh: Z up). The wrong
    side of a ring or of a sphere always has points far from the part; its other points may lie on neighbouring faces."""
    if not fm.verts or not fm.tris:
        return float("inf")
    step = max(1, len(fm.verts) // 200)
    worst = 0.0
    for v in fm.verts[::step]:
        hit = bvh.find_nearest((v[0], -v[2], v[1]))
        worst = max(worst, hit[3] if hit[0] is not None else float("inf"))
    return worst


def _mesh(mesher, f, bvh):
    """Mesh of one face (on a ring or a sphere, the side closest to the display mesh)."""
    fm = mesher.mesh_face(f)
    if bvh is not None and _two_sided(mesher, f):
        alt = mesher.mesh_face(f, other_side=True)
        if mesher.stale(fm):   # the second try cut an edge of the face: the first one again, with the same points
            fm = mesher.mesh_face(f)
        if _gap(bvh, alt) < _gap(bvh, fm):
            fm = alt
    return fm


def _boxes(model, k):
    """{(k, body): (face ids, box of its vertices or None)} of one source, box in metres, Z up."""
    N, out = model.N, {}
    for f in model.faces():
        ids, pts = out.setdefault((k, model.body_of(f)), (set(), []))
        ids.add(f.f[0])
        for lp in model.loops(f):
            for fn in lp:
                vx = N.get(N[fn].f[4])
                pt = vx and N.get(vx.f[5])
                if pt is not None and pt.type == 29:
                    x, y, z = pt.f[5]
                    pts.append((x, -z, y))
    return {key: (ids, ([min(p[c] for p in pts) for c in range(3)], [max(p[c] for p in pts) for c in range(3)]) if pts else None)
            for key, (ids, pts) in out.items()}


HIDDEN = (("COARSE", "Coarse", "Hidden bodies imported hidden, meshed at five times the tolerance (faster)"),
          ("FULL", "Full", "Hidden bodies imported hidden, meshed at the tolerance"),
          ("SKIP", "Skip", "Hidden bodies left out"))


def _body_object(context, root, label, body, lst, hidden, face_ids):
    """Object of one body from its face meshes (None if it has no triangle)."""
    index, verts, tris, nrm, owner = {}, [], [], [], []
    for fm in lst:
        ids = []
        for v in fm.verts:
            key = (round(v[0] * 1e8), round(v[1] * 1e8), round(v[2] * 1e8))
            if key not in index:
                index[key] = len(verts)
                verts.append((v[0], -v[2], v[1]))   # Y vers le haut -> Z vers le haut
            ids.append(index[key])
        for a, b, c in fm.tris:
            ia, ib, ic = ids[a], ids[b], ids[c]
            if ia == ib or ib == ic or ia == ic:
                continue
            tris.append((ia, ib, ic))
            nrm += [(n[0], -n[2], n[1]) for n in (fm.normals[a], fm.normals[b], fm.normals[c])]
            owner.append(fm.face)
    if not tris:
        return None
    if callable(hidden):
        hidden = hidden(verts)
    me = bpy.data.meshes.new(label)
    me.from_pydata(verts, [], tris)
    me.shade_smooth()
    me.normals_split_custom_set(nrm)   # normales de la vraie surface : ombrage exact, arêtes nettes entre les faces
    if face_ids:
        attr = me.attributes.new("cad_face", 'INT', 'FACE')
        attr.data.foreach_set("value", owner)
    ob = bpy.data.objects.new(label + (" (hidden)" if hidden else ""), me)
    context.collection.objects.link(ob)
    ob.parent = root
    if hidden:
        ob.hide_render = True
        ob.hide_viewport = True
    return ob


def exact_steps(context, path, tol=0.0002, face_ids=False, hidden_bodies="COARSE"):
    """Import of the exact geometry in small steps: yields the share done (0 to 1), returns (root, stats).

    Each body becomes an object as soon as its faces are meshed, so that a modal import shows the part growing."""
    from mathutils.bvhtree import BVHTree
    from . import brep, parasolid, shown, tessellate
    try:
        sources = [("", parasolid.parse(parasolid.partition_stream(path)))]
    except reader.FormatError:
        sources = []
    yield 0.0
    if not sources or not brep.Model(sources[0][1]).faces():   # bodies kept per feature (imported bodies, newer files)
        sources = []
        for nm, st in parasolid.body_streams(path):
            sources.append((nm, parasolid.parse(st)))
            yield 0.0
    # bodies shown by the original software: those of the display mesh, matched by the ids of their faces
    models = [(k, nm, brep.Model(N)) for k, (nm, N) in enumerate(sources)]
    info = {}
    for k, _nm, model in models:
        info.update(_boxes(model, k))
    try:
        visible = shown.shown(reader.bodies(path), info)
    except Exception:   # noqa: BLE001 — affichage illisible : pas de verdict
        visible = None
    yield 0.0
    jobs = []
    for k, nm, model in models:
        faces = model.faces()
        if visible is not None:
            hid = [f for f in faces if (k, model.body_of(f)) not in visible]
            faces = [f for f in faces if (k, model.body_of(f)) in visible]
            if hid and hidden_bodies != "SKIP":
                jobs.append((k, nm, tessellate.Mesher(model, tol=tol * (5 if hidden_bodies == "COARSE" else 1)), hid))
        jobs.append((k, nm, tessellate.Mesher(model, tol=tol), faces))
    total = sum(len(f) for _k, _n, _m, f in jobs)
    if not total:
        raise reader.FormatError("no surface found in the exact geometry")
    # maillage d'affichage : départage les anneaux à deux côtés possibles (et, sans verdict ci-dessus, repère les corps
    # cachés : leurs sommets sont loin du maillage d'affichage)
    try:
        rv, rt, _o = reader.read(path)
        bvh = BVHTree.FromPolygons(rv, rt, all_triangles=True)   # même repère que les points ci-dessous (Z vers le haut)
    except reader.FormatError:
        bvh = None
    yield 0.0

    def far(verts):
        # the display mesh is coarse (chords of a few mm on a large curved skin): a body is hidden only when most of its
        # points are far from it compared with its own size, not with the meshing tolerance
        ds = sorted(bvh.find_nearest(verts[k])[3] for k in range(0, len(verts), max(1, len(verts) // 200)))
        lo = [min(v[c] for v in verts) for c in range(3)]
        hi = [max(v[c] for v in verts) for c in range(3)]
        size = math.sqrt(sum((hi[c] - lo[c]) ** 2 for c in range(3)))
        return ds[len(ds) // 2] > max(0.002, 5 * tol, 0.005 * size)

    name = os.path.splitext(os.path.basename(path))[0]
    root = bpy.data.objects.new(name, None)
    root.empty_display_size = 0.05
    context.collection.objects.link(root)
    root.location = context.scene.cursor.location
    made, nfaces, nflat, done = [], 0, 0, 0
    try:
        for k, nm, mesher, faces in jobs:
            per_body = {}
            for f in faces:
                per_body.setdefault(mesher.M.body_of(f), []).append(f)
            for body, flist in per_body.items():   # bodies share no edge: each one is done, then made, on its own
                fms = []
                for f in flist:
                    fms.append(_mesh(mesher, f, bvh))
                    done += 1
                    yield done / total
                for _round in range(3):   # faces meshed before a neighbour cut one of their edges: again, with the new points
                    redo = [i for i, fm in enumerate(fms) if mesher.stale(fm)]
                    if not redo:
                        break
                    for i in redo:
                        fms[i] = _mesh(mesher, flist[i], bvh)
                        yield done / total
                nfaces += len(fms)
                nflat += sum(1 for fm in fms if fm.fallback)
                key = (k, body)
                hidden = (key not in visible) if visible is not None else (far if bvh is not None else False)
                ob = _body_object(context, root, nm or f"{name} {len(made) + 1}", key, fms, hidden, face_ids)
                if ob is not None:
                    made.append(ob)
    except GeneratorExit:   # import cancelled: what is made stays
        if not made:
            bpy.data.objects.remove(root)
        raise
    except Exception:   # the caller falls back to the display mesh: no half part left beside it
        for ob in made:
            me = ob.data
            bpy.data.objects.remove(ob)
            bpy.data.meshes.remove(me)
        bpy.data.objects.remove(root)
        raise
    if not made:   # nothing meshed: no empty left behind, the caller falls back to the display mesh
        bpy.data.objects.remove(root)
        raise reader.FormatError("no surface could be meshed")
    for o in context.selected_objects:
        o.select_set(False)
    for o in made:
        if not o.hide_viewport:
            o.select_set(True)
    root.select_set(True)
    context.view_layer.objects.active = root
    return root, {"faces": nfaces, "fallback": nflat}


def load_exact(context, path, tol=0.0002, face_ids=False, hidden_bodies="COARSE"):
    """Part meshed from its exact geometry (Parasolid): one object per body under an empty named after the file."""
    steps = exact_steps(context, path, tol, face_ids, hidden_bodies)
    wm = context.window_manager
    wm.progress_begin(0, 1000)
    try:
        while True:
            try:
                wm.progress_update(int(next(steps) * 1000))
            except StopIteration as end:
                return end.value
    finally:
        wm.progress_end()


SOURCES = (("EXACT", "Exact Geometry", "Mesh the exact surfaces of the part (Parasolid) at the tolerance below"),
           ("DISPLAY", "Display Mesh", "The mesh stored in the file for display (fast, fixed detail)"))


class IMPORT_OT_sldprt(bpy.types.Operator, ImportHelper):
    """Import .sldprt part files: exact geometry meshed at a tolerance, or the display mesh stored in the file"""
    bl_idname = "import_mesh.sldprt"
    bl_label = "Import SLDPRT"
    bl_options = {'REGISTER', 'UNDO', 'PRESET'}

    filename_ext = ".sldprt"
    filter_glob: bpy.props.StringProperty(default="*.sldprt;*.SLDPRT", options={'HIDDEN'})
    files: bpy.props.CollectionProperty(type=bpy.types.OperatorFileListElement, options={'HIDDEN', 'SKIP_SAVE'})
    directory: bpy.props.StringProperty(subtype='DIR_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    source: bpy.props.EnumProperty(name="Source", items=SOURCES, default="EXACT")
    tolerance: bpy.props.FloatProperty(name="Tolerance", description="Largest gap between the mesh and the exact surface (smaller: more triangles)",
                                       default=0.0002, min=0.00001, max=0.01, subtype='DISTANCE', unit='LENGTH', precision=3)
    detail: bpy.props.EnumProperty(name="Detail", items=LEVELS, default="FILE",
                                   description="Number of triangles: the file stores one mesh; Fine and Very Fine bend new points onto the curved surfaces")
    merge: bpy.props.BoolProperty(name="Merge Faces", description="Join the faces of the part into one closed mesh", default=True)
    smooth_angle: bpy.props.FloatProperty(name="Sharp Edges Above", description="Edges sharper than this stay sharp, the others are smoothed",
                                          subtype='ANGLE', default=math.radians(30), min=0, max=math.pi)
    face_ids: bpy.props.BoolProperty(name="Keep CAD Faces", description="Store the CAD face of each triangle in the attribute « cad_face »", default=False)
    hidden_bodies: bpy.props.EnumProperty(name="Hidden Bodies", items=HIDDEN, default="COARSE",
                                          description="Bodies hidden in the part (not in its display mesh)")

    def _paths(self):
        return [os.path.join(self.directory, f.name) for f in self.files if f.name] or [self.filepath]

    def _display(self, context, p):
        load(context, p, self.merge, self.smooth_angle, self.face_ids, self.detail)

    def _ended(self, p, st):
        if st["fallback"]:
            self.report({'WARNING'}, f"{os.path.basename(p)}: {st['fallback']} face(s) of {st['faces']} meshed flat (surface not handled)")

    def _failed(self, context, p, e):
        """Exact geometry not readable: the display mesh instead (True if something was imported)."""
        msg = str(e) if isinstance(e, reader.FormatError) else f"{type(e).__name__}: {e}"
        self.report({'WARNING'}, f"{os.path.basename(p)}: exact geometry not readable ({msg}), display mesh used")
        try:
            self._display(context, p)
            return True
        except (reader.FormatError, OSError) as e2:
            self.report({'WARNING'}, f"{os.path.basename(p)}: {e2}")
            return False

    def execute(self, context):
        if not bpy.app.background and context.window is not None:   # with a window: step by step, Blender stays usable
            self._queue, self._steps, self._cur, self._done, self._share = self._paths(), None, None, 0, 0.0
            self._count = len(self._queue)
            wm = context.window_manager
            self._timer = wm.event_timer_add(0.01, window=context.window)
            wm.modal_handler_add(self)
            wm.progress_begin(0, 1000)
            self._status(context)
            return {'RUNNING_MODAL'}
        done = 0   # scripts, background: in one go
        for p in self._paths():
            try:
                if self.source == "EXACT":
                    try:
                        _root, st = load_exact(context, p, self.tolerance, self.face_ids, self.hidden_bodies)
                        self._ended(p, st)
                        done += 1
                    except Exception as e:   # noqa: BLE001 — géométrie exacte illisible : maillage d'affichage
                        done += self._failed(context, p, e)
                else:
                    self._display(context, p)
                    done += 1
            except (reader.FormatError, OSError) as e:
                self.report({'WARNING'}, f"{os.path.basename(p)}: {e}")
        if not done:
            return {'CANCELLED'}
        self.report({'INFO'}, f"{done} part(s) imported")
        return {'FINISHED'}

    def _status(self, context):
        n = self._count - len(self._queue)
        name = os.path.basename(self._cur) if self._cur else ""
        part = f" ({n}/{self._count})" if self._count > 1 else ""
        context.workspace.status_text_set(f"Importing {name}{part}: {int(self._share * 100)} %   Esc to stop")
        context.window_manager.progress_update(int(((n - 1 + self._share) / max(1, self._count)) * 1000))

    def _step(self, context):
        """One small step of the import; False when every file is done."""
        if self._steps is None:
            if not self._queue:
                return False
            self._cur, self._share = self._queue.pop(0), 0.0
            if self.source != "EXACT":
                try:
                    self._display(context, self._cur)
                    self._done += 1
                except (reader.FormatError, OSError) as e:
                    self.report({'WARNING'}, f"{os.path.basename(self._cur)}: {e}")
                return True
            self._steps = exact_steps(context, self._cur, self.tolerance, self.face_ids, self.hidden_bodies)
        try:
            self._share = next(self._steps)
        except StopIteration as end:
            self._steps = None
            self._ended(self._cur, end.value[1])
            self._done += 1
        except Exception as e:   # noqa: BLE001 — géométrie exacte illisible : maillage d'affichage
            self._steps = None
            self._done += self._failed(context, self._cur, e)
        return True

    def _finish(self, context, stopped=False):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
        wm.progress_end()
        context.workspace.status_text_set(None)
        if stopped:
            self.report({'WARNING'}, "Import stopped: the bodies already meshed are kept")
            return {'FINISHED'}
        if not self._done:
            return {'CANCELLED'}
        self.report({'INFO'}, f"{self._done} part(s) imported")
        return {'FINISHED'}

    def modal(self, context, event):
        if event.type == 'ESC' and event.value == 'PRESS':
            if self._steps is not None:
                self._steps.close()   # the bodies made stay
                self._steps = None
            return self._finish(context, stopped=True)
        if event.type == 'Z' and (event.ctrl or event.oskey):
            return {'RUNNING_MODAL'}   # no undo while objects are being made
        if event.type != 'TIMER':
            return {'PASS_THROUGH'}
        t0 = time.perf_counter()
        try:
            while time.perf_counter() - t0 < 0.05:   # about 20 redraws a second
                if not self._step(context):
                    return self._finish(context)
        except ReferenceError:   # the objects being made were deleted meanwhile
            self._steps = None
            self.report({'WARNING'}, "Import stopped: its objects were deleted")
            return self._finish(context)
        self._status(context)
        return {'PASS_THROUGH'}

    def cancel(self, context):   # Blender closes the window or the file during the import
        if getattr(self, "_steps", None) is not None:
            self._steps.close()
        context.window_manager.event_timer_remove(self._timer)
        context.window_manager.progress_end()
        context.workspace.status_text_set(None)

    def draw(self, context):
        lay = self.layout
        lay.use_property_split = True
        lay.prop(self, "source")
        if self.source == "EXACT":
            lay.prop(self, "tolerance")
            lay.prop(self, "hidden_bodies")
        else:
            lay.prop(self, "detail")
            lay.prop(self, "merge")
            lay.prop(self, "smooth_angle")
        lay.prop(self, "face_ids")

    def invoke(self, context, event):
        if self.directory and self.files:
            return self.execute(context)   # glisser-déposer
        return ImportHelper.invoke(self, context, event)


class IO_FH_sldprt(bpy.types.FileHandler):
    bl_idname = "IO_FH_sldprt"
    bl_label = "SLDPRT"
    bl_import_operator = "import_mesh.sldprt"
    bl_file_extensions = ".sldprt;.SLDPRT"

    @classmethod
    def poll_drop(cls, context):
        return context.area is not None and context.area.type == 'VIEW_3D'


def menu(self, context):
    self.layout.operator(IMPORT_OT_sldprt.bl_idname, text="SLDPRT (.sldprt)")


CLASSES = (IMPORT_OT_sldprt, IO_FH_sldprt)


def register():
    for c in CLASSES:
        bpy.utils.register_class(c)
    bpy.types.TOPBAR_MT_file_import.append(menu)


def unregister():
    bpy.types.TOPBAR_MT_file_import.remove(menu)
    for c in reversed(CLASSES):
        bpy.utils.unregister_class(c)
