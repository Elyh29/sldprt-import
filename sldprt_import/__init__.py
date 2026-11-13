# SPDX-License-Identifier: GPL-3.0-or-later
"""SLDPRT Import: File > Import > SLDPRT (.sldprt), or drop the file in the 3D View."""
import math
import os

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


def load_exact(context, path, tol=0.0002, face_ids=False):
    """Part meshed from its exact geometry (Parasolid): one object per body under an empty named after the file."""
    from mathutils.bvhtree import BVHTree
    from . import brep, parasolid, tessellate
    N = parasolid.parse(parasolid.partition_stream(path))
    model = brep.Model(N)
    mesher = tessellate.Mesher(model, tol=tol)
    faces = model.faces()
    wm = context.window_manager
    wm.progress_begin(0, len(faces))
    fms = []
    for k, f in enumerate(faces):
        fms.append(mesher.mesh_face(f))
        wm.progress_update(k)
    wm.progress_end()
    # corps affichés par le logiciel d'origine : leurs sommets sont sur le maillage d'affichage (les autres sont cachés dans le fichier)
    try:
        rv, rt, _o = reader.read(path)
        bvh = BVHTree.FromPolygons(rv, rt, all_triangles=True)   # même repère que les points ci-dessous (Z vers le haut)
    except reader.FormatError:
        bvh = None
    bodies = {}
    for fm in fms:
        bodies.setdefault(fm.body, []).append(fm)
    name = os.path.splitext(os.path.basename(path))[0]
    root = bpy.data.objects.new(name, None)
    root.empty_display_size = 0.05
    context.collection.objects.link(root)
    root.location = context.scene.cursor.location
    made = []
    for i, (body, lst) in enumerate(sorted(bodies.items(), key=lambda kv: -sum(len(f.tris) for f in kv[1]))):
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
            continue
        hidden = False
        if bvh is not None:
            ds = sorted(bvh.find_nearest(verts[k])[3] for k in range(0, len(verts), max(1, len(verts) // 50)))
            hidden = ds[len(ds) // 2] > max(0.001, 5 * tol)
        me = bpy.data.meshes.new(f"{name} {i + 1}")
        me.from_pydata(verts, [], tris)
        me.shade_smooth()
        me.normals_split_custom_set(nrm)   # normales de la vraie surface : ombrage exact, arêtes nettes entre les faces
        if face_ids:
            attr = me.attributes.new("cad_face", 'INT', 'FACE')
            attr.data.foreach_set("value", owner)
        ob = bpy.data.objects.new(f"{name} {i + 1}" + (" (hidden)" if hidden else ""), me)
        context.collection.objects.link(ob)
        ob.parent = root
        if hidden:
            ob.hide_render = True
            ob.hide_viewport = True
        made.append(ob)
    for o in context.selected_objects:
        o.select_set(False)
    for o in made:
        if not o.hide_viewport:
            o.select_set(True)
    root.select_set(True)
    context.view_layer.objects.active = root
    return root, mesher.stats


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

    def execute(self, context):
        paths = [os.path.join(self.directory, f.name) for f in self.files if f.name] or [self.filepath]
        done = 0
        for p in paths:
            try:
                if self.source == "EXACT":
                    try:
                        _root, st = load_exact(context, p, self.tolerance, self.face_ids)
                        if st["fallback"]:
                            self.report({'WARNING'}, f"{os.path.basename(p)}: {st['fallback']} face(s) of {st['faces']} meshed flat (surface not handled)")
                    except Exception as e:   # noqa: BLE001 — géométrie exacte illisible : maillage d'affichage
                        msg = str(e) if isinstance(e, reader.FormatError) else f"{type(e).__name__}: {e}"
                        self.report({'WARNING'}, f"{os.path.basename(p)}: exact geometry not readable ({msg}), display mesh used")
                        load(context, p, self.merge, self.smooth_angle, self.face_ids, self.detail)
                else:
                    load(context, p, self.merge, self.smooth_angle, self.face_ids, self.detail)
                done += 1
            except (reader.FormatError, OSError) as e:
                self.report({'WARNING'}, f"{os.path.basename(p)}: {e}")
        if not done:
            return {'CANCELLED'}
        self.report({'INFO'}, f"{done} part(s) imported")
        return {'FINISHED'}

    def draw(self, context):
        lay = self.layout
        lay.use_property_split = True
        lay.prop(self, "source")
        if self.source == "EXACT":
            lay.prop(self, "tolerance")
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
