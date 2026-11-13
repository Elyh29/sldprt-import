# SLDPRT Import

Open `.sldprt` part files in Blender, without the CAD software that made them.

**File ‣ Import ‣ SLDPRT (.sldprt)**, or drag the file into the 3D Viewport.

## What you get

- **Exact Geometry** (default): the part's real surfaces (planes, cylinders, cones, spheres, tori, B-splines, fillets, swept and revolved surfaces) are read from the file and meshed at the tolerance you choose. Smaller tolerance, more triangles.
  - One object per body, under an empty named after the file.
  - Normals taken from the true surface: round faces shade smoothly, edges between faces stay crisp.
  - Bodies hidden in the original file are imported hidden (`(hidden)` in their name).
- **Display Mesh**: the light mesh stored in the file for previews. Fast, with four detail levels (Coarse, As Saved, Fine, Very Fine).
- **Keep CAD Faces**: stores the CAD face of each triangle in the `cad_face` attribute, handy to select a whole face.
- Several files at once; nothing to install besides the add-on.

If the exact geometry of a part cannot be read, the add-on falls back to the display mesh and tells you so: the import never fails silently.

## Install

Blender 4.3 or newer. **Edit ‣ Preferences ‣ Get Extensions ‣ ⌄ ‣ Install from Disk…** and pick `sldprt_import-x.y.z.zip`.

## Limits (early version)

- Parts only (`.sldprt`); assemblies (`.sldasm`) and drawings are not read.
- Some large parts still have a few open edges where faces meet; a few rare fillet types may be off.
- No colours, materials or feature history: geometry only.

Bugs and sample files that don't open well: [issues](https://github.com/Elyh29/sldprt-import/issues). Please only attach files you are allowed to share.

## How it works

A `.sldprt` file stores, among other things, a preview mesh and the exact geometry of the part in the Parasolid transmit format. The add-on reads both directly from the file, in pure Python: no CAD software, no third-party library, no network access. The reading was written from the analysis of the files themselves and from public documentation of the Parasolid format.

## Trademarks

Not affiliated with or endorsed by Dassault Systèmes or Siemens. SOLIDWORKS is a registered trademark of Dassault Systèmes; Parasolid is a trademark of Siemens Digital Industries Software. They are named only to describe file compatibility.

## Licence

GPL-3.0-or-later. © 2026 Ely
