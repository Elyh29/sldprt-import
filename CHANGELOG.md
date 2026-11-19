# Changelog

## 0.1.9
- Blender no longer freezes during the import: the part is meshed step by step, each body appears as soon as it is done, the progress shows in the status bar, Esc stops the import and keeps the bodies already made.
- Undo is held back while the import runs.

## 0.1.8
- Hidden bodies read from the file: a body is hidden when SolidWorks does not show it (its faces are not in the display mesh of the part), instead of being guessed from distances.
- New option Hidden Bodies: meshed coarsely (default, faster), at the full tolerance, or left out.

## 0.1.7
- Closed meshes on large parts: faces far from the origin lost points along their edges (thin slits between faces), and edges meeting at a sharp angle could cut across each other.
- Offset surfaces (skins with a thickness) meshed up to 80 times faster, with even triangles.
- Large bodies stored in several blocks are read whole.
- A large curved panel is no longer taken for a hidden body.

## 0.1.6
- Parts saved by recent versions (Parasolid 30 and later): the exact geometry is now read; before, the display mesh was used instead.
- Meshes follow the surfaces within the tolerance: curved edges and faces are refined where needed (large skins and thin strips were faceted).
- No more faces meshed flat on the test parts: cone tips, spheres stored as revolved circles, closed edges without vertex, fillets along circles, ellipses and splines.
- Far fewer open edges: faces that wrap around (cylinders, cones, tori) share the points of their neighbours.
- Sphere caps and rings meshed on the right side.
- Faster on parts with fillets and splines.

## 0.1.5
- Parts that keep their bodies per feature (imported bodies, recent files): the exact geometry is now read; each body is named after its feature.
- Rings and torus-like surfaces of revolution meshed on the right side (they could stray by several centimetres).
- When the exact geometry gives nothing, the display mesh is used with a warning, instead of an empty object.
- A cut file no longer shows a Python error.

## 0.1.4
- First public release.

## 0.1.3
- Large parts (thousands of faces) read correctly.
- A face that cannot be meshed no longer stops the import: it is meshed flat, and the display mesh is used if the whole geometry is unreadable.

## 0.1.2
- Exact geometry: surfaces read from the file and meshed at a tolerance, one object per body, normals from the true surface, hidden bodies imported hidden.

## 0.1.1
- Detail levels for the display mesh (Coarse, As Saved, Fine, Very Fine).

## 0.1.0
- Display mesh import, drag and drop, several files at once.
