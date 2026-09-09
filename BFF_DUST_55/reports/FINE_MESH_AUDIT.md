# FINE mesh audit

Classification: **production diagnostic mesh with residual convergence
uncertainty**. Production status: **blocked**.

The reconstructed FINE discretization contains 4284 quadrilateral surface
panels, 4445 points, and 59 coupling nodes. Coordinates and panel normals are
finite; no zero-normal or exactly duplicated panels were found. Left/right
point symmetry is exact at the audit tolerance. The six planform-region
boundaries are present exactly, and no non-manifold edges were found.

Panel area ranges from 0.3930 to 17.1508 in². Aspect ratio ranges from 1.0253 to
24.8568, with mean 3.4456. About 98.60% of panels pass the local upper/lower
outward-normal test; the remaining fraction requires visual inspection around
leading/trailing-edge transitions before a production claim.

The audit does **not** pass the scientific gate:

- the available geometry has no winglet, so no wing/winglet junction exists;
- the original mesh-study input files were not present and the levels were
  reconstructed from the available no-winglet topology;
- most requested hinge endpoints are not exact spanwise mesh lines (nearest
  errors reach 1.3378 in);
- general geometric intersection testing has not certified the absence of all
  possible non-exact overlaps;
- the supplied MEDIUM→FINE differences remain material, especially pitching
  moment and spanwise distribution.

Open `fine_mesh_audit.vtu` together with `fine_coupling_nodes.vtu` in ParaView.
Color the mesh by `panel_aspect_ratio`, display surface edges, and render the
coupling-node file as points or glyphs.
