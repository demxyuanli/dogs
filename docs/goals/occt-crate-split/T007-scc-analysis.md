# T007 — `occt-topo` SCC and back-edge inventory

Scope: the 449 module directories/files left in `crates/occt-topo/src` after T005
plus the 17 extracted by T009. Edge definition: `crate::<module>` occurrences
inside a module's `.rs` files, with `//` comment lines and `#[cfg(test)]` blocks
excluded (the raw-with-comments counts are in `brep-back-edges.txt`).

## 1. The knot

`SCC(shape)` — the modules mutually reachable with `shape` — has **204 members**
(list in `scc-shape.txt`). It spans every layer the SRS wanted to separate:

- shape representation: `shape`, `tshape`, `tgeometry`, `iterator`, `transform`, `topexp`
- BRep access: `brep_tool`, `brep_tools`, `brep_surface`, `pcurve`, `pcurve_full`
- boolean engine: `bop_*`, `pave_*`, `bopds`, `boptools_2d`
- intersection: `int_face_face`, `intpatch_*`, `intwalk_*`, `int_tools_*`
- blending: `chfi3d`, `chfi_ds`, `chfi_k_part`, `fillet_*`, `blend`, `draft`
- meshing: `meshing`, `fclass2d`, `shape_mesh`, `mesh`
- healing: `shhealing`, `shape_fix_compose_shell`, `shell_splitter`, `wire_splitter`

Because the SCC is irreducible at module granularity, any split of it must
invert (or accept) a set of back-edges. "Algorithm crates may not depend on each
other" (§4.1.3) is therefore not achievable inside this SCC; what is achievable
is a set of extracted crates whose cross-crate edges follow a DAG, paid for by a
finite, enumerated number of inverted edges.

## 2. Candidate boundaries for the `occt-brep` / `occt-topology` split

`A` is the lower crate; a back-edge is `A -> outside(A)`, i.e. code in the lower
crate that reaches up into an algorithm.

| id | modules | members | back-edges |
|----|---------|---------|-----------|
| B1  | TKBRep without `tgeometry` / `topo_tools` | 15 | 9 |
| B1b | B1 + `tgeometry` + `topo_tools` | 17 | **1** |
| B2  | B1b + `brep_tools{,_uv,_modifier,_wire_explorer}`, `brep_surface`, `brep_bnd_lib`, `brep_uv_bounds`, `brep_faces`, `brep_measure` | 26 | 17 |
| B3  | B2 + `pcurve`, `pcurve_full`, `brep_check`, `shape_analysis{,_edge}`, `locope`, `bnd_tools`, `brep_adaptor_*`, `shape_naming` | 36 | 61 |

B1's 9 edges are all `tgeometry`/`topo_tools` references, so B1b is the honest
minimal TKBRep boundary. B3's 61 edges are dominated by `locope` (26) — that
module is a consumer of the whole algorithm stack and does not belong in the
lower crate at all.

Full lists: `cut-B1.txt`, `cut-B1b.txt`, `cut-B2.txt`, `cut-B3.txt`.

## 3. Decision (T008)

Approved boundary: **B1b**, extracted as `occt-topology` (T009). It is the only
candidate whose inversion cost is a single edge, and it is exactly the TKBRep
scope the SRS names for `occt-brep`:

```
abs, shape, tshape, iterator, transform, topexp, tools, validate,
debug_flags, builder, brep_tool, brep_point_representation,
brep_builder_update_vertex, brep_lprop, brep_tools_substitution,
tgeometry, topo_tools
```

The one back-edge and its inversion:

```
tgeometry -> pcurve_full::classify_surface_kind
```

`classify_surface_kind` + `SurfaceKind` are an `Adaptor3d_Surface::GetType()`
wrapper with no algorithm content. They moved down to
`occt_geom::surface_kind` and both `pcurve_full` and `tgeometry` now import from
there, which is a genuine dependency inversion rather than a moved predicate.

Deferred boundary: **B2** (17 edges) is the next tier. Its edges are port
defects rather than architectural necessities — `brep_tools -> mesh`,
`brep_tools -> meshing` exist because `BRep_Tool::Triangulation` reads the
triangulation through `meshing`, whereas OCCT's `BRep_Tool` only touches
`Poly_Triangulation`. Those are fixable one at a time with the `.cxx` branch as
the reference.

Not pursued: **B3**'s `locope` (26 edges) and the rest of the 204-SCC. `locope`
depends on the entire algorithm stack by construction, so it is an algorithm
module that happens to look like a `BRepTools` sibling.

## 4. Why the remaining 204-SCC is not split further

The remaining SCC cannot be partitioned without inventing dependency-inversion
layers that OCCT does not have. OCCT itself keeps TKBRep and TKTopAlgo as
toolkits of the same library, linked as one unit; the CALL-DAG inside them is
cyclic in the same way. Reproducing a DAG here would mean rewriting call
direction (`shape_mesh` reaching into `brep_tool`, `fclass2d` reaching into
`meshing`, `tshape` reaching into `brep_tool`), which changes control flow, and
§2 of the porting gates forbids that.
