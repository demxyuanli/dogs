//! `BRepOffset_Tool::ExtentFace` (`BRepOffset_Tool.cxx:3793-4333`), the pass
//! `BRepOffset_MakeOffset` runs on every face it has to grow before the
//! offset result is rebuilt.
//!
//! Control flow, one wire at a time (`cxx:3877-4331`):
//!
//! 1. `EnLargeFace(F, EF, true)` (`cxx:3818`) builds the enlarged face `EF`
//!    that the offset face `NF` will be cut from.
//! 2. When the *surface itself* changed (`SurfaceChange`), every pcurve of `F`
//!    is re-attached to `EF` (`cxx:3826-3873`), including the patch edge the
//!    caller keeps in `ConstShapes`, and the range is re-written on every
//!    representation (`BRep_Builder::Range` 3-arg).
//! 3. `MapVertexEdges(W, MVE)` (`cxx:3886`) and `BuildNeighbour(W, F, NOnV1,
//!    NOnV2)` (`cxx:3887`) give the wire's local topology.
//! 4. First edge pass (`cxx:3895-3918`): an edge already lying on `EF`
//!    (`TryProject`) is dropped from `ToBuild` and needs no extension.
//! 5. Second edge pass (`cxx:3920-3968`): the remaining `ToBuild` edges are
//!    intersected with `EF` through `Inter3D`, keeping one section
//!    (`SelectEdge` when there are several), and the two neighbour edges of
//!    the wire are extended with `ExtentEdge`.
//! 6. Vertex pass (`cxx:3973-4136`): each endpoint of a rebuilt edge becomes
//!    the intersection of the two rebuilt neighbours, or the projection of the
//!    already known vertex, or the constant vertex re-parametrised on the new
//!    edge.
//! 7. Reconstruction (`cxx:4149-4328`): every edge of `W` is copied into a new
//!    wire, re-oriented and re-bounded by the new vertices, and the wire is
//!    added to `NF` (`cxx:4329`).
//!
//! Port guards (each one replaces an OCCT exception, never a branch):
//!
//! - `NOnV1(E)` / `NOnV2(E)` are `NCollection_DataMap::operator()` calls
//!   (`cxx:3955`, `:3961`, `:3986`, `:4061`), which raise
//!   `Standard_NoSuchObject` on an unbound key - an open wire with a single
//!   edge never reaches the closing pair of `BuildNeighbour`
//!   (`cxx:3769-3788`). The port skips the edge (and, in the vertex pass,
//!   aborts the extension) instead of throwing.
//! - `MVE(V1)` (`cxx:4023`, `:4101`) is the same `operator()`: an unbound
//!   vertex counts as zero incident edges, so the projection is skipped.
//! - `ConstShapes(V1)` / `ConstShapes(V2)` in the reconstruction (`cxx:4166`)
//!   can only be unbound after the abort above; the port drops the edge.
//! - `BRepTools::Update(NF)` (`cxx:4331`) only drops cached UV points and
//!   triangulation; the port has neither, so it is a no-op.
//!
//! `Extent` (`cxx:3806`) and the commented-out context box (`cxx:3809-3815`)
//! are dead in OCCT and are not reproduced.

use std::collections::HashMap;

use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType, State};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::shape_fix_compose_shell::empty_copied_edge;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same_occurrence, shape_map_key, ShapeMapKey};

use super::enlarge_face::EnLargeFaceOptions;
use super::pcurve_builder::build_pcurves;
use super::BRepOffsetTool;

/// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>`:
/// `ConstShapes`, `ToBuild` and the local `Build` of `ExtentFace`
/// (`cxx:3801-3804`).
pub type ShapeMap = HashMap<ShapeMapKey, TopoShape>;

impl BRepOffsetTool {
    /// `BRepOffset_Tool::ExtentFace(F, ConstShapes, ToBuild, Side, TolConf,
    /// NF)` (`cxx:3793-4333`).
    ///
    /// `const_shapes` and `to_build` are the caller's maps and are updated in
    /// place (`ToBuild` entries are unbound once handled, `ConstShapes`
    /// receives every vertex/edge the extension produced). `None` is a null
    /// `BRep_Tool::Surface(F)` (`cxx:3461` inside `EnLargeFace`) or one of the
    /// `return` statements of `cxx:4020` / `:4095`.
    pub fn extent_face(
        f: &Face,
        const_shapes: &mut ShapeMap,
        to_build: &mut ShapeMap,
        side: State,
        tol_conf: f64,
    ) -> Option<Face> {
        let reg = GeometryRegistry::global();
        let builder = TopoBuilder::new();

        // `cxx:3818-3819`: `SurfaceChange = EnLargeFace(F, EF, true)`. The
        // three-argument overload keeps OCCT's defaults, i.e.
        // `UpDatePCurve = false`, `enlargeU/Vfirst/Vlast = true`,
        // `theExtensionMode = 1` and the four lengths at `-1.`
        // (`BRepOffset_Tool.hxx:154-165`).
        let opts = EnLargeFaceOptions {
            can_extent_surface: true,
            ..Default::default()
        };
        let (ef, surface_change) = Self::enlarge_face(f, &opts)?;

        // `cxx:3822-3825`: `NF = TopoDS::Face(EF.EmptyCopied())`,
        // `NF.Orientation(FORWARD)`. `BRep_TFace::EmptyCopy` carries the
        // surface, the tolerance and the natural-restriction flag over to the
        // new TShape, so they are re-registered here.
        let mut nf = Face(ef.0.empty_copied());
        if let Some(geom) = reg.face_geom(&ef.0) {
            reg.set_face(&nf.0, geom);
        }
        nf.0.set_orientation(Orientation::Forward);

        // `cxx:3826-3873`: the pcurves of the enlarged surface. `Emap` is a
        // `TopExp::MapShapes` map, so every edge of `F` appears once.
        if surface_change {
            let fforward = Face(f.0.oriented(Orientation::Forward));
            let ef_key = GeometryRegistry::shape_key(&ef.0);
            for ce in crate::topo_tools::edges_of(&fforward.0) {
                let ce = Edge(ce.0.oriented(Orientation::Forward));
                // `BRep_Tool::CurveOnSurface(CE, Fforward, f, l)`: the COS
                // range, and the null test of `cxx:3843`.
                let Some((c2, first, last)) = boptools_2d::curve_on_surface_range(&ce, &fforward)
                else {
                    continue;
                };
                let (mut first, mut last) = (first, last);
                // `if (ConstShapes.IsBound(CE)) { Ecs = ConstShapes(CE);
                // BRep_Tool::Range(Ecs, f, l); }` (`cxx:3845-3849`): the patch
                // edge's *3D* range replaces the COS range.
                let ecs = const_shapes
                    .get(&shape_map_key(&ce.0))
                    .and_then(|s| Edge::wrap(s.clone()));
                if let Some(e) = &ecs {
                    let (a, b) = BRepTool::edge_parameters(e);
                    first = a;
                    last = b;
                }
                // `BRep_Builder::UpdateEdge(CE, C2[, C2R], EF,
                // BRep_Tool::Tolerance(CE))` (`cxx:3854-3868`) only ever
                // *raises* the edge tolerance; the pcurve registry has no
                // per-representation tolerance, so that argument is not used.
                // `cxx:3850-3871`: a closed edge gets both pcurves, an open one
                // the single `C2`.
                let c2r = if BRepTool::is_closed_edge_face(&ce, &fforward) {
                    let rev = Edge(ce.0.oriented(Orientation::Reversed));
                    boptools_2d::curve_on_surface(&rev, &fforward)
                } else {
                    None
                };
                let curves = match c2r {
                    Some(c2r) => vec![c2, c2r],
                    None => vec![c2],
                };
                reg.set_edge_pcurves(&ce.0, ef_key, curves.clone());
                if let Some(e) = &ecs {
                    reg.set_edge_pcurves(&e.0, ef_key, curves);
                }
                // `B.Range(CE, f, l)` (`cxx:3870`): `BRep_Builder::Range` with
                // `Only3d = false` re-writes every representation's range
                // (`BRep_Builder.cxx:1091-1117`).
                reg.set_ranges_all(&ce.0, first, last);
                if let Some(e) = &ecs {
                    reg.set_ranges_all(&e.0, first, last);
                }
            }
        }

        // `Build` is declared once for the whole call (`cxx:3804`) and each
        // rebuilt edge is unbound from it as soon as it is consumed
        // (`cxx:4307`).
        let mut build: ShapeMap = HashMap::new();

        // `cxx:3877-4330`.
        for w_shape in explorer_shapes(&f.0.oriented(Orientation::Forward), ShapeType::Wire) {
            let Some(w) = Wire::wrap(w_shape) else {
                continue;
            };
            // `cxx:3883-3888`.
            let mve = Self::map_vertex_edges(&w.0);
            let mut n_on_v1: ShapeMap = HashMap::new();
            let mut n_on_v2: ShapeMap = HashMap::new();
            let _ = Self::build_neighbour(&w, f, &mut n_on_v1, &mut n_on_v2);

            // `cxx:3890`: shared by both edge passes; `TryProject` and
            // `Inter3D` clear them on entry.
            let mut lint1: Vec<TopoShape> = Vec::new();
            let mut lint2: Vec<TopoShape> = Vec::new();
            let edges = explorer_shapes(&w.0.oriented(Orientation::Forward), ShapeType::Edge);

            //------------------------------------------------------------
            // `cxx:3895-3918`: an edge lying on the enlarged face needs no
            // extension.
            //------------------------------------------------------------
            for e_shape in &edges {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                let ek = shape_map_key(&e.0);
                if const_shapes.contains_key(&ek) {
                    to_build.remove(&ek);
                }
                let Some(tb) = to_build.get(&ek).and_then(|s| Face::wrap(s.clone())) else {
                    continue;
                };
                // `TryProject(ToBuild(E), EF, LOE, LInt2, LInt1, Side,
                // TolConf)` (`cxx:3906-3912`): OCCT fills its first output with
                // the edges oriented on `ToBuild(E)` and its second with the
                // same edges oriented on `EF`, so the fourth and fifth
                // arguments are swapped with respect to the local lists. The
                // test below is on `LInt1`, i.e. on the `EF` orientation.
                let mut on_tb: Vec<Edge> = Vec::new();
                let mut on_ef: Vec<Edge> = Vec::new();
                let ok = Self::try_project(
                    &tb,
                    &ef,
                    std::slice::from_ref(&e),
                    &mut on_tb,
                    &mut on_ef,
                    side,
                    tol_conf,
                );
                lint2 = on_tb.into_iter().map(|x| x.0).collect();
                lint1 = on_ef.into_iter().map(|x| x.0).collect();
                if ok && !lint1.is_empty() {
                    to_build.remove(&ek);
                }
            }

            //------------------------------------------------------------
            // `cxx:3920-3968`: intersect the remaining edges with `EF` and
            // extend the wire's neighbours.
            //------------------------------------------------------------
            for e_shape in &edges {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                let ek = shape_map_key(&e.0);
                if const_shapes.contains_key(&ek) {
                    to_build.remove(&ek);
                }
                let Some(tb) = to_build.get(&ek).and_then(|s| Face::wrap(s.clone())) else {
                    continue;
                };
                // `EnLargeFace(TopoDS::Face(ToBuild(E)), StopFace, false)`
                // (`cxx:3929`).
                let stop_opts = EnLargeFaceOptions {
                    can_extent_surface: false,
                    ..Default::default()
                };
                let Some((stop_face, _)) = Self::enlarge_face(&tb, &stop_opts) else {
                    // Port guard: `EnLargeFace` dereferences
                    // `BRep_Tool::Surface(F)`.
                    continue;
                };
                // `Inter3D(EF, StopFace, LInt1, LInt2, Side, E, NullFace,
                // NullFace)` (`cxx:3931`): the section edges oriented on `EF`
                // land in `LInt1`.
                Self::inter3d(
                    &ef,
                    &stop_face,
                    &mut lint1,
                    &mut lint2,
                    side,
                    Some(&e),
                    None,
                    None,
                );
                // `cxx:3936-3939`: no intersection, e.g. a chosen planar face
                // with an offsetted cylindrical neighbour whose radius
                // shrinks.
                if lint1.is_empty() {
                    continue;
                }
                // `cxx:3939-3943`: several edges means a seam crossing.
                if lint1.len() > 1 {
                    let mut selected: Vec<Edge> = lint1
                        .iter()
                        .filter_map(|s| Edge::wrap(s.clone()))
                        .collect();
                    Self::select_edge(&e, &mut selected);
                    lint1 = selected.into_iter().map(|x| x.0).collect();
                    if lint1.is_empty() {
                        // Port guard: `LInt1.First()` (`cxx:3945`) raises on
                        // the empty list `SelectEdge` can leave behind.
                        continue;
                    }
                }
                let Some(ne) = lint1.first().and_then(|s| Edge::wrap(s.clone())) else {
                    continue;
                };
                // `TE->Tolerance(TE->Tolerance() * 10.)` (`cxx:3946-3947`: the
                // `BRep_TEdge` is written directly, so this is not
                // `UpdateTolerance` - it may lower a larger value).
                reg.set_edge_tolerance(&ne.0, BRepTool::edge_tolerance(&ne) * 10.0);
                // `cxx:3948-3953`.
                let ne_oriented = if ne.0.orientation() == e.0.orientation() {
                    ne.0.oriented(Orientation::Forward)
                } else {
                    ne.0.oriented(Orientation::Reversed)
                };
                build.insert(ek.clone(), ne_oriented);

                // `cxx:3954-3967`: extend both wire neighbours on `EF`. The
                // `NE` of the previous block is overwritten by `ExtentEdge`
                // (`cxx:2081`), so each neighbour gets its own edge.
                for n_map in [&n_on_v1, &n_on_v2] {
                    let Some(e_on_v) = n_map.get(&ek).and_then(|s| Edge::wrap(s.clone())) else {
                        // Port guard: `TopoDS::Edge(NOnV1(E))` raises in OCCT.
                        continue;
                    };
                    let nk = shape_map_key(&e_on_v.0);
                    if to_build.contains_key(&nk)
                        || const_shapes.contains_key(&nk)
                        || build.contains_key(&nk)
                    {
                        continue;
                    }
                    if let Some(ext) = Self::extent_edge(f, &ef, &e_on_v) {
                        build.insert(nk, ext.0.oriented(Orientation::Forward));
                    }
                }
            }

            //------------------------------------------------------------
            // `cxx:3973-4136`: the endpoints of every rebuilt edge.
            //------------------------------------------------------------
            for e_shape in &edges {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                let (v1, v2) = edge_vertices(&e);
                let (Some(v1), Some(v2)) = (v1, v2) else {
                    continue;
                };
                let ek = shape_map_key(&e.0);
                let Some(be) = build.get(&ek).and_then(|s| Edge::wrap(s.clone())) else {
                    continue;
                };
                let ne_on_v1 = n_on_v1.get(&ek).and_then(|s| Edge::wrap(s.clone()));
                let ne_on_v2 = n_on_v2.get(&ek).and_then(|s| Edge::wrap(s.clone()));
                // `cxx:3987-4057` then `cxx:4062-4134`: the first endpoint uses
                // `LV.First()` / `TopExp::FirstVertex` / `FORWARD`, the last the
                // mirrored choice.
                let v_for_v1 = Self::extent_face_vertex(
                    true,
                    &e,
                    &v1,
                    &be,
                    ne_on_v1.as_ref(),
                    &mve,
                    const_shapes,
                    to_build,
                    &build,
                    &ef,
                    tol_conf,
                )?;
                const_shapes.insert(shape_map_key(&v1.0), v_for_v1.0.clone());
                build.insert(shape_map_key(&v1.0), v_for_v1.0);
                let v_for_v2 = Self::extent_face_vertex(
                    false,
                    &e,
                    &v2,
                    &be,
                    ne_on_v2.as_ref(),
                    &mve,
                    const_shapes,
                    to_build,
                    &build,
                    &ef,
                    tol_conf,
                )?;
                const_shapes.insert(shape_map_key(&v2.0), v_for_v2.0.clone());
                build.insert(shape_map_key(&v2.0), v_for_v2.0);
            }

            //------------------------------------------------------------
            // `cxx:4140-4329`: rebuild the wire.
            //------------------------------------------------------------
            // `B.MakeWire(NW)` (`cxx:4147`): a fresh empty wire, orientation
            // FORWARD.
            let mut nw = Wire(TopoShape::new(ShapeType::Wire));
            for e_shape in &edges {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                let (v1, v2) = edge_vertices(&e);
                let (Some(v1), Some(v2)) = (v1, v2) else {
                    continue;
                };
                let ek = shape_map_key(&e.0);
                let mut ne;
                if let Some(be) = build.get(&ek).and_then(|s| Edge::wrap(s.clone())) {
                    // `cxx:4155-4160`: `BRep_Tool::Range(NE, f, l)`.
                    let (ne_first, ne_last) = BRepTool::edge_parameters(&be);
                    let or = be.0.orientation();
                    let (Some(nv1), Some(nv2)) = (
                        const_shapes
                            .get(&shape_map_key(&v1.0))
                            .and_then(|s| Vertex::wrap(s.clone())),
                        const_shapes
                            .get(&shape_map_key(&v2.0))
                            .and_then(|s| Vertex::wrap(s.clone())),
                    ) else {
                        // Port guard: `ConstShapes(V1)` raises in OCCT.
                        continue;
                    };
                    // `cxx:4163-4179`: `U1` on `NE.Oriented(INTERNAL)`, `U2` on
                    // `NE.Oriented(FORWARD)`.
                    let mut u1 = BRepTool::parameter_of_vertex(
                        &be.0.oriented(Orientation::Internal),
                        &nv1.0.oriented(Orientation::Internal),
                    );
                    let mut u2 = BRepTool::parameter_of_vertex(
                        &be.0.oriented(Orientation::Forward),
                        &nv2.0.oriented(Orientation::Internal),
                    );
                    // `cxx:4181-4184`: `NE = NE.EmptyCopied()`,
                    // `NE.Orientation(FORWARD)`; `BRep_TEdge::EmptyCopy`
                    // (`BRep_TEdge.cxx:104-129`) keeps the curve
                    // representations, so the pcurves of the built edge survive.
                    let mut ne_e = Edge(empty_copied_edge(&be.0));
                    ne_e.0.set_orientation(Orientation::Forward);
                    // `eps` of `cxx:4142`.
                    let eps = CONFUSION;
                    if is_same_occurrence(&nv1.0, &nv2.0) {
                        // `cxx:4186-4248`: a closed edge; the two bounds are
                        // taken from the range written at `cxx:4156`.
                        if or == Orientation::Forward {
                            u1 = ne_first;
                            u2 = ne_last;
                        } else {
                            u1 = ne_last;
                            u2 = ne_first;
                        }
                        if or == Orientation::Forward {
                            if u1 > u2 {
                                if (u1 - ne_last).abs() < eps {
                                    u1 = ne_first;
                                }
                                if (u2 - ne_first).abs() < eps {
                                    u2 = ne_last;
                                }
                            }
                            builder.add(&mut ne_e.0, &nv1.0.oriented(Orientation::Forward));
                            builder.add(&mut ne_e.0, &nv2.0.oriented(Orientation::Reversed));
                            reg.set_ranges_all(&ne_e.0, u1, u2);
                            const_shapes.insert(ek.clone(), ne_e.0.clone());
                            ne_e.0.set_orientation(e.0.orientation());
                        } else {
                            if u2 > u1 {
                                if (u2 - ne_last).abs() < eps {
                                    u2 = ne_first;
                                }
                                if (u1 - ne_first).abs() < eps {
                                    u1 = ne_last;
                                }
                            }
                            builder.add(&mut ne_e.0, &nv2.0.oriented(Orientation::Forward));
                            builder.add(&mut ne_e.0, &nv1.0.oriented(Orientation::Reversed));
                            reg.set_ranges_all(&ne_e.0, u2, u1);
                            const_shapes.insert(ek.clone(), ne_e.0.oriented(Orientation::Reversed));
                            ne_e.0.set_orientation(e.0.orientation().reversed());
                        }
                    } else if or == Orientation::Forward {
                        // `cxx:4250-4279`.
                        if u1 > u2 {
                            builder.add(&mut ne_e.0, &nv2.0.oriented(Orientation::Forward));
                            builder.add(&mut ne_e.0, &nv1.0.oriented(Orientation::Reversed));
                            reg.set_ranges_all(&ne_e.0, u2, u1);
                        } else {
                            builder.add(&mut ne_e.0, &nv1.0.oriented(Orientation::Forward));
                            builder.add(&mut ne_e.0, &nv2.0.oriented(Orientation::Reversed));
                            reg.set_ranges_all(&ne_e.0, u1, u2);
                        }
                        const_shapes.insert(ek.clone(), ne_e.0.clone());
                        ne_e.0.set_orientation(e.0.orientation());
                    } else {
                        // `cxx:4279-4304`.
                        if u2 > u1 {
                            builder.add(&mut ne_e.0, &nv1.0.oriented(Orientation::Forward));
                            builder.add(&mut ne_e.0, &nv2.0.oriented(Orientation::Reversed));
                            reg.set_ranges_all(&ne_e.0, u1, u2);
                            const_shapes.insert(ek.clone(), ne_e.0.clone());
                            ne_e.0.set_orientation(e.0.orientation());
                        } else {
                            builder.add(&mut ne_e.0, &nv2.0.oriented(Orientation::Forward));
                            builder.add(&mut ne_e.0, &nv1.0.oriented(Orientation::Reversed));
                            reg.set_ranges_all(&ne_e.0, u2, u1);
                            const_shapes.insert(ek.clone(), ne_e.0.oriented(Orientation::Reversed));
                            ne_e.0.set_orientation(e.0.orientation().reversed());
                        }
                    }
                    ne = ne_e;
                    // `cxx:4307`: the rebuilt edge is consumed.
                    build.remove(&ek);
                } else if let Some(cs) = const_shapes.get(&ek).and_then(|s| Edge::wrap(s.clone()))
                {
                    // `cxx:4309-4321`: a constant edge only needs its pcurves
                    // on `NF`.
                    ne = cs;
                    let _ = build_pcurves(&ne, &nf);
                    let or = ne.0.orientation();
                    let new_or = if or == Orientation::Reversed {
                        e.0.orientation().reversed()
                    } else {
                        e.0.orientation()
                    };
                    ne.0.set_orientation(new_or);
                } else {
                    // `cxx:4322-4326`: an untouched edge is copied as is and
                    // registered forward.
                    ne = e.clone();
                    const_shapes.insert(ek, e.0.oriented(Orientation::Forward));
                }
                // `B.Add(NW, NE)` (`cxx:4327`).
                builder.add(&mut nw.0, &ne.0);
            }
            // `B.Add(NF, NW.Oriented(W.Orientation()))` (`cxx:4329`).
            builder.add(&mut nf.0, &nw.0.oriented(w.0.orientation()));
        }

        // `cxx:4330-4331`.
        nf.0.set_orientation(f.0.orientation());
        // `BRepTools::Update(NF)`: only the cached UV points / triangulation
        // change; the port keeps neither.
        Some(nf)
    }

    /// The first (`cxx:3983-4057`) or last (`cxx:4059-4134`) endpoint of a
    /// rebuilt edge: the shared vertex of the two rebuilt neighbours, the
    /// projection of an already known vertex, or the constant vertex
    /// re-parametrised on the new edge.
    ///
    /// `None` groups the two aborts OCCT spells differently but that both leave
    /// `ExtentFace`: the bare `return` of `cxx:4020` / `:4095` (`Inter2d` found
    /// no intersection) and the `NCollection_DataMap::operator()` exception on
    /// an unbound neighbour / vertex.
    #[allow(clippy::too_many_arguments)]
    fn extent_face_vertex(
        at_first: bool,
        e: &Edge,
        the_v: &Vertex,
        be: &Edge,
        neighbour: Option<&Edge>,
        mve: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        const_shapes: &ShapeMap,
        to_build: &ShapeMap,
        build: &ShapeMap,
        ef: &Face,
        tol_conf: f64,
    ) -> Option<Vertex> {
        let ek = shape_map_key(&e.0);
        let vk = shape_map_key(&the_v.0);
        let ne = neighbour?;
        let nk = shape_map_key(&ne.0);
        // `Build.IsBound(NEOnV1) && (ToBuild.IsBound(E) ||
        // ToBuild.IsBound(NEOnV1))` (`cxx:3988`, `:4062`).
        let neighbour_built =
            build.contains_key(&nk) && (to_build.contains_key(&ek) || to_build.contains_key(&nk));
        if neighbour_built {
            // `E.IsSame(NEOnV1)` (`cxx:3989`, `:4064`): a closed edge is its
            // own neighbour on both endpoints.
            if is_same_occurrence(&e.0, &ne.0) {
                let (first, last) = edge_vertices(be);
                return if at_first { first } else { last };
            }
            if !build.contains_key(&vk) {
                // `cxx:3995-4022` / `:4074-4097`: the endpoint is the
                // intersection of the two rebuilt edges.
                let bne = Edge(build.get(&nk)?.clone());
                let mut lv: Vec<Vertex> = Vec::new();
                Self::inter2d(ef, be, &bne, &mut lv, CONFUSION);
                if lv.is_empty() {
                    return None;
                }
                let forward = be.0.orientation() == Orientation::Forward;
                let pick_first = if at_first { forward } else { !forward };
                return if pick_first {
                    lv.first().cloned()
                } else {
                    lv.last().cloned()
                };
            }
            // `cxx:4025-4036` / `:4100-4110`: the vertex is already known and
            // is projected onto the new edge when more than two edges meet
            // there.
            let mut v = Vertex(build.get(&vk)?.clone());
            let incident = mve.get(&vk).map_or(0, |v| v.len());
            if incident > 2 {
                let mut o = if at_first {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                if be.0.orientation() == Orientation::Reversed {
                    o = o.reversed();
                }
                v.0.set_orientation(o);
                Self::project_vertex_on_edge(&v, be, tol_conf);
            }
            return Some(v);
        }
        // `cxx:4038-4057` / `:4113-4133`: plain projection of the constant
        // vertex (or of the one already put in `ConstShapes`).
        let mut v = const_shapes
            .get(&vk)
            .and_then(|s| Vertex::wrap(s.clone()))
            .unwrap_or_else(|| the_v.clone());
        let mut o = if at_first {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        if be.0.orientation() == Orientation::Reversed {
            o = o.reversed();
        }
        v.0.set_orientation(o);
        if !Self::try_parameter(e, &v, be, tol_conf) {
            Self::project_vertex_on_edge(&v, be, tol_conf);
        }
        Some(v)
    }
}
