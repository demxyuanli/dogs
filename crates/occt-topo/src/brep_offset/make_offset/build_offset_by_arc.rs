//! `BRepOffset_MakeOffset::BuildOffsetByArc` (`BRepOffset_MakeOffset.cxx:
//! 1893-2152`) and `ToContext` (`:2171-2340`), the arc-join construction.
//!
//! `BuildOffsetByArc` runs in five steps, each one the `.cxx` block:
//!
//! 1. `:1911-1918` `MakeOffsetFaces(MapSF)`: one `BRepOffset_Offset` per face,
//!    which fills `MapSF` with the face offsets.
//! 2. `:1920-2006` tubes on the edges of `myFaceComp`: for every edge whose two
//!    ancestors have the concavity `OT` (`Convex` for a positive offset,
//!    `Concave` for a negative one), take the two offset edges the faces
//!    generated and build `BRepOffset_Offset(E, EOn1, EOn2, CurOffset,
//!    E1f, E1l)` where `E1f` / `E1l` are the end edges the adjacent tubes
//!    already generated. An edge with a single ancestor is a free border: its
//!    offset edge goes straight into `myInitOffsetEdge`.
//! 3. `:2008-2073` spheres on the vertices of `myFaceComp`: a vertex with three
//!    or more incident edges, all of them `OT` and all of them tangential,
//!    becomes `BRepOffset_Offset(V, LOE, CurOffset)`; the free-border block
//!    then rebinds every `LOE` image into `myInitOffsetEdge`.
//! 4. `:2079-2082` `ToContext(MapSF)` when `myFaces` is not empty.
//! 5. `:2087-2147` the `myInitOffsetFace` / `myImageOffset` / `myAsDes` update
//!    of every offset kept in `MapSF`, skipping the `Reversed` and
//!    `Degenerated` ones.
//!
//! Port guards, each replacing an OCCT exception or a null-handle
//! dereference, never a branch:
//!
//! * `MapSF(Anc.First())` / `MapSF(E)` (`:1947-1950`, `:1997`, `:2034`,
//!   `:2065`) raise `Standard_NoSuchObject` on a missing key, and
//!   `TopoDS::Edge` / `TopoDS::Face` on the resulting null shape then runs into
//!   an exception inside `BRepOffset_Offset`. The port records
//!   [`Error::UnknownError`] and returns, like `MakeOffsetFaces` does for the
//!   constructor exceptions.
//! * `ExtentFace` (`:2260`) answering `None` is a null `NF`; the port skips the
//!   face and leaves it out of `myInitOffsetFace`, instead of OCCT's reuse of
//!   the previous iteration's `NF`.
//! * `Message_ProgressScope` (`:1905`, `:1926-1932`, `:2011-2017`, `:2094-
//!   2100`): no progress machinery in the port, same as the rest of
//!   `BRepOffset`.

use std::collections::{HashMap, HashSet};

use occt_core::kernel::geomabs::Shape;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_offset::{
    BRepOffsetOffset, BRepOffsetStatus, BRepOffsetTool, Error, MapSf, ShapeMap,
};
use crate::chfi_ds::TypeOfConcavity;
use crate::shape::{TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key, ShapeMapKey};

use super::make_offset_faces::make_offset_faces;
use super::BRepOffsetMakeOffset;

/// `BuildOffsetByArc(theRange)` (`cxx:1893-2152`).
pub fn build_offset_by_arc(this: &mut BRepOffsetMakeOffset) {
    // Key order of `MapSF`: the `NCollection_DataMap` the `.cxx` iterates at
    // `:2094` is a hash map too, so this only has to stay stable within a run.
    let mut order: Vec<TopoShape> = Vec::new();
    let mut map_sf: MapSf = HashMap::new();

    // `cxx:1911-1918`.
    make_offset_faces(this, &mut map_sf);
    if this.my_error != Error::NoError {
        return;
    }
    // Collect the keys `MakeOffsetFaces` bound. It is not observable from the
    // outside, so the tubes step below records them instead: every face of
    // `myFaceComp` that `MakeOffsetFaces` walked has an entry.
    for f in explorer_shapes(&this.my_face_comp.0, ShapeType::Face) {
        if map_sf.contains_key(&shape_map_key(&f)) {
            order.push(f);
        }
    }

    // `cxx:1920-1924`.
    let ot = if this.my_offset < 0.0 {
        TypeOfConcavity::Concave
    } else {
        TypeOfConcavity::Convex
    };

    // `cxx:1926-2006`: tubes on the edges.
    let mut done: HashSet<ShapeMapKey> = HashSet::new();
    for s in explorer_shapes(&this.my_face_comp.0, ShapeType::Edge) {
        let Some(e) = crate::shape::Edge::wrap(s) else {
            continue;
        };
        if !done.insert(shape_map_key(&e.0)) {
            continue;
        }
        let anc: Vec<TopoShape> = this.my_analyse.ancestors(&e.0).to_vec();
        if anc.len() == 2 {
            let l: Vec<crate::brep_offset::Interval> = this.my_analyse.type_of(&e).to_vec();
            if l.is_empty() || l[0].kind() != ot {
                continue;
            }
            let cur_offset = this.my_face_offset.get(&anc[0]).unwrap_or(this.my_offset);
            let (Some(gen1), Some(gen2)) = (
                map_sf
                    .get(&shape_map_key(&anc[0]))
                    .and_then(|of| of.generated(&e.0)),
                map_sf
                    .get(&shape_map_key(&anc[1]))
                    .and_then(|of| of.generated(&e.0)),
            ) else {
                this.my_error = Error::UnknownError;
                return;
            };
            let (Some(mut e_on1), Some(mut e_on2)) = (
                crate::shape::Edge::wrap(gen1),
                crate::shape::Edge::wrap(gen2),
            ) else {
                this.my_error = Error::UnknownError;
                return;
            };

            // `cxx:1952-1959`: the tangent edges at both endpoints.
            let (v1f, v1l) = edge_vertices(&e);
            let mut e1f: Option<crate::shape::Edge> = None;
            let mut e1l: Option<crate::shape::Edge> = None;
            if let Some(v) = v1f.as_ref() {
                let mut tang: Vec<TopoShape> = Vec::new();
                this.my_analyse.tangent_edges(&e, v, &mut tang);
                for t in &tang {
                    if let Some(of) = map_sf.get(&shape_map_key(t)) {
                        e1f = of.generated(&v.0).and_then(crate::shape::Edge::wrap);
                        break;
                    }
                }
            }
            if let Some(v) = v1l.as_ref() {
                let mut tang: Vec<TopoShape> = Vec::new();
                this.my_analyse.tangent_edges(&e, v, &mut tang);
                for t in &tang {
                    if let Some(of) = map_sf.get(&shape_map_key(t)) {
                        e1l = of.generated(&v.0).and_then(crate::shape::Edge::wrap);
                        break;
                    }
                }
            }

            // `cxx:1988-1989`.
            let mut of = BRepOffsetOffset::new();
            if of
                .init_pipe_first_last(
                    &e,
                    &mut e_on1,
                    &mut e_on2,
                    cur_offset,
                    e1f.as_ref(),
                    e1l.as_ref(),
                    false,
                    1.0e-4,
                    Shape::C1,
                )
                .is_err()
            {
                this.my_error = Error::UnknownError;
                return;
            }
            map_sf.insert(shape_map_key(&e.0), of);
            order.push(e.0.clone());
        } else {
            // `cxx:1993-2002`: a free border.
            let Some(a0) = anc.first() else {
                this.my_error = Error::UnknownError;
                return;
            };
            let Some(gen) = map_sf
                .get(&shape_map_key(a0))
                .and_then(|of| of.generated(&e.0))
            else {
                this.my_error = Error::UnknownError;
                return;
            };
            this.my_init_offset_edge.set_root(&e.0);
            let _ = this.my_init_offset_edge.bind(&e.0, &gen);
        }
    }

    // `cxx:2008-2073`: spheres on the vertices.
    let mut done: HashSet<ShapeMapKey> = HashSet::new();
    for s in explorer_shapes(&this.my_face_comp.0, ShapeType::Vertex) {
        let Some(v) = Vertex::wrap(s) else {
            continue;
        };
        if !done.insert(shape_map_key(&v.0)) {
            continue;
        }
        let la: Vec<TopoShape> = this.my_analyse.ancestors(&v.0).to_vec();
        let le = this.my_analyse.edges_of_vertex(&v, ot);

        if le.len() >= 3 && le.len() == la.len() {
            // `cxx:2032-2037`: the tubes' images of `V`, reversed.
            let mut loe: Vec<crate::shape::Edge> = Vec::with_capacity(le.len());
            for it in &le {
                let Some(gen) = map_sf
                    .get(&shape_map_key(&it.0))
                    .and_then(|of| of.generated(&v.0))
                else {
                    this.my_error = Error::UnknownError;
                    return;
                };
                let mut g = gen;
                // `MapSF(...).Generated(V).Reversed()` (`cxx:2034`): an
                // unconditional orientation flip, not "make it Reversed".
                g.reverse();
                let Some(ge) = crate::shape::Edge::wrap(g) else {
                    this.my_error = Error::UnknownError;
                    return;
                };
                loe.push(ge);
            }

            // `cxx:2038-2045`: the offset of the first ancestor.
            let Some(la0) = la.first() else {
                this.my_error = Error::UnknownError;
                return;
            };
            let lla: Vec<TopoShape> = this.my_analyse.ancestors(la0).to_vec();
            let Some(ff) = lla.first() else {
                this.my_error = Error::UnknownError;
                return;
            };
            let cur_offset = this.my_face_offset.get(ff).unwrap_or(this.my_offset);

            // `cxx:2046-2048`.
            let mut of = BRepOffsetOffset::new();
            if of
                .init_vertex(&v, &loe, cur_offset, false, 1.0e-4, Shape::C1)
                .is_err()
            {
                this.my_error = Error::UnknownError;
                return;
            }
            map_sf.insert(shape_map_key(&v.0), of);
            order.push(v.0.clone());
        }

        // `cxx:2052-2071`: a vertex on at least one free border.
        let lbf = this
            .my_analyse
            .edges_of_vertex(&v, TypeOfConcavity::FreeBound);
        if !lbf.is_empty() {
            let mut first = true;
            for it in &le {
                let Some(gen) = map_sf
                    .get(&shape_map_key(&it.0))
                    .and_then(|of| of.generated(&v.0))
                else {
                    this.my_error = Error::UnknownError;
                    return;
                };
                if first {
                    this.my_init_offset_edge.set_root(&v.0);
                    let _ = this.my_init_offset_edge.bind(&v.0, &gen);
                    first = false;
                } else {
                    let _ = this.my_init_offset_edge.add(&v.0, &gen);
                }
            }
        }
    }

    // `cxx:2079-2082`.
    if this.my_faces.extent() != 0 {
        to_context(this, &mut map_sf);
    }

    // `cxx:2087-2147`: update `myInitOffsetFace`, `myImageOffset` and
    // `myAsDes`.
    let rt = if this.my_offset < 0.0 {
        TypeOfConcavity::Convex
    } else {
        TypeOfConcavity::Concave
    };
    for si in order {
        let Some(sf) = map_sf.get(&shape_map_key(&si)) else {
            continue;
        };
        // `cxx:2105-2109`: `Reversed` / `Degenerated` offsets are not stored.
        if matches!(
            sf.status(),
            BRepOffsetStatus::Reversed | BRepOffsetStatus::Degenerated
        ) {
            continue;
        }
        let Some(of) = sf.face().cloned() else {
            continue;
        };

        // `cxx:2112-2115`.
        let _ = this.my_init_offset_face.bind(&si, &of.0);
        this.my_init_offset_face.set_root(&si);
        this.my_image_offset.set_root(&of.0);

        if si.shape_type() == ShapeType::Face {
            // `cxx:2118-2136`: the edges of the face whose interval type is
            // not `RT`.
            for es in explorer_shapes(&si.oriented(Orientation::Forward), ShapeType::Edge) {
                let Some(e) = crate::shape::Edge::wrap(es) else {
                    continue;
                };
                let l: Vec<crate::brep_offset::Interval> = this.my_analyse.type_of(&e).to_vec();
                if l.is_empty() || l[0].kind() == rt {
                    continue;
                }
                let oo = e.orientation();
                let Some(oe) = sf.generated(&e.0) else {
                    continue;
                };
                this.my_as_des.add(&of.0, &oe.oriented(oo));
            }
        } else {
            // `cxx:2138-2145`: a tube.
            for es in explorer_shapes(&of.0.oriented(Orientation::Forward), ShapeType::Edge) {
                this.my_as_des.add(&of.0, &es);
            }
        }
    }
}

/// `ToContext(MapSF)` (`cxx:2171-2340`): extend the offset faces of the
/// context faces (`myFaces`) to the shape and rebuild their edges.
pub fn to_context(this: &mut BRepOffsetMakeOffset, map_sf: &mut MapSf) {
    let mut created: ShapeMap = HashMap::new();
    let mut mef: ShapeMap = HashMap::new();
    // `NCollection_IndexedMap<TopoDS_Shape> FacesToBuild`.
    let mut faces_to_build: Vec<TopoShape> = Vec::new();
    let side = State::Out;

    // `cxx:2207-2240`: the edges and vertices the context faces reconstruct.
    for j in 1..=this.my_faces.extent() {
        let Some(cf) = this.my_faces.value(j).cloned() else {
            continue;
        };
        for es in explorer_shapes(&cf.oriented(Orientation::Forward), ShapeType::Edge) {
            let Some(e) = crate::shape::Edge::wrap(es) else {
                continue;
            };
            if !this.my_analyse.has_ancestor(&e.0) {
                continue;
            }
            let lea: Vec<TopoShape> = this.my_analyse.ancestors(&e.0).to_vec();
            for it in &lea {
                let Some(of) = map_sf.get(&shape_map_key(it)) else {
                    continue;
                };
                if !faces_to_build.iter().any(|s| shape_map_key(s) == shape_map_key(it)) {
                    faces_to_build.push(it.clone());
                }
                if let Some(g) = of.generated(&e.0) {
                    mef.insert(shape_map_key(&g), cf.clone());
                }
            }
            let (vv1, vv2) = edge_vertices(&e);
            for vi in [vv1, vv2].into_iter().flatten() {
                let lva: Vec<TopoShape> = this.my_analyse.ancestors(&vi.0).to_vec();
                for it in &lva {
                    let Some(ev) = crate::shape::Edge::wrap(it.clone()) else {
                        continue;
                    };
                    let Some(of) = map_sf.get(&shape_map_key(&ev.0)) else {
                        continue;
                    };
                    if !faces_to_build
                        .iter()
                        .any(|s| shape_map_key(s) == shape_map_key(&ev.0))
                    {
                        faces_to_build.push(ev.0.clone());
                    }
                    if let Some(g) = of.generated(&vi.0) {
                        mef.insert(shape_map_key(&g), cf.clone());
                    }
                }
            }
        }
    }

    // `cxx:2244-2250`.
    let rt = if this.my_offset < 0.0 {
        TypeOfConcavity::Convex
    } else {
        TypeOfConcavity::Concave
    };

    // `cxx:2254-2311`: rebuild the faces.
    // `Created` is a plain `DataMap<Shape, Shape>`, so its keys are recovered
    // from the offset faces: `ExtentFace` binds each rebuilt edge with the
    // edge of `F` it came from (`cxx:4156-4270`), so every key is an edge of
    // one of the `BOF.Face()` shapes below. This is a container detail, not a
    // predicate.
    let mut created_keys: Vec<TopoShape> = Vec::new();
    let mut seen_keys: HashSet<ShapeMapKey> = HashSet::new();
    for s in &faces_to_build {
        let Some(bof) = map_sf.get(&shape_map_key(s)).cloned() else {
            continue;
        };
        let Some(f) = bof.face().cloned() else {
            continue;
        };

        // `cxx:2259-2260`.
        let Some(nf) = BRepOffsetTool::extent_face(&f, &mut created, &mut mef, side, this.my_tol)
        else {
            // Port guard: no extended face, see the module notes.
            map_sf.remove(&shape_map_key(s));
            continue;
        };
        map_sf.remove(&shape_map_key(s));

        // `cxx:2313`: the keys of `Created`, see the note above.
        for es in explorer_shapes(&f.0.oriented(Orientation::Forward), ShapeType::Edge) {
            let k = shape_map_key(&es);
            if created.contains_key(&k) && seen_keys.insert(k) {
                created_keys.push(es);
            }
        }

        // `cxx:2264-2267`.
        let _ = this.my_init_offset_face.bind(s, &nf.0);
        this.my_init_offset_face.set_root(s);
        this.my_image_offset.set_root(&nf.0);

        if s.shape_type() == ShapeType::Face {
            for es in explorer_shapes(&s.oriented(Orientation::Forward), ShapeType::Edge) {
                let Some(e) = crate::shape::Edge::wrap(es) else {
                    continue;
                };
                let l: Vec<crate::brep_offset::Interval> = this.my_analyse.type_of(&e).to_vec();
                let Some(mut oe) = bof.generated(&e.0) else {
                    continue;
                };
                let or = e.orientation();
                oe.set_orientation(or);
                if l.is_empty() || l[0].kind() == rt {
                    continue;
                }
                match created.get(&shape_map_key(&oe)) {
                    Some(ne) => {
                        let mut ne = ne.clone();
                        let o = if ne.orientation() == Orientation::Reversed {
                            or.reversed()
                        } else {
                            or
                        };
                        ne.set_orientation(o);
                        this.my_as_des.add(&nf.0, &ne);
                    }
                    None => this.my_as_des.add(&nf.0, &oe),
                }
            }
        } else {
            // `cxx:2301-2309`: a tube.
            for es in explorer_shapes(&nf.0.oriented(Orientation::Forward), ShapeType::Edge) {
                this.my_as_des.add(&nf.0, &es);
            }
        }
    }

    // `cxx:2313-2337`: the free borders.
    for oe in &created_keys {
        if !this.my_init_offset_edge.is_image(oe) {
            continue;
        }
        let Ok(e) = this.my_init_offset_edge.image_from(oe) else {
            continue;
        };
        let Some(mut ne) = created.get(&shape_map_key(oe)).cloned() else {
            continue;
        };
        let first = this.my_init_offset_edge.image(&e);
        let Some(or) = first.first().map(|s| s.orientation()) else {
            continue;
        };
        let o = if ne.orientation() == Orientation::Reversed {
            or.reversed()
        } else {
            or
        };
        ne.set_orientation(o);
        let _ = this.my_init_offset_edge.remove(oe);
        let _ = this.my_init_offset_edge.bind(&e, &ne);
    }
}
