//! `BRepFill_Evolved::PrepareSpine` (`BRepFill_Evolved.cxx:1716-1779`) with
//! its statics `CutEdge` (`cxx:3277-3380`) and `EdgeVertices` (`cxx:172-186`).
//!
//! `PrepareSpine` rebuilds the spine as a face over the *same* surface and
//! location whose wires are the input wires with each edge cut at the extrema
//! of curvature and at the inflexion points (`MAT2d_CutCurve`). A closed
//! circular edge is instead split into three arcs (`cxx:3293-3335`). The
//! resulting 3D curves are built afterwards by `BRepLib::BuildCurves3d`.
//!
//! `myMapSpine` records which original edge — and which original vertices —
//! each piece came from, exactly like `BRepFill_OffsetWire::PrepareSpine`.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_lib_build_curves3d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::mat2d::Mat2dCutCurve;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, edges_of, wires_of_face};

use super::shape_map::ShapeMap;
use super::side::edge_vertices_oriented;

/// `BRep_Builder::Add` on a freshly `E.EmptyCopied()` edge: carry the curve
/// representations (3D curve and pcurve slots) over to the copy, which is what
/// `B.Range` then re-trims (`BRepFill_OffsetWire`'s spine cut does the same).
fn copied_edge(e: &Edge) -> TopoShape {
    let registry = GeometryRegistry::global();
    let ne = e.0.empty_copied();
    if let Some(geom) = registry.edge_geom(e) {
        registry.set_edge(&ne, geom);
    }
    registry.copy_edge_pcurve_slots(&ne, &e.0);
    ne
}

/// `BRep_Builder::MakeVertex(V)` for `BRepLib_MakeVertex(P)` (`cxx:3302-3303`,
/// `cxx:3366-3367`): the 2D point is lifted onto the `z = 0` plane.
fn make_vertex_at(p: occt_core::gp::GpPnt2d) -> Vertex {
    // `BRepLib_MakeVertex` leaves the `BRep_TVertex` tolerance at its
    // constructor value (`RealEpsilon`, `BRep_TVertex.cxx:27`).
    TopoBuilder::new().make_vertex(
        GpPnt::new(p.x(), p.y(), 0.0),
        occt_core::precision::REAL_EPSILON,
    )
}

/// `CutEdge(E, F, Cuts)` (`BRepFill_Evolved.cxx:3277-3380`).
///
/// `Cuts` stays empty when the edge is not modified.
fn cut_edge(e: &Edge, f: &Face) -> Vec<TopoShape> {
    let mut cuts: Vec<TopoShape> = Vec::new();
    let builder = TopoBuilder::new();
    let registry = GeometryRegistry::global();

    // `cxx:3287-3288`: `TopExp::Vertices(E, V1, V2)`.
    let (v1, v2) = edge_vertices(e);

    // `cxx:3289-3290`: `C2d = BRep_Tool::CurveOnSurface(E, F, f, l)`,
    // `CT2d = new Geom2d_TrimmedCurve(C2d, f, l)`.
    let Some((c2d, f_par, l_par)) = curve_on_surface_range(e, f) else {
        return cuts;
    };
    let ct2d = Geom2dTrimmedCurve::new(c2d, f_par, l_par);

    // `cxx:3292-3336`: cut a closed circle into three arcs.
    if ct2d.basis_curve().is_geom2d_circle() && BRepTool::is_closed(&e.0) {
        let m1 = (2.0 * f_par + l_par) / 3.0;
        let m2 = (f_par + 2.0 * l_par) / 3.0;
        let p1 = ct2d.d0(m1);
        let p2 = ct2d.d0(m2);

        let vl1 = make_vertex_at(p1);
        let vl2 = make_vertex_at(p2);

        let (Some(v1), Some(v2)) = (v1, v2) else {
            // `TopoDS::Vertex` of a null handle raises in OCCT.
            return cuts;
        };

        // `cxx:3306-3318`.
        let mut fe = copied_edge(e);
        let mut me = copied_edge(e);
        let mut le = copied_edge(e);
        fe.set_orientation(Orientation::Forward);
        me.set_orientation(Orientation::Forward);
        le.set_orientation(Orientation::Forward);

        // `cxx:3319-3328`.
        builder.add(&mut fe, &v1.0.oriented(Orientation::Forward));
        builder.add(&mut fe, &vl1.0.oriented(Orientation::Reversed));
        registry.set_ranges_all(&fe, f_par, m1);

        builder.add(&mut me, &vl1.0.oriented(Orientation::Forward));
        builder.add(&mut me, &vl2.0.oriented(Orientation::Reversed));
        registry.set_ranges_all(&me, m1, m2);

        builder.add(&mut le, &vl2.0.oriented(Orientation::Forward));
        builder.add(&mut le, &v2.0.oriented(Orientation::Forward));
        registry.set_ranges_all(&le, m2, l_par);

        // `cxx:3330-3332`.
        for piece in [fe, me, le] {
            cuts.push(piece.oriented(e.0.orientation()));
        }
        return cuts;
    }

    // `cxx:3343-3345`: `Cuter.Perform(CT2d)`.
    let mut cuter = Mat2dCutCurve::new();
    cuter.perform(&ct2d);

    if cuter.un_modified() {
        // `cxx:3347-3351`: edge not modified.
        return cuts;
    }

    // `cxx:3355-3379`: one edge per piece; only the last one keeps the
    // original end vertex.
    let Some(mut vf) = v1 else {
        return cuts;
    };
    let Some(v2) = v2 else {
        return cuts;
    };
    let nb = cuter.nb_curves();
    for k in 1..=nb {
        let cc = cuter.value(k as i32);
        let vl = if k == nb {
            v2.clone()
        } else {
            make_vertex_at(cc.d0(cc.last_parameter()))
        };

        let mut ne = copied_edge(e);
        ne.set_orientation(Orientation::Forward);
        builder.add(&mut ne, &vf.0.oriented(Orientation::Forward));
        builder.add(&mut ne, &vl.0.oriented(Orientation::Reversed));
        registry.set_ranges_all(&ne, cc.first_parameter(), cc.last_parameter());
        cuts.push(ne.oriented(e.0.orientation()));

        vf = vl;
    }
    cuts
}

/// `BRepFill_Evolved::PrepareSpine(WorkSpine, MapSpine)`
/// (`BRepFill_Evolved.cxx:1716-1779`).
pub fn prepare_spine(
    spine: &Face,
    work_spine: &mut Face,
    map_spine: &mut ShapeMap<TopoShape>,
) {
    let builder = TopoBuilder::new();

    // `cxx:1726-1728`: the face keeps the spine's surface, location and
    // tolerance.
    let Some(surface) = BRepTool::face_surface(spine) else {
        return;
    };
    let tol_f = BRepTool::face_tolerance(spine);
    let mut work = builder.make_face_with_tolerance(surface, &[], tol_f);
    work.0.set_location(spine.0.location());

    // `cxx:1730`: `TopoDS_Iterator IteF(mySpine)`.
    for wire in wires_of_face(spine) {
        // `cxx:1734-1735`.
        let mut new_wire = builder.make_wire(&[]);
        let is_closed = wire.0.is_closed();

        // `cxx:1737`: `TopoDS_Iterator IteW(IteF.Value())`.
        for edge in edges_of(&wire.0) {
            // `cxx:1741-1744`: bind both vertices to themselves.
            let (v1, v2) = edge_vertices_oriented(&edge);
            if let Some(v) = &v1 {
                map_spine.bind(&v.0, v.0.clone());
            }
            if let Some(v) = &v2 {
                map_spine.bind(&v.0, v.0.clone());
            }

            // `cxx:1748`: `CutEdge(E, mySpine, Cuts)`.
            let cuts = cut_edge(&edge, spine);

            if cuts.is_empty() {
                // `cxx:1749-1753`.
                builder.add(&mut new_wire.0, &edge.0);
                map_spine.bind(&edge.0, edge.0.clone());
            } else {
                // `cxx:1755-1771`.
                for cut in cuts {
                    map_spine.bind(&cut, edge.0.clone());
                    builder.add(&mut new_wire.0, &cut);
                    let (c1, c2) = edge_vertices(&Edge(cut));
                    if let Some(v) = c1 {
                        if !map_spine.is_bound(&v.0) {
                            map_spine.bind(&v.0, edge.0.clone());
                        }
                    }
                    if let Some(v) = c2 {
                        if !map_spine.is_bound(&v.0) {
                            map_spine.bind(&v.0, edge.0.clone());
                        }
                    }
                }
            }
        }

        // `cxx:1774-1776`.
        new_wire.0.set_closed(is_closed);
        builder.add_wire(&mut work, &new_wire);
    }

    // `cxx:1778-1779`: `BRepLib::BuildCurves3d(WorkSpine)`.
    brep_lib_build_curves3d::build_curves3d(&work.0);

    *work_spine = work;
}

/// The 2D curve of an edge on a face, for the callers that need the same
/// `Geom2dAdaptor_Curve` OCCT's `CutEdge` works with.
pub fn edge_pcurve(edge: &Edge, face: &Face) -> Option<Arc<dyn Curve2d>> {
    curve_on_surface_range(edge, face).map(|(c, _, _)| c)
}
