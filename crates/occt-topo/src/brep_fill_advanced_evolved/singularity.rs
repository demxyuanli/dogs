//! `RebuildFaces` (`BRepFill_AdvancedEvolved.cxx:1127-1155`),
//! `MakeEdgeDegenerated` (`cxx:1157-1214`), `InsertEDegenerated`
//! (`cxx:1216-1428`) and `CheckSingularityAndAdd` (`cxx:1430-1610`).
//!
//! `CheckSingularityAndAdd` splits a face when its own edges interfere (a
//! planar face is always checked, the other surfaces only when they can carry a
//! singularity: cone, sphere, Bezier / BSpline, surface of revolution), plugs
//! the singular vertex of a cone-like face with a degenerated edge
//! (`InsertEDegenerated`), and rebuilds the face from the split edges
//! (`RebuildFaces`).
//!
//! Translation notes:
//! - `BRep_Builder::UpdateEdge(anEdegen, aL, theFace, Confusion())` stores a
//!   `BRep_CurveOnSurface` whose location is identity; the port keeps the
//!   pcurve in [`GeometryRegistry`] under the face's key
//!   (`set_edge_pcurve`), which is what `BRep_Tool::CurveOnSurface(E, F)`
//!   reads back.
//! - `BRepAdaptor_Surface::UResolution` / `VResolution` (`cxx:1170`) are
//!   [`SurfaceTool::u_resolution`] / [`SurfaceTool::v_resolution`]
//!   (`BRepApprox_SurfaceTool.lxx:206-214`).
//! - `BRepAdaptor_Surface(theFace, false)` is
//!   [`SurfaceAdaptor::of_face`](crate::topopebrep::tool::adaptor::SurfaceAdaptor::of_face)
//!   (`BRepAdaptor_Surface.cxx:77-80`).
//! - `TopoDS_Iterator(aFace)` is
//!   [`ShapeIterator::of_shape`](crate::iterator::ShapeIterator::of_shape).
//! - `CurveOnSurface(E, F)` answering a null handle (an edge without a pcurve
//!   on that face) makes the reference dereference null; the port skips that
//!   edge (`anE1 = anE2`) instead.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::{Orientation, ShapeType};
use crate::bop_builder2::BopBuilder;
use crate::bopalgo_tools_wires::{edges_to_wires, wires_to_faces};
use crate::brep_approx::surface_tool::SurfaceTool;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::builder_area::AreaBuilder;
use crate::builder_face::FaceBuilder;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::iterator::ShapeIterator;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same, map_shapes_and_unique_ancestors, shape_map_key,
    ShapeMapKey,
};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::AdvancedEvolved;

/// `BOPAlgo_Tools::EdgesToWires(theEdges, theWires, theShared, theAngTol)`'s
/// and `WiresToFaces`'s default `theAngTol` (`BOPAlgo_Tools.hxx:154`, `:171`).
const TOOLS_ANG_TOL: f64 = 1.0e-8;

/// `BRep_Tool::CurveOnSurface(E, F, f, l)`: the pcurve of `edge` on `face`.
fn pcurve_on_face(edge: &TopoShape, face: &TopoShape) -> Option<Arc<dyn Curve2d>> {
    let face_key = GeometryRegistry::shape_key(face);
    GeometryRegistry::global().edge_pcurve(edge, face_key)
}

/// `BRep_Tool::Parameter(V, E)` (`BRep_Tool.cxx:1542-1580`).
fn parameter_of_vertex(edge: &TopoShape, vertex: &TopoShape) -> f64 {
    BRepTool::parameter_of_vertex(edge, vertex)
}

/// `RebuildFaces(theLE, theSourceFace, theList)` (`cxx:1127-1155`).
fn rebuild_faces(the_le: &[TopoShape], the_source_face: &TopoShape, the_list: &mut Vec<TopoShape>) {
    // `cxx:1130-1141`: `BOPAlgo_BuilderFace` over the edges, on the FORWARD
    // source face.
    let mut a_bf = FaceBuilder::new();
    let a_f = the_source_face.oriented(Orientation::Forward);
    a_bf.set_face(&Face(a_f));
    a_bf.set_shapes(the_le);
    let _ = a_bf.perform();

    // `cxx:1143-1153`.
    let a_lfr = a_bf.areas().to_vec();
    if a_lfr.is_empty() {
        the_list.push(the_source_face.clone());
        return;
    }
    for f in a_lfr {
        the_list.push(f);
    }
}

/// `MakeEdgeDegenerated(theV, theFace, thePf, thePl, theLEdges)`
/// (`cxx:1157-1214`): returns `true` when a degenerated edge has been created
/// and appends it to `the_edges` **twice** (the caller needs both orientations
/// to feed a loop).
fn make_edge_degenerated(
    the_v: &TopoShape,
    the_face: &TopoShape,
    the_pf: &GpPnt2d,
    the_pl: &GpPnt2d,
    the_edges: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:1169-1176`.
    let an_as = SurfaceAdaptor::of_face(&Face(the_face.clone()));
    let a_tol = 2.0 * BRepTool::vertex_tolerance(&Vertex(the_v.clone()));
    let a_tol_u = SurfaceTool::u_resolution(&an_as, a_tol);
    let a_tol_v = SurfaceTool::v_resolution(&an_as, a_tol);

    if (the_pf.x() - the_pl.x()).abs() < a_tol_u && (the_pf.y() - the_pl.y()).abs() < a_tol_v {
        return false;
    }

    // `cxx:1179-1184`.
    let a_vf = the_v.oriented(Orientation::Forward);
    let a_vl = the_v.oriented(Orientation::Reversed);

    let a_v = GpVec2d::new(the_pl.x() - the_pf.x(), the_pl.y() - the_pf.y());
    let Ok(a_dir1) = GpDir2d::from_vec2d(&a_v) else {
        // `gp_Dir2d(aV)` raises on a null vector.
        return false;
    };
    let a_l1 = Geom2dLine::from_pnt_dir(*the_pf, a_dir1);
    let Ok(a_dir2) = GpDir2d::from_vec2d(&a_v.reversed()) else {
        return false;
    };
    let a_l2 = Geom2dLine::from_pnt_dir(*the_pl, a_dir2);

    // `cxx:1187-1212`: `MakeEdge` + `UpdateEdge(edge, aL, theFace, Confusion())`
    // + `Degenerated(edge, true)` + `Range(edge, 0, aLPar)` + the two vertices.
    let builder = TopoBuilder::new();
    let registry = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(the_face);
    let a_l_par = a_v.magnitude();
    let a_point = BRepTool::vertex_point(&Vertex(the_v.clone()));

    let mut an_edegen1 = builder.make_degenerated_edge(&a_point);
    registry.set_edge_range(&an_edegen1.0, 0.0, a_l_par);
    registry.set_edge_pcurve(&an_edegen1.0, face_key, Arc::new(a_l1));
    builder.add_edge_vertices(&mut an_edegen1, &Vertex(a_vf.clone()), &Vertex(a_vl.clone()));

    let mut an_edegen2 = builder.make_degenerated_edge(&a_point);
    registry.set_edge_range(&an_edegen2.0, 0.0, a_l_par);
    registry.set_edge_pcurve(&an_edegen2.0, face_key, Arc::new(a_l2));
    builder.add_edge_vertices(&mut an_edegen2, &Vertex(a_vf), &Vertex(a_vl));

    // `cxx:1209-1211`.
    the_edges.push(an_edegen1.0);
    the_edges.push(an_edegen2.0);

    true
}

/// `InsertEDegenerated(theFace, theLEdges)` (`cxx:1216-1428`).
fn insert_e_degenerated(the_face: &TopoShape, the_edges: &mut Vec<TopoShape>) {
    let builder = TopoBuilder::new();

    // `cxx:1220-1230`: the edges go into one wired loop.
    let mut a_wir = builder.make_wire(&[]);
    for e in the_edges.iter() {
        builder.add_edge(&mut a_wir, &Edge(e.clone()));
    }

    // `cxx:1232-1234`.
    let a_map_ve = map_shapes_and_unique_ancestors(
        &a_wir.0,
        ShapeType::Vertex,
        ShapeType::Edge,
        false,
    );

    // `cxx:1236-1247`: the explorer positions on the first edge; an empty walk
    // (only INTERNAL / EXTERNAL edges) returns.
    let mut an_exp = WireExplorer::of_wire_on_face(&a_wir, &Face(the_face.clone()));
    let Some(a_first_edge) = an_exp.current_edge().map(|e| e.0) else {
        return;
    };

    let mut an_e1 = a_first_edge.clone();
    let mut a_last_edge: Option<TopoShape> = None;
    an_exp.next();

    // `cxx:1249-1271`: the vertices of every degenerated edge of the wire.
    let mut a_map_v_of_de: Vec<ShapeMapKey> = Vec::new();
    for e in explorer_shapes(&a_wir.0, ShapeType::Edge) {
        if !BRepTool::is_degenerated(&Edge(e.clone())) {
            continue;
        }
        if let (Some(v1), _) = edge_vertices(&Edge(e)) {
            let key = shape_map_key(&v1.0);
            if !a_map_v_of_de.contains(&key) {
                a_map_v_of_de.push(key);
            }
        }
    }

    // `cxx:1272-1358`.
    while an_exp.more() {
        let Some(an_e2) = an_exp.current_edge().map(|e| e.0) else {
            break;
        };
        a_last_edge = Some(an_e2.clone());

        let Some(a_vert_curr) = an_exp.current_vertex().map(|v| v.0.clone()) else {
            break;
        };

        // `cxx:1279-1285`: "necessary degenerated edge has already been
        // created".
        let a_vert_key = shape_map_key(&a_vert_curr);
        if a_map_v_of_de.contains(&a_vert_key) {
            an_e1 = an_e2;
            an_exp.next();
            continue;
        }

        // `cxx:1287-1297`.
        let (a_c1, a_c2) = match (
            pcurve_on_face(&an_e1, the_face),
            pcurve_on_face(&an_e2, the_face),
        ) {
            (Some(c1), Some(c2)) => (c1, c2),
            // The reference dereferences the null handle here; the port skips
            // the pair instead.
            _ => {
                an_e1 = an_e2;
                an_exp.next();
                continue;
            }
        };
        let a_f = parameter_of_vertex(&an_e1, &a_vert_curr);
        let a_l = parameter_of_vertex(&an_e2, &a_vert_curr);
        let a_pf = a_c1.d0(a_f);
        let a_pl = a_c2.d0(a_l);

        // `cxx:1298-1303`.
        if make_edge_degenerated(&a_vert_curr, the_face, &a_pf, &a_pl, the_edges) {
            a_map_v_of_de.push(a_vert_key);
            an_e1 = an_e2;
            an_exp.next();
            continue;
        }

        // `cxx:1300-1308`: `aMapVE.Seek(aVertCurr)`.
        let an_e_list = a_map_ve.get(&a_vert_key).cloned().unwrap_or_default();
        if an_e_list.len() <= 2 {
            an_e1 = an_e2;
            an_exp.next();
            continue;
        }

        // `cxx:1310-1327`: "case like cone with apex".
        let mut has_degenerated = false;
        for edge in &an_e_list {
            if BRepTool::is_degenerated(&Edge(edge.clone())) {
                has_degenerated = true;
                break;
            }
        }
        if has_degenerated {
            an_e1 = an_e2;
            an_exp.next();
            continue;
        }

        // `cxx:1329-1356`: "look for the pair for anE1 and anE2 edges".
        let mut pair_found = false;
        for i in 0..2 {
            let a_point = if i == 1 { a_pl } else { a_pf };
            for edge in &an_e_list {
                if is_same(edge, &an_e1) || is_same(edge, &an_e2) {
                    continue;
                }
                let Some(a_c) = pcurve_on_face(edge, the_face) else {
                    continue;
                };
                let a_f = parameter_of_vertex(edge, &a_vert_curr);
                let a_p = a_c.d0(a_f);
                if make_edge_degenerated(&a_vert_curr, the_face, &a_point, &a_p, the_edges) {
                    a_map_v_of_de.push(a_vert_key.clone());
                    pair_found = true;
                    break;
                }
            }
            if pair_found {
                break;
            }
        }

        an_e1 = an_e2;
        an_exp.next();
    }

    // `cxx:1360-1363`.
    let Some(a_last_edge) = a_last_edge else {
        return;
    };

    // `cxx:1365-1383`: `TopExp::Vertices(aFirstEdge, aV[0], aV[1])` and the
    // "not-closed edge covered by a single vertex" bail-outs.
    let (a_v0, a_v1) = edge_vertices(&Edge(a_first_edge.clone()));
    match (&a_v0, &a_v1) {
        (Some(v0), Some(v1)) if is_same(&v0.0, &v1.0) => return,
        _ => {}
    }
    let (a_v2, a_v3) = edge_vertices(&Edge(a_last_edge.clone()));
    match (&a_v2, &a_v3) {
        (Some(v2), Some(v3)) if is_same(&v2.0, &v3.0) => return,
        _ => {}
    }

    let a_v = [a_v0, a_v1, a_v2, a_v3];

    // `cxx:1386-1427`.
    for an_id_fe in 0..2 {
        for an_id_le in 2..4 {
            let (Some(a_v_fe), Some(a_v_le)) = (&a_v[an_id_fe], &a_v[an_id_le]) else {
                continue;
            };
            if !is_same(&a_v_fe.0, &a_v_le.0) {
                continue;
            }

            let a_vert_key = shape_map_key(&a_v_fe.0);
            let an_e_list = a_map_ve.get(&a_vert_key).cloned().unwrap_or_default();
            if an_e_list.len() > 2 {
                // Non-manifold vertex, or the seam case of `cxx:1397-1416`.
                continue;
            }

            let (Some(a_c1), Some(a_c2)) = (
                pcurve_on_face(&a_first_edge, the_face),
                pcurve_on_face(&a_last_edge, the_face),
            ) else {
                continue;
            };
            let a_f = parameter_of_vertex(&a_first_edge, &a_v_fe.0);
            let a_l = parameter_of_vertex(&a_last_edge, &a_v_le.0);
            let a_pf = a_c1.d0(a_f);
            let a_pl = a_c2.d0(a_l);

            make_edge_degenerated(&a_v_fe.0, the_face, &a_pf, &a_pl, the_edges);
        }
    }
}

impl AdvancedEvolved {
    /// `CheckSingularityAndAdd(theF, theFuzzyToler, theListOfFaces,
    /// theListOfSplits)` (`cxx:1430-1610`): answers whether `the_f` has been
    /// split.
    pub(super) fn check_singularity_and_add(
        &self,
        the_f: &TopoShape,
        the_fuzzy_toler: f64,
        the_list_of_faces: &mut Vec<TopoShape>,
        the_list_of_splits: &mut Vec<TopoShape>,
    ) -> bool {
        // `cxx:1440-1444`.
        let an_as = SurfaceAdaptor::of_face(&Face(the_f.clone()));
        let mut a_s_type = an_as.get_type();
        if a_s_type == GeomAbsSurfaceType::OffsetSurface {
            if let Some(basis) = an_as.surface().and_then(|s| s.offset_basis_surface()) {
                a_s_type = surface_kind(basis.as_ref());
            }
        }

        // `cxx:1447-1510`: a planar face always gets its edges checked for
        // interferences.
        if a_s_type == GeomAbsSurfaceType::Plane {
            let mut a_me: Vec<ShapeMapKey> = Vec::new();
            let mut a_le: Vec<TopoShape> = Vec::new();
            for e in explorer_shapes(the_f, ShapeType::Edge) {
                let key = shape_map_key(&e);
                if !a_me.contains(&key) {
                    a_me.push(key);
                    a_le.push(e);
                }
            }

            // `cxx:1463-1473`: "split interfered edges".
            let mut a_pf = PaveFiller::new();
            a_pf.set_arguments(&a_le);
            let _ = a_pf.perform();
            if a_pf.has_errors() {
                the_list_of_faces.push(the_f.clone());
                return false;
            }

            // `cxx:1475-1482`: "interfered edges have not been detected".
            if a_pf.ds().nb_shapes() == a_pf.ds().nb_source_shapes() {
                the_list_of_faces.push(the_f.clone());
                return false;
            }

            // `cxx:1484-1492`.
            let mut a_builder = BopBuilder::new();
            a_builder.set_arguments(&a_le);
            a_builder.set_filler(a_pf);
            let Ok(an_edges) = a_builder.perform_with_filler() else {
                the_list_of_faces.push(the_f.clone());
                return false;
            };
            if a_builder.has_errors() {
                the_list_of_faces.push(the_f.clone());
                return false;
            }

            // `cxx:1494-1508`.
            let (_nb, a_comp_w) = edges_to_wires(&an_edges, true, TOOLS_ANG_TOL);
            let (_ok, a_comp_f) = wires_to_faces(&a_comp_w, TOOLS_ANG_TOL);
            for f in explorer_shapes(&a_comp_f, ShapeType::Face) {
                the_list_of_splits.push(f);
            }
            return true;
        }

        // `cxx:1512-1518`: the surfaces that can carry a singularity.
        if a_s_type != GeomAbsSurfaceType::Cone
            && a_s_type != GeomAbsSurfaceType::Sphere
            && a_s_type != GeomAbsSurfaceType::BezierSurface
            && a_s_type != GeomAbsSurfaceType::BSplineSurface
            && a_s_type != GeomAbsSurfaceType::SurfaceOfRevolution
        {
            the_list_of_faces.push(the_f.clone());
            return false;
        }

        // `cxx:1522-1528`.
        let mut is_split = false;
        let mut a_list_edges: Vec<TopoShape> = Vec::new();

        let a_face = the_f.oriented(Orientation::Forward);

        // `cxx:1530-1598`.
        for wire in ShapeIterator::of_shape(&a_face) {
            let a_wir = Wire(wire);

            // `cxx:1535-1540`.
            let a_lgf: Vec<TopoShape> = crate::topo_tools::edges_of_wire(&a_wir)
                .into_iter()
                .map(|e| e.0)
                .collect();

            // `cxx:1542-1550`.
            let mut a_pf = PaveFiller::new();
            a_pf.set_arguments(&a_lgf);
            a_pf.set_fuzzy_value(the_fuzzy_toler);
            // `cxx:1563`: `aBuilder.SetNonDestructive(true)` -- the port keeps
            // the option on the filler the builder reads it from
            // (`BopBuilder::non_destructive`).
            a_pf.set_non_destructive(true);
            let _ = a_pf.perform();
            if a_pf.has_errors() {
                continue;
            }

            // `cxx:1552-1558`: "no new shapes have been created".
            if a_pf.ds().nb_shapes() == a_pf.ds().nb_source_shapes() {
                continue;
            }

            // `cxx:1560-1567`.
            let mut a_builder = BopBuilder::new();
            a_builder.set_arguments(&a_lgf);
            a_builder.set_non_destructive(true);
            a_builder.set_filler(a_pf);
            if a_builder.perform_with_filler().is_err() || a_builder.has_errors() {
                continue;
            }

            // `cxx:1569-1588`: keep the original edge when it was not modified
            // (or is degenerated), else its images.
            let mut a_le: Vec<TopoShape> = Vec::new();
            for a_sh in &a_lgf {
                let a_lm = a_builder.history().modified(a_sh);
                if a_lm.map_or(true, |l| l.is_empty())
                    || BRepTool::is_degenerated(&Edge(a_sh.clone()))
                {
                    a_le.push(a_sh.clone());
                    continue;
                }
                for em in a_lm.unwrap_or(&[]) {
                    a_le.push(em.clone());
                }
            }

            // `cxx:1589-1592`.
            is_split = true;
            insert_e_degenerated(&a_face, &mut a_le);
            a_list_edges.extend(a_le);
        }

        // `cxx:1593-1597`.
        if !is_split {
            the_list_of_faces.push(the_f.clone());
            return false;
        }

        // `cxx:1600-1608`.
        rebuild_faces(&a_list_edges, the_f, the_list_of_splits);
        let a_splits = the_list_of_splits.clone();
        for f in &a_splits {
            the_list_of_faces.push(f.oriented(the_f.orientation()));
        }

        true
    }
}
