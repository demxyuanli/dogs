//! `Build` phase 3: `SectionsEdges` -> `FinalWire`.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MiddlePath.cxx`
//! (OCCT 8.0.0): section faces `cxx:845-890`, middle edges `cxx:892-1017`,
//! the interpolation of the missed edges `cxx:1020-1103` and the final wire
//! `cxx:1105-1147`.
//!
//! Faithfully kept:
//!
//! * every section is closed with `BRepLib_MakeEdge(V2, V1)` on an open
//!   section (`cxx:858-863`), then turned into a face with
//!   `BRepLib_MakeFace(aWire, true)` ("Only plane"); when the face build fails
//!   the wire itself is the section shape (`cxx:866-872`);
//! * the centre of mass is `BRepGProp::SurfaceProperties` for a face and
//!   `BRepGProp::LinearProperties` for a wire (`cxx:878-889`);
//! * a middle edge is a straight line between the centres when every path
//!   segment of that section is a line (`cxx:921-924`), a circle when they are
//!   all circles with a common axis and the same arc (`cxx:925-1016`), and
//!   otherwise leaves `MidEdges(i)` null so the interpolation pass fills it;
//! * the interpolation pass averages the path tangents through
//!   `GeomLib::Inertia` (`cxx:1041-1089`) and runs `GeomAPI_Interpolate` with
//!   those tangents and `theFlags` all true (`cxx:1091-1100`).
//!
//! Arms that raise in OCCT, and answer `Err` here:
//!
//! * `TopoDS::Edge` / `TopoDS::Vertex` on a null or mistyped shape;
//! * `BRepLib_MakeEdge(P1, P2)` and `BRepLib_MakeEdge(V1, V2)` with coincident
//!   points (`BRepLib_LineThroughIdenticPoints`);
//! * `gp_Dir(gp_Vec)` on a null vector and `gp_Ax2(P, N, Vx)` with a `Vx` that
//!   is not normal to `N` (`Standard_ConstructionError`), both of which the
//!   `GC_MakeCircle` arm of `cxx:1001-1011` would raise;
//! * `GeomAPI_Interpolate::Curve()` on a not-done interpolation
//!   (`StdFail_NotDone`), and `BRepLib_MakeWire::Wire()` on a wire that never
//!   reached `Done` (`StdFail_NotDone`).
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * `std::cout << "Interpolation failed"` (`cxx:1096-1098`).
//! * The commented-out "Temporary" compound dump (`cxx:1117-1145`) and the
//!   commented-out `VecProd` axis fix-up (`cxx:986-992`).
//! * `FinalWire` (`cxx:1115`) is only read through `myShape` (`cxx:1116`), so
//!   the port stores the single wire.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx2, GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::REAL_SMALL;
use occt_geom::geom_lib::geom_lib_inertia;
use occt_geom::{Curve, GeomCircle};

use crate::abs::{Orientation, ShapeType};
use crate::brep_gprop_full::{linear_properties, surface_properties};
use crate::brep_lib_make_face::make_face_from_wire;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_offset_api_middle_path::statics::{
    make_edge_over_curve, make_edge_range, tangent_of_edge, type_of_edge,
};
use crate::brep_offset_api_middle_path::{BRepOffsetApiMiddlePath, Ctx};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::wire_vertices;

/// `BRepLib_MakeEdge(V1, V2)` (`BRepLib_MakeEdge.cxx:185-198`): `None` is
/// `BRepLib_LineThroughIdenticPoints`.
fn make_edge_vv(v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    if p1.distance(&p2) <= REAL_SMALL {
        return None;
    }
    Some(TopoBuilder::new().make_edge_segment_with_vertices(&p1, &p2, v1, v2))
}

/// `BRepLib_MakeEdge(P1, P2)` (`BRepLib_MakeEdge.cxx:201-215`): `None` is
/// `BRepLib_LineThroughIdenticPoints`.
fn make_edge_pp(p1: &GpPnt, p2: &GpPnt) -> Option<Edge> {
    if p1.distance(p2) <= REAL_SMALL {
        return None;
    }
    Some(TopoBuilder::new().make_edge_segment(p1, p2))
}

/// `myPaths(indp)(k)` with the reference's 1-based indices.
fn path_at<'a>(
    mp: &'a BRepOffsetApiMiddlePath,
    indp: usize,
    k: usize,
) -> Result<&'a TopoShape, String> {
    mp.my_paths
        .get(indp - 1)
        .and_then(|p| p.get(k - 1))
        .ok_or_else(|| "Standard_OutOfRange: myPaths".to_string())
}

/// `TopoDS::Edge(aShape)`.
fn as_edge(s: &TopoShape) -> Result<Edge, String> {
    Edge::wrap(s.clone()).ok_or_else(|| "Standard_TypeMismatch: TopoDS::Edge".to_string())
}

/// `BRep_Tool::Curve(anEdge, fpar, lpar)` with the `Geom_TrimmedCurve` unwrap
/// of `cxx:938-941`, plus `BRep_Tool::Curve` of the edge.
fn edge_curve_untrimmed(an_edge: &Edge) -> Result<(Arc<dyn Curve>, f64, f64), String> {
    let curve = BRepTool::edge_curve(an_edge)
        .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
    let (fpar, lpar) = BRepTool::edge_parameters(an_edge);
    let curve = match curve.untrimmed_basis() {
        Some((basis, _, _)) => basis,
        None => curve,
    };
    Ok((curve, fpar, lpar))
}

/// `cxx:925-1016`: `TypeOfMidEdge == GeomAbs_Circle`.
fn mid_edge_circle(
    mp: &BRepOffsetApiMiddlePath,
    centers: &[GpPnt],
    i: usize,
    lin_tol: f64,
    ang_tol: f64,
) -> Result<Option<Edge>, String> {
    // `cxx:926-931`.
    let mut the_axis = GpAx1::new(GpPnt::zero(), GpDir::default_dir());
    let mut the_dir1 = GpDir::default_dir();
    let mut the_dir2 = GpDir::default_dir();
    let mut the_angle = 0.0f64;
    let mut the_tangent = GpVec::zero();
    let mut similar_arcs = true;

    // `cxx:934-977`.
    for j in 1..=mp.my_paths.len() {
        let shape = path_at(mp, j, i)?;
        let an_edge = as_edge(shape)?;
        let (a_curve, fpar, lpar) = edge_curve_untrimmed(&an_edge)?;
        // `cxx:941-943`.
        let mut pnt1 = a_curve.value(fpar);
        let pnt2 = a_curve.value(lpar);
        let a_circ = a_curve
            .gp_circ()
            .ok_or_else(|| "Standard_TypeMismatch: Geom_Circle".to_string())?;

        if j == 1 {
            // `cxx:946-958`.
            the_axis = a_circ.axis();
            the_dir1 = GpDir::from_vec(&GpVec::from_pnts(&a_circ.location(), &pnt1))
                .map_err(|e| format!("Standard_ConstructionError: gp_Dir({e})"))?;
            the_dir2 = GpDir::from_vec(&GpVec::from_pnts(&a_circ.location(), &pnt2))
                .map_err(|e| format!("Standard_ConstructionError: gp_Dir({e})"))?;
            the_angle = lpar - fpar;
            let the_param = if an_edge.0.orientation() == Orientation::Forward {
                fpar
            } else {
                lpar
            };
            let (p, mut tangent) = a_curve.d1(the_param);
            pnt1 = p;
            if an_edge.0.orientation() == Orientation::Reversed {
                tangent.reverse();
            }
            the_tangent = tangent;
        } else {
            // `cxx:959-977`.
            let an_axis = a_circ.axis();
            let a_lin = GpLin::new(an_axis);
            if !a_lin.contains(the_axis.location(), lin_tol)
                || !an_axis
                    .direction()
                    .is_parallel_tol(the_axis.direction(), ang_tol)
            {
                similar_arcs = false;
                break;
            }
            let a_dir1 = GpDir::from_vec(&GpVec::from_pnts(&a_circ.location(), &pnt1))
                .map_err(|e| format!("Standard_ConstructionError: gp_Dir({e})"))?;
            let a_dir2 = GpDir::from_vec(&GpVec::from_pnts(&a_circ.location(), &pnt2))
                .map_err(|e| format!("Standard_ConstructionError: gp_Dir({e})"))?;
            if !((a_dir1.is_equal_tol(&the_dir1, ang_tol)
                && a_dir2.is_equal_tol(&the_dir2, ang_tol))
                || (a_dir1.is_equal_tol(&the_dir2, ang_tol)
                    && a_dir2.is_equal_tol(&the_dir1, ang_tol)))
            {
                similar_arcs = false;
                break;
            }
        }
    }

    if !similar_arcs {
        // `cxx:979`: `MidEdges(i)` keeps its null value.
        return Ok(None);
    }

    // `cxx:980-984`.
    let axis_loc = *the_axis.location().xyz();
    let axis_dir = *the_axis.direction().xyz();
    let parameter = centers[i - 1].xyz().subtracted(&axis_loc).dot(&axis_dir);
    let the_center_of_circ = GpPnt::from_xyz(&axis_loc.added(&axis_dir.multiplied(parameter)));

    // `cxx:986-993`.
    let vec1 = GpVec::from_pnts(&the_center_of_circ, &centers[i - 1]);
    let vec2 = GpVec::from_pnts(&the_center_of_circ, &centers[i]);

    let mut an_angle = vec1.angle_with_ref(&vec2, &GpVec::from_xyz(the_axis.direction().xyz()));
    if an_angle < 0.0 {
        an_angle += 2.0 * std::f64::consts::PI;
    }
    // `cxx:996-1001`.
    if (an_angle - the_angle).abs() > ang_tol {
        the_axis.reverse();
    }

    // `cxx:1002-1011`.
    let mut ax2 = make_ax2(&the_center_of_circ, &the_axis, &vec1)?;
    let mut the_circle = make_circle(&ax2, vec1.magnitude());
    let a_tangent = the_circle.d1(0.0).1;
    if a_tangent.dot(&the_tangent) < 0.0 {
        the_axis.reverse();
        ax2 = make_ax2(&the_center_of_circ, &the_axis, &vec1)?;
        the_circle = make_circle(&ax2, vec1.magnitude());
    }

    // `cxx:1012-1015`: `BRepLib_MakeEdge(aME(theCircle, 0., theAngle))`.
    Ok(make_edge_range(&the_circle, 0.0, the_angle))
}

/// `gp_Ax2(theCenterOfCirc, theAxis.Direction(), Vec1)` (`cxx:1002`, `:1009`).
fn make_ax2(center: &GpPnt, axis: &GpAx1, vx: &GpVec) -> Result<GpAx2, String> {
    let x_dir =
        GpDir::from_vec(vx).map_err(|e| format!("Standard_ConstructionError: gp_Dir({e})"))?;
    GpAx2::new(*center, *axis.direction(), x_dir)
        .map_err(|e| format!("Standard_ConstructionError: gp_Ax2({e})"))
}

/// `GC_MakeCircle(theAx2, R)` (`cxx:1003`, `:1010`): `Geom_Circle(A2, R)`.
fn make_circle(ax2: &GpAx2, radius: f64) -> Arc<dyn Curve> {
    Arc::new(GeomCircle::new(occt_core::gp::GpCirc::new(*ax2, radius))) as Arc<dyn Curve>
}

/// `mypaths` tangents of the section `k` (`cxx:1041-1089`).
fn averaged_tangent(
    mp: &BRepOffsetApiMiddlePath,
    k: usize,
    i: usize,
    j: usize,
) -> Result<GpVec, String> {
    let mut pnts: Vec<GpPnt> = Vec::new();
    for indp in 1..=mp.my_paths.len() {
        let mut a_tangent;
        if k == i {
            // `cxx:1047-1055`: the tangent at the beginning of the first
            // segment.
            let a_shape = path_at(mp, indp, k)?;
            if a_shape.shape_type() == ShapeType::Vertex {
                continue;
            }
            a_tangent = tangent_of_edge(a_shape, true);
        } else if k == j {
            // `cxx:1056-1064`: the tangent at the end of the last segment.
            let a_shape = path_at(mp, indp, k - 1)?;
            if a_shape.shape_type() == ShapeType::Vertex {
                continue;
            }
            a_tangent = tangent_of_edge(a_shape, false);
        } else {
            // `cxx:1065-1075`: the average of the two neighbouring tangents.
            let prev = path_at(mp, indp, k - 1)?;
            let cur = path_at(mp, indp, k)?;
            if prev.shape_type() == ShapeType::Vertex || cur.shape_type() == ShapeType::Vertex {
                continue;
            }
            let tangent1 = tangent_of_edge(prev, false);
            let tangent2 = tangent_of_edge(cur, true);
            a_tangent = tangent1.added(&tangent2);
        }
        // `cxx:1076-1077`.
        a_tangent.normalize();
        pnts.push(GpPnt::from_xyz(a_tangent.xyz()));
    }
    // `cxx:1082-1089`: `GeomLib::Inertia(PntArray, theBary, ...)`, then
    // `gp_Vec theTangent(theBary.XYZ())`.
    let inertia = geom_lib_inertia(&pnts);
    Ok(GpVec::from_xyz(inertia.bary.xyz()))
}

/// `Build` after the sections (`cxx:845-1147`).
pub(super) fn build_middle_path(
    mp: &mut BRepOffsetApiMiddlePath,
    ctx: &mut Ctx,
) -> Result<(), String> {
    // `cxx:845-847`.
    let nb_sec_faces = ctx.sections_edges.len();
    let nb_e = ctx.edge_seq.len();

    // `cxx:848-873`: one face (or wire) per section.
    let mut sec_faces: Vec<TopoShape> = Vec::with_capacity(nb_sec_faces);
    for i in 0..nb_sec_faces {
        let mut mw = MakeWire::new();
        for j in 0..nb_e {
            let an_edge = as_edge(
                ctx.sections_edges[i]
                    .get(j)
                    .and_then(|s| s.as_ref())
                    .ok_or_else(|| "Standard_OutOfRange: SectionsEdges".to_string())?,
            )?;
            // `cxx:856`: `MW.Add(anEdge)`; a failing add leaves the wire not
            // done, which `Wire()` tests below.
            let _ = mw.add(&an_edge);
        }
        if !mp.my_closed_section {
            // `cxx:859-863`.
            let a_wire = mw.wire();
            let (v1, v2) = wire_vertices(&a_wire);
            let (v1, v2) = (
                v1.ok_or_else(|| "Standard_NullObject: TopExp::Vertices".to_string())?,
                v2.ok_or_else(|| "Standard_NullObject: TopExp::Vertices".to_string())?,
            );
            let closing = make_edge_vv(&v2, &v1)
                .ok_or_else(|| "BRepLib_LineThroughIdenticPoints: BRepLib_MakeEdge".to_string())?;
            let _ = mw.add(&closing);
        }
        // `cxx:864-872`.
        let a_wire = mw.wire();
        match make_face_from_wire(&a_wire, true) {
            Some(face) => sec_faces.push(face.0),
            None => sec_faces.push(a_wire.0),
        }
    }

    // `cxx:876-890`: the centre of mass of every section.
    let mut centers: Vec<GpPnt> = Vec::with_capacity(nb_sec_faces);
    for shape in &sec_faces {
        // `cxx:881-886`.
        let props = if shape.shape_type() == ShapeType::Face {
            surface_properties(shape)?.0
        } else {
            linear_properties(shape)?
        };
        centers.push(props.center());
    }

    // `cxx:892-1017`.
    let lin_tol = 1.0e-5;
    let ang_tol = 1.0e-7;
    let nb_mid = nb_sec_faces.saturating_sub(1);
    let mut mid_edges: Vec<Option<Edge>> = vec![None; nb_mid];
    for i in 1..nb_sec_faces {
        // `cxx:898-920`: the common type of every path segment of `i`.
        let mut type_of_mid_edge = GeomAbsCurveType::OtherCurve;
        for j in 1..=mp.my_paths.len() {
            let a_shape = path_at(mp, j, i)?;
            if a_shape.shape_type() == ShapeType::Vertex {
                type_of_mid_edge = GeomAbsCurveType::OtherCurve;
                break;
            }
            let an_edge = as_edge(a_shape)?;
            let a_type = type_of_edge(&an_edge);
            if j == 1 {
                type_of_mid_edge = a_type;
            } else if a_type != type_of_mid_edge {
                type_of_mid_edge = GeomAbsCurveType::OtherCurve;
                break;
            }
        }

        if type_of_mid_edge == GeomAbsCurveType::Line {
            // `cxx:921-924`.
            mid_edges[i - 1] =
                Some(make_edge_pp(&centers[i - 1], &centers[i]).ok_or_else(|| {
                    "BRepLib_LineThroughIdenticPoints: BRepLib_MakeEdge".to_string()
                })?);
        } else if type_of_mid_edge == GeomAbsCurveType::Circle {
            mid_edges[i - 1] = mid_edge_circle(mp, &centers, i, lin_tol, ang_tol)?;
        }
    }

    // `cxx:1020-1103`: build the missed edges.
    let mut i = 1usize;
    while i < nb_sec_faces {
        if mid_edges[i - 1].is_none() {
            // `cxx:1025-1030`: the first non-null edge after `i`.
            let mut j = i + 1;
            while j < nb_sec_faces {
                if mid_edges[j - 1].is_some() {
                    break;
                }
                j += 1;
            }
            // From `i` to `j-1` every edge is null.
            let mut the_points: Vec<GpPnt> = Vec::with_capacity(j - i + 1);
            let mut the_tangents: Vec<GpVec> = Vec::with_capacity(j - i + 1);
            // `cxx:1037-1040`.
            for k in i..=j {
                the_points.push(centers[k - 1]);
            }
            // `cxx:1041-1089`.
            for k in i..=j {
                the_tangents.push(averaged_tangent(mp, k, i, j)?);
            }
            // `cxx:1090`: `theFlags->Init(true)`.
            let the_flags = vec![true; the_points.len()];

            // `cxx:1092-1101`.
            let mut interpol = occt_geom::geom_api_interpolate::GeomAPIInterpolate::new(
                the_points, false, lin_tol,
            )
            .map_err(|e| format!("Standard_ConstructionError: GeomAPI_Interpolate({e})"))?;
            interpol
                .load_tangents(the_tangents, the_flags, true)
                .map_err(|e| {
                    format!("Standard_ConstructionError: GeomAPI_Interpolate::Load({e})")
                })?;
            interpol.perform();
            // `cxx:1094-1098`: the "Interpolation failed" console write is
            // UNPORTED; `Curve()` raises when the interpolation is not done.
            let inter_curve = interpol
                .curve()
                .ok_or_else(|| "StdFail_NotDone: GeomAPI_Interpolate::Curve".to_string())?;
            let inter_curve: Arc<dyn Curve> = Arc::new(inter_curve.clone());
            // `cxx:1100`.
            mid_edges[i - 1] = make_edge_over_curve(&inter_curve);
            // `cxx:1101`.
            i = j;
        }
        i += 1;
    }

    // `cxx:1105-1116`.
    let mut make_final_wire = MakeWire::new();
    for mid in &mid_edges {
        if let Some(edge) = mid {
            let _ = make_final_wire.add(edge);
        }
    }
    if !make_final_wire.is_done() {
        return Err("StdFail_NotDone: BRepLib_MakeWire".to_string());
    }
    mp.my_shape = Some(make_final_wire.wire().0);

    Ok(())
}
