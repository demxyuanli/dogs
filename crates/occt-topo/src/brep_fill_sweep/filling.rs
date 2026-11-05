//! `Filling` (`BRepFill_Sweep.cxx:867-1338`): the corner face of a `Round`
//! transition, a surface of revolution swept between the two adjacent corner
//! profiles.
//!
//! The reference has a `GeomPlate` refinement block at `cxx:1003-1179`; it is
//! entirely inside a `/* ... */` comment in 7.9.3, so the live control flow has
//! no plate surface (and no `pointu_f` / `pointu_l`), and neither has the port.
//!
//! Translation notes:
//! - `BRep_Tool::Curve(E, f, l)` and `BRep_Tool::Pnt(V)` apply the shape
//!   location, hence [`BRepTool::edge_curve_world`] /
//!   [`BRepTool::vertex_point_world`].
//! - `BRep_Builder::UpdateEdge(E, C2d, S, L, Tol)` is
//!   [`GeometryRegistry::update_edge_pcurve_on_surface`]; the seam overload
//!   `UpdateEdge(E, C1, C2, S, L, Tol)` is
//!   [`GeometryRegistry::set_edge_pcurves_on_surface`].
//! - `BRep_Builder::Range(E, f, l)` is [`GeometryRegistry::set_ranges_all`];
//!   `BRep_TVertex::Tolerance(T)` is
//!   [`GeometryRegistry::set_vertex_tolerance`].
//! - `BRepAdaptor_Curve(E)::LastParameter()` is the 3D range of `E`, and
//!   `BRepAdaptor_Surface(F)::D1` differs from the surface's own `D1` only by
//!   the face location, which leaves the normal's direction unchanged.
//! - `Surf->VIso(V)` is [`Surface::v_iso_curve`] with the generic
//!   [`iso_u`] fallback, and `Surf->DN(0, V, 1, 0)` is
//!   [`GeomSurfaceOfRevolution::eval_dn_rev`].

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx2, GpAx2d, GpDir2d, GpPnt2d, GpTrsf, GpVec, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION, SQUARE_CONFUSION};
use occt_geom::geom_lib_same_range_2d::same_range_2d;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom::surface::Surface;
use occt_geom::surface_of_revolution::GeomSurfaceOfRevolution;
use occt_geom::surface_ops::iso_u;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_lib_continuity::set_stored_continuity;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, ShapeMapKey};

use super::build_edge::build_edge_3d;
use super::build_face::build_face;

/// `gp::DX2d()`.
fn dx2d() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).expect("DX2d")
}

/// `gp::DY2d()`.
fn dy2d() -> GpDir2d {
    GpDir2d::new(0.0, 1.0).expect("DY2d")
}

/// What `Filling` returns: the reference's `bool` plus its three out
/// parameters (`Aux1`, `Aux2`, `Result`).
pub struct FillingResult {
    /// The reference's return value.
    pub done: bool,
    /// `Aux1` after the call (`cxx:1262-1270`).
    pub aux1: Option<Edge>,
    /// `Aux2` after the call (`cxx:1272-1280`).
    pub aux2: Option<Edge>,
    /// `Result`, null when `BuildFace` failed.
    pub result: Option<Face>,
}

/// `Filling(EF, F1, EL, F2, EEmap, Tol, Axe, TangentOnPart1, Aux1, Aux2,
/// Result)` (`cxx:867-1338`).
///
/// `Aux1` / `Aux2` are the reference's in/out edges: `None` is a null
/// `TopoDS_Edge`. `Err` is reserved for the branches where the reference
/// throws (a missing 3D curve, a degenerate rectangular trim) or where the port
/// cannot reproduce the reference value (`BuildEdge` answering null); the
/// reference's own `return false` arms come back as `done = false`.
#[allow(clippy::too_many_arguments)]
pub fn filling(
    ef: &Edge,
    f1: &Face,
    el: &Edge,
    f2: &Face,
    ee_map: &mut HashMap<ShapeMapKey, crate::shape::TopoShape>,
    tol: f64,
    axe: &GpAx2,
    tangent_on_part1: &GpVec,
    aux1: Option<Edge>,
    aux2: Option<Edge>,
) -> Result<FillingResult, String> {
    let reg = GeometryRegistry::global();

    // `cxx:888-898` — `TopExp::Vertices` with `CumOri = false`.
    let e1 = ef;
    let e2 = el;
    let (vf, vl) = crate::topo_tools::edge_vertices(e1);
    let mut vf = vf.ok_or_else(|| "Filling: E1 has no first vertex (cxx:892)".to_string())?;
    let mut vl = vl.ok_or_else(|| "Filling: E1 has no last vertex (cxx:892)".to_string())?;
    vf.0.set_orientation(Orientation::Forward);
    vl.0.set_orientation(Orientation::Forward);

    let (a, b) = crate::topo_tools::edge_vertices(e2);
    let mut v1 = a.ok_or_else(|| "Filling: E2 has no first vertex (cxx:896)".to_string())?;
    let mut v2 = b.ok_or_else(|| "Filling: E2 has no last vertex (cxx:896)".to_string())?;
    v1.0.set_orientation(Orientation::Reversed);
    v2.0.set_orientation(Orientation::Reversed);

    // `cxx:900-904` — `B.MakeEdge(E3)` / `B.MakeEdge(E4)` make empty edges,
    // not null ones, so two of them are distinct TShapes.
    let mut e3 = Edge::new();
    let mut e4 = Edge::new();
    let mut with_e3 = false;
    let mut with_e4 = false;

    // `cxx:906-911`.
    if let Some(a) = aux1.as_ref() {
        if !is_same(&vf.0, &v1.0) {
            e3 = a.clone();
            with_e3 = true;
        }
    }

    // `cxx:913-923`.
    if is_same(&vf.0, &vl.0) {
        e4 = e3.clone();
        e4.0.reverse();
        with_e4 = with_e3;
    } else if let Some(a) = aux2.as_ref() {
        if !is_same(&vl.0, &v2.0) {
            e4 = a.clone();
            with_e4 = true;
        }
    }

    // `cxx:929-933`.
    let prof1 = BRepTool::edge_curve_world(e1)
        .ok_or_else(|| "Filling: BRep_Tool::Curve(E1) is null (cxx:930)".to_string())?;
    let (e1_first, e1_last) = BRepTool::edge_parameters(e1);
    let prof2 = BRepTool::edge_curve_world(e2)
        .ok_or_else(|| "Filling: BRep_Tool::Curve(E2) is null (cxx:933)".to_string())?;
    let (e2_first, e2_last) = BRepTool::edge_parameters(e2);

    // `cxx:936`.
    let same_curve_domain = ef.0.orientation() != el.0.orientation();

    // `cxx:941-942` — `aTf.SetTransformation(Axe)`.
    let mut a_tf = GpTrsf::identity();
    a_tf.set_transformation(&axe.to_ax3());

    // `cxx:947-963` — the sample furthest from the axis of revolution.
    let a_prm = [e1_first, 0.5 * (e1_first + e1_last), e1_last];
    let a_p1: Vec<occt_core::gp::GpPnt> =
        a_prm.iter().map(|p| prof1.d0(*p).transformed(&a_tf)).collect();
    let mut a_max_idx: i32 = -1;
    let mut a_max_dist = f64::MIN; // `RealFirst()`
    for (i, p) in a_p1.iter().enumerate() {
        let d = p.x() * p.x() + p.z() * p.z();
        if d > a_max_dist {
            a_max_dist = d;
            a_max_idx = i as i32;
        }
    }

    // `cxx:965-970`.
    let a_prm2 = [e2_first, 0.5 * (e2_first + e2_last), e2_last];
    let idx2 = if same_curve_domain { a_max_idx } else { 2 - a_max_idx };
    let a_p2 = prof2.d0(a_prm2[idx2 as usize]).transformed(&a_tf);
    let a_v1 = GpVec2d::new(a_p1[a_max_idx as usize].z(), a_p1[a_max_idx as usize].x());
    let a_v2 = GpVec2d::new(a_p2.z(), a_p2.x());
    if a_v1.square_magnitude() <= RESOLUTION || a_v2.square_magnitude() <= RESOLUTION {
        return Ok(FillingResult { done: false, aux1, aux2, result: None });
    }

    // `cxx:972-981`.
    let mut angle = a_v1.angle(&a_v2);
    let mut axe_rot = GpAx1::new(axe.location(), *axe.y_direction());
    if angle < 0.0 {
        angle = -angle;
        axe_rot.reverse();
    }

    // `cxx:983-985` — `Geom_SurfaceOfRevolution` trimmed to
    // `[0, Angle] x [f1, l1]`.
    let rev = Arc::new(GeomSurfaceOfRevolution::new(prof1.clone(), axe_rot));
    let surf: Arc<dyn Surface> = Arc::new(GeomRectangularTrimmedSurface::uv(
        rev.clone(),
        0.0,
        angle,
        e1_first,
        e1_last,
    )?);

    // `cxx:987-1000`: control the direction of the rotation. `Surf->DN` is the
    // basis revolution's derivative — the rectangular trim does not change
    // parameter values.
    let mut to_reverse_result = false;
    let d1u = rev.eval_dn_rev(0.0, a_prm[a_max_idx as usize], 1, 0);
    if d1u.angle(tangent_on_part1) > std::f64::consts::PI / 2.0 {
        to_reverse_result = true;
    }

    // `cxx:1011-1019`.
    let l = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 0.0), dy2d()));
    let c1: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(
        Arc::new(l),
        e1_first,
        e1_last,
    ));
    let l = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(angle, 0.0), dy2d()));
    let mut c2_trimmed = Geom2dTrimmedCurve::new(Arc::new(l), e1_first, e1_last);

    // `cxx:1021-1037`: control the direction and the range. The 2D curve is
    // reversed when `E2` is traversed opposite to `E1`.
    let p2d = c2_trimmed.d0(e1_first);
    let p_at_first = surf.d0(p2d.x(), p2d.y());
    let p2d = c2_trimmed.d0(e1_last);
    let p_at_last = surf.d0(p2d.x(), p2d.y());
    let p = BRepTool::vertex_point_world(&v1);
    if p.distance(&p_at_last) + tol < p.distance(&p_at_first) {
        <Geom2dTrimmedCurve as Curve2d>::reverse(&mut c2_trimmed);
        std::mem::swap(&mut v1, &mut v2);
    }

    // `cxx:1038-1045` — `GeomLib::SameRange` onto `[f2, l2]`.
    let c2: Arc<dyn Curve2d> = Arc::new(c2_trimmed);
    let c3_same = same_range_2d(
        PCONFUSION,
        &c2,
        c2.first_parameter(),
        c2.last_parameter(),
        e2_first,
        e2_last,
    )?;
    let c2 = c3_same;

    // `cxx:1058-1063`: the two seam pcurves across the rotation angle.
    let l = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, e1_first), dx2d()));
    let c3: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(Arc::new(l), 0.0, angle));
    let l = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, e1_last), dx2d()));
    let c4: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(Arc::new(l), 0.0, angle));

    // `cxx:1181-1187` — `B.UpdateEdge(E1/E2, C, Surf, Loc, Precision::Confusion())`.
    reg.update_edge_pcurve_on_surface(&e1.0, &surf, Some(c1.clone()), CONFUSION);
    reg.update_edge_pcurve_on_surface(&e2.0, &surf, Some(c2.clone()), CONFUSION);

    // `cxx:1188-1211`: the two borders coincide (`E3.IsSame(E4)`).
    if is_same(&e3.0, &e4.0) {
        if !with_e3 {
            let c3d = surf
                .v_iso_curve(e1_first)
                .unwrap_or_else(|| iso_u(&*surf, e1_first));
            e3 = build_edge_3d(&c3d, &c3, &surf, &vf, &v1, 0.0, angle, CONFUSION).ok_or_else(
                || "Filling: BuildEdge(C3d, C3, Surf, Vf, V1, 0, Angle) failed (cxx:1192)".to_string(),
            )?;
        } else {
            // `cxx:1197-1204`: `BRepAdaptor_Curve(E3).LastParameter()` is the
            // 3D range of `E3`.
            let angle_old = BRepTool::edge_parameters(&e3).1;
            if angle > angle_old {
                reg.set_ranges_all(&e3.0, 0.0, angle);
                let v = super::build_shell::oriented_edge_vertices(&e3)
                    .1
                    .ok_or_else(|| "Filling: TopExp::LastVertex(E3) is null (cxx:1201)".to_string())?;
                reg.set_vertex_tolerance(&v.0, CONFUSION);
            }
        }

        // `cxx:1206-1209`.
        reg.set_edge_pcurves_on_surface(&e3.0, &surf, vec![c3.clone(), c4.clone()], CONFUSION);
        e4 = e3.clone();
        e4.0.reverse();
    } else {
        // `cxx:1212-1232`.
        if !with_e3 {
            let c3d = surf
                .v_iso_curve(e1_first)
                .unwrap_or_else(|| iso_u(&*surf, e1_first));
            e3 = build_edge_3d(&c3d, &c3, &surf, &vf, &v1, 0.0, angle, CONFUSION).ok_or_else(
                || "Filling: BuildEdge(C3d, C3, Surf, Vf, V1, 0, Angle) failed (cxx:1216)".to_string(),
            )?;
        } else {
            let angle_old = BRepTool::edge_parameters(&e3).1;
            if angle > angle_old {
                reg.set_ranges_all(&e3.0, 0.0, angle);
                let v = super::build_shell::oriented_edge_vertices(&e3)
                    .1
                    .ok_or_else(|| "Filling: TopExp::LastVertex(E3) is null (cxx:1225)".to_string())?;
                reg.set_vertex_tolerance(&v.0, CONFUSION);
            }
            reg.update_edge_pcurve_on_surface(&e3.0, &surf, Some(c3.clone()), CONFUSION);
        }

        // `cxx:1234-1251`.
        if !with_e4 {
            let c3d = surf
                .v_iso_curve(e1_last)
                .unwrap_or_else(|| iso_u(&*surf, e1_last));
            e4 = build_edge_3d(&c3d, &c4, &surf, &vl, &v2, 0.0, angle, CONFUSION).ok_or_else(
                || "Filling: BuildEdge(C3d, C4, Surf, Vl, V2, 0, Angle) failed (cxx:1237)".to_string(),
            )?;
        } else {
            let angle_old = BRepTool::edge_parameters(&e4).1;
            if angle > angle_old {
                reg.set_ranges_all(&e4.0, 0.0, angle);
                let v = super::build_shell::oriented_edge_vertices(&e4)
                    .1
                    .ok_or_else(|| "Filling: TopExp::LastVertex(E4) is null (cxx:1246)".to_string())?;
                reg.set_vertex_tolerance(&v.0, CONFUSION);
            }
            reg.update_edge_pcurve_on_surface(&e4.0, &surf, Some(c4.clone()), CONFUSION);
        }
    }

    // `cxx:1254-1255` — `BuildFace(Surf, E1, E3, E2, E4, EEmap, false, false, Result)`.
    let mut result = build_face(&surf, e1, &e3, e2, &e4, ee_map, false, false);

    // `cxx:1258-1259` — `B.Continuity(E, F, Result, GeomAbs_G1)`.
    if let Some(r) = result.as_ref() {
        set_stored_continuity(e1, f1, r, 1);
        set_stored_continuity(e2, f2, r, 1);
    }

    // `cxx:1262-1281` — render the calculated borders; a degenerate border is
    // nullified (`B.MakeEdge(Aux)` on an edge, i.e. no usable border).
    let out_aux1 = if !BRepTool::is_degenerated(&e3) { Some(e3.clone()) } else { None };
    let out_aux2 = if !BRepTool::is_degenerated(&e4) { Some(e4.clone()) } else { None };

    // `cxx:1286-1316` — provide a correct normal computation: the normal is
    // not evaluated at the singularity point.
    let face1_key = GeometryRegistry::shape_key(&f1.0);
    let mut c1_var = c1.clone();
    let mut angle_out = f64::MIN; // `RealFirst()`
    for i in 0..3 {
        let p2d = c1_var.d0(a_prm[i]);
        let (_, du, dv) = surf.d1(p2d.x(), p2d.y());
        let n1 = du.crossed(&dv);
        if n1.square_magnitude() < SQUARE_CONFUSION {
            continue;
        }

        // `cxx:1302-1305` — `C1` becomes the pcurve of `E1` on `F1`, and stays
        // so for the following iterations.
        let cos = reg.edge_pcurve(&e1.0, face1_key).ok_or_else(|| {
            "Filling: BRep_Tool::CurveOnSurface(E1, F1) is null (cxx:1302)".to_string()
        })?;
        c1_var = cos;
        let p2d = c1_var.d0(a_prm[i]);

        // `cxx:1306-1308` — `BRepAdaptor_Surface(F1)->D1`.
        let surf1 = BRepTool::face_surface(f1)
            .ok_or_else(|| "Filling: BRepAdaptor_Surface(F1) has no surface (cxx:1306)".to_string())?;
        let (_, du, dv) = surf1.d1(p2d.x(), p2d.y());
        let trsf = f1.0.location().transformation();
        let n2 = du.transformed(&trsf).crossed(&dv.transformed(&trsf));
        if n2.square_magnitude() < SQUARE_CONFUSION {
            continue;
        }

        angle_out = n1.angle(&n2);
        break;
    }

    // `cxx:1318-1320`.
    if angle_out == f64::MIN {
        return Ok(FillingResult { done: false, aux1, aux2, result });
    }

    // `cxx:1322-1334`.
    if let Some(r) = result.as_mut() {
        if (f1.0.orientation() == Orientation::Reversed) ^ (angle_out > std::f64::consts::PI / 2.0)
        {
            r.0.set_orientation(Orientation::Reversed);
        } else {
            r.0.set_orientation(Orientation::Forward);
        }
        if to_reverse_result {
            r.0.reverse();
        }
    }

    Ok(FillingResult { done: true, aux1: out_aux1, aux2: out_aux2, result })
}
