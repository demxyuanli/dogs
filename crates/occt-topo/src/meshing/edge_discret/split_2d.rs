//! 2D pcurve splitting by deflection (`BRepMesh_EdgeTessellator` / k-part helpers).

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

use super::prelude::*;
use super::*;

/// Distance from `p` to the segment `a..b`.
/// `BRepMesh_CurveTessellator::splitByDeflection2d` (`cxx:157-188`).
///
/// `aNodesNb` is captured before the pcurve loop. Each pcurve rebuilds
/// `aParamArray` from the current discretizer's first `aNodesNb` parameters
/// (`cxx:175-186`), so splits from an earlier pcurve shift which original
/// spans the later pcurves still see.
pub fn split_by_deflection2d(
    curve: &dyn Curve,
    params: &mut Vec<f64>,
    pcurves: &[(Arc<dyn Curve2d>, Arc<dyn Surface>)],
    lin_def: f64,
    min_size: f64,
) {
    if params.len() < 2 {
        return;
    }
    let sq_def = lin_def * lin_def;
    let sq_min = sq_def.max(min_size * min_size);
    let nodes_nb = params.len();
    for (c2d, surf) in pcurves {
        // `BRepMesh_CurveTessellator.cxx:168-171`: skip only
        // `STANDARD_TYPE(Geom_Plane)`, not a sampled is_planar classify.
        // Offset / BSpline / Rev stay in the loop (unlike Value cxx:241-246).
        if surf.gp_pln().is_some() {
            continue;
        }
        params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        params.dedup_by(|a, b| (*a - *b).abs() < PCONFUSION);
        let n = nodes_nb.min(params.len());
        let arr = params[..n].to_vec();
        for w in arr.windows(2) {
            split_segment_2d(
                curve,
                c2d.as_ref(),
                surf.as_ref(),
                w[0],
                w[1],
                1,
                sq_def,
                sq_min,
                params,
            );
        }
    }
    params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    params.dedup_by(|a, b| (*a - *b).abs() < PCONFUSION);
}

/// `BRepMesh_CurveTessellator::splitSegment` (`cxx:274-340`).
pub(super) fn split_segment_2d(
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    surf: &dyn Surface,
    first: f64,
    last: f64,
    iter: i32,
    sq_def: f64,
    sq_min: f64,
    params: &mut Vec<f64>,
) {
    if iter > 10 {
        return;
    }
    if (last - first).abs() < 2.0 * PCONFUSION {
        return;
    }
    let (cf, cl) = (c2d.first_parameter(), c2d.last_parameter());
    if cf.is_finite() && first - cf < -PCONFUSION {
        return;
    }
    if cl.is_finite() && last - cl > PCONFUSION {
        return;
    }
    let uvf = c2d.d0(first);
    let uvl = c2d.d0(last);
    let p3df = surf.d0(uvf.x(), uvf.y());
    let p3dl = surf.d0(uvl.x(), uvl.y());
    if p3df.square_distance(&p3dl) < sq_min {
        return;
    }
    let uvm = GpPnt2d::new(0.5 * (uvf.x() + uvl.x()), 0.5 * (uvf.y() + uvl.y()));
    let mid_surf = surf.d0(uvm.x(), uvm.y());
    let v1 = GpVec::from_pnts(&p3df, &mid_surf);
    if v1.square_magnitude() < sq_min {
        return;
    }
    let mut vec = GpVec::from_pnts(&p3df, &p3dl);
    let mag = vec.magnitude();
    if mag < 1e-30 {
        return;
    }
    vec = vec.multiplied_scalar(1.0 / mag);
    let proj = vec.multiplied_scalar(v1.dot(&vec));
    let dist = v1.subtracted(&proj);
    if dist.square_magnitude() < sq_def {
        return;
    }
    let midpar = 0.5 * (first + last);
    // `cxx:335-336` AddPoint(myCurve.D0(midpar), midpar, false); params-only
    // path records midpar — Tessellate3d evaluates the SameParam 3D curve later.
    params.push(midpar);
    split_segment_2d(
        c3d, c2d, surf, first, midpar, iter + 1, sq_def, sq_min, params,
    );
    split_segment_2d(
        c3d, c2d, surf, midpar, last, iter + 1, sq_def, sq_min, params,
    );
}

/// `BRepMesh_CurveTessellator::Value` (`cxx:223-270`) for CurveOnSurface:
/// drop an interior sample whose pcurve UV is outside the surface range
/// (padded by `U/VResolution(Confusion)`) when the 3D point is farther than
/// the edge tolerance from `surface(UV)`. Analytic and periodic surfaces keep
/// every sample.
///
/// Type gate matches `cxx:241-246`: only `BSpline` / `Bezier` / `OtherSurface`
/// run the UV out-of-range check. `OffsetSurface` is excluded even when the
/// basis exposes poles (`GetType() == GeomAbs_OffsetSurface`).
pub fn curve_tessellator_value_ok(
    pcurve: &dyn Curve2d,
    surface: &dyn Surface,
    parameter: f64,
    point: &GpPnt,
    edge_tol: f64,
) -> bool {
    use super::super::range_splitter::MeshSurfaceType;
    let ty = super::super::range_splitter::mesh_surface_type(surface);
    match ty {
        MeshSurfaceType::BSplineSurface
        | MeshSurfaceType::BezierSurface
        | MeshSurfaceType::OtherSurface => {}
        _ => return true,
    }
    if surface.is_u_periodic() || surface.is_v_periodic() {
        return true;
    }
    let uv = pcurve.d0(parameter);
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();
    let (du, dv) = surface
        .uv_resolution(CONFUSION)
        .unwrap_or((CONFUSION, CONFUSION));
    if uv.x() > u0 - du && uv.x() < u1 + du && uv.y() > v0 - dv && uv.y() < v1 + dv {
        return true;
    }
    let on_surf = surface.d0(uv.x(), uv.y());
    point.square_distance(&on_surf) < edge_tol * edge_tol
}

/// Parameter of an INTERNAL vertex on `edge`. PointOnCurve is unported;
/// `GeomAPI_ProjectPointOnCurve` on the 3D curve stands in for
/// `BRep_Tool::Parameter(V,E)` (`BRep_Tool.cxx:1503-1510`).
pub(super) fn vertex_parameter_on_edge(vertex: &Vertex, edge: &Edge) -> Option<f64> {
    let pnt = BRepTool::vertex_point_world(vertex);
    let curve = BRepTool::edge_curve_adaptor_world(edge)?;
    let tol = BRepTool::vertex_tolerance(vertex);
    geom_api::project_point_on_curve(curve.as_ref(), &pnt, tol).map(|p| p.parameter)
}

/// `BRepMesh_CurveTessellator.cxx:100-114` min-point threshold from `GetType`.
/// Circle / Ellipse / Parabola / Hyperbola use 4. Periodic non-NURBS (e.g.
/// trimmed conics) also use 4. BSpline / Bezier stay at 2.
pub fn tessellator_min_points(curve: &dyn Curve) -> usize {
    if curve.circle_radius().is_some()
        || curve.gp_ellipse().is_some()
        || curve.gp_parabola().is_some()
        || curve.gp_hyperbola().is_some()
    {
        return 4;
    }
    if curve.is_periodic() && curve.bspline_poles().is_none() && curve.bezier_poles().is_none() {
        return 4;
    }
    2
}

/// `Adaptor3d_CurveOnSurface::EvalKPart` reverse of a classified circle
/// (`cxx:1597-1602` and the same SetDirection pattern on every iso Circle).
pub(super) fn reverse_kpart_circ(circ: &mut GpCirc) {
    let mut ax = circ.position();
    ax.set_direction(ax.direction().reversed());
    circ.set_position(&ax);
}

/// `to3d(Pl, Circ2d)` (`Adaptor3d_CurveOnSurface.cxx:57-83`).
pub(super) fn circ2d_to3d_on_plane(pl: &occt_core::gp::GpPln, c: &occt_core::gp::GpCirc2d) -> Option<GpCirc> {
    let loc = c.position().location();
    let p = slib::plane_value(pl, loc.x(), loc.y());
    let vx2 = *c.position().x_direction();
    let vy2 = *c.position().y_direction();
    let vx = GpVec::from_xyz(
        &pl.pos
            .x_direction()
            .xyz()
            .multiplied(vx2.x)
            .added(&pl.pos.y_direction().xyz().multiplied(vx2.y)),
    );
    let vy = GpVec::from_xyz(
        &pl.pos
            .x_direction()
            .xyz()
            .multiplied(vy2.x)
            .added(&pl.pos.y_direction().xyz().multiplied(vy2.y)),
    );
    let n = vx.crossed(&vy);
    let xd = GpDir::from_vec(&vx).ok()?;
    let nd = GpDir::from_vec(&n).ok()?;
    let ax = GpAx2::new(p, nd, xd).ok()?;
    Some(GpCirc::new(ax, c.radius()))
}

/// Rotate a V-iso circle by the 2d line U (`cxx:1594-1596`).
pub(super) fn rotate_iso_v(circ: &mut GpCirc, axis: &occt_core::gp::GpAx3, u: f64) {
    let drev = axis
        .x_direction()
        .crossed(axis.y_direction())
        .unwrap_or(*axis.x_direction());
    let axe = GpAx1::new(axis.location(), drev);
    circ.rotate(&axe, u);
}

/// Kparts of `Adaptor3d_CurveOnSurface` (`cxx:1552-1733`): `myCirc` / `myLin`
/// are cached by `EvalKPart`, the three conics are the on-demand queries
/// `Ellipse` / `Hyperbola` / `Parabola` (`cxx:1399-1416`).
#[derive(Default)]
pub(super) struct KPart {
    pub(super) circ: Option<GpCirc>,
    pub(super) lin: Option<GpLin>,
    pub(super) elips: Option<GpElips>,
    pub(super) parab: Option<GpParab>,
    pub(super) hypr: Option<GpHypr>,
    /// `Bezier()` (`cxx:1458-1485`) and `BSpline()` (`cxx:1489-1520`): the
    /// pole class lifted onto the plane with `to3d(mySurface->Plane(), Pole2d)`.
    /// Both cxx bodies are `GeomLib::To3d`'s Bezier (`GeomLib.cxx:582-602`) and
    /// B-spline (`GeomLib.cxx:603-632`) arms verbatim, so the port reuses that
    /// single faithful implementation instead of writing a second one.
    pub(super) poles: Option<Arc<dyn Curve>>,
}

/// `EvalKPart` plane arm (`cxx:1558-1577`): `myType = CTy` (`cxx:1560`), so the
/// pcurve class alone decides, and each elementary class answers its 3D query
/// through `to3d(mySurface->Plane(), myCurve->...)` (`cxx:1381-1416`).
///
/// `GeomAbs_OffsetCurve` / `GeomAbs_OtherCurve` pcurves leave every field empty,
/// which is the cxx behaviour for both (`Offset()` is not overridden on this
/// adaptor). A pole-classified pcurve fills `poles` instead.
pub(super) fn plane_k_part(pc: &dyn Curve2d, surf: &dyn Surface, pl: &GpPln) -> KPart {
    let mut k = KPart::default();
    if let Some(c2) = pc.gp_circ2d() {
        k.circ = circ2d_to3d_on_plane(pl, &c2);
    } else if pc.is_line() {
        let (uv, duv) = pc.d1(0.0);
        let (p, d1u, d1v) = surf.d1(uv.x(), uv.y());
        let v = d1u
            .multiplied_scalar(duv.x())
            .added(&d1v.multiplied_scalar(duv.y()));
        if let Ok(dir) = GpDir::from_vec(&v) {
            k.lin = Some(GpLin::from_pnt_dir(p, dir));
        }
    } else if let Some(e2) = pc.gp_elips2d() {
        // `Ellipse()` (`cxx:1399-1402`).
        k.elips = clib::to_3d_elips(&pl.pos.ax2(), &e2).ok();
    } else if let Some(p2) = pc.gp_parab2d() {
        // `Parabola()` (`cxx:1413-1416`).
        k.parab = clib::to_3d_parab(&pl.pos.ax2(), &p2).ok();
    } else if let Some(h2) = pc.gp_hypr2d() {
        // `Hyperbola()` (`cxx:1406-1409`).
        k.hypr = clib::to_3d_hypr(&pl.pos.ax2(), &h2).ok();
    } else if pc.is_bezier2d() || pc.is_bspline2d() {
        // `Bezier()` (`cxx:1458-1485`) / `BSpline()` (`cxx:1489-1520`): both lift
        // the pcurve poles with `to3d(mySurface->Plane(), Pole2d)` and branch on
        // `myCurve->IsRational()`, which is `GeomLib::To3d`'s Bezier
        // (`GeomLib.cxx:582-602`) and B-spline (`GeomLib.cxx:603-632`) arms.
        k.poles = occt_geom::geom_lib::to_3d(&pl.pos.ax2(), pc);
    }
    k
}

/// `Adaptor3d_CurveOnSurface::EvalKPart` (`cxx:1552-1733`).
///
/// The pole-based plane classes `Bezier()` (`cxx:1458-1485`) and `BSpline()`
/// (`cxx:1489-1520`) are reached through `plane_k_part`: a plane surface with a
/// pole-classified pcurve leaves `myType` at `GeomAbs_BezierCurve` /
/// `GeomAbs_BSplineCurve` (`cxx:1560`), and the lifted 3D curve is cached in
/// `kpart_poles` so `Curve::bezier_poles` / `bspline_poles` / `bspline_knots` /
/// `bspline_weights` / `nurbs_degree` answer for that cos exactly as
/// `Extrema_CurveTool` forwards them.
///
/// UNPORTED: `GeomAbs_OffsetCurve`. `Offset()` is not overridden on this
/// adaptor, and the class is the `default:` arm in all three consumers
/// (`Extrema_GGExtPC.hxx:390-502`, `Extrema_GFuncExtPC.hxx:115-128`,
/// `Extrema_GLocateExtPC.hxx:110-140`), so dropping it is unobservable.
///
/// `GeomAdaptor_Surface` of a rectangular trim reports the basis GetType.
pub(super) fn eval_k_part(pc: &dyn Curve2d, surf: &dyn Surface) -> KPart {
    if let Some(basis) = surf.rectangular_trimmed_basis() {
        return eval_k_part(pc, basis.as_ref());
    }
    if let Some(pl) = surf.gp_pln() {
        return plane_k_part(pc, surf, &pl);
    }
    let Some(l2) = pc.gp_lin2d() else {
        return KPart::default();
    };
    let mut circ = None;
    let mut lin = None;
    let d = *l2.direction();
    let loc = l2.location();
    const DX2D: GpDir2d = GpDir2d { x: 1.0, y: 0.0 };
    const DY2D: GpDir2d = GpDir2d { x: 0.0, y: 1.0 };
    if d.is_parallel(&DX2D, ANGULAR) {
        if let Some(sph) = surf.gp_sphere() {
            if (loc.y().abs() - std::f64::consts::FRAC_PI_2).abs() >= PCONFUSION {
                let axis = sph.position();
                let mut c = slib::sphere_v_iso(axis, sph.radius(), loc.y());
                rotate_iso_v(&mut c, axis, loc.x());
                if d.is_opposite(&DX2D, ANGULAR) {
                    reverse_kpart_circ(&mut c);
                }
                circ = Some(c);
            }
        } else if let Some(cyl) = surf.gp_cylinder() {
            let axis = cyl.position();
            let mut c = slib::cylinder_v_iso(&axis, cyl.radius(), loc.y());
            rotate_iso_v(&mut c, &axis, loc.x());
            if d.is_opposite(&DX2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        } else if let Some(cone) = surf.gp_cone() {
            let axis = cone.position();
            let mut c = slib::cone_v_iso(&axis, cone.radius(), cone.semi_angle(), loc.y());
            rotate_iso_v(&mut c, &axis, loc.x());
            if d.is_opposite(&DX2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        } else if let Some(tor) = surf.gp_torus() {
            let axis = tor.position();
            let mut c = slib::torus_v_iso(axis, tor.major_radius(), tor.minor_radius(), loc.y());
            rotate_iso_v(&mut c, axis, loc.x());
            if d.is_opposite(&DX2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        }
    } else if d.is_parallel(&DY2D, ANGULAR) {
        if let Some(sph) = surf.gp_sphere() {
            let axis = sph.position();
            let mut c = slib::sphere_u_iso(axis, sph.radius(), 0.0);
            let drev = axis
                .x_direction()
                .crossed(&axis.direction())
                .unwrap_or(*axis.x_direction());
            let axe_y = GpAx1::new(axis.location(), drev);
            c.rotate(&axe_y, loc.y());
            rotate_iso_v(&mut c, axis, loc.x());
            if d.is_opposite(&DY2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        } else if let Some(cyl) = surf.gp_cylinder() {
            let mut l = slib::cylinder_u_iso(&cyl.position(), cyl.radius(), loc.x());
            let tr = GpVec::from_xyz(l.direction().xyz()).multiplied_scalar(loc.y());
            l = l.translated_vec(&tr);
            if d.is_opposite(&DY2D, ANGULAR) {
                l.set_direction(l.direction().reversed());
            }
            lin = Some(l);
        } else if let Some(cone) = surf.gp_cone() {
            let mut l = slib::cone_u_iso(
                &cone.position(),
                cone.radius(),
                cone.semi_angle(),
                loc.x(),
            );
            let tr = GpVec::from_xyz(l.direction().xyz()).multiplied_scalar(loc.y());
            l = l.translated_vec(&tr);
            if d.is_opposite(&DY2D, ANGULAR) {
                l.set_direction(l.direction().reversed());
            }
            lin = Some(l);
        } else if let Some(tor) = surf.gp_torus() {
            let mut c = slib::torus_u_iso(
                tor.position(),
                tor.major_radius(),
                tor.minor_radius(),
                loc.x(),
            );
            let axe = *c.position().axis();
            c.rotate(&axe, loc.y());
            if d.is_opposite(&DY2D, ANGULAR) {
                reverse_kpart_circ(&mut c);
            }
            circ = Some(c);
        }
    }
    KPart {
        circ,
        lin,
        ..KPart::default()
    }
}
