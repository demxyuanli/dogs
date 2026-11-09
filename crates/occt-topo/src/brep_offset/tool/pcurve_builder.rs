//! File-local `BuildPCurves` (`BRepOffset_Tool.cxx:327-483`) with its two
//! helpers `FindPeriod` (`:155-190`) and `PutInBounds` (`:197-316`), the pcurve
//! generator `TryProject` (`:1945-2006`) and `ExtentFace` (`:4311`) call.
//!
//! Three arms, in the reference's order:
//! 1. an existing `CurveOnSurface(E, F)` returns immediately (`cxx:330-334`);
//! 2. a Bezier / B-spline (unwrapped through `Geom_OffsetSurface`) basis tries
//!    to reuse one of the boundary isolation curves: the edge endpoints are
//!    projected on every edge of `BRepLib_MakeFace(theSurf, Confusion)` and the
//!    first boundary edge whose both projections are within `TolProj = 1.e-4`
//!    yields a `Geom2d_TrimmedCurve` (`cxx:347-428`);
//! 3. otherwise `ProjLib_ProjectedCurve(AS, AC, Confusion)` and the elementary
//!    curve type it reports (`cxx:430-466`).
//!
//! A `Standard_ConstructionError` (`cxx:480`) is the `Err` of the port.
//!
//! Deviations: (a) `C2d->Translate(...)` inside `PutInBounds` mutates a
//! `handle(Geom2d_Curve)&`, which the port answers as a returned
//! `Arc<dyn Curve2d>` (the caller reassigns); (b) `BRep_Builder::UpdateEdge(E,
//! C2d, F, tol)` writes only the pcurve slot here, because OCCT's overload
//! merges into the existing representation list with the *same* tolerances
//! (`BRep_Builder.cxx:1078-1130`), which the registry keeps per (edge, face)
//! key anyway.

use std::sync::Arc;

use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision::CONFUSION;
use occt_geom::projlib::ProjectedCurve;
use occt_geom::{Curve, GeomTrimmedCurveBasis, Surface};
use occt_geom2d::{
    trimmed::Geom2dTrimmedCurve, Curve2d, Geom2dCircle, Geom2dEllipse, Geom2dHyperbola, Geom2dLine,
    Geom2dParabola,
};

use crate::abs::ShapeType;
use crate::boptools_2d;
use crate::brep_lib_make_face::make_face_from_surface;
use crate::brep_lib_same_range::same_range;
use crate::brep_tool::BRepTool;
use crate::geom_abs_curve_type::{curve2d_kind, GeomAbsCurveType};
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::tool::adaptor::{GeomAdaptorCurve, SurfaceAdaptor};

/// `BRepAdaptor_Curve(E)`: the edge's 3D curve bounded by the edge range.
fn edge_adaptor(e: &Edge) -> Option<GeomAdaptorCurve> {
    let c = BRepTool::edge_curve_world(e)?;
    let (f, l) = BRepTool::edge_parameters(e);
    Some(GeomAdaptorCurve::with_range(c, f, l))
}

/// `BuildPCurves(E, F)` (`cxx:327-483`).
pub(super) fn build_pcurves(e: &Edge, f: &Face) -> Result<(), String> {
    // `cxx:330-334`.
    if boptools_2d::curve_on_surface_range(e, f).is_some() {
        return Ok(());
    }

    // `cxx:337`: `Tolerance = Precision::Confusion()` (the edge tolerance is
    // commented out in the reference).
    const TOLERANCE: f64 = CONFUSION;

    // `cxx:339-340`.
    let a_s = SurfaceAdaptor::of_face(f);
    let a_c = edge_adaptor(e)
        .ok_or_else(|| "BRepOffset_Tool::BuildPCurves: edge without a 3D curve".to_string())?;

    // `cxx:343-349`: `Geom_OffsetSurface` is unwrapped before the type test.
    let the_surf = BRepTool::face_surface(f)
        .ok_or_else(|| "BRepOffset_Tool::BuildPCurves: face without a surface".to_string())?;
    let basis_of_offset = the_surf.offset_basis_surface();
    let typ_s: &dyn Surface = basis_of_offset.as_deref().unwrap_or(the_surf.as_ref());

    if typ_s.is_bezier_surface() || typ_s.is_bspline_surface() {
        // `cxx:350-428`.
        let fpoint = a_c.d0(a_c.first_parameter());
        let lpoint = a_c.d0(a_c.last_parameter());
        let the_face = make_face_from_surface(the_surf.clone(), CONFUSION);
        let mut u1 = 0.0;
        let mut u2 = 0.0;
        const TOL_PROJ: f64 = 1.0e-4;
        let mut the_edge: Option<Edge> = None;
        for shape in explorer_shapes(&the_face.0, ShapeType::Edge) {
            let an_edge = Edge(shape);
            let Some(a_curve) = edge_adaptor(&an_edge) else {
                continue;
            };
            // `Extrema_ExtPC fextr(fpoint, aCurve)` (`cxx:363`).
            let fextr = occt_geom::extrema_pc::extrema_ext_pc_range(
                &a_curve,
                &fpoint,
                a_curve.first_parameter(),
                a_curve.last_parameter(),
            );
            if fextr.is_empty() {
                continue;
            }
            let mut dist2min = f64::MAX;
            for sol in &fextr {
                let d2 = sol.point.square_distance(&fpoint);
                if d2 < dist2min {
                    dist2min = d2;
                    u1 = sol.u;
                }
            }
            if dist2min > TOL_PROJ * TOL_PROJ {
                continue;
            }
            // `Extrema_ExtPC lextr(lpoint, aCurve)` (`cxx:381`).
            let lextr = occt_geom::extrema_pc::extrema_ext_pc_range(
                &a_curve,
                &lpoint,
                a_curve.first_parameter(),
                a_curve.last_parameter(),
            );
            if lextr.is_empty() {
                continue;
            }
            dist2min = f64::MAX;
            for sol in &lextr {
                let d2 = sol.point.square_distance(&lpoint);
                if d2 < dist2min {
                    dist2min = d2;
                    u2 = sol.u;
                }
            }
            if dist2min <= TOL_PROJ * TOL_PROJ {
                the_edge = Some(an_edge);
                break;
            }
        }

        if let Some(the_edge) = the_edge {
            // `cxx:406-427`.
            if u2 < u1 {
                std::mem::swap(&mut u1, &mut u2);
            }
            if let Some((pc, _cf, _cl)) = boptools_2d::curve_on_surface_range(&the_edge, &the_face)
            {
                let mut c2d: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(pc, u1, u2));
                if the_surf.is_u_periodic() || the_surf.is_v_periodic() {
                    c2d = put_in_bounds(f, e, c2d);
                }
                let reg = GeometryRegistry::global();
                reg.set_edge_pcurve(&e.0, GeometryRegistry::shape_key(&f.0), c2d);
                same_range(&e.0, CONFUSION);
            }
            return Ok(());
        }
    }

    // `cxx:430-436`: `ProjLib_ProjectedCurve Proj(HS, HC, Tolerance)`.
    let hc = a_c
        .curve()
        .cloned()
        .ok_or_else(|| "BRepOffset_Tool::BuildPCurves: empty curve adaptor".to_string())?;
    let hs = a_s
        .surface()
        .cloned()
        .ok_or_else(|| "BRepOffset_Tool::BuildPCurves: empty surface adaptor".to_string())?;
    let proj = ProjectedCurve::with_surface_curve_tol(hs, hc, TOLERANCE);

    // `cxx:437-466`.
    let mut c2d: Option<Arc<dyn Curve2d>> = match proj.get_type() {
        CurveType::Line => proj.line().ok().map(|l| {
            Arc::new(Geom2dLine::from_pnt_dir(l.location(), *l.direction())) as Arc<dyn Curve2d>
        }),
        CurveType::Circle => proj
            .circle()
            .ok()
            .map(|c| Arc::new(Geom2dCircle::new(c)) as _),
        CurveType::Ellipse => proj
            .ellipse()
            .ok()
            .map(|c| Arc::new(Geom2dEllipse::new(c)) as _),
        CurveType::Parabola => proj
            .parabola()
            .ok()
            .map(|c| Arc::new(Geom2dParabola::new(c)) as _),
        CurveType::Hyperbola => proj
            .hyperbola()
            .ok()
            .map(|c| Arc::new(Geom2dHyperbola::new(c)) as _),
        CurveType::BezierCurve => proj.bezier().map(|c| c as Arc<dyn Curve2d>),
        CurveType::BSplineCurve => proj.bspline().map(|c| c as Arc<dyn Curve2d>),
        // `cxx:465-466` `default: break` leaves `C2d` null.
        _ => None,
    };

    // `cxx:469-471`: `AS.IsUPeriodic() || AS.IsVPeriodic()`.
    let periodic = a_s
        .surface()
        .map(|s| s.is_u_periodic() || s.is_v_periodic())
        .unwrap_or(false);
    if periodic {
        if let Some(c) = c2d.take() {
            c2d = Some(put_in_bounds(f, e, c));
        }
    }

    // `cxx:472-481`.
    match c2d {
        Some(c2d) => {
            GeometryRegistry::global().set_edge_pcurve(
                &e.0,
                GeometryRegistry::shape_key(&f.0),
                c2d,
            );
            Ok(())
        }
        None => Err("Standard_ConstructionError: BRepOffset_Tool::BuildPCurves".to_string()),
    }
}

/// `FindPeriod(F, umin, umax, vmin, vmax)` (`cxx:155-190`): the 2D box of every
/// pcurve of `F`, sampled at 2 parameters for lines and 20 otherwise.
///
/// An edge without a pcurve returns `None`, which is OCCT's early `return`
/// (`cxx:167-170`) - the caller's four out-parameters then keep whatever they
/// held, so the port keeps them unset.
fn find_period(f: &Face) -> Option<(f64, f64, f64, f64)> {
    let mut b = occt_core::bnd::BndBox2d::new();
    for shape in explorer_shapes(&f.0, ShapeType::Edge) {
        let e = Edge(shape);
        let (c, pf, pl) = boptools_2d::curve_on_surface_range(&e, f)?;
        let mut nbp = 20.0;
        if curve2d_kind(c.as_ref()) == GeomAbsCurveType::Line {
            nbp = 2.0;
        }
        let step = (pl - pf) / nbp;
        // `PC.D0(pf, P); B.Add(P)` then `for (i = 2; i < nbp; i++)` and the
        // final `PC.D0(pl, P)` (`cxx:177-187`).
        let mut u = pf;
        b.add_point(&c.d0(u));
        let mut i = 2.0;
        while i < nbp {
            u += step;
            b.add_point(&c.d0(u));
            i += 1.0;
        }
        b.add_point(&c.d0(pl));
    }
    b.get()
}

/// `PutInBounds(F, E, C2d)` (`cxx:197-...`): translate the pcurve so that it
/// starts inside the UV box of `F`, when the surface is periodic in U or V.
///
/// Returns the translated curve: OCCT mutates the `handle(Geom2d_Curve)&` in
/// place (`cxx:242`, `:271`), which an `Arc` cannot do.
pub(super) fn put_in_bounds(f: &Face, e: &Edge, c2d: Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {    let (fpar, lpar) = BRepTool::edge_parameters(e);
    let mut s = match BRepTool::face_surface(f) {
        Some(s) => s,
        None => return c2d,
    };
    // `cxx:206-209`.
    if let Some(basis) = s.rectangular_trimmed_basis() {
        s = basis;
    }
    // `cxx:213-216`.
    if !s.is_u_periodic() && !s.is_v_periodic() {
        return c2d;
    }
    // `cxx:218`: `FindPeriod` leaves the box unset when it returns early; the
    // port answers the surface bounds in that case, which is the only reading
    // that keeps the two branches below defined.
    let (umin, umax, vmin, vmax) = find_period(f).unwrap_or_else(|| {
        let (u1, u2) = s.u_range();
        let (v1, v2) = s.v_range();
        (u1, u2, v1, v2)
    });

    let mut c2d = c2d;
    // `cxx:220-264`: the U pass.
    if s.is_u_periodic() {
        let period = s.u_period();
        let eps = period * 1.0e-6;
        let pf = c2d.d0(fpar);
        let pl = c2d.d0(lpar);
        let pm = c2d.d0(0.34 * fpar + 0.66 * lpar);
        let mut min_c = pf.x().min(pl.x()).min(pm.x());
        let mut max_c = pf.x().max(pl.x()).max(pm.x());
        let mut du = 0.0;
        if min_c < umin - eps {
            du = (((umin - min_c) / period).trunc() + 1.0) * period;
            // `int((umin - minC) / period) + 1` truncates toward zero; the
            // expression is positive here so `trunc` is the C cast.
        }
        if min_c > umax + eps {
            du = -(((min_c - umax) / period).trunc() + 1.0) * period;
        }
        if du != 0.0 {
            let mut t = GpTrsf2d::identity();
            t.set_translation_vec(&GpVec2d::new(du, 0.0));
            c2d = c2d.transformed(&t).into();
            min_c += du;
            max_c += du;
        }
        // `cxx:247-262`: `du` is the *same* variable as above, so when
        // `d2 >= d1` it keeps the first clamp's value and the curve is
        // translated a second time. That is the reference behaviour.
        if max_c > umax + 100.0 * eps {
            let d1 = max_c - umax;
            let d2 = umin - min_c + period;
            if d2 < d1 {
                du = -period;
            }
            if du != 0.0 {
                let mut t = GpTrsf2d::identity();
                t.set_translation_vec(&GpVec2d::new(du, 0.0));
                c2d = c2d.transformed(&t).into();
            }
        }
    }

    // `cxx:264-308`: the V pass, symmetric.
    if s.is_v_periodic() {
        let period = s.v_period();
        let eps = period * 1.0e-6;
        let pf = c2d.d0(fpar);
        let pl = c2d.d0(lpar);
        let pm = c2d.d0(0.34 * fpar + 0.66 * lpar);
        let mut min_c = pf.y().min(pl.y()).min(pm.y());
        let mut max_c = pf.y().max(pl.y()).max(pm.y());
        let mut dv = 0.0;
        if min_c < vmin - eps {
            dv = (((vmin - min_c) / period).trunc() + 1.0) * period;
        }
        if min_c > vmax + eps {
            dv = -(((min_c - vmax) / period).trunc() + 1.0) * period;
        }
        if dv != 0.0 {
            let mut t = GpTrsf2d::identity();
            t.set_translation_vec(&GpVec2d::new(0.0, dv));
            c2d = c2d.transformed(&t).into();
            min_c += dv;
            max_c += dv;
        }
        if max_c > vmax + 100.0 * eps {
            let d1 = max_c - vmax;
            let d2 = vmin - min_c + period;
            if d2 < d1 {
                dv = -period;
            }
            if dv != 0.0 {
                let mut t = GpTrsf2d::identity();
                t.set_translation_vec(&GpVec2d::new(0.0, dv));
                c2d = c2d.transformed(&t).into();
            }
        }
    }
    c2d
}

/// `BRepAdaptor_Curve` over a `[f, l]` window of `c` - the `Geom_TrimmedCurve`
/// view `TryProject` builds at `cxx:1979` (`C = new Geom_TrimmedCurve(C, f, l)`).
pub(super) fn trimmed_view(c: Arc<dyn Curve>, f: f64, l: f64) -> Arc<dyn Curve> {
    Arc::new(GeomTrimmedCurveBasis::new(c, f, l))
}
