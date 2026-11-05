//! The static helpers of `BRepFill_TrimSurfaceTool::IntersectWith`.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_TrimSurfaceTool.cxx:100-378`.
//!
//! * `Bubble` (`cxx:100-120`) — order the sequence of intersection points by
//!   increasing bisector parameter.
//! * `EvalPhase` (`cxx:122-152`) — the `M_PI` phase of a cone / surface of
//!   revolution whose `(0, 0)` point sits on the opposite side of the axis.
//! * `EvalParameters` (`cxx:155-378`) — intersect the projection of one edge
//!   on the `xOy` plane with the bisector (the regular branch), or, for a
//!   degenerated edge, answer the single point of its vertex on the bisector
//!   with the angle around the revolved surface's axis (the degenerated
//!   branch).
//!
//! The `OCCT_DEBUG` printouts of `cxx:254-258` are not ported.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPnt, GpPnt2d, GpVec};
use occt_core::kernel::geomabs::{CurveType, SurfaceType};
use occt_core::precision;
use occt_core::gp::GpPln;
use occt_geom::curve::Curve;
use occt_geom::geom_proj_lib;
use occt_geom::plane::GeomPlane;
use occt_geom::projlib;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom::Surface;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Vertex};
use crate::topopebrep::build::tools2d::first_vertex;

use super::{Error, PI, TOL};

/// `Bubble(Seq)` (`cxx:102-120`): order the sequence by increasing `X`
/// (ascending bisector parameter), exchanging neighbours until no swap
/// happens. The `for` bound is the **original** length (`cxx:107`), as in
/// OCCT.
fn bubble(seq: &mut [GpPnt]) {
    let nb_points = seq.len();
    let mut invert = true;
    while invert {
        invert = false;
        for i in 0..nb_points.saturating_sub(1) {
            let p1 = seq[i];
            let p2 = seq[i + 1];
            if p2.x() < p1.x() {
                seq.swap(i, i + 1);
                invert = true;
            }
        }
    }
}

/// `EvalPhase(Edge, Face, GAS, Axis)` (`cxx:124-152`).
///
/// `VDeg` is the `V` of the edge's first UV point; the first *other* edge of
/// the face supplies `V` (the coordinate which differs most from `VDeg`), and
/// the phase is `M_PI` when that `(0, V)` point projects negatively on the
/// axis' `XDirection`.
fn eval_phase(edge: &Edge, face: &Face, gas: &dyn Surface, axis: &GpAx3) -> f64 {
    let (pe1, _pe2) = BRepTool::uv_points(edge, face);
    let v_deg = pe1.y();
    let mut v = 0.0;
    for e in crate::topo_tools::explorer_shapes(&face.0, crate::abs::ShapeType::Edge) {
        let Some(other) = Edge::wrap(e) else {
            continue;
        };
        if crate::topo_tools::is_same(&other.0, &edge.0) {
            continue;
        }
        let (pf1, pf2) = BRepTool::uv_points(&other, face);
        v = if (pf1.y() - v_deg).abs() > (pf2.y() - v_deg).abs() {
            pf1.y()
        } else {
            pf2.y()
        };
        break;
    }
    let p = gas.d0(0.0, v);

    // `gp_Vec(Axis.Location(), P).Dot(Axis.XDirection()) < 0.` (`cxx:145`).
    let loc = axis.location();
    let xd = axis.x_direction();
    let dir = GpVec::from_pnts(&loc, &p).dot(&GpVec::new(xd.x(), xd.y(), xd.z()));
    if dir < 0.0 {
        PI
    } else {
        0.0
    }
}

/// `EvalParameters(Edge, Face, Bis, Seq)` (`cxx:155-378`).
///
/// Appends `gp_Pnt(bisector parameter, edge parameter, 0)` for every
/// intersection of the edge's `xOy` projection with the bisector. For a
/// degenerated edge it appends the single `gp_Pnt(bisector first parameter,
/// angle, 0)` of the edge's vertex.
pub fn eval_parameters(
    edge: &Edge,
    face: &Face,
    bis: &Arc<dyn Curve2d>,
    seq: &mut Vec<GpPnt>,
) -> Result<(), Error> {
    let degener = BRepTool::is_degenerated(edge);

    // `double Tol = 1.e-6; double TolC = 0.;` (`cxx:178-179`).
    let tol = TOL;
    let tol_c = 0.0;

    if !degener {
        // `cxx:184-188`: `CT = Geom_TrimmedCurve(C, f, l)` transformed by `L`,
        // then projected on `Geom_Plane(0, 0, 1, 0)`.
        let plane: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));
        let curve = BRepTool::edge_curve_world(edge).ok_or(Error::NullCurve)?;
        let (f, l) = BRepTool::edge_parameters(edge);
        let ct: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(curve, f, l));
        let c2d = geom_proj_lib::curve_2d(&ct, &plane).ok_or(Error::NullProjection)?;

        // `Geom2dAdaptor_Curve AC(C2d); Geom2dAdaptor_Curve ABis(Bis);`
        // (`cxx:189-190`).
        let mut ac = Geom2dAdaptorCurve::from_curve(c2d.clone());
        let mut abis = Geom2dAdaptorCurve::from_curve(bis.clone());

        // `Geom2dInt_GInter(ABis, AC, TolC, Tol)` (`cxx:192`).
        let mut inter = Geom2dIntGInter::new();
        perform(&mut inter, &abis, &ac, tol_c, tol);

        // `cxx:194-197`.
        if !inter.is_done() {
            return Err(Error::IntersectionNotDone);
        }

        let mut nb_points = inter.nb_points();

        if nb_points < 1 {
            // `cxx:203-235`: elongate both curves and enlarge the tolerance.
            let c_type = ac.get_type().unwrap_or(CurveType::OtherCurve);
            let bis_type = abis.get_type().unwrap_or(CurveType::OtherCurve);
            let can_elongate_c = !matches!(
                c_type,
                CurveType::BezierCurve
                    | CurveType::BSplineCurve
                    | CurveType::OffsetCurve
                    | CurveType::OtherCurve
            );
            let can_elongate_bis = !matches!(
                bis_type,
                CurveType::BezierCurve
                    | CurveType::BSplineCurve
                    | CurveType::OffsetCurve
                    | CurveType::OtherCurve
            );

            if can_elongate_c {
                if let Some(mut t_c2d) = c2d.trimmed_copy2d() {
                    t_c2d.set_trim_range(
                        t_c2d.first_parameter() - tol,
                        t_c2d.last_parameter() + tol,
                    );
                    ac.load(Arc::new(t_c2d));
                }
            }
            if can_elongate_bis {
                if let Some(mut t_bis) = bis.trimmed_copy2d() {
                    t_bis.set_trim_range(
                        t_bis.first_parameter() - tol,
                        t_bis.last_parameter() + tol,
                    );
                    abis.load(Arc::new(t_bis));
                }
            }
            // `Intersector = Geom2dInt_GInter(ABis, AC, TolC, Tol * 10);`
            // (`cxx:228`).
            perform(&mut inter, &abis, &ac, tol_c, tol * 10.0);

            // `cxx:230-233`.
            if !inter.is_done() {
                return Err(Error::IntersectionNotDone);
            }
            nb_points = inter.nb_points();
        }

        // `cxx:236-246`: one point per intersection.
        if nb_points > 0 {
            for i in 1..=nb_points {
                let u1 = inter.point(i).param_on_first();
                let u2 = inter.point(i).param_on_second();
                seq.push(GpPnt::new(u1, u2, 0.0));
            }
        }

        // `cxx:248-270`: one point at the middle of each overlapping segment.
        let nb_segments = inter.nb_segments();
        if nb_segments > 0 {
            for i in 1..=nb_segments {
                let seg = inter.segment(i);
                let mut u1 = seg.first_point().param_on_first();
                u1 += seg.last_point().param_on_first();
                u1 /= 2.0;
                let mut u2 = seg.first_point().param_on_second();
                u2 += seg.last_point().param_on_second();
                u2 /= 2.0;
                seq.push(GpPnt::new(u1, u2, 0.0));
            }
        }

        // `cxx:270-286`: order by increasing bisector parameter, then drop the
        // duplicates. The bound of the removal loop is `NbPoints` — the count
        // of the intersector, **not** the sequence length that grew with the
        // segments above (`cxx:276`).
        bubble(seq);
        let mut nb_points = nb_points as i64;
        let mut i: i64 = 1;
        while i < nb_points {
            let p1 = seq[(i - 1) as usize];
            let p2 = seq[i as usize];
            if p2.x() - p1.x() < tol {
                seq.remove((i - 1) as usize);
                i -= 1;
                nb_points -= 1;
            }
            i += 1;
        }
        return Ok(());
    }

    // `cxx:287-377`: the degenerated branch. The single vertex of the edge,
    // as `(x, y)` of its 3D point.
    let p3d = first_vertex(&edge.0)
        .map(|s| BRepTool::vertex_point_world(&Vertex(s)))
        .ok_or(Error::NullCurve)?;
    let p2d = GpPnt2d::new(p3d.x(), p3d.y());

    let mut ubis = bis.first_parameter();
    let mut pbis = bis.d0(ubis);

    // `cxx:300-319`: the first / last bisector parameter whose point sits on
    // the vertex, `Precision::IsPositiveInfinite` guarding the `Distance`.
    if super::is_positive_infinite(pbis.x())
        || super::is_positive_infinite(pbis.y())
        || pbis.distance(&p2d) > tol
    {
        ubis = bis.last_parameter();
        if ubis >= precision::INFINITE {
            return Ok(());
        }
        pbis = bis.d0(ubis);
        if pbis.distance(&p2d) > tol {
            return Ok(());
        }
    }

    // `cxx:321-326`: the surface and its adaptor type.
    let gs = BRepTool::face_surface(face).ok_or(Error::NullCurve)?;

    let axis: GpAx3;
    let mut phase = 0.0;

    match projlib::surface_type(gs.as_ref()) {
        // `cxx:330-332`.
        SurfaceType::Sphere => axis = gs.gp_sphere().ok_or(Error::NullCurve)?.pos,
        // `cxx:334-343`: the phase flips when the face is on the other side of
        // the apex than the `(0, 0)` point of the cone.
        SurfaceType::Cone => {
            axis = gs.gp_cone().ok_or(Error::NullCurve)?.pos;
            phase = eval_phase(edge, face, gs.as_ref(), &axis);
        }
        // `cxx:344-346`.
        SurfaceType::Torus => axis = gs.gp_torus().ok_or(Error::NullCurve)?.pos,
        // `cxx:347-349`.
        SurfaceType::Cylinder => axis = gs.gp_cylinder().ok_or(Error::NullCurve)?.pos,
        // `cxx:350-361`.
        SurfaceType::SurfaceOfRevolution => {
            // `GeomAdaptor_SurfaceOfRevolution ASRev(HC, GAS.AxeOfRevolution());
            // Axis = ASRev.Axis();` (`cxx:356-358`).
            let basis = gs.revolution_basis_curve().ok_or(Error::NullCurve)?;
            let axe = gs.revolution_axis().ok_or(Error::NullCurve)?;
            let (bf, bl) = (basis.first_parameter(), basis.last_parameter());
            axis = crate::brep_sweep::revolution_ax3(basis.as_ref(), bf, bl, &axe);
            phase = eval_phase(edge, face, gs.as_ref(), &axis);
        }
        // `cxx:361-362`.
        _ => return Err(Error::NotImplemented),
    }

    // `cxx:364-377`: `U = Axis.XDirection().AngleWithRef(D1, XDirection ^
    // YDirection) + Phase`, brought into `[0, 2 * M_PI)`.
    let d12d = bis.eval_dn(ubis, 1);
    let d1 = GpVec::new(d12d.x(), d12d.y(), 0.0);
    let xdir = *axis.x_direction();
    let ydir = *axis.y_direction();
    let xd = GpVec::new(xdir.x(), xdir.y(), xdir.z());
    let yd = GpVec::new(ydir.x(), ydir.y(), ydir.z());
    let vref = GpDir::from_vec(&xd.crossed(&yd)).map_err(|_| Error::ConstructionError)?;
    let other = GpDir::from_vec(&d1).map_err(|_| Error::ConstructionError)?;
    let mut u = xdir.angle_with_ref(&other, &vref);
    u += phase;
    if u < 0.0 {
        u += 2.0 * PI;
    }

    // `P = gp_Pnt(Bis->FirstParameter(), U, 0.)` (`cxx:374`): note the first
    // component is the bisector's **first** parameter, not `UBis`.
    seq.push(GpPnt::new(bis.first_parameter(), u, 0.0));
    Ok(())
}

/// `Geom2dInt_GInter(A1, A2, TolConf, Tol)` (`Geom2dInt_GInter_0.cxx:29-89`):
/// both domains are `ComputeDomain(adaptor, max(TolConf, Tol))` over the
/// adaptor's own `myFirst` / `myLast`, which `Load` leaves untouched
/// (`Geom2dAdaptor_Curve.cxx:250-258`).
fn perform(
    inter: &mut Geom2dIntGInter,
    abis: &Geom2dAdaptorCurve,
    ac: &Geom2dAdaptorCurve,
    tol_conf: f64,
    tol: f64,
) {
    let tol_domain = tol.max(tol_conf);
    let c1 = abis.curve().expect("loaded bisector");
    let c2 = ac.curve().expect("loaded edge curve");
    let d1 = Geom2dIntGInter::domain_over_range(
        c1.as_ref(),
        abis.first_parameter(),
        abis.last_parameter(),
        tol_domain,
    );
    let d2 = Geom2dIntGInter::domain_over_range(
        c2.as_ref(),
        ac.first_parameter(),
        ac.last_parameter(),
        tol_domain,
    );
    inter.perform(c1.as_ref(), &d1, c2.as_ref(), &d2, tol_conf, tol);
}
