//! The `BRepFill_MultiLine(Face1, Face2, Edge1, Edge2, Inv1, Inv2, Bissec)`
//! constructor.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_MultiLine.cxx:78-391`.
//!
//! The body collects the (u, v) box of each face from the p-curve points of
//! its edges, picks the two iso lines of the surface that carry the bisector
//! (`UU1` / `VV1`, `UU2` / `VV2`), projects them to 2D so that a parameter of
//! the bisector can be matched against them (`cxx:227-235`), and finally
//! decides whether the configuration is a particular case (`cxx:358-390`).
//!
//! `BRep_Tool::Surface(F, L)` followed by `S->Transformed(L)` (`cxx:146-149`,
//! `cxx:263-267`) is [`BRepTool::face_surface_world`]: the location of a face
//! only moves its surface, so the parameter values of `BRep_Tool::UVPoints`
//! stay valid.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln};
use occt_core::kernel::geomabs::{CurveType, SurfaceType};
use occt_core::{elib, precision};
use occt_geom::curve::Curve;
use occt_geom::geom_proj_lib;
use occt_geom::plane::GeomPlane;
use occt_geom::projlib;
use occt_geom::surface::Surface;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::topo_tools::explorer_shapes;

use super::value_on_face::value_on_face;
use super::{Error, MultiLine};

/// `static bool isIsoU(const TopoDS_Face&, const TopoDS_Edge&)`
/// (`cxx:52-68`): the p-curve leaves its start point mostly along V, so the
/// edge is a `U = const` iso of the face.
fn is_iso_u(face: &Face, edge: &Edge) -> Result<bool, Error> {
    let (curve, first, _last) = crate::boptools_2d::curve_on_surface_range(edge, face)
        .ok_or(Error::EdgeWithoutPCurve)?;
    // `gp_Dir2d D = C->DN(f, 1);` (`cxx:65`) — `gp_Dir2d` from a `gp_Vec2d`
    // normalizes, which does not change the comparison below.
    let d = curve.d1(first).1;
    Ok(d.x().abs() < d.y().abs())
}

/// The (u, v) box of a face, taken from `BRep_Tool::UVPoints` of every edge
/// (`cxx:114-141`, `cxx:236-260`).
fn face_uv_bounds(face: &Face) -> Option<(f64, f64, f64, f64)> {
    let mut first = true;
    let (mut umin, mut umax, mut vmin, mut vmax) = (0.0, 0.0, 0.0, 0.0);
    for e in explorer_shapes(&face.0, ShapeType::Edge) {
        let edge = Edge::wrap(e)?;
        let (p1, p2) = BRepTool::uv_points(&edge, face);
        if first {
            first = false;
            umin = p1.x().min(p2.x());
            umax = p1.x().max(p2.x());
            vmin = p1.y().min(p2.y());
            vmax = p1.y().max(p2.y());
        } else {
            umin = umin.min(p1.x().min(p2.x()));
            umax = umax.max(p1.x().max(p2.x()));
            vmin = vmin.min(p1.y().min(p2.y()));
            vmax = vmax.max(p1.y().max(p2.y()));
        }
    }
    if first {
        None
    } else {
        Some((umin, umax, vmin, vmax))
    }
}

/// `GeomAdaptor_Curve Dummy(UU1); Dummy.GetType() == GeomAbs_Circle &&
/// Dummy.Circle().Radius() < eps` (`cxx:156-159` and the seven other copies):
/// the iso collapses to the degenerate pole of a cone / sphere, so the other
/// iso bound is used instead.
fn is_small_circle(iso: &Arc<dyn Curve>, eps: f64) -> bool {
    projlib::curve_type(iso.as_ref()) == CurveType::Circle
        && iso.gp_circ().map(|c| c.radius() < eps).unwrap_or(false)
}

/// The `UU` / `VV` pair of one face (`cxx:151-173` when `is_iso_u`, the
/// `cxx:174-193` mirror otherwise): `UU` is the iso that runs along the edge
/// (`UIso` for an iso-U edge) and `VV` the crossing one.
fn select_isos(
    s: &Arc<dyn Surface>,
    is_iso_u: bool,
    degenerated: bool,
    inv: bool,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    eps: f64,
) -> Result<(Arc<dyn Curve>, Arc<dyn Curve>), Error> {
    let lower = !degenerated && !inv;
    let (uu, vv) = if is_iso_u {
        // `cxx:151-173`.
        if lower {
            let mut uu = s.u_iso_curve(umin).ok_or(Error::NullCurve)?;
            if is_small_circle(&uu, eps) {
                uu = s.u_iso_curve(umax).ok_or(Error::NullCurve)?;
            }
            (uu, s.v_iso_curve(vmin).ok_or(Error::NullCurve)?)
        } else {
            let mut uu = s.u_iso_curve(umax).ok_or(Error::NullCurve)?;
            if is_small_circle(&uu, eps) {
                uu = s.u_iso_curve(umin).ok_or(Error::NullCurve)?;
            }
            (uu, s.v_iso_curve(vmin).ok_or(Error::NullCurve)?)
        }
    } else {
        // `cxx:174-193`: the roles are swapped — `UU` is a `VIso` and `VV` a
        // `UIso`.
        if lower {
            let mut uu = s.v_iso_curve(vmin).ok_or(Error::NullCurve)?;
            if is_small_circle(&uu, eps) {
                uu = s.v_iso_curve(vmax).ok_or(Error::NullCurve)?;
            }
            (uu, s.u_iso_curve(umin).ok_or(Error::NullCurve)?)
        } else {
            let mut uu = s.v_iso_curve(vmax).ok_or(Error::NullCurve)?;
            if is_small_circle(&uu, eps) {
                uu = s.v_iso_curve(vmin).ok_or(Error::NullCurve)?;
            }
            (uu, s.u_iso_curve(umin).ok_or(Error::NullCurve)?)
        }
    };
    Ok((uu, vv))
}

/// `BRepFill_MultiLine(Face1, Face2, Edge1, Edge2, Inv1, Inv2, Bissec)`
/// (`cxx:78-391`).
pub(super) fn build(
    face1: Face,
    face2: Face,
    edge1: &Edge,
    edge2: &Edge,
    inv1: bool,
    inv2: bool,
    bissec: Arc<dyn Curve2d>,
) -> Result<MultiLine, Error> {
    // `constexpr double mult = 5.; constexpr double eps = mult *
    // Precision::Confusion();` (`cxx:90-92`).
    const MULT: f64 = 5.0;
    let eps = MULT * precision::CONFUSION;

    let iso_u1 = is_iso_u(&face1, edge1)?;
    let iso_u2 = is_iso_u(&face2, edge2)?;

    // `BasisPlane = new Geom_Plane(0., 0., 1., 0.)` (`cxx:105`): the XY plane.
    let basis_plane: Arc<dyn Surface> =
        Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));

    // ---- Face1 (`cxx:113-235`) -----------------------------------------
    let (mut umin, mut umax, mut vmin, mut vmax) =
        face_uv_bounds(&face1).ok_or(Error::NullSurface)?;
    let s1 = BRepTool::face_surface_world(&face1).ok_or(Error::NullSurface)?;
    let degenerated1 = BRepTool::is_degenerated(edge1);
    let (uu1, vv1) =
        select_isos(&s1, iso_u1, degenerated1, inv1, umin, umax, vmin, vmax, eps)?;

    // `cxx:196-203`: for an iso-U edge the roles of the two bounds swap.
    if iso_u1 {
        std::mem::swap(&mut umin, &mut vmin);
        std::mem::swap(&mut umax, &mut vmax);
    }
    // `cxx:208-225`: bring the projected range into the iso's period.
    if uu1.is_periodic() {
        elib::clib::adjust_periodic(
            uu1.first_parameter(),
            uu1.last_parameter(),
            precision::PCONFUSION,
            &mut umin,
            &mut umax,
        );
    }
    if vv1.is_periodic() {
        elib::clib::adjust_periodic(
            vv1.first_parameter(),
            vv1.last_parameter(),
            precision::PCONFUSION,
            &mut vmin,
            &mut vmax,
        );
    }

    // `myU1 = Geom2dAdaptor_Curve(GeomProjLib::Curve2d(UU1, BasisPlane), Umin,
    // Umax)` (`cxx:227`).
    let u1 = Geom2dAdaptorCurve::with_range(
        geom_proj_lib::curve_2d(&uu1, &basis_plane).ok_or(Error::NullProjection)?,
        umin,
        umax,
    );
    // `UU1->D1(Umin, P, DZ); RefPlane = new Geom_Plane(gp_Ax3(P, DZ, gp::DZ()));`
    // (`cxx:229-230`).
    let (p, dz) = uu1.d1(umin);
    let dz = GpDir::from_vec(&dz).map_err(|_| Error::ConstructionError)?;
    let x_dir = GpDir::new(0.0, 0.0, 1.0).expect("gp::DZ");
    let ref_plane: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(
        GpAx3::from_point_dir_x(p, dz, &x_dir).map_err(|_| Error::ConstructionError)?,
    )));
    // `myV1 = Geom2dAdaptor_Curve(GeomProjLib::Curve2d(VV1, RefPlane), Vmin,
    // Vmax)` (`cxx:232`).
    let v1 = Geom2dAdaptorCurve::with_range(
        geom_proj_lib::curve_2d(&vv1, &ref_plane).ok_or(Error::NullProjection)?,
        vmin,
        vmax,
    );

    // ---- Face2 (`cxx:236-350`) -----------------------------------------
    let (mut umin, mut umax, mut vmin, mut vmax) =
        face_uv_bounds(&face2).ok_or(Error::NullSurface)?;
    let s2 = BRepTool::face_surface_world(&face2).ok_or(Error::NullSurface)?;
    let degenerated2 = BRepTool::is_degenerated(edge2);
    let (uu2, vv2) =
        select_isos(&s2, iso_u2, degenerated2, inv2, umin, umax, vmin, vmax, eps)?;

    if iso_u2 {
        std::mem::swap(&mut umin, &mut vmin);
        std::mem::swap(&mut umax, &mut vmax);
    }
    if uu2.is_periodic() {
        elib::clib::adjust_periodic(
            uu2.first_parameter(),
            uu2.last_parameter(),
            precision::PCONFUSION,
            &mut umin,
            &mut umax,
        );
    }
    if vv2.is_periodic() {
        elib::clib::adjust_periodic(
            vv2.first_parameter(),
            vv2.last_parameter(),
            precision::PCONFUSION,
            &mut vmin,
            &mut vmax,
        );
    }

    let u2 = Geom2dAdaptorCurve::with_range(
        geom_proj_lib::curve_2d(&uu2, &basis_plane).ok_or(Error::NullProjection)?,
        umin,
        umax,
    );
    let (p, dz) = uu2.d1(umin);
    let dz = GpDir::from_vec(&dz).map_err(|_| Error::ConstructionError)?;
    let ref_plane: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(
        GpAx3::from_point_dir_x(p, dz, &x_dir).map_err(|_| Error::ConstructionError)?,
    )));
    let v2 = Geom2dAdaptorCurve::with_range(
        geom_proj_lib::curve_2d(&vv2, &ref_plane).ok_or(Error::NullProjection)?,
        vmin,
        vmax,
    );

    let mut m = MultiLine {
        face1,
        face2,
        u1,
        v1,
        u2,
        v2,
        iso_u1,
        iso_u2,
        bis: Geom2dAdaptorCurve::from_curve(bissec),
        k_part: 0,
        cont: occt_core::kernel::geomabs::Shape::C0,
    };

    // `cxx:352-390`: a straight bisector that stays on an iso of `Face1`
    // (`myKPart = 1`, continuity `G1`), and the two-plane case (`myKPart = 2`).
    if m.bis.get_type() == Some(CurveType::Line) {
        let delta_u = m.bis.last_parameter() - m.bis.first_parameter();
        let p1 = value_on_face(
            m.bis.first_parameter() + 0.1 * delta_u,
            &m.bis,
            &m.u1,
            &m.v1,
            m.iso_u1,
        )?;
        let p2 = value_on_face(
            m.bis.first_parameter() + 0.9 * delta_u,
            &m.bis,
            &m.u1,
            &m.v1,
            m.iso_u1,
        )?;
        if m.iso_u1 {
            if (p1.y() - p2.y()).abs() < eps {
                m.k_part = 1;
            }
        } else if (p1.x() - p2.x()).abs() < eps {
            m.k_part = 1;
        }
        if m.k_part == 1 {
            m.cont = occt_core::kernel::geomabs::Shape::G1;
        }
        // `GeomAdaptor_Surface GAS1(S) / GAS2(S)` (`cxx:205-206`, `cxx:323-324`).
        let type1 = projlib::surface_type(s1.as_ref());
        let type2 = projlib::surface_type(s2.as_ref());
        if type1 == SurfaceType::Plane && type2 == SurfaceType::Plane {
            m.k_part = 2;
        }
    }
    Ok(m)
}
