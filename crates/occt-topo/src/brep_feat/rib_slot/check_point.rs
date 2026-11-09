//! `BRepFeat_RibSlot::CheckPoint` / `Normal` / `IntPar`
//! (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:443-576`).
//!
//! The three helpers answer the material side of the profile (`CheckPoint`,
//! `Normal`) and the parameter of a point on a profile curve (`IntPar`);
//! `MakeLinearForm` / `MakeRevolutionForm` call them while they build the tool.

use occt_core::cslib::{normal_d1, CSLibDerivativeStatus};
use occt_core::elib::{clib, slib};
use occt_core::gp::{GpDir, GpPln, GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};

/// `CheckPoint(e, bnd, Pln)` (`cxx:443-476`): a point pulled off the middle of
/// `e` towards the material side of `pln`.
///
/// `D = -(T x N) / 10` where `T` is the edge tangent at the mid parameter
/// (reversed for a reversed edge) and `N` the plane normal; `bnd` is unused by
/// the reference.
pub fn check_point(e: &Edge, _bound: f64, pln: &GpPln) -> GpPnt {
    let (f, l) = BRepTool::edge_parameters(e);
    let par = (f + l) / 2.0;
    // `BRep_Tool::Curve(e, f, l)` (`cxx:458`). A null handle is a broken edge
    // in the reference, which dereferences it at `cxx:462`.
    let Some(cc) = BRepTool::edge_curve_world(e) else {
        return GpPnt::default();
    };
    let (mut pp, mut tgt) = cc.d1(par);
    if e.0.orientation() == Orientation::Reversed {
        tgt.reverse();
    }
    let n = pln.position().direction();
    let d = tgt
        .crossed(&GpVec::new(n.x(), n.y(), n.z()))
        .divided(-10.0);
    pp.translate_vec(&d);
    pp
}

/// `Normal(F, P)` (`cxx:480-529`): the normal of `F` at `P`.
///
/// The reference asks `ElSLib::Parameters` for the `(U, V)` of `P` on the
/// elementary surface, then `CSLib::Normal` on the first derivatives of the
/// adaptor. A non-elementary surface answers `gp_Dir::X` (`cxx:525-527`).
pub fn normal(f: &Face, p: &GpPnt) -> GpDir {
    let Some(surf) = BRepTool::face_surface_world(f) else {
        return GpDir::default_dir();
    };
    // `switch (AS.GetType())` (`cxx:494-517`).
    let (u, v) = if let Some(pln) = surf.gp_pln() {
        slib::plane_parameters(&pln.position(), p)
    } else if let Some(cyl) = surf.gp_cylinder() {
        slib::cylinder_parameters(&cyl.pos, p)
    } else if let Some(cone) = surf.gp_cone() {
        slib::cone_parameters(&cone.pos, cone.radius, cone.semi_angle, p)
    } else if let Some(torus) = surf.gp_torus() {
        slib::torus_parameters(&torus.pos, torus.major_radius, torus.minor_radius, p)
    } else {
        return GpDir::default_dir();
    };

    let (_, d1u, d1v) = surf.d1(u, v);
    let (status, dir) = normal_d1(&d1u, &d1v, CONFUSION);
    // `CSLib::Normal` leaves `N` untouched when `St != CSLib_Done`; the port
    // answers the same `X` default as the non-elementary arm above.
    let mut n = match (status, dir) {
        (CSLibDerivativeStatus::Done, Some(d)) => d,
        _ => GpDir::default_dir(),
    };
    if f.0.orientation() == Orientation::Forward {
        n.reverse();
    }
    n
}

/// `IntPar(C, P)` (`cxx:536-575`): the parameter of `P` on the elementary curve
/// `C`, `0` for a null handle or any other curve type.
pub fn int_par(c: Option<&dyn Curve>, p: &GpPnt) -> f64 {
    let Some(c) = c else {
        return 0.0;
    };
    // `switch (AC.GetType())` (`cxx:548-572`).
    if let Some(l) = c.gp_line() {
        clib::line_parameter(&l, p)
    } else if let Some(ci) = c.gp_circ() {
        clib::circle_parameter(&ci.position(), p)
    } else if let Some(e) = c.gp_ellipse() {
        clib::ellipse_parameter(e.position(), e.major_radius(), e.minor_radius(), p)
    } else if let Some(h) = c.gp_hyperbola() {
        clib::hyperbola_parameter(h.position(), h.major_radius, h.minor_radius, p)
    } else if let Some(pa) = c.gp_parabola() {
        clib::parabola_parameter(pa.position(), p)
    } else {
        0.0
    }
}
