//! `LocOpe::TgtFaces` (`LocOpe.cxx:130-193`).
//!
//! Two faces are tangent along an edge when the sampled angle between their
//! unit normals never exceeds `0.0001` rad. Both p-curves are evaluated at the
//! *3D* curve parameters of the edge, which the reference takes for granted
//! (`cxx:149-175`) - the same assumption holds here because every p-curve of an
//! edge shares that parameterisation.
//!
//! Not ported from the reference: `BRepAdaptor_Surface bs(F1, false)`
//! (`cxx:132`) is constructed and never read.

use std::f64::consts::PI;

use crate::abs::Orientation;
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `LocOpe::TgtFaces(E, F1, F2)` (`cxx:130-193`).
///
/// A missing p-curve on either face is the null `myCurve` of
/// `BRepAdaptor_Curve2d` (`BRepAdaptor_Curve2d.cxx:44-58`), which the reference
/// dereferences at `cxx:166` / `cxx:172`; the port answers `false`.
pub fn tgt_faces(e: &Edge, f1: &Face, f2: &Face) -> bool {
    // `cxx:138-146`.
    let surf1 = SurfaceAdaptor::of_face(f1);
    let surf2 = SurfaceAdaptor::of_face(f2);
    let Some((curve1, _f1, _l1)) = boptools_2d::curve_on_surface_range(e, f1) else {
        return false;
    };
    let Some((curve2, _f2, _l2)) = boptools_2d::curve_on_surface_range(e, f2) else {
        return false;
    };

    // `cxx:140-147`.
    let mut e = e.clone();
    e.0.set_orientation(Orientation::Forward);
    let rev1 = f1.0.orientation() == Orientation::Reversed;
    let rev2 = f2.0.orientation() == Orientation::Reversed;

    // `cxx:148-153`: `ta` and the 1% end trim that keeps the samples off sharp
    // patch corners.
    let ta = 0.0001;
    let (mut f, mut l) = BRepTool::edge_parameters(&e);
    let eps = (l - f) / 100.0;
    f += eps;
    l -= eps;

    // `cxx:159-189`.
    let mut angmin = PI;
    let mut angmax = -PI;
    for i in 0..=20 {
        let u = f + (l - f) * i as f64 / 20.0;

        let p = curve1.d0(u);
        let Some((_, du, dv)) = surf1.d1(p.x(), p.y()) else {
            return false;
        };
        let mut d1 = du.crossed(&dv).normalized();
        if rev1 {
            d1.reverse();
        }

        let p = curve2.d0(u);
        let Some((_, du, dv)) = surf2.d1(p.x(), p.y()) else {
            return false;
        };
        let mut d2 = du.crossed(&dv).normalized();
        if rev2 {
            d2.reverse();
        }

        // `cxx:180-187`.
        let ang = d1.angle(&d2);
        if ang <= angmin {
            angmin = ang;
        }
        if ang >= angmax {
            angmax = ang;
        }
    }
    // `cxx:190`: only the maximum is tested.
    let _ = angmin;
    angmax <= ta
}
