//! `BRepCheck::PrecCurve` (`BRepCheck.cxx:70-99`) and
//! `BRepCheck::PrecSurface` (`BRepCheck.cxx:103-130`).
//!
//! Both answer `RealEpsilon()` for every adaptor type but the ellipse (curve)
//! and the cone (surface), where they answer the largest `Epsilon()`
//! (`Standard_Real.hxx:242-246`) of the characteristic quantities instead.
//! `BRepLib_ValidateEdge::correctTolerance` (`BRepLib_ValidateEdge.cxx:68-77`)
//! adds `max(PrecCurve, PrecSurface)` to the tolerance it checks.

use occt_core::precision::{epsilon, REAL_EPSILON};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;

/// `BRepCheck::PrecCurve` (`BRepCheck.cxx:70-99`).
pub fn prec_curve(a_ac3d: &dyn Curve) -> f64 {
    let mut x_emax = REAL_EPSILON;
    // `if (aCT == GeomAbs_Ellipse)` (`cxx:74-96`).
    if let Some(el) = a_ac3d.gp_ellipse() {
        // `aEL3D.Location().Coord(aX[0], aX[1], aX[2])` (`cxx:80`).
        let loc = el.location();
        // `aX[3] = MajorRadius(); aX[4] = MinorRadius();` (`cxx:81-82`).
        let a_x = [
            loc.x(),
            loc.y(),
            loc.z(),
            el.major_radius(),
            el.minor_radius(),
        ];
        x_emax = -1.0;
        for v in a_x {
            // `if (aX[i] < 0.) aX[i] = -aX[i];` (`cxx:86-89`).
            let a = v.abs();
            // `double aXE = Epsilon(aX[i]); if (aXE > aXEmax) aXEmax = aXE;`
            // (`cxx:90-95`).
            let a_xe = epsilon(a);
            if a_xe > x_emax {
                x_emax = a_xe;
            }
        }
    }
    x_emax
}

/// `BRepCheck::PrecSurface` (`BRepCheck.cxx:103-130`).
pub fn prec_surface(a_ah_surf: &dyn Surface) -> f64 {
    let mut x_emax = REAL_EPSILON;
    // `if (aST == GeomAbs_Cone)` (`cxx:108-128`).
    if let Some(cone) = a_ah_surf.gp_cone() {
        // `aCone.Location().Coord(aX[0], aX[1], aX[2])` (`cxx:113`).
        let loc = cone.location();
        // `aX[3] = aCone.RefRadius();` (`cxx:114`).
        let a_x = [loc.x(), loc.y(), loc.z(), cone.radius()];
        x_emax = -1.0;
        for v in a_x {
            let a = v.abs();
            let a_xe = epsilon(a);
            if a_xe > x_emax {
                x_emax = a_xe;
            }
        }
    }
    x_emax
}
