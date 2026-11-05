//! The two `Reduce` statics of `ChFi3d_Builder_C2.cxx:63-99`.
//!
//! Both cut a parameter range in two halves at the 49 % / 51 % split:
//!
//! - the surface overload splits the `v` range and gives the low half to the
//!   first adaptor and the high half to the second (swapped when `p1 >= p2`),
//! - the curve overload keeps the half of `[First, Last]` that `p1` falls in.
//!
//! The split is a fixed 49/51 of the *whole* basis range, not of the current
//! trim: `surf->Bounds` / `c->FirstParameter()` are read from the basis
//! surface / curve, exactly like OCCT.

use std::sync::Arc;

use crate::topopebrep::tool::adaptor::{GeomAdaptorCurve, SurfaceAdaptor};

/// `Reduce(p1, p2, hs1, hs2)` (`cxx:63-84`).
///
/// Both adaptors must carry a surface; a null one is the OCCT null dereference
/// and leaves both adaptors untouched.
pub fn reduce_surface(p1: f64, p2: f64, hs1: &mut SurfaceAdaptor, hs2: &mut SurfaceAdaptor) {
    // `const occ::handle<Geom_Surface>& surf = s1.Surface(); double ud, uf,
    // vd, vf; surf->Bounds(ud, uf, vd, vf);` (`cxx:69-71`).
    let Some(surf) = hs1.surface().cloned() else {
        return;
    };
    let (ud, uf) = surf.u_range();
    let (vd, vf) = surf.v_range();
    // `double milmoins = 0.51 * vd + 0.49 * vf, milplus = 0.49 * vd + 0.51 *
    // vf;` (`cxx:72`).
    let milmoins = 0.51 * vd + 0.49 * vf;
    let milplus = 0.49 * vd + 0.51 * vf;
    // `if (p1 < p2) { s1.Load(surf, ud, uf, vd, milmoins); s2.Load(surf, ud,
    // uf, milplus, vf); } else { ... }` (`cxx:73-83`).
    if p1 < p2 {
        *hs1 = SurfaceAdaptor::of_surface_with_range(Arc::clone(&surf), ud, uf, vd, milmoins);
        *hs2 = SurfaceAdaptor::of_surface_with_range(Arc::clone(&surf), ud, uf, milplus, vf);
    } else {
        *hs1 = SurfaceAdaptor::of_surface_with_range(Arc::clone(&surf), ud, uf, milplus, vf);
        *hs2 = SurfaceAdaptor::of_surface_with_range(surf, ud, uf, vd, milmoins);
    }
}

/// `Reduce(p1, p2, hc)` (`cxx:86-99`).
///
/// A null curve is the OCCT null dereference and leaves the adaptor untouched.
pub fn reduce_curve(p1: f64, p2: f64, hc: &mut GeomAdaptorCurve) {
    // `GeomAdaptor_Curve& c = *hc; double f = c.FirstParameter(); double l =
    // c.LastParameter();` (`cxx:87-89`).
    let Some(curve) = hc.curve().cloned() else {
        return;
    };
    let f = hc.first_parameter();
    let l = hc.last_parameter();
    // `double milmoins = 0.51 * f + 0.49 * l, milplus = 0.49 * f + 0.51 * l;`
    // (`cxx:90`).
    let milmoins = 0.51 * f + 0.49 * l;
    let milplus = 0.49 * f + 0.51 * l;
    // `if (p1 < p2) c.Load(c.Curve(), f, milmoins); else c.Load(c.Curve(),
    // milplus, l);` (`cxx:91-98`).
    let (first, last) = if p1 < p2 {
        (f, milmoins)
    } else {
        (milplus, l)
    };
    *hc = GeomAdaptorCurve::with_range(curve, first, last);
}
