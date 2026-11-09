//! `Recadre` (`Contap_Contour.cxx:229-281`).

use std::f64::consts::TAU;

use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `Recadre(myHS1, u1, v1)` (`cxx:229-281`).
///
/// The period is `M_PI + M_PI` - the source comment keeps `UPeriod()` out
/// (`cxx:256`, `cxx:269`), i.e. the period is assumed to be `2*pi` for every
/// periodic elementary surface.
pub(crate) fn recadre(surf: &SurfaceAdaptor, u1: &mut f64, v1: &mut f64) {
    let (is_u_periodic, is_v_periodic) = match surf.get_type() {
        GeomAbsSurfaceType::Cylinder
        | GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere => (true, false),
        GeomAbsSurfaceType::Torus => (true, true),
        _ => (false, false),
    };

    if is_u_periodic {
        let lmf = TAU;
        let f = surf.first_u_parameter();
        let l = surf.last_u_parameter();
        while *u1 < f {
            *u1 += lmf;
        }
        while *u1 > l {
            *u1 -= lmf;
        }
    }
    if is_v_periodic {
        let lmf = TAU;
        let f = surf.first_v_parameter();
        let l = surf.last_v_parameter();
        while *v1 < f {
            *v1 += lmf;
        }
        while *v1 > l {
            *v1 -= lmf;
        }
    }
}
