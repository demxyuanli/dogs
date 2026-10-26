//! `ChFi3d_Recale`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:487-551`
//! (OCCT 8.0.0).
//!
//! When two `gp_Pnt2d` are more than half a period apart on a periodic
//! surface, one of the two coordinates is shifted by a whole period so that the
//! pair follows the same branch. `refon1` selects which of the two is moved
//! (`p2` when set, `p1` otherwise), and the direction depends on which of the
//! two coordinates is the smaller - the four-arm chain of the source.
//!
//! Divergence that is not a heuristic: `Bs.GeomSurfaceOriginal()` followed by
//! `down_cast<Geom_RectangularTrimmedSurface>` is
//! [`occt_geom::surface::Surface::rectangular_trimmed_basis`]; the periodicity
//! and the period are read on that basis.

use occt_core::gp::GpPnt2d;

use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFi3d_Recale(Bs, p1, p2, refon1)` (`cxx:487-551`): folds `p1` / `p2` onto
/// one period of the surface of `Bs`.
///
/// A `Bs` with no surface loaded has no periodicity to fold against, so the two
/// points are left untouched (OCCT would fault on the null handle).
pub fn ch_fi3d_recale(bs: &SurfaceAdaptor, p1: &mut GpPnt2d, p2: &mut GpPnt2d, refon1: bool) {
    let Some(original) = bs.surface() else {
        return;
    };
    // `occ::handle<Geom_Surface> surf = Bs.GeomSurfaceOriginal(); ts =
    //  down_cast<Geom_RectangularTrimmedSurface>(surf);
    //  if (!ts.IsNull()) surf = ts->BasisSurface();`
    let surf = original
        .rectangular_trimmed_basis()
        .unwrap_or_else(|| std::sync::Arc::clone(original));

    if surf.is_u_periodic() {
        let mut u1 = p1.x();
        let mut u2 = p2.x();
        let uper = surf.u_period();
        if (u2 - u1).abs() > 0.5 * uper {
            if u2 < u1 && refon1 {
                u2 += uper;
            } else if u2 < u1 && !refon1 {
                u1 -= uper;
            } else if u1 < u2 && refon1 {
                u2 -= uper;
            } else if u1 < u2 && !refon1 {
                u1 += uper;
            }
        }
        p1.set_x(u1);
        p2.set_x(u2);
    }

    if surf.is_v_periodic() {
        let mut v1 = p1.y();
        let mut v2 = p2.y();
        let vper = surf.v_period();
        if (v2 - v1).abs() > 0.5 * vper {
            if v2 < v1 && refon1 {
                v2 += vper;
            } else if v2 < v1 && !refon1 {
                v1 -= vper;
            } else if v1 < v2 && refon1 {
                v2 -= vper;
            } else if v1 < v2 && !refon1 {
                v1 += vper;
            }
        }
        p1.set_y(v1);
        p2.set_y(v2);
    }
}
