use occt_core::bnd::BndBox;
use occt_geom::Curve;

use super::boxes::point_box_distance;
use super::curve_kind::GeomCurveKind;
use super::resolution::resolution_of;

// ---------------------------------------------------------------------------
// `FindParameters`
// ---------------------------------------------------------------------------

/// `IntTools_EdgeEdge::FindParameters` (`IntTools_EdgeEdge.cxx:553-671`).
///
/// Walks from each end of `[t1, a_t2]` towards the box, growing the step by
/// the local curve resolution, then bisects (golden section) for the boundary
/// parameter. `tb1`/`tb2` are the caller's *current* range values, which
/// OCCT passes by reference and only overwrites when a bounding point is
/// found. `None` reproduces OCCT's `false` return (the edge is out of the box).
#[allow(clippy::too_many_arguments)]
pub(in crate::edge_edge) fn find_parameters(
    c: &dyn Curve,
    kind: GeomCurveKind,
    t1: f64,
    t2: f64,
    tol: f64,
    res: f64,
    ptol: f64,
    res_coeff: f64,
    cbox: &BndBox,
    tb1_in: f64,
    tb2_in: f64,
) -> Option<(f64, f64)> {
    let cf = 0.618_033_988_749_894_8; // = 0.5*(1+sqrt(5))/2
    let mut cbx = *cbox;
    cbx.set_gap(cbx.gap() + tol);
    let max_dt = (t2 - t1) * 0.01;

    let mut tb1 = tb1_in;
    let mut tb2 = tb2_in;
    let mut b_ret = false;

    for i in 0..2 {
        let mut tb = if i == 0 { t1 } else { t2 };
        let mut t = if i == 0 { t2 } else { tb1 };
        let c_val = if i == 0 { 1.0 } else { -1.0 };
        let mut dt = res;
        let mut distp = 0.0;
        b_ret = false;
        let mut k = 1.0;
        // Looking for the point on the edge which is in the box.
        while c_val * (t - tb) >= 0.0 {
            let p = c.d0(tb);
            let dist = point_box_distance(cbox, &p);
            if dist > tol {
                if distp > 0.0 {
                    let mut to_grow = false;
                    if (distp - dist).abs() / distp < 0.1 {
                        dt = resolution_of(c, kind, res_coeff, k * dist);
                        if dt < max_dt {
                            to_grow = true;
                            k *= 2.0;
                        }
                    }
                    if !to_grow {
                        k = 1.0;
                        dt = resolution_of(c, kind, res_coeff, dist);
                    }
                }
                tb += c_val * dt;
            } else {
                b_ret = true;
                break;
            }
            distp = dist;
        }
        if !b_ret {
            if i == 0 {
                // Edge is out of the box.
                return None;
            }
            b_ret = true;
            tb = tb1;
            dt = t2 - tb1;
        }
        t = if i == 0 { t1 } else { t2 };
        if tb != t {
            // One point IN, one point OUT: look for the bounding point.
            let mut tin = tb;
            let mut tout = tb - c_val * dt;
            let mut diff = tin - tout;
            while diff.abs() > ptol {
                tb = tout + diff * cf;
                let p = c.d0(tb);
                if cbx.is_out(&p) {
                    tout = tb;
                } else {
                    tin = tb;
                }
                diff = tin - tout;
            }
            if i == 0 {
                tb1 = tb;
            } else {
                tb2 = tb;
            }
        }
    }
    if b_ret {
        Some((tb1, tb2))
    } else {
        None
    }
}
