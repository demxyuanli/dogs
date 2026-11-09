//! `GetOrientation` (`LocOpe_Gluer.cxx:354-451`): the file-local static
//! `LocOpe_Gluer::Bind` uses to decide `FUSE` from `CUT`.
//!
//! It samples the p-curve of every edge of the glued face at 21 parameters,
//! takes the surface normal there (`Sn->D1` frame, flipped when the face is
//! reversed), projects the sample point onto the basis surface with
//! `Extrema_ExtPS`, and compares the two normals. The first sample that yields
//! a usable pair decides: normals pointing the same way mean the glued face was
//! given `REVERSED` (`cxx:440-444`). No sample means `TopAbs_INTERNAL`
//! (`cxx:451`), which `LocOpe_Gluer::Bind` reads as "no decision".
//!
//! Deviations, all for a null shape the reference would dereference:
//! - a face without a surface (`Sn` / `Sb` null) answers `INTERNAL`;
//! - an edge without a p-curve on the glued face (`C2d` null, `cxx:386`)
//!   answers `INTERNAL` as well instead of crashing.
//!
//! `GeomAdaptor_Surface`'s parameter bounds are what `cxx:374-377` feeds
//! `Extrema_ExtPS::Initialize`: an unbounded direction is
//! `Precision::Infinite()`, not an IEEE infinity ([`adaptor_bound`]).

use occt_core::precision::{Precision, CONFUSION, INFINITE};
use occt_geom::approx_same_parameter;
use occt_geom::extrema_surf::{ExtPs, ExtremaExtFlag};

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::topo_tools::explorer_shapes;

/// `RealLast()` (`Standard_Real.hxx`): the `dist2min` seed of `cxx:417`.
const REAL_LAST: f64 = f64::MAX;

/// `GeomAdaptor_Surface::FirstUParameter()` and friends: an unbounded direction
/// is reported as `Precision::Infinite()` (`GeomAdaptor_Surface.cxx:402-431`
/// copies `Geom_Surface::FirstUParameter()`, which is
/// `-Precision::Infinite()` for an unbounded surface).
fn adaptor_bound(x: f64) -> f64 {
    if x.is_finite() {
        x
    } else if x > 0.0 {
        INFINITE
    } else {
        -INFINITE
    }
}

/// `GetOrientation(const TopoDS_Face& Fn, const TopoDS_Face& Fb)`
/// (`cxx:354-451`).
pub(super) fn get_orientation(f_new: &Face, f_base: &Face) -> Orientation {
    // `cxx:356-358`.
    let (Some(s_new), Some(s_base)) = (
        BRepTool::face_surface(f_new),
        BRepTool::face_surface(f_base),
    ) else {
        return Orientation::Internal;
    };

    // `cxx:369-380`: the projector is initialised once for the basis surface.
    let tol_u = approx_same_parameter::u_resolution(s_base.as_ref(), CONFUSION);
    let tol_v = approx_same_parameter::v_resolution(s_base.as_ref(), CONFUSION);
    let (u_first, u_last) = s_base.u_range();
    let (v_first, v_last) = s_base.v_range();
    let mut proj = ExtPs::new();
    proj.set_flag(ExtremaExtFlag::Min);
    proj.initialize(
        s_base.as_ref(),
        adaptor_bound(u_first),
        adaptor_bound(u_last),
        adaptor_bound(v_first),
        adaptor_bound(v_last),
        tol_u,
        tol_v,
    );

    // `cxx:382-450`.
    for raw in explorer_shapes(&f_new.0, ShapeType::Edge) {
        // `cxx:385-386`.
        let Some(edge) = Edge::wrap(raw) else {
            continue;
        };
        let Some((curve_2d, mut f, mut l)) = boptools_2d::curve_on_surface_range(&edge, f_new)
        else {
            return Orientation::Internal;
        };
        // `cxx:387-397`: a p-curve range may be infinite on an unbounded
        // surface; the reference samples a 200-wide window around the finite end.
        if Precision::is_negative_infinite(f) && Precision::is_positive_infinite(l) {
            f = -100.0;
            l = 100.0;
        } else if Precision::is_negative_infinite(f) {
            f = l - 200.0;
        } else if Precision::is_positive_infinite(l) {
            l = f + 200.0;
        }
        // `cxx:399-448`.
        let delta_u = (l - f) / 20.0;
        for i in 1..=21 {
            let pt_2d = curve_2d.d0(f + (i - 1) as f64 * delta_u);
            let (pt, d1u, d1v) = s_new.d1(pt_2d.x(), pt_2d.y());
            let mut n1 = d1u.crossed(&d1v);
            if n1.magnitude() <= CONFUSION {
                continue;
            }
            n1 = n1.normalized();
            // `cxx:409-412`.
            if f_new.0.orientation() == Orientation::Reversed {
                n1.reverse();
            }

            // `cxx:414-446`.
            proj.perform(&pt);
            if !proj.is_done() {
                continue;
            }
            let mut dist_2min = REAL_LAST;
            let mut jmin = 0usize;
            for j in 1..=proj.nb_ext() {
                let d = proj.square_distance(j);
                if d < dist_2min {
                    jmin = j;
                    dist_2min = d;
                }
            }
            if jmin == 0 {
                continue;
            }
            let (uu, vv, _p) = proj.point(jmin);
            let (_pt, d1u, d1v) = s_base.d1(uu, vv);
            let mut n2 = d1u.crossed(&d1v);
            if n2.magnitude() <= CONFUSION {
                continue;
            }
            n2 = n2.normalized();
            // `cxx:437-439`.
            if f_base.0.orientation() == Orientation::Reversed {
                n2.reverse();
            }
            // `cxx:440-445`.
            if n1.dot(&n2) > 0.0 {
                return Orientation::Reversed;
            }
            return Orientation::Forward;
        }
    }
    // `cxx:451`.
    Orientation::Internal
}
