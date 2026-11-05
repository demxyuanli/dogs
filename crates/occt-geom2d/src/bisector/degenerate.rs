//! `Degenerate` (`Bisector_BisecAna.cxx:1736-1802`): replaces a degenerated
//! conic bisector by the straight line it collapses to.

use occt_core::gp::GpLin2d;

use crate::gcc_int::GccIntBisec;

/// `static bool Degenerate(GccInt_Bisec&, double)` (`cxx:1736-1802`).
///
/// A hyperbola with a vanishing major radius degenerates to its Y axis, one
/// with a vanishing minor radius to its X axis; a parabola with a vanishing
/// focal length to its symmetry axis, and an ellipse with a vanishing minor
/// radius to its major axis. The two hyperbola tests are independent `if`s in
/// OCCT, so a doubly degenerated hyperbola ends on the X axis.
pub fn degenerate(a_bisector: &mut GccIntBisec, tolerance: f64) -> bool {
    let replacement = match a_bisector {
        GccIntBisec::Hyperbola(h) => {
            // Both tests run, the second one wins (`cxx:1753-1771`).
            let mut line = None;
            if h.major_radius < tolerance {
                line = Some(GpLin2d::new(h.y_axis()));
            }
            if h.minor_radius < tolerance {
                line = Some(GpLin2d::new(h.x_axis()));
            }
            line
        }
        GccIntBisec::Parabola(p) => {
            if p.focal < tolerance {
                Some(GpLin2d::new(p.mirror_axis()))
            } else {
                None
            }
        }
        GccIntBisec::Ellipse(e) => {
            if e.minor_radius < tolerance {
                Some(GpLin2d::new(e.x_axis()))
            } else {
                None
            }
        }
        _ => None,
    };
    match replacement {
        Some(line) => {
            *a_bisector = GccIntBisec::Line(line);
            true
        }
        None => false,
    }
}
