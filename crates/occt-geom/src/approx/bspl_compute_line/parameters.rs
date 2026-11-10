//! `Approx_BSplComputeLine::Parameters` (`Approx_BSplComputeLine.gxx:811-894`):
//! the initial parametrisation of a line slice.

use crate::approx::line_tool::LineTool;
use crate::approx::bspl_compute_line::BsplComputeLine;
use occt_core::approx::ApproxParametrizationType;

/// `Parameters(Line, firstP, LastP, TheParameters)` (`gxx:811-894`).
///
/// OCCT fills the caller's `math_Vector` at the absolute point indices; the
/// port returns the same values indexed from `first_p` (slot `i - first_p`).
pub(super) fn compute_parameters<L: LineTool>(
    this: &BsplComputeLine,
    line: &L,
    first_p: i32,
    last_p: i32,
) -> Vec<f64> {
    let nbp = (last_p - first_p + 1).max(1) as usize;
    let mut the_parameters = vec![0.0f64; nbp];

    // The first parameter is always zero; scaling it would divide by zero when
    // the slice holds a single point (`gxx:824-826`).
    the_parameters[0] = 0.0;

    let a_nbp = last_p - first_p + 1;
    if a_nbp == 2 {
        the_parameters[nbp - 1] = 1.0;
    } else if this.par == ApproxParametrizationType::ChordLength
        || this.par == ApproxParametrizationType::Centripetal
    {
        let nb_p3d = line.nb_p3d();
        let nb_p2d = line.nb_p2d();

        let mut dist;
        for i in first_p + 1..=last_p {
            let (prev, prev2d) = line.value(i - 1);
            let (cur, cur2d) = line.value(i);

            dist = 0.0;
            for j in 0..nb_p3d {
                let a1 = prev[j as usize];
                let a2 = cur[j as usize];
                dist += square_distance(&a1, &a2);
            }
            for j in 0..nb_p2d {
                let a12d = prev2d[j as usize];
                let a22d = cur2d[j as usize];
                dist += square_distance2d(&a12d, &a22d);
            }

            dist = dist.sqrt();
            let slot = (i - first_p) as usize;
            if this.par == ApproxParametrizationType::ChordLength {
                the_parameters[slot] = the_parameters[slot - 1] + dist;
            } else {
                // Par == Approx_Centripetal (`gxx:879-882`).
                the_parameters[slot] = the_parameters[slot - 1] + dist.sqrt();
            }
        }
        let last = the_parameters[nbp - 1];
        for i in first_p + 1..=last_p {
            the_parameters[(i - first_p) as usize] /= last;
        }
    } else {
        // `gxx:886-893`.
        for i in first_p + 1..=last_p {
            the_parameters[(i - first_p) as usize] =
                (i as f64 - first_p as f64) / (last_p as f64 - first_p as f64);
        }
    }

    the_parameters
}

/// `gp_Pnt::SquareDistance` (`gp_Pnt.hxx:145-155`).
fn square_distance(a: &occt_core::gp::GpPnt, b: &occt_core::gp::GpPnt) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    let dz = a.z() - b.z();
    dx * dx + dy * dy + dz * dz
}

/// `gp_Pnt2d::SquareDistance` (`gp_Pnt2d.hxx:117-127`).
fn square_distance2d(
    a: &occt_core::gp::GpPnt2d,
    b: &occt_core::gp::GpPnt2d,
) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    dx * dx + dy * dy
}



