//! `Approx_ComputeLine::Parameters` (`gxx:1249-1322`): chord-length,
//! centripetal or iso parametrisation of the points `firstP..lastP`.

use occt_core::approx::ApproxParametrizationType;

use crate::approx::line_tool::LineTool;
use crate::approx::compute_line::ComputeLine;
use occt_math::MathVector;

impl ComputeLine {
    /// `Parameters(Line, firstP, lastP, TheParameters)` (`gxx:1249-1322`).
    pub(crate) fn compute_parameters<L: LineTool>(
        &self,
        line: &L,
        first_p: i32,
        last_p: i32,
    ) -> MathVector {
        let mut the_parameters =
            MathVector::new(first_p.max(0) as usize, last_p.max(0) as usize);

        if self.par == ApproxParametrizationType::ChordLength
            || self.par == ApproxParametrizationType::Centripetal
        {
            let nb_p3d = line.nb_p3d();
            let nb_p2d = line.nb_p2d();

            the_parameters.set_value(first_p as usize, 0.0);
            let mut dist: f64;

            for i in first_p + 1..=last_p {
                let (tab_p, tab_p2d) = if nb_p3d != 0 && nb_p2d != 0 {
                    line.value(i - 1)
                } else if nb_p2d != 0 {
                    (Vec::new(), line.value_2d(i - 1))
                } else {
                    (line.value_3d(i - 1), Vec::new())
                };
                let (tab_pp, tab_pp2d) = if nb_p3d != 0 && nb_p2d != 0 {
                    line.value(i)
                } else if nb_p2d != 0 {
                    (Vec::new(), line.value_2d(i))
                } else {
                    (line.value_3d(i), Vec::new())
                };

                dist = 0.0;
                for j in 0..nb_p3d as usize {
                    dist += square_distance(&tab_p[j], &tab_pp[j]);
                }
                for j in 0..nb_p2d as usize {
                    dist += square_distance2d(&tab_p2d[j], &tab_pp2d[j]);
                }

                dist = dist.sqrt();
                if self.par == ApproxParametrizationType::ChordLength {
                    let prev = the_parameters.value((i - 1) as usize);
                    the_parameters.set_value(i as usize, prev + dist);
                } else {
                    // Par == Approx_Centripetal
                    let prev = the_parameters.value((i - 1) as usize);
                    the_parameters.set_value(i as usize, prev + dist.sqrt());
                }
            }
            let last = the_parameters.value(last_p as usize);
            for i in first_p..=last_p {
                let v = the_parameters.value(i as usize);
                the_parameters.set_value(i as usize, v / last);
            }
        } else {
            for i in first_p..=last_p {
                let v = (i - first_p) as f64 / (last_p - first_p) as f64;
                the_parameters.set_value(i as usize, v);
            }
        }
        the_parameters
    }
}

/// `gp_Pnt::SquareDistance` (`gp_Pnt.hxx`).
fn square_distance(a: &occt_core::gp::GpPnt, b: &occt_core::gp::GpPnt) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    let dz = a.z() - b.z();
    dx * dx + dy * dy + dz * dz
}

/// `gp_Pnt2d::SquareDistance` (`gp_Pnt2d.hxx`).
fn square_distance2d(a: &occt_core::gp::GpPnt2d, b: &occt_core::gp::GpPnt2d) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    dx * dx + dy * dy
}





