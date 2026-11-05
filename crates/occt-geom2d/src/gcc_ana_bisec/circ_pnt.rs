//! `GccAna_CircPnt2dBisec` (`GccAna_CircPnt2dBisec.cxx:36-196`): the bisector
//! of a circle and a point - a circle when the point is the centre, an ellipse
//! inside, a straight line on, and one or two hyperbola branches outside.

use occt_core::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpPnt2d};

use super::{check_index, GccAnaBisecError};
use crate::gcc_int::GccIntBisec;

/// `class GccAna_CircPnt2dBisec`.
#[derive(Debug, Clone, Copy)]
pub struct GccAnaCircPnt2dBisec {
    circle: GpCirc2d,
    point: GpPnt2d,
    my_tolerance: f64,
    /// `-1` inside, `0` on, `1` outside the circle.
    theposition: i32,
    nbr_sol: i32,
    well_done: bool,
}

impl GccAnaCircPnt2dBisec {
    /// `GccAna_CircPnt2dBisec::GccAna_CircPnt2dBisec(const gp_Circ2d&, const gp_Pnt2d&)`
    /// (`cxx:36-42`): the tolerance is fixed to `1.0e-10`.
    pub fn new(circle: &GpCirc2d, point: &GpPnt2d) -> Self {
        Self::with_tolerance(circle, point, 0.0)
    }

    /// `GccAna_CircPnt2dBisec::GccAna_CircPnt2dBisec(const gp_Circ2d&, const gp_Pnt2d&, double)`
    /// (`cxx:44-56`): `myTolerance` starts at `1.0e-10` and is only *raised* by
    /// the argument, never lowered.
    pub fn with_tolerance(circle: &GpCirc2d, point: &GpPnt2d, tolerance: f64) -> Self {
        let mut my_tolerance: f64 = 1.0e-10;
        if my_tolerance < tolerance {
            my_tolerance = tolerance;
        }
        let mut me = Self {
            circle: *circle,
            point: *point,
            my_tolerance,
            theposition: 0,
            nbr_sol: 0,
            well_done: false,
        };
        me.define_solutions();
        me
    }

    /// `GccAna_CircPnt2dBisec::DefineSolutions` (`cxx:58-80`).
    fn define_solutions(&mut self) {
        let dist = self.circle.radius() - self.point.distance(&self.circle.location());
        if dist.abs() < self.my_tolerance {
            self.theposition = 0;
            self.nbr_sol = 1;
        } else if dist > 0.0 {
            self.theposition = -1;
            self.nbr_sol = 1;
        } else {
            self.theposition = 1;
            self.nbr_sol = 2;
        }
        self.well_done = true;
    }

    /// `GccAna_CircPnt2dBisec::IsDone` (`cxx:187-190`).
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// `GccAna_CircPnt2dBisec::NbSolutions` (`cxx:192-196`).
    pub fn nb_solutions(&self) -> i32 {
        self.nbr_sol
    }

    /// `GccAna_CircPnt2dBisec::ThisSolution` (`cxx:94-184`).
    pub fn this_solution(&self, index: i32) -> Result<GccIntBisec, GccAnaBisecError> {
        check_index(self.well_done, self.nbr_sol, index)?;

        let xpoint = self.point.x();
        let ypoint = self.point.y();
        let xcencir = self.circle.location().x();
        let ycencir = self.circle.location().y();
        let r1 = self.circle.radius();
        let dist = self.point.distance(&self.circle.location());

        if dist < self.my_tolerance {
            // Concentred point: the bisector is the circle of half the radius.
            let circ = GpCirc2d::new(
                GpAx22d::from_xdir(self.point, GpDir2d::default()),
                r1 / 2.0,
            );
            return Ok(GccIntBisec::Circle(circ));
        }

        let center = GpPnt2d::new((xpoint + xcencir) / 2.0, (ypoint + ycencir) / 2.0);
        let majax = GpAx2d::new(
            center,
            GpDir2d::new(xpoint - xcencir, ypoint - ycencir).unwrap_or(GpDir2d::default()),
        );

        if self.theposition == -1 {
            // Point inside the circle: an ellipse of half the circle radius and
            // half the half-chord.
            let elips = GpElips2d::new(
                GpAx22d::from_xdir(majax.loc, majax.vdir),
                r1 / 2.0,
                (r1 * r1 - dist * dist).sqrt() / 2.0,
            );
            Ok(GccIntBisec::Ellipse(elips))
        } else if self.theposition == 0 {
            // Point on the circle: a straight line through point and centre.
            let dirsol = if self.circle.is_direct() {
                GpDir2d::new(xcencir - xpoint, ycencir - ypoint)
            } else {
                GpDir2d::new(xpoint - xcencir, ypoint - ycencir)
            };
            let line = GpLin2d::from_pnt_dir(self.point, dirsol.unwrap_or(GpDir2d::default()));
            Ok(GccIntBisec::Line(line))
        } else {
            // Point outside: the two main branches of a hyperbola.
            let d1 = (dist * dist - r1 * r1).sqrt() / 2.0;
            let d2 = r1 / 2.0;
            let hyp = GpHypr2d::new(GpAx22d::from_xdir(majax.loc, majax.vdir), d2, d1);
            if index == 1 {
                Ok(GccIntBisec::Hyperbola(hyp))
            } else {
                Ok(GccIntBisec::Hyperbola(hyp.other_branch()))
            }
        }
    }
}
