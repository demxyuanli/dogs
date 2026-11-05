//! `Bisector_PolyBis` (`Bisector_PolyBis.cxx:18-131`, `.hxx:30-57`): the
//! polygon of `Bisector_PointOnBis` vertices that gives the initial parameter
//! estimate used by `Bisector_BisecCC`.
//!
//! `thePoints` is a fixed `Bisector_PointOnBis[30]` array (`.hxx:53`) with
//! 1-based `nbPoints`. `Append` is called at most `NbPnts = 21` times by
//! `Bisector_BisecCC::Perform` (`cxx:116`, `:204`), so the array is never
//! exhausted there; the Rust port panics on overflow where OCCT writes out of
//! bounds.

use occt_core::gp::trsf2d::GpTrsf2d;
use occt_core::precision::REAL_SMALL;

use super::point_on_bis::BisectorPointOnBis;

/// `Bisector_PointOnBis thePoints[30]` (`.hxx:53`).
const MAX_POINTS: usize = 30;

/// `class Bisector_PolyBis` (`Bisector_PolyBis.hxx:30-57`).
#[derive(Debug, Clone)]
pub struct BisectorPolyBis {
    /// `thePoints` (`.hxx:53`), 1-based like OCCT: index `k` is OCCT's
    /// `thePoints[k]`.
    the_points: [BisectorPointOnBis; MAX_POINTS],
    /// `nbPoints` (`.hxx:54`).
    nb_points: usize,
}

impl Default for BisectorPolyBis {
    fn default() -> Self {
        Self::new()
    }
}

impl BisectorPolyBis {
    /// `Bisector_PolyBis()` (`cxx:20-23`).
    pub fn new() -> Self {
        Self {
            the_points: [BisectorPointOnBis::new(); MAX_POINTS],
            nb_points: 0,
        }
    }

    /// `Append(P)` (`cxx:27-31`).
    pub fn append(&mut self, p: BisectorPointOnBis) {
        // OCCT writes `thePoints[++nbPoints]` without a bounds check (`.hxx:53`
        // is a 30-slot array); the Rust port reports the overflow instead of
        // writing past the end.
        assert!(
            self.nb_points < MAX_POINTS,
            "Bisector_PolyBis: more than {MAX_POINTS} points"
        );
        self.nb_points += 1;
        self.the_points[self.nb_points - 1] = p;
    }

    /// `Length()` (`cxx:35-38`).
    pub fn length(&self) -> usize {
        self.nb_points
    }

    /// `IsEmpty()` (`cxx:42-45`).
    pub fn is_empty(&self) -> bool {
        self.nb_points == 0
    }

    /// `Value(Index)` (`cxx:49-52`); 1-based, panics where OCCT reads
    /// `thePoints[Index]` out of range.
    pub fn value(&self, index: usize) -> BisectorPointOnBis {
        assert!(index >= 1 && index <= self.nb_points, "Bisector_PolyBis: OutOfRange");
        self.the_points[index - 1]
    }

    /// `First()` (`cxx:57-60`).
    pub fn first(&self) -> BisectorPointOnBis {
        assert!(self.nb_points >= 1, "Bisector_PolyBis: empty");
        self.the_points[0]
    }

    /// `Last()` (`cxx:64-67`).
    pub fn last(&self) -> BisectorPointOnBis {
        assert!(self.nb_points >= 1, "Bisector_PolyBis: empty");
        self.the_points[self.nb_points - 1]
    }

    /// `Interval(U)` (`cxx:81-121`): the 1-based index of the polygon segment
    /// that brackets `U`.
    pub fn interval(&self, u: f64) -> usize {
        if self.last().param_on_bis() - u < REAL_SMALL {
            return self.nb_points - 1;
        }
        let d_u =
            (self.last().param_on_bis() - self.first().param_on_bis()) / ((self.nb_points - 1) as f64);
        if d_u <= REAL_SMALL {
            return 1;
        }

        let mut int_u = ((u - self.first().param_on_bis()).abs() / d_u) as usize;
        int_u += 1;

        if self.the_points[int_u - 1].param_on_bis() >= u {
            let mut i = int_u as isize;
            while i >= 1 {
                if self.the_points[(i as usize) - 1].param_on_bis() <= u {
                    int_u = i as usize;
                    break;
                }
                i -= 1;
            }
        } else {
            let mut i = int_u;
            while i <= self.nb_points - 1 {
                if self.the_points[i - 1].param_on_bis() >= u {
                    int_u = i - 1;
                    break;
                }
                i += 1;
            }
        }
        int_u
    }

    /// `Transform(T)` (`cxx:125-131`).
    pub fn transform(&mut self, t: &GpTrsf2d) {
        for i in 0..self.nb_points {
            let mut p = self.the_points[i].point();
            p.transform(t);
            self.the_points[i].set_point(p);
        }
    }
}
