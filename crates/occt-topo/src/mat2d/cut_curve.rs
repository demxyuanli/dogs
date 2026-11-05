//! `MAT2d_CutCurve` (`MAT2d_CutCurve.cxx:27-118`, `.hxx:...`): splits a 2D curve
//! at its curvature extrema and inflection points.
//!
//! `Perform` (`cxx:39-85`) runs `GeomLProp_CurAndInf2d` on the curve and cuts it
//! wherever the characteristic parameter is far enough from the running start
//! (`PTol = PConfusion * 10`) and the point far enough from the running start
//! point (`Tol = Confusion * 10`). A cut is kept only when at least one interior
//! parameter survived; the trailing piece `[UF, UL]` is appended then.
//!
//! `NbCurves` / `Value` raise `Standard_OutOfRange` while the curve was left
//! untouched, which the port reports as a `panic!`.

use std::sync::Arc;

use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::cur_and_inf::GeomLPropCurAndInf2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::Curve2d;

/// `class MAT2d_CutCurve` (`MAT2d_CutCurve.hxx`).
#[derive(Clone, Default)]
pub struct Mat2dCutCurve {
    /// `theCurves` (`.hxx`), a `NCollection_Sequence<Handle(Geom2d_TrimmedCurve)>`.
    the_curves: Vec<Geom2dTrimmedCurve>,
}

impl Mat2dCutCurve {
    /// `MAT2d_CutCurve()` (`cxx:27`).
    pub fn new() -> Self {
        Self { the_curves: Vec::new() }
    }

    /// `MAT2d_CutCurve(C)` (`cxx:33-36`).
    pub fn with_curve(c: &dyn Curve2d) -> Self {
        let mut cut = Self::new();
        cut.perform(c);
        cut
    }

    /// `Perform(C)` (`cxx:39-85`).
    pub fn perform(&mut self, c: &dyn Curve2d) {
        self.the_curves.clear();

        let mut sommets = GeomLPropCurAndInf2d::new();
        // `cxx:45-46`.
        let p_tol = PCONFUSION * 10.0;
        let tol = CONFUSION * 10.0;
        let mut ya_cut = false;

        sommets.perform(c);

        // The curve handle is needed by every `Geom2d_TrimmedCurve` this
        // method builds; `c` is borrowed for the whole call, so each trimmed
        // curve owns its own clone.
        let curve: Arc<dyn Curve2d> = Arc::from(c.clone_dyn());

        if sommets.is_done() && !sommets.points().is_empty() {
            let mut uf = c.first_parameter();
            let ul = c.last_parameter();
            let mut pf = c.d0(uf);
            let pl = c.d0(ul);

            for i in 1..=sommets.points().nb_points() as i32 {
                let uc = sommets.points().parameter(i);
                let pc = c.d0(uc);
                if uc - uf > p_tol && pc.distance(&pf) > tol {
                    if ul - uc < p_tol || pl.distance(&pc) < tol {
                        break;
                    }
                    self.the_curves
                        .push(Geom2dTrimmedCurve::new(curve.clone(), uf, uc));
                    uf = uc;
                    pf = pc;
                    ya_cut = true;
                }
            }
            if ya_cut {
                // `cxx:79-80`.
                self.the_curves
                    .push(Geom2dTrimmedCurve::new(curve.clone(), uf, ul));
            }
        }
    }

    /// `UnModified()` (`cxx:89-92`).
    pub fn un_modified(&self) -> bool {
        self.the_curves.is_empty()
    }

    /// `NbCurves()` (`cxx:94-101`); panics where OCCT raises
    /// `Standard_OutOfRange` on an unmodified curve.
    pub fn nb_curves(&self) -> usize {
        if self.un_modified() {
            panic!("MAT2d_CutCurve::NbCurves: Standard_OutOfRange");
        }
        self.the_curves.len()
    }

    /// `Value(Index)` (`cxx:103-118`); panics where OCCT raises
    /// `Standard_OutOfRange`.
    pub fn value(&self, index: i32) -> &Geom2dTrimmedCurve {
        if self.un_modified() {
            panic!("MAT2d_CutCurve::Value: Standard_OutOfRange");
        }
        if index < 1 || (index as usize) > self.the_curves.len() {
            panic!("MAT2d_CutCurve::Value: Standard_OutOfRange");
        }
        &self.the_curves[(index as usize) - 1]
    }

    /// The pieces produced by [`Mat2dCutCurve::perform`] (empty while
    /// unmodified).
    pub fn curves(&self) -> &[Geom2dTrimmedCurve] {
        &self.the_curves
    }
}
