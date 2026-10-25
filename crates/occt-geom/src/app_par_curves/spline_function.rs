//! `AppParCurves::SplineFunction` — the B-spline basis matrix and its
//! derivative that `AppParCurves_LeastSquare` fills for the non-Bezier case.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves.hxx:20-48` and
//! `AppParCurves.cxx:159-256` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The OCCT body copies `flatknots` into a local `NCollection_Array1`
//!   `Aflatknots` with the same bounds and reads it 1-based (`cxx:177-182`);
//!   the port copies the `MathVector` into a `Vec` the same way and reads it
//!   through [`flat_knot`], which panics with `Standard_OutOfRange` where
//!   OCCT would raise the same exception.
//! - `BSplCLib::LocateParameter(deg, Aflatknots, U, false, deg1, nbpoles + 1,
//!   kindex, NewU)` is `occt_core::bspl::locate::locate_parameter_flat`
//!   (`BSplCLib.cxx:189-214`); `NewU` is only an out-parameter that the body
//!   never reads, so it is dropped.
//! - `locpoles` / `locdpoles` are `math_Vector(1, deg1)`; the port keeps them
//!   as `deg1`-sized vectors behind the same 1-based helpers.

use occt_core::bspl::locate::locate_parameter_flat;
use occt_math::{MathIntVector, MathMatrix, MathVector};

use crate::app_par_curves::bernstein::AppParCurves;

/// 1-based read of the copied `Aflatknots` array.
fn flat_knot(aflatknots: &[f64], i: i32) -> f64 {
    match i.checked_sub(1).and_then(|k| aflatknots.get(k as usize)) {
        Some(v) => *v,
        None => panic!("Standard_OutOfRange: AppParCurves::SplineFunction"),
    }
}

/// 1-based read of a local pole vector.
fn pole(v: &[f64], i: i32) -> f64 {
    match i.checked_sub(1).and_then(|k| v.get(k as usize)) {
        Some(x) => *x,
        None => panic!("Standard_OutOfRange: AppParCurves::SplineFunction"),
    }
}

/// 1-based write of a local pole vector.
fn set_pole(v: &mut [f64], i: i32, x: f64) {
    match i.checked_sub(1).and_then(|k| v.get_mut(k as usize)) {
        Some(slot) => *slot = x,
        None => panic!("Standard_OutOfRange: AppParCurves::SplineFunction"),
    }
}

impl AppParCurves {
    /// `AppParCurves::SplineFunction(nbpoles, deg, Parameters, flatknots, A,
    /// DA, index)` (`cxx:159-256`).
    ///
    /// `index` receives `kindex - deg - 1` for every parameter row
    /// (`cxx:198`).
    pub fn spline_function(
        nbpoles: i32,
        deg: i32,
        parameters: &MathVector,
        flatknots: &MathVector,
        a: &mut MathMatrix,
        da: &mut MathMatrix,
        index: &mut MathIntVector,
    ) {
        let deg1 = deg + 1;
        let mut locpoles = vec![0.0; deg1 as usize];
        let mut locdpoles = vec![0.0; deg1 as usize];
        let firstp = parameters.lower();
        let lastp = parameters.upper();

        // `cxx:177-182`: the local copy of the flat knots.
        let aflatknots: Vec<f64> = (flatknots.lower()..=flatknots.upper())
            .map(|i| flatknots.value(i))
            .collect();

        // OCCT seeds `oldkindex = 1` (`cxx:185`) and passes it as the locate
        // guess; `BSplCLib::LocateParameter`'s flat-knot overload ignores the
        // incoming guess, so the seed is written before it can be read.
        let mut oldkindex: i32;

        for i in firstp..=lastp {
            let u = parameters.value(i);
            let (kindex, _new_u) = locate_parameter_flat(
                deg,
                &aflatknots,
                u,
                false,
                deg1,
                nbpoles + 1,
            );
            oldkindex = kindex;

            index.set_value(i, kindex - deg - 1);

            set_pole(&mut locpoles, 1, 1.0);

            let mut qq = 2;
            while qq <= deg {
                set_pole(&mut locpoles, qq, 0.0);
                let mut pp = 1;
                while pp <= qq - 1 {
                    let inverse = 1.0
                        / (flat_knot(&aflatknots, kindex + pp)
                            - flat_knot(&aflatknots, kindex - qq + pp + 1));
                    let saved = (u - flat_knot(&aflatknots, kindex - qq + pp + 1))
                        * inverse
                        * pole(&locpoles, pp);
                    let cur = pole(&locpoles, pp);
                    set_pole(
                        &mut locpoles,
                        pp,
                        cur * (flat_knot(&aflatknots, kindex + pp) - u) * inverse,
                    );
                    let cur = pole(&locpoles, pp) + pole(&locpoles, qq);
                    set_pole(&mut locpoles, pp, cur);
                    set_pole(&mut locpoles, qq, saved);
                    pp += 1;
                }
                qq += 1;
            }

            qq = deg + 1;
            for pp in 1..=deg {
                let v = pole(&locpoles, pp);
                set_pole(&mut locdpoles, pp, v);
            }

            let mut locqq = 0.0;
            let mut locdqq = 0.0;
            for pp in 1..=deg {
                let inverse = 1.0
                    / (flat_knot(&aflatknots, kindex + pp)
                        - flat_knot(&aflatknots, kindex - qq + pp + 1));
                let saved = (u - flat_knot(&aflatknots, kindex - qq + pp + 1))
                    * inverse
                    * pole(&locpoles, pp);
                let cur = pole(&locpoles, pp);
                set_pole(
                    &mut locpoles,
                    pp,
                    cur * (flat_knot(&aflatknots, kindex + pp) - u) * inverse,
                );
                let cur = pole(&locpoles, pp) + locqq;
                set_pole(&mut locpoles, pp, cur);
                locqq = saved;
                let local_inverse = (deg as f64) * inverse;
                let saved = local_inverse * pole(&locdpoles, pp);
                let cur = pole(&locdpoles, pp) * -local_inverse;
                set_pole(&mut locdpoles, pp, cur);
                let cur = pole(&locdpoles, pp) + locdqq;
                set_pole(&mut locdpoles, pp, cur);
                locdqq = saved;
            }

            set_pole(&mut locpoles, qq, locqq);
            set_pole(&mut locdpoles, qq, locdqq);

            for j in 1..=deg1 {
                let val = pole(&locpoles, j);
                let theindex = j + oldkindex - deg1;
                a.set_value(i, theindex as usize, val);
                da.set_value(i, theindex as usize, pole(&locdpoles, j));
            }

            for j in 1..oldkindex - deg {
                a.set_value(i, j as usize, 0.0);
                da.set_value(i, j as usize, 0.0);
            }
            for j in oldkindex + 1..=nbpoles {
                a.set_value(i, j as usize, 0.0);
                da.set_value(i, j as usize, 0.0);
            }
        }
    }
}
