//! `math_FunctionAllRoots`.
//! Source: `math_FunctionAllRoots.hxx:34-99`, `.lxx:18-64`,
//! `math_FunctionAllRoots.cxx:26-238`.
//!
//! This algorithm uses a sample of the function to find all intervals on which
//! the function is null, and afterwards uses the `math_FunctionRoots` algorithm
//! to find the points where the function is null outside the null intervals.
//! Knowledge of the derivative is required.

use crate::math_fn::MathFunctionWithDerivative;
use crate::math_function_roots::FunctionRoots;
use crate::math_function_sample::FunctionSample;

/// `math_FunctionAllRoots` (`math_FunctionAllRoots.hxx:34-99`).
#[derive(Debug, Clone)]
pub struct FunctionAllRoots {
    done: bool,
    pdeb: Vec<f64>,
    pfin: Vec<f64>,
    piso: Vec<f64>,
    ideb: Vec<i32>,
    ifin: Vec<i32>,
    iiso: Vec<i32>,
}

impl FunctionAllRoots {
    /// `math_FunctionAllRoots(F, S, EpsX, EpsF, EpsNul)`
    /// (`math_FunctionAllRoots.cxx:26-238`).
    pub fn new<F: MathFunctionWithDerivative>(
        f: &mut F,
        s: &FunctionSample,
        eps_x: f64,
        eps_f: f64,
        eps_nul: f64,
    ) -> Self {
        let mut pdeb: Vec<f64> = Vec::new();
        let mut pfin: Vec<f64> = Vec::new();
        let mut piso: Vec<f64> = Vec::new();
        let mut ideb: Vec<i32> = Vec::new();
        let mut ifin: Vec<i32> = Vec::new();
        let mut iiso: Vec<i32> = Vec::new();

        let nbp = s.nb_points();
        let mut val = 0.0;
        f.value(s.get_parameter(1), &mut val);
        let mut p_nul = val.abs() <= eps_nul;
        let mut val_sav = 0.0;
        if !p_nul {
            val_sav = val;
        }
        let mut inter_nul = false;
        let mut nul_d = false;
        let mut nul_f = false;
        let mut deb_nul = 0.0;
        let mut fin_nul;
        let mut ind_d = 0;
        let mut ind_f;

        let mut i = 2;
        let mut fini = i > nbp;

        while !fini {
            f.value(s.get_parameter(i), &mut val);
            let nul = val.abs() <= eps_nul;
            if !nul {
                val_sav = val;
            }
            if inter_nul && !nul {
                inter_nul = false;
                pdeb.push(deb_nul);
                ideb.push(ind_d);
                let mut cst = if val > 0.0 { eps_nul } else { -eps_nul };
                let res1 = FunctionRoots::new(
                    f,
                    s.get_parameter(i - 1),
                    s.get_parameter(i),
                    10,
                    eps_x,
                    eps_f,
                    0.0,
                    cst,
                );
                assert!(
                    res1.is_done() && !res1.is_all_null() && res1.nb_solutions() != 0,
                    "Standard_NumericError in math_FunctionAllRoots"
                );

                fin_nul = res1.value(1);
                ind_f = res1.state_number(1);

                cst = -cst;
                let res2 = FunctionRoots::new(
                    f,
                    s.get_parameter(i - 1),
                    s.get_parameter(i),
                    10,
                    eps_x,
                    eps_f,
                    0.0,
                    cst,
                );
                assert!(
                    res2.is_done() && !res2.is_all_null(),
                    "Standard_NumericError in math_FunctionAllRoots"
                );

                if res2.nb_solutions() != 0 {
                    if res2.value(1) < fin_nul {
                        fin_nul = res2.value(1);
                        ind_f = res2.state_number(1);
                    }
                }
                pfin.push(fin_nul);
                ifin.push(ind_f);
            } else if !inter_nul && p_nul && nul {
                inter_nul = true;
                if i == 2 {
                    deb_nul = s.get_parameter(1);
                    let mut val_bid = 0.0;
                    f.value(deb_nul, &mut val_bid);
                    ind_d = f.get_state_number();
                    nul_d = true;
                } else {
                    let mut cst = if val_sav > 0.0 { eps_nul } else { -eps_nul };
                    let res1 = FunctionRoots::new(
                        f,
                        s.get_parameter(i - 2),
                        s.get_parameter(i - 1),
                        10,
                        eps_x,
                        eps_f,
                        0.0,
                        cst,
                    );
                    assert!(
                        res1.is_done() && !res1.is_all_null() && res1.nb_solutions() != 0,
                        "Standard_NumericError in math_FunctionAllRoots"
                    );
                    deb_nul = res1.value(res1.nb_solutions());
                    ind_d = res1.state_number(res1.nb_solutions());

                    cst = -cst;
                    let res3 = FunctionRoots::new(
                        f,
                        s.get_parameter(i - 2),
                        s.get_parameter(i - 1),
                        10,
                        eps_x,
                        eps_f,
                        0.0,
                        cst,
                    );
                    assert!(
                        res3.is_done() && !res3.is_all_null(),
                        "Standard_NumericError in math_FunctionAllRoots"
                    );

                    if res3.nb_solutions() != 0 {
                        if res3.value(res3.nb_solutions()) > deb_nul {
                            deb_nul = res3.value(res3.nb_solutions());
                            ind_d = res3.state_number(res3.nb_solutions());
                        }
                    }
                }
            }
            i += 1;
            p_nul = nul;
            fini = i > nbp;
        }

        if inter_nul {
            pdeb.push(deb_nul);
            ideb.push(ind_d);
            fin_nul = s.get_parameter(nbp);
            let mut val_bid = 0.0;
            f.value(fin_nul, &mut val_bid);
            ind_f = f.get_state_number();
            pfin.push(fin_nul);
            ifin.push(ind_f);
            nul_f = true;
        }

        if pdeb.is_empty() {
            let res = FunctionRoots::new(
                f,
                s.get_parameter(1),
                s.get_parameter(nbp),
                nbp,
                eps_x,
                eps_f,
                0.0,
                0.0,
            );
            assert!(
                res.is_done() && !res.is_all_null(),
                "Standard_NumericError in math_FunctionAllRoots"
            );

            for j in 1..=res.nb_solutions() {
                piso.push(res.value(j));
                iiso.push(res.state_number(j));
            }
        } else {
            let nbp_min = 3;
            if !nul_d {
                let nbrpt = (((pdeb[0] - s.get_parameter(1))
                    / (s.get_parameter(nbp) - s.get_parameter(1)))
                .abs()
                    * nbp as f64)
                    .trunc() as i32;
                let res = FunctionRoots::new(
                    f,
                    s.get_parameter(1),
                    pdeb[0],
                    if nbrpt > nbp_min { nbrpt } else { nbp_min },
                    eps_x,
                    eps_f,
                    0.0,
                    0.0,
                );
                assert!(
                    res.is_done() && !res.is_all_null(),
                    "Standard_NumericError in math_FunctionAllRoots"
                );

                for j in 1..=res.nb_solutions() {
                    piso.push(res.value(j));
                    iiso.push(res.state_number(j));
                }
            }
            for k in 2..=pdeb.len() {
                let nbrpt = (((pdeb[k - 1] - pfin[k - 2])
                    / (s.get_parameter(nbp) - s.get_parameter(1)))
                .abs()
                    * nbp as f64)
                    .trunc() as i32;
                let res = FunctionRoots::new(
                    f,
                    pfin[k - 2],
                    pdeb[k - 1],
                    if nbrpt > nbp_min { nbrpt } else { nbp_min },
                    eps_x,
                    eps_f,
                    0.0,
                    0.0,
                );
                assert!(
                    res.is_done() && !res.is_all_null(),
                    "Standard_NumericError in math_FunctionAllRoots"
                );

                for j in 1..=res.nb_solutions() {
                    piso.push(res.value(j));
                    iiso.push(res.state_number(j));
                }
            }
            if !nul_f {
                let last = pdeb.len();
                let nbrpt = (((s.get_parameter(nbp) - pfin[last - 1])
                    / (s.get_parameter(nbp) - s.get_parameter(1)))
                .abs()
                    * nbp as f64)
                    .trunc() as i32;
                let res = FunctionRoots::new(
                    f,
                    pfin[last - 1],
                    s.get_parameter(nbp),
                    if nbrpt > nbp_min { nbrpt } else { nbp_min },
                    eps_x,
                    eps_f,
                    0.0,
                    0.0,
                );
                assert!(
                    res.is_done() && !res.is_all_null(),
                    "Standard_NumericError in math_FunctionAllRoots"
                );

                for j in 1..=res.nb_solutions() {
                    piso.push(res.value(j));
                    iiso.push(res.state_number(j));
                }
            }
        }

        Self {
            done: true,
            pdeb,
            pfin,
            piso,
            ideb,
            ifin,
            iiso,
        }
    }

    /// `IsDone()` (`math_FunctionAllRoots.lxx:18-21`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbIntervals()` (`math_FunctionAllRoots.lxx:30-34`).
    pub fn nb_intervals(&self) -> i32 {
        assert!(self.done, "StdFail_NotDone in math_FunctionAllRoots::NbIntervals");
        self.pdeb.len() as i32
    }

    /// `GetInterval(Index, A, B)` (`math_FunctionAllRoots.lxx:36-41`), 1-based.
    pub fn get_interval(&self, index: i32) -> (f64, f64) {
        assert!(self.done, "StdFail_NotDone in math_FunctionAllRoots::GetInterval");
        (self.pdeb[(index - 1) as usize], self.pfin[(index - 1) as usize])
    }

    /// `GetIntervalState(Index, IFirst, ILast)` (`math_FunctionAllRoots.lxx:43-49`),
    /// 1-based.
    pub fn get_interval_state(&self, index: i32) -> (i32, i32) {
        assert!(
            self.done,
            "StdFail_NotDone in math_FunctionAllRoots::GetIntervalState"
        );
        (self.ideb[(index - 1) as usize], self.ifin[(index - 1) as usize])
    }

    /// `NbPoints()` (`math_FunctionAllRoots.lxx:51-55`).
    pub fn nb_points(&self) -> i32 {
        assert!(self.done, "StdFail_NotDone in math_FunctionAllRoots::NbPoints");
        self.piso.len() as i32
    }

    /// `GetPoint(Index)` (`math_FunctionAllRoots.lxx:57-61`), 1-based.
    pub fn get_point(&self, index: i32) -> f64 {
        assert!(self.done, "StdFail_NotDone in math_FunctionAllRoots::GetPoint");
        self.piso[(index - 1) as usize]
    }

    /// `GetPointState(Index)` (`math_FunctionAllRoots.lxx:63-68`), 1-based.
    pub fn get_point_state(&self, index: i32) -> i32 {
        assert!(self.done, "StdFail_NotDone in math_FunctionAllRoots::GetPointState");
        self.iiso[(index - 1) as usize]
    }
}
