//! `math_Uzawa` — resolution of `C * X = B` around a starting point `X0` by
//! the Uzawa algorithm, minimising `Norm(X - X0)`.
//!
//! Source: `FoundationClasses/TKMath/math/math_Uzawa.cxx:39-285` and
//! `math_Uzawa.hxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - Two constructors (`cxx:39-54` and `cxx:57-73`); both funnel into
//!   `Perform` (`cxx:77-263`). The second one passes `Nce` then `Nci`
//!   (`cxx:73`), matching the body's parameter order (`cxx:81-82`); the
//!   header spells the names the other way round (`hxx:120`) but the types
//!   are identical, so the behaviour is the first-`Nce`-rows-are-equal.
//! - When `Nci == 0` (`cxx:108-157`) the resolution is direct: `C * T(C)` is
//!   assembled in the **lower** triangle only (`cxx:113-122`, "a utiliser avec
//!   Crout"), inverted by [`MathCrout`], and `Vardua` is `Inv(C*T(C)) *
//!   Errinit` read symmetrically (`cxx:129-141`). `NbIter` is set to 1 and
//!   `Done` to `true`.
//! - Otherwise (`cxx:158-263`) the dual variables are iterated with
//!   `Rho = 1 / (sqrt(2) * Norm(C))` (`cxx:173-186`); the loop leaves
//!   `NbIter == NbIterations + 1` on normal termination (`cxx:190`),
//!   exactly as the C++ `for` does.
//! - `Resul` stays at its zero initialisation when the dual problem diverges
//!   (`cxx:246-251`) or when the iteration budget runs out (`cxx:262`).
//! - The `Xmax`/`ErrMax` initialisation quirk of `cxx:200-235` (assigning at
//!   `i == 1` and then `max`-ing again) is kept verbatim.
//! - `math_Vector(1, 0)` / `math_Matrix(1, 0, ...)` are legal in OCCT but not
//!   in this crate's containers, which need `upper >= lower`; a zero-sized
//!   direction keeps a 1x1 stand-in that the empty loops never touch.
//! - `IsDone`/`Value`/`InitialError`/`Error`/`NbIterations`/`InverseCont`
//!   raise `StdFail_NotDone` when the resolution failed
//!   (`math_Uzawa.lxx`); the port panics.

use crate::math_crout::MathCrout;
use crate::{MathMatrix, MathVector};

/// `math_Uzawa` (`math_Uzawa.hxx`).
#[derive(Debug, Clone)]
pub struct Uzawa {
    /// `Resul` (`hxx`): the solution vector, `Cont.ColNumber()` long.
    resul: MathVector,
    /// `Erruza` (`hxx`): `X - X0`.
    erruza: MathVector,
    /// `Errinit` (`hxx`): `Cont * X0 - Secont`.
    errinit: MathVector,
    /// `Vardua` (`hxx`): the dual variables, `Cont.RowNumber()` long.
    vardua: MathVector,
    /// `CTCinv` (`hxx`): `Inv(Cont * T(Cont))`, lower triangle only.
    ctc_inv: MathMatrix,
    /// `NbIter` (`hxx`).
    nb_iter: usize,
    /// `Done` (`hxx`).
    done: bool,
}

impl Uzawa {
    /// `math_Uzawa(Cont, Secont, StartingPoint, EpsLix = 1e-6,
    /// EpsLic = 1e-6, NbIterations = 500)` (`cxx:39-54`): all the equations of
    /// `Cont` are equalities.
    pub fn new(cont: &MathMatrix, secont: &MathVector, starting_point: &MathVector) -> Self {
        Self::new_with_tolerances(cont, secont, starting_point, 1.0e-6, 1.0e-6, 500)
    }

    /// Same as [`Uzawa::new`] with the caller's tolerances and iteration
    /// budget (`hxx:51-53`).
    pub fn new_with_tolerances(
        cont: &MathMatrix,
        secont: &MathVector,
        starting_point: &MathVector,
        eps_lix: f64,
        eps_lic: f64,
        nb_iterations: usize,
    ) -> Self {
        let nlig = cont.row_count();
        Self::build(cont, secont, starting_point, nlig, 0, eps_lix, eps_lic, nb_iterations)
    }

    /// `math_Uzawa(Cont, Secont, StartingPoint, Nce, Nci, ...)`
    /// (`cxx:57-73`): the `Nce` first equations are equalities, the `Nci`
    /// last ones are inequalities (`<`).
    pub fn new_with_constraints(
        cont: &MathMatrix,
        secont: &MathVector,
        starting_point: &MathVector,
        nce: usize,
        nci: usize,
        eps_lix: f64,
        eps_lic: f64,
        nb_iterations: usize,
    ) -> Self {
        Self::build(cont, secont, starting_point, nce, nci, eps_lix, eps_lic, nb_iterations)
    }

    /// Member initialisation shared by both constructors (`cxx:45-51`,
    /// `cxx:63-69`).
    fn build(
        cont: &MathMatrix,
        secont: &MathVector,
        starting_point: &MathVector,
        nce: usize,
        nci: usize,
        eps_lix: f64,
        eps_lic: f64,
        nb_iterations: usize,
    ) -> Self {
        let ncol = cont.col_count();
        let nlig = cont.row_count();
        let mut this = Self {
            resul: MathVector::with_init(1, ncol.max(1), 0.0),
            erruza: MathVector::with_init(1, ncol.max(1), 0.0),
            errinit: MathVector::with_init(1, ncol.max(1), 0.0),
            vardua: MathVector::with_init(1, nlig.max(1), 0.0),
            ctc_inv: MathMatrix::with_init(1, nlig.max(1), 1, nlig.max(1), 0.0),
            nb_iter: 0,
            done: false,
        };
        this.perform(cont, secont, starting_point, nce, nci, eps_lix, eps_lic, nb_iterations);
        this
    }

    /// `Perform(Cont, Secont, StartingPoint, Nce, Nci, EpsLix, EpsLic,
    /// NbIterations)` (`cxx:77-263`).
    fn perform(
        &mut self,
        cont: &MathMatrix,
        secont: &MathVector,
        starting_point: &MathVector,
        nce: usize,
        nci: usize,
        eps_lix: f64,
        eps_lic: f64,
        nb_iterations: usize,
    ) {
        let nlig = cont.row_count();
        let ncol = cont.col_count();

        // `Standard_DimensionError_Raise_if((Secont.Length() != Nlig) ||
        // ((Nce + Nci) != Nlig), " ")` (`cxx:94`).
        if secont.len() != nlig || (nce + nci) != nlig {
            panic!("Standard_DimensionError: math_Uzawa::Perform");
        }

        // Calcul du vecteur Cont*X0 - D: (erreur initiale) (`cxx:98-106`).
        for i in 1..=nlig {
            let mut acc = cont.value(i, 1) * starting_point.value(1) - secont.value(i);
            for j in 2..=ncol {
                acc += cont.value(i, j) * starting_point.value(j);
            }
            self.errinit.set_value(i, acc);
        }

        // Cas de resolution directe (`cxx:108-157`).
        if nci == 0 {
            self.nb_iter = 1;
            // Calcul de Cont*T(Cont), lower triangle only (`cxx:111-123`).
            for i in 1..=nlig {
                for j in 1..=i {
                    let mut acc = cont.value(i, 1) * cont.value(j, 1);
                    for k in 2..=ncol {
                        acc += cont.value(i, k) * cont.value(j, k);
                    }
                    self.ctc_inv.set_value(i, j, acc);
                }
            }
            // Calcul de l inverse de CTCinv (`cxx:127-128`).
            let inv = MathCrout::new(&self.ctc_inv);
            self.ctc_inv = inv.invert();
            // Vardua = Inv(C*T(C)) * Errinit, read symmetrically
            // (`cxx:129-141`).
            for i in 1..=nlig {
                let mut scale = self.ctc_inv.value(i, 1) * self.errinit.value(1);
                for j in 2..=i {
                    scale += self.ctc_inv.value(i, j) * self.errinit.value(j);
                }
                for j in (i + 1)..=nlig {
                    scale += self.ctc_inv.value(j, i) * self.errinit.value(j);
                }
                self.vardua.set_value(i, scale);
            }
            // Erruza = -T(Cont) * Vardua (`cxx:142-149`).
            for i in 1..=ncol {
                let mut acc = -cont.value(1, i) * self.vardua.value(1);
                for j in 2..=nlig {
                    acc -= cont.value(j, i) * self.vardua.value(j);
                }
                self.erruza.set_value(i, acc);
            }
            // Restitution des valeurs calculees (`cxx:152-155`).
            self.restore(starting_point, ncol);
            self.done = true;
            return;
        }

        // Initialisation des variables duales (`cxx:160-173`).
        for i in 1..=nlig {
            if i <= nce {
                self.vardua.set_value(i, 0.0);
            } else {
                self.vardua.set_value(i, 1.0);
            }
        }

        // Calcul du coefficient Rho (`cxx:175-186`).
        let coef = 1.0 / 2.0_f64.sqrt();
        let mut normat = 0.0;
        for i in 1..=nlig {
            let mut normli = cont.value(i, 1) * cont.value(i, 1);
            for j in 2..=ncol {
                normli += cont.value(i, j) * cont.value(i, j);
            }
            normat += normli;
        }
        let rho = coef / normat;

        // Boucle des iterations de la methode d Uzawa (`cxx:190-261`).
        let mut xmax: f64 = 0.0;
        let mut errmax: f64 = 0.0;
        self.nb_iter = 1;
        while self.nb_iter <= nb_iterations {
            // Erruza = -T(Cont) * Vardua and the step of X (`cxx:192-207`).
            for i in 1..=ncol {
                let xian = self.erruza.value(i);
                let mut acc = -cont.value(1, i) * self.vardua.value(1);
                for j in 2..=nlig {
                    acc -= cont.value(j, i) * self.vardua.value(j);
                }
                self.erruza.set_value(i, acc);
                if self.nb_iter > 1 {
                    if i == 1 {
                        xmax = (self.erruza.value(i) - xian).abs();
                    }
                    xmax = xmax.max((self.erruza.value(i) - xian).abs());
                }
            }

            // Calcul de Xmu a l iteration NbIter et evaluation de l erreur
            // sur la verification des contraintes (`cxx:212-236`).
            for i in 1..=nlig {
                let mut err = cont.value(i, 1) * self.erruza.value(1) + self.errinit.value(i);
                for j in 2..=ncol {
                    err += cont.value(i, j) * self.erruza.value(j);
                }
                let err1;
                if i <= nce {
                    self.vardua.set_value(i, self.vardua.value(i) + rho * err);
                    err1 = (rho * err).abs();
                } else {
                    let xmuian = self.vardua.value(i);
                    let nv = (self.vardua.value(i) + rho * err).max(0.0);
                    self.vardua.set_value(i, nv);
                    err1 = (nv - xmuian).abs();
                }
                if i == 1 {
                    errmax = err1;
                }
                errmax = errmax.max(err1);
            }

            if self.nb_iter > 1 && xmax <= eps_lix {
                if errmax <= eps_lic {
                    self.done = true;
                } else {
                    self.done = false;
                    return;
                }
                // Restitution des valeurs calculees (`cxx:255-257`).
                self.restore(starting_point, ncol);
                self.done = true;
                return;
            }
            self.nb_iter += 1;
        } // fin de la boucle d iterations.
        self.done = false;
    }

    /// `Resul = StartingPoint + Erruza` (`cxx:152`, `cxx:255`).
    fn restore(&mut self, starting_point: &MathVector, ncol: usize) {
        for i in 1..=ncol {
            self.resul.set_value(i, starting_point.value(i) + self.erruza.value(i));
        }
    }

    /// `IsDone()` (`math_Uzawa.lxx`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Value()` (`math_Uzawa.lxx`): the solution of the system.
    pub fn value(&self) -> &MathVector {
        if !self.done {
            panic!("StdFail_NotDone: math_Uzawa::Value");
        }
        &self.resul
    }

    /// `InitialError()` (`math_Uzawa.lxx`): `Cont * StartingPoint - Secont`.
    pub fn initial_error(&self) -> &MathVector {
        if !self.done {
            panic!("StdFail_NotDone: math_Uzawa::InitialError");
        }
        &self.errinit
    }

    /// `Duale(V)` (`cxx:266-268`): the dual variables.
    pub fn duale(&self) -> &MathVector {
        &self.vardua
    }

    /// `Error()` (`math_Uzawa.lxx`): `X - StartingPoint`.
    pub fn error(&self) -> &MathVector {
        if !self.done {
            panic!("StdFail_NotDone: math_Uzawa::Error");
        }
        &self.erruza
    }

    /// `NbIterations()` (`math_Uzawa.lxx`): the number of iterations really
    /// done.
    pub fn nb_iterations(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: math_Uzawa::NbIterations");
        }
        self.nb_iter
    }

    /// `InverseCont()` (`math_Uzawa.lxx`): `Inv(Cont * T(Cont))`, lower
    /// triangle only. Needed for the gradient when approximating a curve.
    pub fn inverse_cont(&self) -> &MathMatrix {
        if !self.done {
            panic!("StdFail_NotDone: math_Uzawa::InverseCont");
        }
        &self.ctc_inv
    }
}
