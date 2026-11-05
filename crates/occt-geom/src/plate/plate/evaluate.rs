//! `Plate_Plate::Evaluate` / `Plate_Plate::EvaluateDerivative`.
//!
//! Source: `Plate_Plate.cxx:1257-1331`.
//!
//! Both functions start with the two `solution == nullptr` / `!OK` guards
//! (`cxx:1259-1266`, `:1298-1305`) and return `gp_XYZ(0, 0, 0)`.

use occt_core::gp::{GpXY, GpXyz};

use super::Plate;

impl Plate {
    /// `Evaluate(point2d)` (`cxx:1257-1291`).
    pub fn evaluate(&self, point2d: GpXY) -> GpXyz {
        if self.solution.is_empty() {
            return GpXyz::zero();
        }
        if !self.ok {
            return GpXyz::zero();
        }

        let mut valeur = GpXyz::zero();

        if !self.polynomial_part_only {
            for i in 0..self.n_el {
                let mut signe = 1.0;
                if (self.deru(i) + self.derv(i)) % 2 == 1 {
                    signe = -1.0;
                }
                valeur = valeur.added(
                    &self
                        .solution_at(i)
                        .multiplied(signe * self.sol_em(point2d.subtracted(&self.point_at(i)), self.deru(i), self.derv(i))),
                );
            }
        }
        let mut i = self.n_el;
        for idu in 0..self.order {
            let mut idv = 0;
            while idu + idv < self.order {
                valeur = valeur.added(&self.solution_at(i).multiplied(self.polm(point2d, idu, idv, 0, 0)));
                i += 1;
                idv += 1;
            }
        }
        valeur
    }

    /// `EvaluateDerivative(point2d, iu, iv)` (`cxx:1295-1331`).
    pub fn evaluate_derivative(&self, point2d: GpXY, iu: i32, iv: i32) -> GpXyz {
        if self.solution.is_empty() {
            return GpXyz::zero();
        }
        if !self.ok {
            return GpXyz::zero();
        }

        let mut valeur = GpXyz::zero();
        if !self.polynomial_part_only {
            for i in 0..self.n_el {
                let mut signe = 1.0;
                if (self.deru(i) + self.derv(i)) % 2 == 1 {
                    signe = -1.0;
                }
                valeur = valeur.added(
                    &self.solution_at(i).multiplied(signe * self.sol_em(point2d.subtracted(&self.point_at(i)), self.deru(i) + iu, self.derv(i) + iv)),
                );
            }
        }
        let mut i = self.n_el;
        for idu in 0..self.order {
            let mut idv = 0;
            while idu + idv < self.order {
                valeur = valeur.added(&self.solution_at(i).multiplied(self.polm(point2d, idu, idv, iu, iv)));
                i += 1;
                idv += 1;
            }
        }
        valeur
    }
}
