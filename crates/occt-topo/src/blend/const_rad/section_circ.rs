//! `BlendFunc_ConstRad::Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstRad.cxx:1141-1199`
//! (OCCT 8.0.0). It is the circle of the simulated fillet section at one
//! blending point, used by `ChFi3d_FilBuilder::SimulSurf`
//! (`ChFi3d_FilBuilder.cxx:648`, `:707`) to fill `ChFiDS_CircSection`.
//!
//! `Pdeb` is always `0.` and `Pfin` is the circle parameter of `pts2`; the
//! two out-parameters become the returned pair.
//!
//! Divergences that are not heuristics:
//! - `C.SetRadius` raises `Standard_ConstructionError` on a negative radius
//!   and `gp_Ax2(Center, np, ns1)` raises when the two directions are parallel;
//!   both answer `Err` here,
//! - `gp_Circ` is immutable through the port's `GpCirc::new`, so the second
//!   `SetPosition` of the `Pfin > 1.5*M_PI` arm rebuilds the circle instead of
//!   mutating it (same `Center` / radius / `ns1`).

use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpCirc, GpDir, GpVec};
use occt_core::precision::PCONFUSION;
use occt_math::MathVector;

use super::constants::EPS;
use super::ConstRad;

impl ConstRad {
    /// `Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)` (`cxx:1141-1199`).
    ///
    /// `Ok((C, Pdeb, Pfin))` is OCCT's out-triple.
    pub fn section_circ(
        &mut self,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> Result<(GpCirc, f64, f64), String> {
        // `math_Vector X(1,4); X(1) = U1; ... ComputeValues(X, 0, true, prm);`
        // (`cxx:1152-1160`).
        let mut x = MathVector::new(1, 4);
        x.set_value(1, u1);
        x.set_value(2, v1);
        x.set_value(3, u2);
        x.set_value(4, v2);
        let prm = param;
        self.compute_values(&x, 0, true, prm);

        // `ns1 = nsurf1; np = nplan;` (`cxx:1162-1163`).
        let mut ns1 = self.nsurf1;
        let mut np = self.nplan;

        // `norm1 = nplan.Crossed(ns1).Magnitude(); if (norm1 < Eps) norm1 = 1;`
        // (`cxx:1166-1170`).
        let mut norm1 = self.nplan.crossed(&ns1).magnitude();
        if norm1 < EPS {
            norm1 = 1.0;
        }
        // `ns1.SetLinearForm(nplan.Dot(ns1)/norm1, nplan, -1./norm1, ns1);`
        // (`cxx:1171`). The argument is the *old* `ns1`, so it is copied out
        // first, exactly like `gp_Vec::SetLinearForm`'s locals.
        let old_ns1 = ns1;
        ns1 = ns1.set_linear_form_2(
            self.nplan.dot(&ns1) / norm1,
            &np,
            -1.0 / norm1,
            &old_ns1,
        );

        // `Center.SetXYZ(pts1.XYZ() + ray1 * ns1.XYZ());` (`cxx:1172`).
        let center = self.pts1.translated_vec(&ns1.multiplied_scalar(self.ray1));

        // `if (ray1 > 0.) ns1.Reverse(); if (choix % 2 != 0) np.Reverse();`
        // (`cxx:1178-1182`).
        if self.ray1 > 0.0 {
            ns1.reverse();
        }
        if self.choix % 2 != 0 {
            np.reverse();
        }

        // `C.SetRadius(abs(ray1)); C.SetPosition(gp_Ax2(Center, np, ns1));`
        // (`cxx:1184-1185`).
        let circ = make_circ(center, &np, &ns1, self.ray1.abs())?;

        // `Pdeb = 0.; Pfin = ElCLib::Parameter(C, pts2);` (`cxx:1186-1187`).
        let mut pfin = clib::circle_parameter(&circ.position(), &self.pts2);

        // `if (Pfin > 1.5*M_PI) { np.Reverse(); C.SetPosition(...);
        // Pfin = ElCLib::Parameter(C, pts2); }` (`cxx:1189-1194`).
        if pfin > 1.5 * std::f64::consts::PI {
            np.reverse();
            let circ2 = make_circ(center, &np, &ns1, self.ray1.abs())?;
            pfin = clib::circle_parameter(&circ2.position(), &self.pts2);
            return Ok(finish(circ2, pfin));
        }

        Ok(finish(circ, pfin))
    }
}

/// `if (Pfin < Precision::PConfusion()) Pfin += Precision::PConfusion();`
/// (`cxx:1195-1198`) and `Pdeb = 0.`.
fn finish(circ: GpCirc, mut pfin: f64) -> (GpCirc, f64, f64) {
    if pfin < PCONFUSION {
        pfin += PCONFUSION;
    }
    (circ, 0.0, pfin)
}

/// `C.SetRadius(R); C.SetPosition(gp_Ax2(Center, np, ns1));` — the two
/// exceptions of `gp_Circ` / `gp_Ax2` doubled into one `Err`.
fn make_circ(center: occt_core::gp::GpPnt, np: &GpVec, ns1: &GpVec, radius: f64) -> Result<GpCirc, String> {
    let z = GpDir::from_vec(np).map_err(|_| {
        "BlendFunc_ConstRad::Section: gp_Ax2 raised on a null normal".to_string()
    })?;
    let x = GpDir::from_vec(ns1)
        .map_err(|_| "BlendFunc_ConstRad::Section: gp_Ax2 raised on a null normal".to_string())?;
    let ax2 = GpAx2::new(center, z, x).map_err(|_| {
        "BlendFunc_ConstRad::Section: gp_Ax2 raised on parallel directions".to_string()
    })?;
    if radius < 0.0 {
        return Err("BlendFunc_ConstRad::Section: gp_Circ::SetRadius raised".to_string());
    }
    Ok(GpCirc::new(ax2, radius))
}
