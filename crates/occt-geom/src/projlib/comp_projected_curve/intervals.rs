//! `ProjLib_CompProjectedCurve::NbIntervals` / `Intervals` / `BuildIntervals`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx:1810-2092`
//! (OCCT 8.0.0).
//!
//! The result is a fused parameter partition of the projected curve: the
//! projection bounds, the curve's own discontinuities, and the parameters at
//! which the projection crosses a `u` or `v` discontinuity of the surface.
//! Those crossings are solved with `ProjLib_PrjResolve` (`myFix` 2 for a `u`
//! discontinuity, 3 for a `v` one).
//!
//! Translation notes:
//! - `GeomLib::FuseIntervals(..., Precision::PConfusion())` (`cxx:2047-2084`)
//!   is [`crate::geom_lib::fuse_intervals`], whose last argument is
//!   `IsAdjustToFirstInterval`; OCCT's 4-argument overload passes `false`.
//! - `NCollection_Sequence::Remove(i--)` on adjacent duplicates
//!   (`cxx:1948-1953`, `cxx:2028-2033`) is a `dedup`-style in-place loop.
//! - `throw Standard_OutOfRange()` for a `GeomAbs_Shape` outside
//!   `C0..CN` (`cxx:1852-1853`) is `Err`.

use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::geom_lib::fuse_intervals;
use crate::projlib::prj_resolve::PrjResolve;

use super::jets::d1;
use super::{CompProjectedCurve, FUNC_TOL};

impl CompProjectedCurve {
    /// `NbIntervals(S)` (`cxx:1811-1818`): nullifies `myTabInt` and rebuilds.
    pub fn nb_intervals(&mut self, shape: u8) -> Result<i32, String> {
        self.tab_int = None;
        let tab = self.build_intervals(shape)?;
        self.tab_int = Some(tab);
        Ok(self.tab_int.as_ref().map(|t| t.len() as i32 - 1).unwrap_or(0))
    }

    /// `Intervals(T, S)` (`cxx:1822-1828`): builds the table on demand.
    pub fn intervals(&mut self, shape: u8) -> Result<&[f64], String> {
        if self.tab_int.is_none() {
            let tab = self.build_intervals(shape)?;
            self.tab_int = Some(tab);
        }
        Ok(self.tab_int.as_deref().unwrap_or(&[]))
    }

    /// `BuildIntervals(S)` (`cxx:1833-2093`).
    fn build_intervals(&self, shape: u8) -> Result<Vec<f64>, String> {
        // `SforS` (`cxx:1835-1856`).
        let sfor_s: u8 = match shape {
            0 => Shape::C1 as u8,
            2 => Shape::C2 as u8,
            4 => Shape::C3 as u8,
            5 => Shape::CN as u8,
            6 => Shape::CN as u8,
            _ => return Err("ProjLib_CompProjectedCurve: unknown GeomAbs_Shape (cxx:1852)".into()),
        };

        let (Some(surface), Some(curve)) = (self.surface.as_ref(), self.curve.as_ref()) else {
            return Ok(Vec::new());
        };
        let s = surface.as_ref();
        let c = curve.as_ref();

        let nb_int_cur = c.nb_intervals(shape);
        let nb_int_sur_u = s.nb_u_intervals(sfor_s);
        let nb_int_sur_v = s.nb_v_intervals(sfor_s);

        let cut_pnts_t = c.parameter_intervals(shape);
        let cut_pnts_u = s.u_intervals(sfor_s);
        let cut_pnts_v = s.v_intervals(sfor_s);

        // `Tl`, `Tr`, `Ul`, `Ur`, `Vl`, `Vr`, `Tol` are assigned per arm below.
        let mut tol;

        // Processing projection bounds (`cxx:1872-1878`).
        let mut barr: Vec<f64> = Vec::with_capacity(2 * self.nb_curves as usize);
        for i in 1..=self.nb_curves {
            if let Some((udeb, ufin)) = self.bounds(i) {
                barr.push(udeb);
                barr.push(ufin);
            }
        }

        // Processing curve discontinuities (`cxx:1881-1888`).
        let carr: Option<Vec<f64>> = if nb_int_cur > 1 {
            Some(
                cut_pnts_t
                    .iter()
                    .skip(1)
                    .take(nb_int_cur as usize - 1)
                    .copied()
                    .collect(),
            )
        } else {
            None
        };

        // Processing U-surface discontinuities (`cxx:1890-1965`).
        let mut tu_disc: Vec<f64> = Vec::new();
        for k in 2..=nb_int_sur_u as usize {
            let cut = cut_pnts_u
                .get(k - 1)
                .copied()
                .ok_or_else(|| "BuildIntervals: UIntervals shorter than NbUIntervals".to_string())?;
            for i in 1..=self.nb_curves {
                let part = self
                    .part(i)
                    .ok_or_else(|| "BuildIntervals: missing part".to_string())?;
                for j in 1..part.len() {
                    let ul = part[j - 1].y();
                    let ur = part[j].y();
                    if (ul - cut).abs() <= self.tol_u {
                        tu_disc.push(part[j - 1].x());
                    } else if (ur - cut).abs() <= self.tol_u {
                        tu_disc.push(part[j].x());
                    } else if (ul < cut && cut < ur) || (ur < cut && cut < ul) {
                        let v = (part[j - 1].z() + part[j].z()) / 2.0;
                        let mut solver =
                            PrjResolve::new(c, s, 2).expect("myFix in 1..=3");

                        let triple = part[j - 1];
                        let d = d1(triple.x(), triple.y(), triple.z(), c, s).expect("d1");
                        tol = if d.x().abs() < CONFUSION {
                            self.tol_u
                        } else {
                            self.tol_u.min(self.tol_u / d.x().abs())
                        };

                        let tl = part[j - 1].x();
                        let tr = part[j].x();

                        solver.perform(
                            (tl + tr) / 2.0,
                            cut,
                            v,
                            (tol, self.tol_v),
                            (tl, s.v_range().0),
                            (tr, s.v_range().1),
                            FUNC_TOL,
                            false,
                        );
                        if solver.is_done() {
                            tu_disc.push(solver.solution().expect("IsDone").0);
                        }
                    }
                }
            }
        }
        // `Remove(i--)` on adjacent duplicates (`cxx:1948-1953`).
        let mut i = 1usize;
        while i < tu_disc.len() {
            if tu_disc[i] - tu_disc[i - 1] < PCONFUSION {
                tu_disc.remove(i);
            } else {
                i += 1;
            }
        }
        let uarr: Option<Vec<f64>> = if tu_disc.is_empty() {
            None
        } else {
            Some(tu_disc)
        };

        // Processing V-surface discontinuities (`cxx:1967-2043`).
        let mut tv_disc: Vec<f64> = Vec::new();
        for k in 2..=nb_int_sur_v as usize {
            let cut = cut_pnts_v
                .get(k - 1)
                .copied()
                .ok_or_else(|| "BuildIntervals: VIntervals shorter than NbVIntervals".to_string())?;
            for i in 1..=self.nb_curves {
                let part = self
                    .part(i)
                    .ok_or_else(|| "BuildIntervals: missing part".to_string())?;
                for j in 1..part.len() {
                    let vl = part[j - 1].z();
                    let vr = part[j].z();
                    if (vl - cut).abs() <= self.tol_v {
                        tv_disc.push(part[j - 1].x());
                    } else if (vr - cut).abs() <= self.tol_v {
                        tv_disc.push(part[j].x());
                    } else if (vl < cut && cut < vr) || (vr < cut && cut < vl) {
                        let u = (part[j - 1].y() + part[j].y()) / 2.0;
                        let mut solver = PrjResolve::new(c, s, 3).expect("myFix in 1..=3");

                        let triple = part[j - 1];
                        let d = d1(triple.x(), triple.y(), triple.z(), c, s).expect("d1");
                        tol = if d.y().abs() < CONFUSION {
                            self.tol_v
                        } else {
                            self.tol_v.min(self.tol_v / d.y().abs())
                        };

                        let tl = part[j - 1].x();
                        let tr = part[j].x();

                        solver.perform(
                            (tl + tr) / 2.0,
                            u,
                            cut,
                            (tol, self.tol_v),
                            (tl, s.u_range().0),
                            (tr, s.u_range().1),
                            FUNC_TOL,
                            false,
                        );
                        if solver.is_done() {
                            tv_disc.push(solver.solution().expect("IsDone").0);
                        }
                    }
                }
            }
        }
        let mut i = 1usize;
        while i < tv_disc.len() {
            if tv_disc[i] - tv_disc[i - 1] < PCONFUSION {
                tv_disc.remove(i);
            } else {
                i += 1;
            }
        }
        let varr: Option<Vec<f64>> = if tv_disc.is_empty() {
            None
        } else {
            Some(tv_disc)
        };

        // Fusion (`cxx:2046-2084`).
        if let Some(carr) = carr {
            barr = fuse_intervals(&barr, &carr, PCONFUSION, false);
        }
        if let Some(uarr) = uarr {
            barr = fuse_intervals(&barr, &uarr, PCONFUSION, false);
        }
        if let Some(varr) = varr {
            barr = fuse_intervals(&barr, &varr, PCONFUSION, false);
        }
        Ok(barr)
    }
}
