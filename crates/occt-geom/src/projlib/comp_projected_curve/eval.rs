//! `ProjLib_CompProjectedCurve` evaluation: `Value`, `D0`, `D1`, `D2`, `DN`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx:1551-1755`
//! (OCCT 8.0.0).
//!
//! `D0` locates the continuous part containing the parameter, interpolates the
//! cached `(u, v)` samples with the cubic through the four neighbouring
//! triples, and hands that guess to `ProjLib_PrjResolve`. When the exact
//! projection fails it retries from an `Extrema_ExtPS` minimum of the curve
//! point, and as a last resort answers the interpolated guess.
//!
//! Translation notes:
//! - `Standard_DomainError("ProjLib_CompProjectedCurve::D0")` (`cxx:1578`) is
//!   `Err` here; an out-of-range index is `Standard_NoSuchObject`.
//! - `Standard_OutOfRange` / `Standard_NotImplemented` in `DN`
//!   (`cxx:1731-1753`) are `Err` too.
//! - `gp_Vec2d` arithmetic is spelled out (`Y1 - Y2` → `subtracted`, `/` →
//!   `divided`).

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::PCONFUSION;

use crate::extrema_surf::{ExtPs, ExtremaExtFlag};
use crate::projlib::prj_resolve::PrjResolve;

use super::jets::{d1, d2};
use super::{CompProjectedCurve, FUNC_TOL};

impl CompProjectedCurve {
    /// `Bounds(Index, Udeb, Ufin)` (`cxx:1503-1513`) as the `(first, last)`
    /// parameter pair of part `index`, `None` for `Standard_NoSuchObject`.
    fn part_bounds(&self, index: i32) -> Option<(f64, f64)> {
        self.bounds(index)
    }

    /// `Value(t)` (`cxx:1553-1558`).
    pub fn value(&self, t: f64) -> Result<GpPnt2d, String> {
        self.d0(t)
    }

    /// `D0(U, P)` (`cxx:1563-1704`).
    pub fn d0(&self, u: f64) -> Result<GpPnt2d, String> {
        let (Some(surface), Some(curve)) = (self.surface.as_ref(), self.curve.as_ref()) else {
            return Err("ProjLib_CompProjectedCurve::D0: no surface or curve".into());
        };
        let s = surface.as_ref();
        let c = curve.as_ref();

        let mut found = None;
        for i in 1..=self.nb_curves {
            if let Some((udeb, ufin)) = self.part_bounds(i) {
                if u >= udeb && u <= ufin {
                    found = Some(i);
                    break;
                }
            }
        }
        let Some(i) = found else {
            return Err("ProjLib_CompProjectedCurve::D0".into());
        };
        let part = self
            .part(i)
            .ok_or_else(|| "ProjLib_CompProjectedCurve::D0".to_string())?;

        let end = part.len();
        let mut j = 1usize;
        while j < end {
            if u >= part[j - 1].x() && u <= part[j].x() {
                break;
            }
            j += 1;
        }

        let mut u0;
        let mut v0;
        // Cubic interpolation (`cxx:1593-1662`).
        if part.len() < 4 || (u - part[j - 1].x()).abs() <= PCONFUSION {
            u0 = part[j - 1].y();
            v0 = part[j - 1].z();
        } else if (u - part[j].x()).abs() <= PCONFUSION {
            u0 = part[j].y();
            v0 = part[j].z();
        } else {
            if j == 1 {
                j = 2;
            }
            if j > part.len() - 2 {
                j = part.len() - 2;
            }

            let x1 = part[j - 2].x();
            let x2 = part[j - 1].x();
            let x3 = part[j].x();
            let x4 = part[j + 1].x();

            let y1 = GpVec2d::new(part[j - 2].y(), part[j - 2].z());
            let y2 = GpVec2d::new(part[j - 1].y(), part[j - 1].z());
            let y3 = GpVec2d::new(part[j].y(), part[j].z());
            let y4 = GpVec2d::new(part[j + 1].y(), part[j + 1].z());

            let i1 = y1.subtracted(&y2).divided(x1 - x2);
            let i2 = y2.subtracted(&y3).divided(x2 - x3);
            let i3 = y3.subtracted(&y4).divided(x3 - x4);

            let i21 = i1.subtracted(&i2).divided(x1 - x3);
            let i22 = i2.subtracted(&i3).divided(x2 - x4);

            let i31 = i21.subtracted(&i22).divided(x1 - x4);

            // `Res = Y1 + (U - X1) * (I1 + (U - X2) * (I21 + (U - X3) * I31))`
            let res = y1.added(
                &i1.added(&i21.added(&i31.multiplied_scalar(u - x3)).multiplied_scalar(u - x2))
                    .multiplied_scalar(u - x1),
            );

            u0 = res.x();
            v0 = res.y();

            if u0 < s.u_range().0 {
                u0 = s.u_range().0;
            } else if u0 > s.u_range().1 {
                u0 = s.u_range().1;
            }
            if v0 < s.v_range().0 {
                v0 = s.v_range().0;
            } else if v0 > s.v_range().1 {
                v0 = s.v_range().1;
            }
        }
        // End of cubic interpolation

        let mut prj = PrjResolve::new(c, s, 1).expect("myFix in 1..=3");
        prj.perform(
            u,
            u0,
            v0,
            (self.tol_u, self.tol_v),
            (s.u_range().0, s.v_range().0),
            (s.u_range().1, s.v_range().1),
            FUNC_TOL,
            false,
        );
        if prj.is_done() {
            let (su, sv) = prj.solution().expect("IsDone");
            return Ok(GpPnt2d::new(su, sv));
        }

        // Fallback: the nearest `Extrema_ExtPS` solution (`cxx:1679-1702`).
        let the_point = c.d0(u);
        let mut ext = ExtPs::new();
        ext.initialize(s, s.u_range().0, s.u_range().1, s.v_range().0, s.v_range().1, self.tol_u, self.tol_v);
        ext.set_flag(ExtremaExtFlag::Min);
        ext.perform(&the_point);
        if ext.is_done() && ext.nb_ext() > 0 {
            let nend = ext.nb_ext();
            let mut imin = 1usize;
            for k in 2..=nend {
                if ext.square_distance(k) < ext.square_distance(imin) {
                    imin = k;
                }
            }
            let (par_u, par_v, _) = ext.point(imin);
            return Ok(GpPnt2d::new(par_u, par_v));
        }
        Ok(GpPnt2d::new(u0, v0))
    }

    /// `D1(t, P, V)` (`cxx:1708-1717`).
    pub fn d1(&self, t: f64) -> Result<(GpPnt2d, GpVec2d), String> {
        let p = self.d0(t)?;
        let (Some(surface), Some(curve)) = (self.surface.as_ref(), self.curve.as_ref()) else {
            return Err("ProjLib_CompProjectedCurve::D1: no surface or curve".into());
        };
        let v = d1(t, p.x(), p.y(), curve.as_ref(), surface.as_ref())?;
        Ok((p, v))
    }

    /// `D2(t, P, V1, V2)` (`cxx:1720-1727`).
    pub fn d2(&self, t: f64) -> Result<(GpPnt2d, GpVec2d, GpVec2d), String> {
        let p = self.d0(t)?;
        let (Some(surface), Some(curve)) = (self.surface.as_ref(), self.curve.as_ref()) else {
            return Err("ProjLib_CompProjectedCurve::D2: no surface or curve".into());
        };
        let (v1, v2) = d2(t, p.x(), p.y(), curve.as_ref(), surface.as_ref())?;
        Ok((p, v1, v2))
    }

    /// `DN(t, N)` (`cxx:1731-1755`): `N < 1` is `Standard_OutOfRange`, `N > 2`
    /// is `Standard_NotImplemented`.
    pub fn dn(&self, t: f64, n: i32) -> Result<GpVec2d, String> {
        if n < 1 {
            return Err("ProjLib_CompProjectedCurve : N must be greater than 0".into());
        } else if n == 1 {
            Ok(self.d1(t)?.1)
        } else if n == 2 {
            Ok(self.d2(t)?.2)
        } else {
            Err("ProjLib_CompProjectedCurve::DN".into())
        }
    }
}
