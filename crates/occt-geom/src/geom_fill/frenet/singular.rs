//! The singular-point machinery of `GeomFill_Frenet`:
//! `IsSingular` / `DoSingular` / `SingularD0` / `SingularD1` / `SingularD2`.
//!
//! Source: `GeomFill_Frenet.cxx:887-1107` (OCCT 8.0.0).
//!
//! `DoSingular` reports its `Delta` out-parameter together with the optional
//! frame: OCCT's `false` return still leaves `Delta = h` for the caller
//! (`cxx:966-969`), and `D0` / `D1` / `D2` resume at `Param + Delta`.

use occt_core::gp::GpVec;
use occt_core::precision::CONFUSION;

use crate::geom_fill::frenet::deriv;
use crate::geom_fill::frenet::Frenet;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw};

/// The frame `DoSingular` derives from the leading non-zero derivatives
/// (`GeomFill_Frenet.cxx:907-985`).
pub(crate) struct DoSingular {
    /// `Tangent` (`cxx:909`).
    pub tangent: GpVec,
    /// `BiNormal` (`cxx:910`).
    pub bi_normal: GpVec,
    /// `n` (`cxx:911`): the derivative order that gives the tangent.
    pub n: i32,
    /// `k` (`cxx:912`): the derivative order that gives the bi-normal.
    pub k: i32,
    /// `TFlag` (`cxx:913`).
    pub t_flag: f64,
    /// `BNFlag` (`cxx:914`).
    pub bn_flag: f64,
}

impl Frenet {
    /// `GeomFill_Frenet::IsSingular` (`cxx:887-903`): the index of the
    /// singular interval containing `u`.
    pub(crate) fn is_singular(&self, u: f64) -> Option<usize> {
        let points = self.singular_points()?;
        for (i, &s) in points.sngl.iter().enumerate() {
            if (u - s).abs() < points.len[i] {
                return Some(i);
            }
        }
        None
    }

    /// `GeomFill_Frenet::DoSingular` (`cxx:907-985`).
    ///
    /// The second element is `Delta`: `h` when the search failed, `0` when the
    /// frame was found.
    pub(crate) fn do_singular(
        &self,
        u: f64,
        index: usize,
    ) -> Result<(Option<DoSingular>, f64), String> {
        let points = self
            .singular_points()
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_Frenet::DoSingular".to_string())?;
        let len = *points
            .len
            .get(index)
            .ok_or_else(|| "Standard_OutOfRange: GeomFill_Frenet::DoSingular".to_string())?;
        let curve = self.trimmed()?;

        const MAX_N: i32 = 20;
        let mut h = 2.0 * len;
        let (a, b) = self.data.interval();
        if u >= (a + b) / 2.0 {
            h = -h;
        }

        let mut i = 1;
        let mut tangent = GpVec::zero();
        while i <= MAX_N {
            tangent = curve.eval_dn(u, i);
            if tangent.magnitude() > CONFUSION {
                break;
            }
            i += 1;
        }
        if i > MAX_N {
            return Ok((None, 0.0));
        }
        tangent.normalize();
        let n = i;

        let mut bi_normal = GpVec::zero();
        let mut k = 0;
        i += 1;
        while i <= MAX_N {
            bi_normal = tangent.crossed(&curve.eval_dn(u, i));
            let magn = bi_normal.magnitude();
            if magn > CONFUSION {
                // modified by jgv, 12.08.03 for OCC605 (`cxx:953-962`)
                let next_bi_normal = tangent.crossed(&curve.eval_dn(u, i + 1));
                if next_bi_normal.magnitude() > magn {
                    i += 1;
                    bi_normal = next_bi_normal;
                }
                k = i;
                break;
            }
            i += 1;
        }
        if k == 0 {
            return Ok((None, h));
        }
        bi_normal.normalize();

        let (t, _n, bn) = self.d0(u + h)?;
        let mut t_flag = 1.0;
        let mut bn_flag = 1.0;
        if tangent.angle(&t) > std::f64::consts::PI / 2.0 {
            t_flag = -1.0;
        }
        if bi_normal.angle(&bn) > std::f64::consts::PI / 2.0 {
            bn_flag = -1.0;
        }

        Ok((
            Some(DoSingular {
                tangent,
                bi_normal,
                n,
                k,
                t_flag,
                bn_flag,
            }),
            0.0,
        ))
    }

    /// `GeomFill_Frenet::SingularD0` (`cxx:987-1007`).
    pub(crate) fn singular_d0(
        &self,
        param: f64,
        index: usize,
    ) -> Result<(Option<(GpVec, GpVec, GpVec)>, f64), String> {
        let (frame, delta) = self.do_singular(param, index)?;
        let Some(ds) = frame else {
            return Ok((None, delta));
        };
        let tangent = ds.tangent.multiplied_scalar(ds.t_flag);
        let bi_normal = ds.bi_normal.multiplied_scalar(ds.bn_flag);
        let normal = bi_normal.crossed(&tangent);
        Ok((Some((tangent, normal, bi_normal)), delta))
    }

    /// `GeomFill_Frenet::SingularD1` (`cxx:1009-1051`).
    pub(crate) fn singular_d1(
        &self,
        param: f64,
        index: usize,
    ) -> Result<(Option<TrihedronD1>, f64), String> {
        let (frame, delta) = self.do_singular(param, index)?;
        let Some(ds) = frame else {
            return Ok((None, delta));
        };
        let curve = self.trimmed()?;
        let f = curve.eval_dn(param, ds.n);
        let df = curve.eval_dn(param, ds.n + 1);
        let mut d_tangent = deriv::fderiv(&f, &df);

        let dtmp = curve.eval_dn(param, ds.k);
        let f = ds.tangent.crossed(&dtmp);
        let df = d_tangent
            .crossed(&dtmp)
            .added(&ds.tangent.crossed(&curve.eval_dn(param, ds.k + 1)));
        let mut d_bi_normal = deriv::fderiv(&f, &df);

        let mut tangent = ds.tangent;
        let mut bi_normal = ds.bi_normal;
        if ds.t_flag < 0.0 {
            tangent = tangent.reversed();
            d_tangent = d_tangent.reversed();
        }
        if ds.bn_flag < 0.0 {
            bi_normal = bi_normal.reversed();
            d_bi_normal = d_bi_normal.reversed();
        }
        let normal = bi_normal.crossed(&tangent);
        let d_normal = d_bi_normal.crossed(&tangent).added(&bi_normal.crossed(&d_tangent));
        Ok((
            Some(TrihedronD1 {
                tangent,
                dtangent: d_tangent,
                normal,
                dnormal: d_normal,
                binormal: bi_normal,
                dbinormal: d_bi_normal,
            }),
            delta,
        ))
    }

    /// `GeomFill_Frenet::SingularD2` (`cxx:1053-1107`).
    pub(crate) fn singular_d2(
        &self,
        param: f64,
        index: usize,
    ) -> Result<(Option<TrihedronD2>, f64), String> {
        let (frame, delta) = self.do_singular(param, index)?;
        let Some(ds) = frame else {
            return Ok((None, delta));
        };
        let curve = self.trimmed()?;
        let f = curve.eval_dn(param, ds.n);
        let df = curve.eval_dn(param, ds.n + 1);
        let d2f = curve.eval_dn(param, ds.n + 2);
        let mut d_tangent = deriv::fderiv(&f, &df);
        let mut d2_tangent = deriv::dderiv(&f, &df, &d2f);

        let dtmp1 = curve.eval_dn(param, ds.k);
        let dtmp2 = curve.eval_dn(param, ds.k + 1);
        let f = ds.tangent.crossed(&dtmp1);
        let df = d_tangent.crossed(&dtmp1).added(&ds.tangent.crossed(&dtmp2));
        let d2f = d2_tangent
            .crossed(&dtmp1)
            .added(&d_tangent.crossed(&dtmp2).multiplied_scalar(2.0))
            .added(&ds.tangent.crossed(&curve.eval_dn(param, ds.k + 2)));
        let mut d_bi_normal = deriv::fderiv(&f, &df);
        let mut d2_bi_normal = deriv::dderiv(&f, &df, &d2f);

        let mut tangent = ds.tangent;
        let mut bi_normal = ds.bi_normal;
        if ds.t_flag < 0.0 {
            tangent = tangent.reversed();
            d_tangent = d_tangent.reversed();
            d2_tangent = d2_tangent.reversed();
        }
        if ds.bn_flag < 0.0 {
            bi_normal = bi_normal.reversed();
            d_bi_normal = d_bi_normal.reversed();
            d2_bi_normal = d2_bi_normal.reversed();
        }
        let normal = bi_normal.crossed(&tangent);
        let d_normal = d_bi_normal.crossed(&tangent).added(&bi_normal.crossed(&d_tangent));
        let d2_normal = d2_bi_normal
            .crossed(&tangent)
            .added(&d_bi_normal.crossed(&d_tangent).multiplied_scalar(2.0))
            .added(&bi_normal.crossed(&d2_tangent));
        Ok((
            Some(TrihedronD2 {
                tangent,
                dtangent: d_tangent,
                d2tangent: d2_tangent,
                normal,
                dnormal: d_normal,
                d2normal: d2_normal,
                binormal: bi_normal,
                dbinormal: d_bi_normal,
                d2binormal: d2_bi_normal,
            }),
            delta,
        ))
    }
}
