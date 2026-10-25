//! `GeomFill_Darboux` — the Darboux trihedron: the frame of a curve lying on
//! a surface, with the surface normal as the binormal.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Darboux.hxx` and
//! `GeomFill_Darboux.cxx:39-561` (OCCT 8.0.0): `FDeriv` `cxx:39-44`,
//! `DDeriv` `cxx:49-58`, `NormalD0` `cxx:62-148`, `NormalD1` `cxx:152-245`,
//! `NormalD2` `cxx:249-349`, `D0` `cxx:365-392`, `D1` `cxx:394-433`,
//! `D2` `cxx:435-511`, `NbIntervals` `cxx:513-516`, `Intervals` `cxx:518-521`,
//! `GetAverageLaw` `cxx:523-551`, `IsConstant` `cxx:553-556`,
//! `IsOnlyBy3dCurve` `cxx:558-561`.
//!
//! Translation notes:
//! - OCCT stores one `Adaptor3d_CurveOnSurface myTrimmed` and recovers the two
//!   halves with `static_cast<Adaptor3d_CurveOnSurface*>(myTrimmed.get())`
//!   (`cxx:369-372`). The port takes the bundle apart the same way
//!   [`crate::geom_fill::BoundWithSurf`] does: the 3D curve is `myCurve` /
//!   `myTrimmed`, and the pcurve + surface are the constructor arguments that
//!   the caller built the adaptor from.
//! - `CSLib::Normal` (`cxx:82`, `cxx:120`, `cxx:174`) is
//!   [`occt_core::cslib::normal_d1_mag`] / [`occt_core::cslib::normal_max_order`];
//!   `CSLib::DNNUV` and `CSLib::DNNormal` are
//!   [`occt_core::cslib::dnnuv`] / [`occt_core::cslib::dn_normal`].
//! - `Surf->DN(U, V, i, j)` is [`Surface::try_eval_dn`], so a surface class
//!   without a ported `EvalDN` makes the law fail with that error instead of
//!   reading the zero default of [`Surface::eval_dn`].
//! - `Geom_UndefinedValue` / `Geom_UndefinedDerivative` from the normal
//!   helpers are `Err` values.
//! - `Adaptor3d_Surface::UContinuity` / `VContinuity` used by the two
//!   continuity guards is [`Surface::continuity`] (the port keeps one value).

use std::sync::Arc;

use occt_core::cslib::{dn_normal, dnnuv, normal_d1_mag, normal_max_order, CSLibNormalStatus};
use occt_core::gp::{GpDir, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_geom2d::curve::Curve2d;

use crate::curve::Curve;
use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::geom_fill::frenet::deriv::{dderiv, fderiv};
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};
use crate::surface::Surface;

/// `MagTol` of the three normal helpers (`cxx:81`, `cxx:173`, `cxx:271`).
const MAG_TOL: f64 = 0.000000001;
/// `MaxOrder` of the singular-normal fallbacks (`cxx:90`, `cxx:183`, `cxx:279`).
const MAX_ORDER: i32 = 3;

/// `NormalD0` (`cxx:62-148`): the surface normal at `(U, V)`, with the
/// `CSLib::Normal` order fallback when `D1` is degenerate.
fn normal_d0(u: f64, v: f64, surf: &dyn Surface) -> Result<GpDir, String> {
    let cont = surf.continuity();
    let (_p, d1u, d1v) = surf.d1(u, v);
    let (status, normal) = normal_d1_mag(&d1u, &d1v, MAG_TOL);
    if status == CSLibNormalStatus::Defined {
        return normal.ok_or_else(|| "Geom_UndefinedValue: NormalD0".to_string());
    }
    if cont == Shape::C0 as u8 || cont == Shape::C1 as u8 {
        return Err("Geom_UndefinedValue: NormalD0".into());
    }

    let n = (MAX_ORDER + 2) as usize;
    let mut der_surf = vec![vec![GpVec::zero(); n]; n];
    for i in 1..=MAX_ORDER + 1 {
        der_surf[i as usize][0] = surf.try_eval_dn(u, v, i, 0)?;
    }
    for i in 0..=MAX_ORDER + 1 {
        for j in 1..=MAX_ORDER + 1 {
            der_surf[i as usize][j as usize] = surf.try_eval_dn(u, v, i, j)?;
        }
    }

    let m = (MAX_ORDER + 1) as usize;
    let mut der_nuv = vec![vec![GpVec::zero(); m]; m];
    for i in 0..=MAX_ORDER {
        for j in 0..=MAX_ORDER {
            der_nuv[i as usize][j as usize] = dnnuv(i, j, &der_surf);
        }
    }

    max_order_normal(u, v, surf, &der_nuv)
}

/// The shared tail of the three helpers: `CSLib::Normal(MaxOrder, ...)` plus
/// the status check (`cxx:120-147`, `cxx:222-239`, `cxx:327-343`).
fn max_order_normal(
    u: f64,
    v: f64,
    surf: &dyn Surface,
    der_nuv: &[Vec<GpVec>],
) -> Result<GpDir, String> {
    let max_order = MAX_ORDER;
    let (u_min, u_max) = surf.u_range();
    let (v_min, v_max) = surf.v_range();
    let (status, normal, _order_u, _order_v) = normal_max_order(
        max_order, der_nuv, MAG_TOL, u, v, u_min, u_max, v_min, v_max,
    );
    if status != CSLibNormalStatus::Defined {
        return Err("Geom_UndefinedValue: CSLib::Normal".into());
    }
    normal.ok_or_else(|| "Geom_UndefinedValue: CSLib::Normal".to_string())
}

/// `NormalD1` (`cxx:152-245`): the normal and its first U/V derivatives.
fn normal_d1(u: f64, v: f64, surf: &dyn Surface) -> Result<(GpDir, GpVec, GpVec), String> {
    let (_p, d1u, d1v, d2u, d2v, d2uv) = surf.d2(u, v);
    // `CSLib::Normal` writes `Normal` here too (`cxx:174`), but the
    // `CSLib::Normal(MaxOrder, ...)` call below (`cxx:222-238`) overwrites it.
    let (status, _normal) = normal_d1_mag(&d1u, &d1v, MAG_TOL);
    let max_order = if status == CSLibNormalStatus::Defined { 0 } else { MAX_ORDER };

    let n = (max_order + 3) as usize;
    let mut der_surf = vec![vec![GpVec::zero(); n]; n];
    der_surf[1][0] = d1u;
    der_surf[0][1] = d1v;
    der_surf[1][1] = d2uv;
    der_surf[2][0] = d2u;
    der_surf[0][2] = d2v;
    for i in 0..=max_order + 1 {
        for j in i..=max_order + 2 {
            if i + j > 2 {
                let d = surf.try_eval_dn(u, v, i, j)?;
                der_surf[i as usize][j as usize] = d;
                if i != j {
                    der_surf[j as usize][i as usize] = surf.try_eval_dn(u, v, j, i)?;
                }
            }
        }
    }

    let m = (max_order + 2) as usize;
    let mut der_nuv = vec![vec![GpVec::zero(); m]; m];
    for i in 0..=max_order + 1 {
        for j in 0..=max_order + 1 {
            der_nuv[i as usize][j as usize] = dnnuv(i, j, &der_surf);
        }
    }

    let (u_min, u_max) = surf.u_range();
    let (v_min, v_max) = surf.v_range();
    let (status, normal, order_u, order_v) = normal_max_order(
        max_order, &der_nuv, MAG_TOL, u, v, u_min, u_max, v_min, v_max,
    );
    if status != CSLibNormalStatus::Defined {
        return Err("Geom_UndefinedValue: NormalD1".into());
    }
    let normal = normal.ok_or_else(|| "Geom_UndefinedValue: NormalD1".to_string())?;
    let d1u_normal = dn_normal(1, 0, &der_nuv, order_u, order_v);
    let d1v_normal = dn_normal(0, 1, &der_nuv, order_u, order_v);
    Ok((normal, d1u_normal, d1v_normal))
}

/// `NormalD2` (`cxx:249-349`): the normal and its first two U/V derivatives.
fn normal_d2(
    u: f64,
    v: f64,
    surf: &dyn Surface,
) -> Result<(GpDir, GpVec, GpVec, GpVec, GpVec, GpVec), String> {
    let (_p, d1u, d1v, d2u, d2v, d2uv) = surf.d2(u, v);
    let d3u = surf.try_eval_dn(u, v, 3, 0)?;
    let d3v = surf.try_eval_dn(u, v, 0, 3)?;
    let d3uuv = surf.try_eval_dn(u, v, 2, 1)?;
    let d3uvv = surf.try_eval_dn(u, v, 1, 2)?;

    // `CSLib::Normal` writes `Normal` here too (`cxx:273`), but the
    // `CSLib::Normal(MaxOrder, ...)` call below (`cxx:329-342`) overwrites it.
    let (status, _normal) = normal_d1_mag(&d1u, &d1v, MAG_TOL);
    let max_order = if status == CSLibNormalStatus::Defined { 0 } else { MAX_ORDER };

    let n = (max_order + 4) as usize;
    let mut der_surf = vec![vec![GpVec::zero(); n]; n];
    der_surf[1][0] = d1u;
    der_surf[0][1] = d1v;
    der_surf[1][1] = d2uv;
    der_surf[2][0] = d2u;
    der_surf[0][2] = d2v;
    der_surf[3][0] = d3u;
    der_surf[2][1] = d3uuv;
    der_surf[1][2] = d3uvv;
    der_surf[0][3] = d3v;
    for i in 0..=max_order + 2 {
        for j in i..=max_order + 3 {
            if i + j > 3 {
                let d = surf.try_eval_dn(u, v, i, j)?;
                der_surf[i as usize][j as usize] = d;
                if i != j {
                    der_surf[j as usize][i as usize] = surf.try_eval_dn(u, v, j, i)?;
                }
            }
        }
    }

    let m = (max_order + 3) as usize;
    let mut der_nuv = vec![vec![GpVec::zero(); m]; m];
    for i in 0..=max_order + 2 {
        for j in 0..=max_order + 2 {
            der_nuv[i as usize][j as usize] = dnnuv(i, j, &der_surf);
        }
    }

    let (u_min, u_max) = surf.u_range();
    let (v_min, v_max) = surf.v_range();
    let (status, normal, order_u, order_v) = normal_max_order(
        max_order, &der_nuv, MAG_TOL, u, v, u_min, u_max, v_min, v_max,
    );
    if status != CSLibNormalStatus::Defined {
        return Err("Geom_UndefinedValue: NormalD2".into());
    }
    let normal = normal.ok_or_else(|| "Geom_UndefinedValue: NormalD2".to_string())?;
    Ok((
        normal,
        dn_normal(1, 0, &der_nuv, order_u, order_v),
        dn_normal(0, 1, &der_nuv, order_u, order_v),
        dn_normal(2, 0, &der_nuv, order_u, order_v),
        dn_normal(0, 2, &der_nuv, order_u, order_v),
        dn_normal(1, 1, &der_nuv, order_u, order_v),
    ))
}

/// `GeomFill_Darboux` (`GeomFill_Darboux.hxx:30-93`).
pub struct Darboux {
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    data: TrihedronLawCurve,
    /// `myTrimmed->GetCurve()` (`cxx:369-370`) of the `CurveOnSurface` bundle.
    pcurve: Option<Arc<dyn Curve2d>>,
    /// `myTrimmed->GetSurface()` (`cxx:371-372`) of the same bundle.
    surface: Option<Arc<dyn Surface>>,
}

impl Darboux {
    /// The `CurveOnSurface` half of the constructor: OCCT's default
    /// `GeomFill_Darboux()` takes no argument and reads the bundle back out of
    /// `myTrimmed` after `SetCurve`; the port keeps the bundle apart, so the
    /// pcurve and the surface are supplied here and the 3D curve through
    /// [`TrihedronLaw::set_curve`].
    pub fn new(pcurve: Arc<dyn Curve2d>, surface: Arc<dyn Surface>) -> Self {
        Self {
            data: TrihedronLawCurve::new(),
            pcurve: Some(pcurve),
            surface: Some(surface),
        }
    }

    fn parts(&self) -> Result<(&Arc<dyn Curve2d>, &Arc<dyn Surface>), String> {
        match (&self.pcurve, &self.surface) {
            (Some(c), Some(s)) => Ok((c, s)),
            _ => Err("Standard_NullObject: GeomFill_Darboux needs a CurveOnSurface".into()),
        }
    }

    /// The `Copy` construction: same bundle, curve installed by `SetCurve`.
    fn from_parts(
        pcurve: Option<Arc<dyn Curve2d>>,
        surface: Option<Arc<dyn Surface>>,
    ) -> Self {
        Self {
            data: TrihedronLawCurve::new(),
            pcurve,
            surface,
        }
    }
}

impl TrihedronLaw for Darboux {
    /// `GeomFill_TrihedronLaw::SetCurve` (`GeomFill_TrihedronLaw.cxx:25-30`):
    /// `myCurve = myTrimmed = C`.
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve)
    }

    /// `GeomFill_Darboux::Copy` (`cxx:354-361`).
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        let opy = Darboux::from_parts(self.pcurve.clone(), self.surface.clone());
        if let Some(curve) = self.data.curve() {
            TrihedronLaw::set_curve(&opy, curve);
        }
        Arc::new(opy)
    }

    /// `GeomFill_Darboux::D0` (`cxx:365-392`).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        let (pcurve, surface) = self.parts()?;
        let (c2d, d2d) = pcurve.d1(param);

        let normal_dir = normal_d0(c2d.x(), c2d.y(), surface.as_ref())?;
        let bi_normal = GpVec::from_xyz(normal_dir.xyz());

        let (_s, ds_du, ds_dv) = surface.d1(c2d.x(), c2d.y());
        let tangent = ds_du
            .multiplied_scalar(d2d.x())
            .added(&ds_dv.multiplied_scalar(d2d.y()))
            .normalized();

        let normal = bi_normal.crossed(&tangent);
        Ok((tangent, normal, bi_normal))
    }

    /// `GeomFill_Darboux::D1` (`cxx:394-433`).
    fn d1(&self, param: f64) -> Result<TrihedronD1, String> {
        let (pcurve, surface) = self.parts()?;
        let (c2d, d2d, d2_2d) = pcurve.d2(param);
        let (_s, ds_du, ds_dv, d2s_du, d2s_dv, d2s_duv) = surface.d2(c2d.x(), c2d.y());

        let f = ds_du
            .multiplied_scalar(d2d.x())
            .added(&ds_dv.multiplied_scalar(d2d.y()));
        let tangent = f.normalized();
        let df = ds_du
            .multiplied_scalar(d2_2d.x())
            .added(&ds_dv.multiplied_scalar(d2_2d.y()))
            .added(&d2s_du.multiplied_scalar(d2d.x() * d2d.x()))
            .added(&d2s_duv.multiplied_scalar(2.0 * d2d.x() * d2d.y()))
            .added(&d2s_dv.multiplied_scalar(d2d.y() * d2d.y()));
        let d_tangent = fderiv(&f, &df);

        let (normal_dir, d1u_normal, d1v_normal) = normal_d1(c2d.x(), c2d.y(), surface.as_ref())?;
        let bi_normal = GpVec::from_xyz(normal_dir.xyz());
        let d_bi_normal = d1u_normal
            .multiplied_scalar(d2d.x())
            .added(&d1v_normal.multiplied_scalar(d2d.y()));

        let normal = bi_normal.crossed(&tangent);
        let d_normal = bi_normal
            .crossed(&d_tangent)
            .added(&d_bi_normal.crossed(&tangent));

        Ok(TrihedronD1 {
            tangent,
            dtangent: d_tangent,
            normal,
            dnormal: d_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
        })
    }

    /// `GeomFill_Darboux::D2` (`cxx:435-511`).
    fn d2(&self, param: f64) -> Result<TrihedronD2, String> {
        let (pcurve, surface) = self.parts()?;
        let (c2d, d2d, d2_2d, d3_2d) = pcurve.d3(param);
        let (_s, ds_du, ds_dv, d2s_du, d2s_dv, d2s_duv) = surface.d2(c2d.x(), c2d.y());
        let (u, v) = (c2d.x(), c2d.y());
        let d3s_du = surface.try_eval_dn(u, v, 3, 0)?;
        let d3s_dv = surface.try_eval_dn(u, v, 0, 3)?;
        let d3s_duuv = surface.try_eval_dn(u, v, 2, 1)?;
        let d3s_duvv = surface.try_eval_dn(u, v, 1, 2)?;

        let (dx, dy) = (d2d.x(), d2d.y());
        let (d2x, d2y) = (d2_2d.x(), d2_2d.y());
        let (d3x, d3y) = (d3_2d.x(), d3_2d.y());

        let f = ds_du
            .multiplied_scalar(dx)
            .added(&ds_dv.multiplied_scalar(dy));
        let tangent = f.normalized();
        let df = ds_du
            .multiplied_scalar(d2x)
            .added(&ds_dv.multiplied_scalar(d2y))
            .added(&d2s_du.multiplied_scalar(dx * dx))
            .added(&d2s_duv.multiplied_scalar(2.0 * dx * dy))
            .added(&d2s_dv.multiplied_scalar(dy * dy));

        let d2f = ds_du
            .multiplied_scalar(d3x)
            .added(&ds_dv.multiplied_scalar(d3y))
            .added(
                &d2s_du
                    .multiplied_scalar(dx)
                    .added(&d2s_duv.multiplied_scalar(dy))
                    .multiplied_scalar(d2x),
            )
            .added(
                &d2s_duv
                    .multiplied_scalar(dx)
                    .added(&d2s_dv.multiplied_scalar(dy))
                    .multiplied_scalar(d2y),
            )
            .added(
                &d3s_du
                    .multiplied_scalar(dx)
                    .added(&d3s_duuv.multiplied_scalar(dy))
                    .multiplied_scalar(dx * dx),
            )
            .added(
                &d3s_duvv
                    .multiplied_scalar(dx)
                    .added(&d3s_dv.multiplied_scalar(dy))
                    .multiplied_scalar(dy * dy),
            )
            .added(
                &d2s_du
                    .multiplied_scalar(dx * d2x)
                    .added(&d2s_dv.multiplied_scalar(dy * d2y))
                    .added(
                        &d3s_duuv
                            .multiplied_scalar(dx * dy * dx)
                            .added(&d3s_duvv.multiplied_scalar(dx * dy * dy)),
                    )
                    .added(&d2s_duv.multiplied_scalar(d2x * dy + dx * d2y))
                    .multiplied_scalar(2.0),
            );

        let d_tangent = fderiv(&f, &df);
        let d2_tangent = dderiv(&f, &df, &d2f);

        let (normal_dir, d1u_normal, d1v_normal, d2u_normal, d2v_normal, d2uv_normal) =
            normal_d2(u, v, surface.as_ref())?;
        let bi_normal = GpVec::from_xyz(normal_dir.xyz());
        let d_bi_normal = d1u_normal
            .multiplied_scalar(dx)
            .added(&d1v_normal.multiplied_scalar(dy));
        let d2_bi_normal = d1u_normal
            .multiplied_scalar(d2x)
            .added(&d1v_normal.multiplied_scalar(d2y))
            .added(&d2u_normal.multiplied_scalar(dx * dx))
            .added(&d2uv_normal.multiplied_scalar(2.0 * dx * dy))
            .added(&d2v_normal.multiplied_scalar(dy * dy));

        let normal = bi_normal.crossed(&tangent);
        let d_normal = bi_normal
            .crossed(&d_tangent)
            .added(&d_bi_normal.crossed(&tangent));
        let d2_normal = bi_normal
            .crossed(&d2_tangent)
            .added(&d_bi_normal.crossed(&d_tangent).multiplied_scalar(2.0))
            .added(&d2_bi_normal.crossed(&tangent));

        Ok(TrihedronD2 {
            tangent,
            dtangent: d_tangent,
            d2tangent: d2_tangent,
            normal,
            dnormal: d_normal,
            d2normal: d2_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
            d2binormal: d2_bi_normal,
        })
    }

    /// `GeomFill_Darboux::NbIntervals` (`cxx:513-516`): delegates to
    /// `myCurve` with the same shape.
    fn nb_intervals(&self, s: u8) -> i32 {
        match self.data.curve() {
            Some(c) => c.nb_intervals(s),
            None => 1,
        }
    }

    /// `GeomFill_Darboux::Intervals` (`cxx:518-521`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        match self.data.curve() {
            Some(c) => c.parameter_intervals(s),
            None => vec![0.0, 0.0],
        }
    }

    /// `GeomFill_TrihedronLaw::SetInterval` (`GeomFill_TrihedronLaw.cxx:57-60`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)
    }

    /// `GeomFill_TrihedronLaw::GetInterval` (`GeomFill_TrihedronLaw.cxx:62-67`).
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_Darboux::GetAverageLaw` (`cxx:523-551`).
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec) {
        const NUM: usize = 20;
        let mut a_tangent = GpVec::zero();
        let mut a_normal = GpVec::zero();
        let mut a_bi_normal = GpVec::zero();
        let (first, last) = self.data.interval();
        let step = (last - first) / NUM as f64;
        for i in 0..=NUM {
            let mut param = first + i as f64 * step;
            if param > last {
                param = last;
            }
            if let Ok((t, n, bn)) = self.d0(param) {
                a_tangent = a_tangent.added(&t);
                a_normal = a_normal.added(&n);
                a_bi_normal = a_bi_normal.added(&bn);
            }
        }
        a_tangent = a_tangent.divided((NUM + 1) as f64);
        a_normal = a_normal.divided((NUM + 1) as f64);
        a_tangent = a_tangent.normalized();
        a_bi_normal = a_tangent.crossed(&a_normal).normalized();
        a_normal = a_bi_normal.crossed(&a_tangent);
        (a_tangent, a_normal, a_bi_normal)
    }

    /// `GeomFill_Darboux::IsConstant` (`cxx:553-556`).
    fn is_constant(&self) -> bool {
        match self.data.curve() {
            Some(c) => matches!(extrema_curve_type(c.as_ref()), GeomAbsCurveType::Line),
            None => false,
        }
    }

    /// `GeomFill_Darboux::IsOnlyBy3dCurve` (`cxx:558-561`): false.
    fn is_only_by_3d_curve(&self) -> bool {
        false
    }
}
