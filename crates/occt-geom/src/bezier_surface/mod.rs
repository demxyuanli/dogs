//! 3D Bezier surface.

use crate::bezier_curve::GeomBezierCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::surface::Surface;
use occt_core::gp::{GpPnt, GpVec};

/// Tensor-product (non-rational) Bezier surface in 3D.
///
/// Poles are stored row-major: `poles[i * n_v + j]` is the control point for
/// the `i`-th u-index and `j`-th v-index.
#[derive(Clone)]
pub struct GeomBezierSurface {
    pub poles: Vec<GpPnt>,
    pub n_u: usize,
    pub n_v: usize,
    /// `Geom_BezierSurface::myWeights`. `None` is the polynomial surface.
    pub weights: Option<Vec<f64>>,
    /// `myURational`: a row of weights changes along V (`cxx:459-468`).
    u_rational: bool,
    /// `myVRational`: a column of weights changes along U (`cxx:447-456`).
    v_rational: bool,
}

fn lerp(a: &GpPnt, b: &GpPnt, t: f64) -> GpPnt {
    GpPnt::new(
        a.x() + (b.x() - a.x()) * t,
        a.y() + (b.y() - a.y()) * t,
        a.z() + (b.z() - a.z()) * t,
    )
}

fn lerp_vec(a: &GpVec, b: &GpVec, t: f64) -> GpVec {
    GpVec::new(
        a.x() + (b.x() - a.x()) * t,
        a.y() + (b.y() - a.y()) * t,
        a.z() + (b.z() - a.z()) * t,
    )
}

fn de_casteljau(poles: &[GpPnt], u: f64) -> GpPnt {
    let mut pts: Vec<GpPnt> = poles.to_vec();
    while pts.len() > 1 {
        for i in 0..(pts.len() - 1) {
            pts[i] = lerp(&pts[i], &pts[i + 1], u);
        }
        pts.pop();
    }
    pts[0]
}

fn de_casteljau_vec(vs: &[GpVec], u: f64) -> GpVec {
    let mut pts: Vec<GpVec> = vs.to_vec();
    while pts.len() > 1 {
        for i in 0..(pts.len() - 1) {
            pts[i] = lerp_vec(&pts[i], &pts[i + 1], u);
        }
        pts.pop();
    }
    pts[0]
}

impl GeomBezierSurface {
    pub fn new(poles: Vec<GpPnt>, n_u: usize, n_v: usize) -> Result<Self, &'static str> {
        if poles.len() != n_u * n_v {
            return Err("GeomBezierSurface: poles.len() must equal n_u * n_v");
        }
        if poles.len() < 4 {
            return Err("GeomBezierSurface requires at least 4 poles");
        }
        Ok(Self {
            poles,
            n_u,
            n_v,
            weights: None,
            u_rational: false,
            v_rational: false,
        })
    }

    /// Weighted constructor (`Geom_BezierSurface.cxx:416-478`).
    /// Weights that do not vary store a polynomial surface.
    pub fn with_weights(
        poles: Vec<GpPnt>,
        n_u: usize,
        n_v: usize,
        weights: Vec<f64>,
    ) -> Result<Self, &'static str> {
        if poles.len() != n_u * n_v || weights.len() != poles.len() {
            return Err("Geom_BezierSurface");
        }
        let max_poles = (Self::max_degree() as usize) + 1;
        if n_u < 2 || n_u > max_poles || n_v < 2 || n_v > max_poles {
            return Err("Geom_BezierSurface");
        }
        for &weight in &weights {
            if weight <= occt_core::precision::REAL_SMALL {
                return Err("Geom_BezierSurface");
            }
        }
        let (u_rational, v_rational) = rational_flags(&weights, n_u, n_v);
        if u_rational || v_rational {
            Ok(Self {
                poles,
                n_u,
                n_v,
                weights: Some(weights),
                u_rational,
                v_rational,
            })
        } else {
            Self::new(poles, n_u, n_v)
        }
    }

    pub fn nb_u_poles(&self) -> usize {
        self.n_u
    }

    pub fn nb_v_poles(&self) -> usize {
        self.n_v
    }

    pub fn u_degree(&self) -> usize {
        self.n_u - 1
    }

    pub fn v_degree(&self) -> usize {
        self.n_v - 1
    }

    /// `Geom_BezierSurface::Segment` (`Geom_BezierSurface.cxx:941-1042`).
    ///
    /// `weights == None` is the polynomial arm (`:979-996`, `:1036-1040`).
    /// `weights == Some` is the rational arm (`:960-977`, `:1013-1035`):
    /// `BuildCache` keeps the weight component, both tables transpose when
    /// `UDegree() <= VDegree()`, then `UTrimming` / `VTrimming` /
    /// `CoefficientsPoles` update poles and weights together.
    pub fn segment(&mut self, u1: f64, u2: f64, v1: f64, v2: f64) {
        let u_deg = self.u_degree() as i32;
        let v_deg = self.v_degree() as i32;
        let mut rows = vec![vec![GpPnt::new(0.0, 0.0, 0.0); self.n_v]; self.n_u];
        let rational = (self.u_rational || self.v_rational)
            && self
                .weights
                .as_ref()
                .map(|w| w.len() == self.n_u * self.n_v)
                .unwrap_or(false);
        let mut wrows = rational.then(|| vec![vec![0.0; self.n_v]; self.n_u]);
        for i in 0..self.n_u {
            for j in 0..self.n_v {
                rows[i][j] = self.poles[i * self.n_v + j];
                if let (Some(grid), Some(src)) = (wrows.as_mut(), self.weights.as_ref()) {
                    grid[i][j] = src[i * self.n_v + j];
                }
            }
        }
        let u_flat = bezier_flat_knots(u_deg);
        let v_flat = bezier_flat_knots(v_deg);
        let mut wcache = Vec::new();
        let mut coefs = occt_core::bspl::build_cache::build_cache(
            0.5,
            0.5,
            0.5,
            0.5,
            false,
            false,
            u_deg,
            v_deg,
            0,
            0,
            &u_flat,
            &v_flat,
            &rows,
            wrows.as_deref(),
            rational.then_some(&mut wcache),
        );
        let mut wcoefs = rational.then_some(wcache);
        if u_deg <= v_deg {
            let rows_n = (u_deg as usize) + 1;
            let cols_n = (v_deg as usize) + 1;
            let mut transposed = vec![vec![GpPnt::new(0.0, 0.0, 0.0); cols_n]; rows_n];
            for ii in 0..rows_n {
                for jj in 0..cols_n {
                    transposed[ii][jj] = coefs[jj][ii];
                }
            }
            coefs = transposed;
            if let Some(src) = wcoefs.as_ref() {
                let mut transposed_w = vec![vec![0.0; cols_n]; rows_n];
                for ii in 0..rows_n {
                    for jj in 0..cols_n {
                        transposed_w[ii][jj] = src[jj][ii];
                    }
                }
                wcoefs = Some(transposed_w);
            }
        }
        let u_first = 2.0 * (u1 - 0.5);
        let u_last = 2.0 * (u2 - 0.5);
        let v_first = 2.0 * (v1 - 0.5);
        let v_last = 2.0 * (v2 - 0.5);
        if let Some(wcoefs) = wcoefs.as_mut() {
            occt_core::bspl::plib::u_trimming_w(u_first, u_last, &mut coefs, wcoefs);
            occt_core::bspl::plib::v_trimming_w(v_first, v_last, &mut coefs, wcoefs);
        } else {
            occt_core::bspl::plib::u_trimming(u_first, u_last, &mut coefs);
            occt_core::bspl::plib::v_trimming(v_first, v_last, &mut coefs);
        }
        if let Some(wcoefs) = wcoefs.as_ref() {
            let (poles, weights) =
                occt_core::bspl::plib::coefficients_poles_grid_w(&coefs, wcoefs);
            if poles.is_empty() {
                return;
            }
            self.n_u = poles.len();
            self.n_v = poles[0].len();
            self.poles = Vec::with_capacity(self.n_u * self.n_v);
            self.weights = Some(Vec::with_capacity(self.n_u * self.n_v));
            for (row, wrow) in poles.into_iter().zip(weights) {
                self.poles.extend(row);
                if let Some(dst) = self.weights.as_mut() {
                    dst.extend(wrow);
                }
            }
            return;
        }
        let poles = occt_core::bspl::plib::coefficients_poles_grid(&coefs);
        if poles.is_empty() {
            return;
        }
        self.n_u = poles.len();
        self.n_v = poles[0].len();
        self.poles = Vec::with_capacity(self.n_u * self.n_v);
        for row in poles {
            self.poles.extend(row);
        }
    }

    /// Copy `myWeights` and the rational flags without rescanning.
    /// `Segment` does not recompute the flags (`Geom_BezierSurface.cxx:1042`).
    pub fn adopt_weights(&mut self, weights: Vec<f64>, u_rational: bool, v_rational: bool) {
        if weights.len() == self.n_u * self.n_v && (u_rational || v_rational) {
            self.weights = Some(weights);
            self.u_rational = u_rational;
            self.v_rational = v_rational;
        }
    }

    pub fn is_u_rational(&self) -> bool {
        self.u_rational
    }

    pub fn is_v_rational(&self) -> bool {
        self.v_rational
    }

    /// `Geom_BezierSurface::MaxDegree` (`Geom_BezierSurface.cxx:483-486`).
    pub fn max_degree() -> i32 {
        25
    }

    /// `Geom_BezierSurface::Weight` (`Geom_BezierSurface.cxx:1866-1877`).
    /// Indices are 1-based. A polynomial surface reports weight 1.
    pub fn weight(&self, u_index: i32, v_index: i32) -> Result<f64, &'static str> {
        if u_index < 1
            || u_index > self.n_u as i32
            || v_index < 1
            || v_index > self.n_v as i32
        {
            return Err("Geom_BezierSurface::Weight");
        }
        let slot = ((u_index as usize) - 1) * self.n_v + (v_index as usize) - 1;
        Ok(self
            .weights
            .as_ref()
            .map(|w| w[slot])
            .unwrap_or(1.0))
    }

    /// `Geom_BezierSurface::SetWeight` (`Geom_BezierSurface.cxx:1179-1215`).
    /// `UIndex` and `VIndex` are 1-based. `myMaxDerivInvOk` (`:1214`) is
    /// not stored. A polynomial surface keeps an implicit unit-weight grid.
    pub fn set_weight(&mut self, u_index: i32, v_index: i32, weight: f64) -> Result<(), &'static str> {
        if weight <= occt_core::precision::REAL_SMALL {
            return Err("Geom_BezierSurface::SetWeight");
        }
        if u_index < 1
            || u_index > self.n_u as i32
            || v_index < 1
            || v_index > self.n_v as i32
        {
            return Err("Geom_BezierSurface::SetWeight");
        }
        let was_rat = self.u_rational || self.v_rational;
        if !was_rat && (weight - 1.0).abs() <= occt_core::precision::REAL_SMALL {
            return Ok(());
        }
        if self.weights.is_none() {
            self.weights = Some(vec![1.0; self.n_u * self.n_v]);
        }
        let slot = ((u_index as usize) - 1) * self.n_v + (v_index as usize) - 1;
        let current = self.weights.as_ref().map(|w| w[slot]).unwrap_or(1.0);
        if (weight - current).abs() > occt_core::precision::REAL_SMALL {
            if let Some(weights) = self.weights.as_mut() {
                weights[slot] = weight;
            }
            let (u_rat, v_rat) = rational_flags(self.weights.as_deref().unwrap_or(&[]), self.n_u, self.n_v);
            self.u_rational = u_rat;
            self.v_rational = v_rat;
            if !(u_rat || v_rat) {
                self.weights = None;
            }
        }
        Ok(())
    }

    /// `Geom_BezierSurface::SetWeightCol` (`Geom_BezierSurface.cxx:1219-1258`).
    /// `v_index` is 1-based. `lower` is the OCCT array lower, and each array
    /// index is the U index written into that column. `myMaxDerivInvOk`
    /// (`:1257`) is not stored.
    pub fn set_weight_col(
        &mut self,
        v_index: i32,
        lower: i32,
        pole_weights: &[f64],
    ) -> Result<(), &'static str> {
        if v_index < 1 || v_index > self.n_v as i32 {
            return Err("Geom_BezierSurface::SetWeightCol");
        }
        if pole_weights.len() != self.n_u {
            return Err("Geom_BezierSurface::SetWeightCol");
        }
        self.unit_weights_if_polynomial();
        let upper = lower + pole_weights.len() as i32 - 1;
        let mut i = lower;
        let mut src = 0usize;
        while i <= upper {
            let weight = pole_weights[src];
            if weight <= occt_core::precision::REAL_SMALL {
                return Err("Geom_BezierSurface::SetWeightCol");
            }
            if i < 1 || i > self.n_u as i32 {
                return Err("Geom_BezierSurface::SetWeightCol");
            }
            let slot = ((i as usize) - 1) * self.n_v + (v_index as usize) - 1;
            if let Some(weights) = self.weights.as_mut() {
                weights[slot] = weight;
            }
            i += 1;
            src += 1;
        }
        self.store_rational_flags();
        Ok(())
    }

    /// `Geom_BezierSurface::SetWeightRow` (`Geom_BezierSurface.cxx:1262-1301`).
    /// `u_index` is 1-based. `lower`..`upper` must lie in `1..=n_v`.
    /// `myMaxDerivInvOk` (`:1300`) is not stored.
    pub fn set_weight_row(
        &mut self,
        u_index: i32,
        lower: i32,
        pole_weights: &[f64],
    ) -> Result<(), &'static str> {
        if u_index < 1 || u_index > self.n_u as i32 {
            return Err("Geom_BezierSurface::SetWeightRow");
        }
        let upper = if pole_weights.is_empty() {
            lower - 1
        } else {
            lower + pole_weights.len() as i32 - 1
        };
        if lower < 1 || lower > self.n_v as i32 || upper < 1 || upper > self.n_v as i32 {
            return Err("Geom_BezierSurface::SetWeightRow");
        }
        self.unit_weights_if_polynomial();
        let mut i = lower;
        let mut src = 0usize;
        while i <= upper {
            let weight = pole_weights[src];
            if weight <= occt_core::precision::REAL_SMALL {
                return Err("Geom_BezierSurface::SetWeightRow");
            }
            let slot = ((u_index as usize) - 1) * self.n_v + (i as usize) - 1;
            if let Some(weights) = self.weights.as_mut() {
                weights[slot] = weight;
            }
            i += 1;
            src += 1;
        }
        self.store_rational_flags();
        Ok(())
    }

    /// `Geom_BezierSurface::UReverse` (`Geom_BezierSurface.cxx:1305-1341`).
    /// `ClearEvalRepresentation` (`:1339`) is skipped: no eval rep is stored.
    pub fn u_reverse(&mut self) {
        let nu = self.n_u;
        let nv = self.n_v;
        let rational = self.u_rational || self.v_rational;
        for col in 0..nv {
            for row in 0..(nu / 2) {
                let a = row * nv + col;
                let b = (nu - row - 1) * nv + col;
                self.poles.swap(a, b);
                if rational {
                    if let Some(weights) = self.weights.as_mut() {
                        weights.swap(a, b);
                    }
                }
            }
        }
    }

    /// `Geom_BezierSurface::VReverse` (`Geom_BezierSurface.cxx:1352-1388`).
    /// `ClearEvalRepresentation` (`:1386`) is skipped: no eval rep is stored.
    pub fn v_reverse(&mut self) {
        let nu = self.n_u;
        let nv = self.n_v;
        let rational = self.u_rational || self.v_rational;
        for row in 0..nu {
            for col in 0..(nv / 2) {
                let a = row * nv + col;
                let b = row * nv + (nv - col - 1);
                self.poles.swap(a, b);
                if rational {
                    if let Some(weights) = self.weights.as_mut() {
                        weights.swap(a, b);
                    }
                }
            }
        }
    }

    fn unit_weights_if_polynomial(&mut self) {
        if !(self.u_rational || self.v_rational) || self.weights.is_none() {
            self.weights = Some(vec![1.0; self.n_u * self.n_v]);
        }
    }

    fn store_rational_flags(&mut self) {
        let (u_rat, v_rat) = self
            .weights
            .as_deref()
            .map(|weights| rational_flags(weights, self.n_u, self.n_v))
            .unwrap_or((false, false));
        self.u_rational = u_rat;
        self.v_rational = v_rat;
        if !(u_rat || v_rat) {
            self.weights = None;
        }
    }

    /// `GeomConvert` Bezier arm after `Segment` (`GeomConvert_1.cxx:714-729`).
    pub fn converted_bspline(&self) -> Option<GeomBSplineSurface> {
        if self.u_rational || self.v_rational {
            return self.rational_bspline();
        }
        self.osculating_bspline()
    }

    /// `Geom_BezierSurface::UIso` (`Geom_BezierSurface.cxx:1769-1810`).
    pub fn u_iso(&self, u: f64) -> Result<GeomBezierCurve, &'static str> {
        if self.u_rational || self.v_rational {
            if let Some(weights) = self.weights.as_ref() {
                let degree = self.n_u - 1;
                let knots = bezier_flat_knots(degree as i32);
                let mut poles = Vec::with_capacity(self.n_v);
                let mut curve_w = Vec::with_capacity(self.n_v);
                for j in 0..self.n_v {
                    let col: Vec<GpPnt> = (0..self.n_u)
                        .map(|i| self.poles[i * self.n_v + j])
                        .collect();
                    let wcol: Vec<f64> = (0..self.n_u)
                        .map(|i| weights[i * self.n_v + j])
                        .collect();
                    let (point, weight) = iso_sample(&col, &wcol, &knots, degree, u);
                    poles.push(point);
                    curve_w.push(weight);
                }
                if self.u_rational {
                    return GeomBezierCurve::with_weights(poles, curve_w);
                }
                return GeomBezierCurve::new(poles);
            }
        }
        let mut poles = Vec::with_capacity(self.n_v);
        for j in 0..self.n_v {
            let col: Vec<GpPnt> = (0..self.n_u)
                .map(|i| self.poles[i * self.n_v + j])
                .collect();
            poles.push(de_casteljau(&col, u));
        }
        GeomBezierCurve::new(poles)
    }

    /// `Geom_BezierSurface::VIso` (`Geom_BezierSurface.cxx:1821-1862`).
    pub fn v_iso(&self, v: f64) -> Result<GeomBezierCurve, &'static str> {
        if self.v_rational || self.u_rational {
            if let Some(weights) = self.weights.as_ref() {
                let degree = self.n_v - 1;
                let knots = bezier_flat_knots(degree as i32);
                let mut poles = Vec::with_capacity(self.n_u);
                let mut curve_w = Vec::with_capacity(self.n_u);
                for i in 0..self.n_u {
                    let start = i * self.n_v;
                    let row = self.poles[start..start + self.n_v].to_vec();
                    let wrow = weights[start..start + self.n_v].to_vec();
                    let (point, weight) = iso_sample(&row, &wrow, &knots, degree, v);
                    poles.push(point);
                    curve_w.push(weight);
                }
                if self.v_rational {
                    return GeomBezierCurve::with_weights(poles, curve_w);
                }
                return GeomBezierCurve::new(poles);
            }
        }
        let mut poles = Vec::with_capacity(self.n_u);
        for i in 0..self.n_u {
            let start = i * self.n_v;
            poles.push(de_casteljau(&self.poles[start..start + self.n_v], v));
        }
        GeomBezierCurve::new(poles)
    }

    fn du(&self, u: f64, v: f64) -> GpVec {
        if self.n_u < 2 {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        let nu = self.n_u as f64;
        // Differentiate the u-direction Bernstein basis: net is (n_u-1) x n_v.
        let mut net: Vec<GpVec> = Vec::with_capacity((self.n_u - 1) * self.n_v);
        for i in 0..(self.n_u - 1) {
            for j in 0..self.n_v {
                let p0 = &self.poles[i * self.n_v + j];
                let p1 = &self.poles[(i + 1) * self.n_v + j];
                net.push(GpVec::new(
                    (p1.x() - p0.x()) * nu,
                    (p1.y() - p0.y()) * nu,
                    (p1.z() - p0.z()) * nu,
                ));
            }
        }
        // Evaluate each row (along v) then along u.
        let mut pts = Vec::with_capacity(self.n_u - 1);
        for i in 0..(self.n_u - 1) {
            let start = i * self.n_v;
            pts.push(de_casteljau_vec(&net[start..start + self.n_v], v));
        }
        de_casteljau_vec(&pts, u)
    }

    fn dv(&self, u: f64, v: f64) -> GpVec {
        if self.n_v < 2 {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        let nv = self.n_v as f64;
        // Differentiate the v-direction Bernstein basis: net is n_u x (n_v-1).
        let mut net: Vec<GpVec> = Vec::with_capacity(self.n_u * (self.n_v - 1));
        for i in 0..self.n_u {
            for j in 0..(self.n_v - 1) {
                let p0 = &self.poles[i * self.n_v + j];
                let p1 = &self.poles[i * self.n_v + j + 1];
                net.push(GpVec::new(
                    (p1.x() - p0.x()) * nv,
                    (p1.y() - p0.y()) * nv,
                    (p1.z() - p0.z()) * nv,
                ));
            }
        }
        // Evaluate each row (along v) then along u.
        let mut pts = Vec::with_capacity(self.n_u);
        for i in 0..self.n_u {
            let start = i * (self.n_v - 1);
            pts.push(de_casteljau_vec(&net[start..start + self.n_v - 1], v));
        }
        de_casteljau_vec(&pts, u)
    }

    /// Rational Bezier as the equivalent B-spline (`UKnots`/`VKnots` are
    /// `0, 1` with multiplicity `degree+1`). Used by `EvalD0` (`:1425`),
    /// `EvalD1` (`:1481`), `EvalD2` (`:1540`), and `EvalDN` (`:1681`).
    /// `TryEvalSurfaceD*` is not taken: this type stores no eval rep.
    fn rational_bspline(&self) -> Option<GeomBSplineSurface> {
        if !(self.u_rational || self.v_rational) {
            return None;
        }
        let src = self.weights.as_ref()?;
        if src.len() != self.n_u * self.n_v || self.n_u < 2 || self.n_v < 2 {
            return None;
        }
        let mut poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); self.n_v]; self.n_u];
        let mut weights = vec![vec![0.0; self.n_v]; self.n_u];
        for i in 0..self.n_u {
            for j in 0..self.n_v {
                poles[i][j] = self.poles[i * self.n_v + j];
                weights[i][j] = src[i * self.n_v + j];
            }
        }
        let deg_u = self.n_u - 1;
        let deg_v = self.n_v - 1;
        GeomBSplineSurface::rational(
            poles,
            weights,
            bezier_flat_knots(deg_u as i32),
            bezier_flat_knots(deg_v as i32),
            deg_u,
            deg_v,
        )
        .ok()
    }
}

/// `Geom_BezierSurface::Rational` (`Geom_BezierSurface.cxx:66-94`).
/// `myVRational` is a change along U. `myURational` is a change along V.
/// The tolerance is `Epsilon(abs(weight))`.
fn rational_flags(weights: &[f64], n_u: usize, n_v: usize) -> (bool, bool) {
    if weights.len() != n_u * n_v {
        return (false, false);
    }
    let mut v_rational = false;
    'cols: for j in 0..n_v {
        for i in 0..n_u.saturating_sub(1) {
            let a = weights[i * n_v + j];
            let b = weights[(i + 1) * n_v + j];
            if (a - b).abs() > f64::EPSILON * a.abs() {
                v_rational = true;
                break 'cols;
            }
        }
    }
    let mut u_rational = false;
    'rows: for i in 0..n_u {
        for j in 0..n_v.saturating_sub(1) {
            let a = weights[i * n_v + j];
            let b = weights[i * n_v + j + 1];
            if (a - b).abs() > f64::EPSILON * a.abs() {
                u_rational = true;
                break 'rows;
            }
        }
    }
    (u_rational, v_rational)
}

/// One output pole of `BSplSLib::Iso` (`BSplSLib.cxx:1704-1723`):
/// homogeneous de Boor on the Bezier knot sequence, then divide by the
/// blended weight. The weight itself is the blended fourth component.
fn iso_sample(poles: &[GpPnt], weights: &[f64], knots: &[f64], degree: usize, param: f64) -> (GpPnt, f64) {
    let wpts: Vec<GpPnt> = weights.iter().copied().map(|w| GpPnt::new(w, 0.0, 0.0)).collect();
    let point = occt_core::bspl::eval::eval_curve_rational(poles, weights, knots, degree, param);
    let weight = occt_core::bspl::eval::eval_curve(&wpts, knots, degree, param).x();
    (point, weight)
}

mod degree;
mod poles;
mod surface_impl;

/// `Geom_BezierSurface::UKnotSequence`: degree+1 copies of 0 then of 1.
fn bezier_flat_knots(degree: i32) -> Vec<f64> {
    let mult = degree + 1;
    occt_core::bspl::banded_interp::knot_sequence(&[0.0, 1.0], &[mult, mult], degree)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: &GpPnt, b: &GpPnt) -> bool {
        (a.x() - b.x()).abs() < 1e-12
            && (a.y() - b.y()).abs() < 1e-12
            && (a.z() - b.z()).abs() < 1e-12
    }

    #[test]
    fn bilinear_surface() {
        let s = GeomBezierSurface::new(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
                GpPnt::new(1.0, 1.0, 0.0),
            ],
            2,
            2,
        )
        .unwrap();

        assert!(approx(&s.d0(0.5, 0.5), &GpPnt::new(0.5, 0.5, 0.0)));
        assert!(approx(&s.d0(0.0, 0.0), &GpPnt::new(0.0, 0.0, 0.0)));
        assert!(approx(&s.d0(1.0, 1.0), &GpPnt::new(1.0, 1.0, 0.0)));
    }

    #[test]
    fn rejects_mismatched_pole_count() {
        assert!(GeomBezierSurface::new(vec![GpPnt::new(0.0, 0.0, 0.0)], 2, 2).is_err());
    }
}
