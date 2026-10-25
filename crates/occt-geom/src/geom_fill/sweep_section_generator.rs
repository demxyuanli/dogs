//! `GeomFill_SweepSectionGenerator` — the section sequence a sweep is
//! approximated from: a path plus either a radius, a first section, or a first
//! and a last section.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_SweepSectionGenerator.cxx:41-700` and
//! `GeomFill_SweepSectionGenerator.hxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `myType` keeps OCCT's numbering: `1` (radius, non-circular path), `4`
//!   (radius, circular path), `2` / `5` (first section), `3` / `6` (first and
//!   last sections), where `+3` means the path is a `GeomAbs_Circle`, and `0`
//!   (the `Init(Adaptor3d_Curve&, Adaptor3d_Curve&, Adaptor3d_Curve&, double)`
//!   overload, `cxx:228-243`) builds each section from the path point and the
//!   two bounding curves.
//! - `Adaptor3d_Curve` handles (`myAdpPath`, `myAdpFirstSect`, `myAdpLastSect`)
//!   are the curves themselves: the adaptor only forwards `GetType` /
//!   `FirstParameter` / `LastParameter` / `Value` / `Line`, all of which
//!   [`crate::curve::Curve`] answers.
//! - `myPath->NbKnots()` and `myFirstSect->NbKnots()` are the **distinct** knot
//!   counts (`Geom_BSplineCurve::NbKnots() { return myKnots.Length(); }`), i.e.
//!   [`GeomBSplineCurve::distinct_knots`]`().len()`, not the flat knot vector.
//! - `Rot` and `TR` are declared outside the `Perform` loop (`cxx:294-295`) and
//!   therefore keep their previous value on an iteration that does not
//!   reassign them; the port keeps them outside the loop for the same effect.
//! - `GCPnts_QuasiUniformDeflection` is
//!   [`crate::gcpnts::QuasiUniformDeflection`].
//! - `GeomFill_Profiler` is [`crate::geom_fill::profiler::Profiler`];
//!   `Profil.Curve(1)` / `Curve(2)` are [`Profiler::curves`]`[0]` / `[1]`.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::elib::clib;
use occt_core::gp::{GpAx1, GpAx2, GpCirc, GpDir, GpPnt, GpTrsf, GpVec};
use occt_core::precision::{epsilon, CONFUSION};

use crate::bspline_curve::GeomBSplineCurve;
use crate::convert_bspl::{bspline_copy, curve_to_bspline_curve, is_bspline};
use crate::curve::Curve;
use crate::gcpnts::QuasiUniformDeflection;
use crate::geom_fill::profiler::Profiler;
use crate::geom_fill::section_convert::circle_arc_bspline_of;

/// `Precision::Angular()`.
const ANGULAR: f64 = occt_core::precision::ANGULAR;

/// `GeomFill_SweepSectionGenerator` (`GeomFill_SweepSectionGenerator.hxx`).
pub struct SweepSectionGenerator {
    /// `myPath` — the path converted to a B-spline.
    path: Option<GeomBSplineCurve>,
    /// `myCircPathAxis`.
    circ_path_axis: Option<GpAx1>,
    /// `myFirstSect`.
    first_sect: Option<GeomBSplineCurve>,
    /// `myLastSect`.
    last_sect: Option<GeomBSplineCurve>,
    /// `myTrsfs`.
    trsfs: Vec<GpTrsf>,
    /// `myRadius`.
    radius: f64,
    /// `myIsDone`.
    is_done: bool,
    /// `myNbSections`.
    nb_sections: i32,
    /// `myType`.
    kind: i32,
    /// `myPolynomial`.
    polynomial: bool,
    /// `myAdpPath` (`hxx`, used by the `myType == 0` sections only).
    adp_path: Option<Arc<dyn Curve>>,
    /// `myAdpFirstSect`.
    adp_first_sect: Option<Arc<dyn Curve>>,
    /// `myAdpLastSect`.
    adp_last_sect: Option<Arc<dyn Curve>>,
}

impl Default for SweepSectionGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl SweepSectionGenerator {
    /// `GeomFill_SweepSectionGenerator()` (`cxx:41-48`).
    pub fn new() -> Self {
        Self {
            path: None,
            circ_path_axis: None,
            first_sect: None,
            last_sect: None,
            trsfs: Vec::new(),
            radius: 0.0,
            is_done: false,
            nb_sections: 0,
            kind: -1,
            polynomial: false,
            adp_path: None,
            adp_first_sect: None,
            adp_last_sect: None,
        }
    }

    /// `GeomFill_SweepSectionGenerator(Path, Radius)` (`cxx:52-57`).
    pub fn with_radius(path: &dyn Curve, radius: f64) -> Result<Self, String> {
        let mut this = Self::new();
        this.init_radius(path, radius)?;
        Ok(this)
    }

    /// `GeomFill_SweepSectionGenerator(Path, FirstSect)` (`cxx:62-66`).
    pub fn with_first_section(path: &dyn Curve, first_sect: &dyn Curve) -> Result<Self, String> {
        let mut this = Self::new();
        this.init_first_section(path, first_sect)?;
        Ok(this)
    }

    /// `GeomFill_SweepSectionGenerator(Path, FirstSect, LastSect)`
    /// (`cxx:69-76`).
    pub fn with_sections(
        path: &dyn Curve,
        first_sect: &dyn Curve,
        last_sect: &dyn Curve,
    ) -> Result<Self, String> {
        let mut this = Self::new();
        this.init_sections(path, first_sect, last_sect)?;
        Ok(this)
    }

    /// `Init(Path, Radius)` (`cxx:90-113`).
    pub fn init_radius(&mut self, path: &dyn Curve, radius: f64) -> Result<(), String> {
        self.is_done = false;
        self.radius = radius;
        if let Some(circ) = path.gp_circ() {
            self.circ_path_axis = Some(circ.position().axis().clone());
            self.kind = 4;
        } else {
            self.kind = 1;
        }
        self.path = Some(path_to_bspline(path)?);
        Ok(())
    }

    /// `Init(Path, FirstSect)` (`cxx:117-158`).
    pub fn init_first_section(&mut self, path: &dyn Curve, first_sect: &dyn Curve) -> Result<(), String> {
        self.is_done = false;
        self.radius = 0.0;
        if let Some(circ) = path.gp_circ() {
            self.circ_path_axis = Some(circ.position().axis().clone());
            self.kind = 5;
        } else {
            self.kind = 2;
        }
        self.path = Some(path_to_bspline(path)?);
        let mut sect = section_to_bspline(first_sect)?;
        if sect.periodic() {
            sect.set_not_periodic();
        }
        self.first_sect = Some(sect);
        Ok(())
    }

    /// `Init(Path, FirstSect, LastSect)` (`cxx:162-224`).
    pub fn init_sections(
        &mut self,
        path: &dyn Curve,
        first_sect: &dyn Curve,
        last_sect: &dyn Curve,
    ) -> Result<(), String> {
        self.is_done = false;
        self.radius = 0.0;
        if let Some(circ) = path.gp_circ() {
            self.circ_path_axis = Some(circ.position().axis().clone());
            self.kind = 6;
        } else {
            self.kind = 3;
        }
        self.path = Some(path_to_bspline(path)?);

        let mut first = section_to_bspline(first_sect)?;
        let mut last = section_to_bspline(last_sect)?;
        if first.periodic() {
            first.set_not_periodic();
        }
        if last.periodic() {
            last.set_not_periodic();
        }

        let mut profil = Profiler::default();
        profil.add_curve(&first)?;
        profil.add_curve(&last)?;
        profil.perform(CONFUSION)?;
        let curves = profil.curves();
        if curves.len() < 2 {
            return Err("GeomFill_SweepSectionGenerator: profiler produced less than 2 curves".into());
        }
        self.first_sect = Some(curves[0].clone());
        self.last_sect = Some(curves[1].clone());
        Ok(())
    }

    /// `Init(const Adaptor3d_Curve& Path, const Adaptor3d_Curve& Curve1,
    /// const Adaptor3d_Curve& Curve2, const double Radius)` (`cxx:228-243`):
    /// `myType = 0`, the path is converted with the default
    /// `GeomConvert::CurveToBSplineCurve` parameterisation and the two curves
    /// are kept as adaptors.
    pub fn init_adaptors(
        &mut self,
        path: Arc<dyn Curve>,
        curve1: Arc<dyn Curve>,
        curve2: Arc<dyn Curve>,
        radius: f64,
    ) -> Result<(), String> {
        self.is_done = false;
        self.radius = radius;
        self.kind = 0;
        self.path = Some(path_to_bspline(&*path)?);
        self.adp_path = Some(path);
        self.adp_first_sect = Some(curve1);
        self.adp_last_sect = Some(curve2);
        Ok(())
    }

    /// `GeomFill_SweepSectionGenerator(Path, Curve1, Curve2, Radius)`
    /// (`cxx:79-86`).
    pub fn with_adaptors(
        path: Arc<dyn Curve>,
        curve1: Arc<dyn Curve>,
        curve2: Arc<dyn Curve>,
        radius: f64,
    ) -> Result<Self, String> {
        let mut this = Self::new();
        this.init_adaptors(path, curve1, curve2, radius)?;
        Ok(this)
    }

    /// `Perform(Polynomial)` (`cxx:247-373`).
    pub fn perform(&mut self, polynomial: bool) -> Result<(), String> {
        let path = self
            .path
            .as_ref()
            .ok_or_else(|| "GeomFill_SweepSectionGenerator: path is not initialised".to_string())?
            .clone();
        self.polynomial = polynomial;

        // Eval myNbSections.
        let n_spans = (path.distinct_knots().len() as i32) - 1;
        self.nb_sections = 21 * n_spans;

        let u1 = path.first_parameter();
        let u2 = path.last_parameter();

        // Calcul de la longueur approximative de la courbe.
        let p1 = path.d0(u1);
        let p2 = path.d0((u1 + u2) / 2.0);
        let p3 = path.d0(u2);
        let length = p1.distance(&p2) + p2.distance(&p3);
        let fleche = 1.0e-5 * length;
        let mut samp = QuasiUniformDeflection::new();
        // `GCPnts_QuasiUniformDeflection(theC, theDeflection)` uses the default
        // `GeomAbs_C1` continuity (`GCPnts_QuasiUniformDeflection.hxx`), i.e. `2`.
        samp.initialize_whole(&path, fleche, 2);

        if samp.is_done() && samp.nb_points() > self.nb_sections {
            self.nb_sections = samp.nb_points();
        }

        // The transformations are calculated on different points of <myPath>
        // corresponding to the path parameter uniformly reparted.
        let delta_u = (u2 - u1) / (self.nb_sections - 1) as f64;
        let mut parameters = vec![0.0; self.nb_sections as usize];
        parameters[0] = 0.0;
        for i in 2..self.nb_sections as usize {
            parameters[i - 1] = (i as f64 - 1.0) * delta_u;
        }
        parameters[self.nb_sections as usize - 1] = u2 - u1;

        let mut tr = GpTrsf::identity();
        let mut cumul_tr = GpTrsf::identity();
        let mut trans = GpTrsf::identity();

        let (mut p_ref, mut d1_ref) = path.d1(u1);

        if self.kind == 1 || self.kind == 4 {
            // We create a circle with radius <myRadius>. This axis is created
            // with main direction <DRef> (first derivative vector of <myPath>
            // at the first point <PRef>). This circle is, after being converted
            // to a B-spline curve, put in <myFirstSect>.
            let axis = GpAx2::from_axis(
                p_ref.clone(),
                GpDir::from_xyz(&d1_ref.coord)
                    .map_err(|_| "GeomFill_SweepSectionGenerator: null path derivative".to_string())?,
            );
            let circ = GpCirc::new(axis, self.radius);
            let sect = circle_arc_bspline_of(&circ, 0.0, 2.0 * std::f64::consts::PI, ParameterisationType::QuasiAngular)
                .map_err(|e| format!("GeomFill_SweepSectionGenerator: {e:?}"))?;
            self.first_sect = Some(sect);
        }

        if (1..=3).contains(&self.kind) {
            // `Rot` is only reassigned when the derivatives are not parallel,
            // so it keeps its previous value otherwise (see the module note).
            let mut rot = GpTrsf::identity();
            for i in 2..=self.nb_sections {
                let mut u = parameters[i as usize - 1] + u1;
                if i == self.nb_sections {
                    u = u2;
                }
                let (p, d1) = path.d1(u);

                // Eval the translation between the (i-1) section and the i-th.
                trans.set_translation_pnts(&p_ref, &p);

                if !d1_ref.is_parallel_ang(&d1, ANGULAR) {
                    // Eval the rotation between the (i-1) section and the i-th.
                    let cross = d1_ref.crossed(&d1);
                    let axis = GpAx1::new(
                        p.clone(),
                        GpDir::from_xyz(&cross.coord)
                            .map_err(|_| "GeomFill_SweepSectionGenerator: null rotation axis".to_string())?,
                    );
                    let angle = d1_ref.angle_with_ref(&d1, &cross);
                    rot = GpTrsf::identity();
                    rot.set_rotation_ax1(&axis, angle)
                        .map_err(|e| format!("GeomFill_SweepSectionGenerator: {e:?}"))?;
                } else if d1_ref.is_opposite_ang(&d1, ANGULAR) {
                    // TR is the transformation between the (i-1) section and
                    // the i-th.
                    tr = rot.multiplied(&trans);
                }
                // cumulTR is the transformation between <myFirstSect> and the
                // i-th section.
                cumul_tr = tr.multiplied(&cumul_tr);
                self.trsfs.push(cumul_tr.clone());

                p_ref = p;
                d1_ref = d1;
            }
        } else if self.kind != 0 {
            let axis = self
                .circ_path_axis
                .as_ref()
                .ok_or_else(|| "GeomFill_SweepSectionGenerator: missing circular path axis".to_string())?;
            for i in 2..=self.nb_sections {
                cumul_tr = GpTrsf::identity();
                cumul_tr
                    .set_rotation_ax1(axis, parameters[i as usize - 1])
                    .map_err(|e| format!("GeomFill_SweepSectionGenerator: {e:?}"))?;
                self.trsfs.push(cumul_tr.clone());
            }
        }

        self.is_done = true;
        Ok(())
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:375-402`).
    pub fn get_shape(&self, nb_poles: &mut i32, nb_knots: &mut i32, degree: &mut i32, nb_poles2d: &mut i32) {
        if self.kind != 0 {
            if let Some(sect) = &self.first_sect {
                *nb_poles = sect.nb_poles() as i32;
                *nb_knots = sect.distinct_knots().len() as i32;
                *degree = sect.degree() as i32;
            }
        } else {
            *nb_poles = 7;
            *nb_knots = 2;
            *degree = 6;
        }
        *nb_poles2d = 0;
    }

    /// `Knots(TKnots)` (`cxx:406-426`).
    pub fn knots(&self) -> Vec<f64> {
        if self.kind != 0 {
            if let Some(sect) = &self.first_sect {
                return sect.distinct_knots().to_vec();
            }
            return Vec::new();
        }
        vec![0.0, 1.0]
    }

    /// `Mults(TMults)` (`cxx:432-447`).
    pub fn mults(&self) -> Vec<i32> {
        if self.kind != 0 {
            if let Some(sect) = &self.first_sect {
                return sect.multiplicities().to_vec();
            }
            return Vec::new();
        }
        vec![7, 7]
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)`
    /// (`cxx:453-537`). `Poles2d` / `DPoles2d` are never written by OCCT here,
    /// so they are not exposed.
    pub fn section(
        &self,
        p: i32,
        poles: &mut [GpPnt],
        d_poles: &mut [GpVec],
        weights: &mut [f64],
        d_weights: &mut [f64],
    ) -> bool {
        self.section_simple(p, poles, weights);

        // For pipes on edges the derivatives are not computed yet.
        if self.kind == 0 {
            return false;
        }

        // The derivatives on the surface are computed by approximating the path
        // in the neighbourhood of P(u) with the osculating circle of the path.
        let u = if p == 1 {
            self.path.as_ref().map(|path| path.first_parameter())
        } else if p == self.nb_sections {
            self.path.as_ref().map(|path| path.last_parameter())
        } else {
            None
        };
        let Some(u) = u else { return false };
        let Some(path) = &self.path else { return false };

        let (pt, d1, d2) = path.d2(u);
        let l = d1.magnitude();
        if l < epsilon(1.0) {
            return false;
        }

        let t = GpDir::from_xyz(&d1.coord).map(|d| GpVec::new(d.x(), d.y(), d.z()));
        let Ok(t) = t else { return false };
        let m = d2.dot(&t);
        let d = d2.subtracted(&t.multiplied_scalar(m));
        let c = d.magnitude() / (l * l);

        let nb_poles = self.first_sect.as_ref().map(|s| s.nb_poles()).unwrap_or(0);
        if c < epsilon(1.0) {
            // Null curvature: equivalent to a translation of the section.
            for i in 0..nb_poles.min(d_poles.len()) {
                d_poles[i] = d1.clone();
            }
        } else {
            let n = GpDir::from_xyz(&d.coord).map(|dir| GpVec::new(dir.x(), dir.y(), dir.z()));
            let Ok(n) = n else { return false };
            let q = pt.translated_vec(&n.multiplied_scalar(1.0 / c));
            for i in 0..nb_poles {
                if i >= poles.len() || i >= d_poles.len() {
                    break;
                }
                let v = GpVec::from_pnts(&q, &poles[i]);
                let x = v.dot(&t);
                let y = v.dot(&n);
                let mut dp = n.multiplied_scalar(x).subtracted(&t.multiplied_scalar(y));
                if dp.magnitude() > epsilon(1.0) {
                    dp = dp.normalized().multiplied_scalar((x * x + y * y).sqrt());
                }
                d_poles[i] = dp;
            }
        }

        for i in 0..nb_poles.min(d_weights.len()) {
            d_weights[i] = 0.0;
        }
        true
    }

    /// `Section(P, Poles, Poles2d, Weigths)` (`cxx:542-665`).
    pub fn section_simple(&self, p: i32, poles: &mut [GpPnt], weights: &mut [f64]) -> bool {
        if self.kind != 0 {
            let Some(first) = &self.first_sect else { return false };
            let n = poles.len().min(first.nb_poles());
            poles[..n].clone_from_slice(&first.poles[..n]);
            let w = first
                .weights
                .clone()
                .unwrap_or_else(|| vec![1.0; first.nb_poles()]);
            let n = weights.len().min(w.len());
            weights[..n].copy_from_slice(&w[..n]);

            if p > 1 {
                let idx = (p - 2) as usize;
                if idx >= self.trsfs.len() {
                    return false;
                }
                let cumul = &self.trsfs[idx];
                if self.kind == 3 || self.kind == 6 {
                    let Some(last) = &self.last_sect else { return false };
                    for i in 1..=first.nb_poles() {
                        let num = (self.nb_sections - p) as f64;
                        let den = (p - 1) as f64;
                        let xyz = first
                            .pole(i - 1)
                            .coord
                            .multiply(num)
                            .added(&last.pole(i - 1).coord.multiply(den))
                            .multiply(1.0 / (self.nb_sections - 1) as f64);
                        if i - 1 < poles.len() {
                            poles[i - 1] = GpPnt::from_xyz(&xyz);
                        }
                        let wnum = num * first.weight(i as i32) + den * last.weight(i as i32);
                        if i - 1 < weights.len() {
                            weights[i - 1] = wnum / (self.nb_sections - 1) as f64;
                        }
                    }
                }
                for pole in poles.iter_mut() {
                    pole.transform(cumul);
                }
            }
            true
        } else {
            // `myType == 0` (`cxx:579-663`): build the section arc from the path
            // point and the two bounding curves.
            let (Some(path), Some(first), Some(last)) =
                (&self.adp_path, &self.adp_first_sect, &self.adp_last_sect)
            else {
                return false;
            };
            let coef = (p - 1) as f64 / (self.nb_sections - 1) as f64;
            let u = (1.0 - coef) * path.first_parameter() + coef * path.last_parameter();
            let p_path = path.d0(u);

            let alpha =
                (u - path.first_parameter()) / (path.last_parameter() - path.first_parameter());

            let mut u1 = (1.0 - alpha) * first.first_parameter() + alpha * first.last_parameter();
            if let Some(line) = first.gp_line() {
                if first.first_parameter().is_infinite() || first.last_parameter().is_infinite() {
                    u1 = clib::line_parameter(&line, &p_path);
                }
            }
            let p1 = first.d0(u1);

            let mut u2 = (1.0 - alpha) * last.first_parameter() + alpha * last.last_parameter();
            if let Some(line) = last.gp_line() {
                if last.first_parameter().is_infinite() || last.last_parameter().is_infinite() {
                    u2 = clib::line_parameter(&line, &p_path);
                }
            }
            let p2 = last.d0(u2);

            let mut angle = 0.0;
            let mut axis = None;
            if p1.distance(&p2) >= CONFUSION {
                let v1 = GpVec::from_pnts(&p_path, &p1);
                let v2 = GpVec::from_pnts(&p_path, &p2);
                let n = v1.crossed(&v2);
                let (Ok(z), Ok(x)) = (GpDir::from_xyz(&n.coord), GpDir::from_xyz(&v1.coord)) else {
                    return false;
                };
                let Ok(ax2) = GpAx2::new(p_path.clone(), z, x) else {
                    return false;
                };
                angle = clib::circle_parameter(&ax2, &p2);
                axis = Some(ax2);
            }

            if angle < ANGULAR {
                // Degenerate arc: every pole collapses onto P1.
                for pole in poles.iter_mut() {
                    *pole = p1.clone();
                }
                for weight in weights.iter_mut() {
                    *weight = 1.0;
                }
            } else {
                let Some(ax2) = axis else { return false };
                let circ = GpCirc::new(ax2, self.radius);
                let parameterisation = if self.polynomial {
                    ParameterisationType::Polynomial
                } else {
                    ParameterisationType::QuasiAngular
                };
                let Ok(bs) = circle_arc_bspline_of(&circ, 0.0, angle, parameterisation) else {
                    return false;
                };
                let n = poles.len().min(bs.nb_poles());
                for i in 0..n {
                    poles[i] = bs.pole(i).clone();
                }
                let nw = weights.len().min(bs.nb_poles());
                for i in 0..nw {
                    weights[i] = bs.weight((i + 1) as i32);
                }
            }
            true
        }
    }

    /// `Transformation(Index)` (`cxx:667-676`), 1-based. OCCT throws
    /// `Standard_RangeError` past the end; the port returns `None`.
    pub fn transformation(&self, index: i32) -> Option<&GpTrsf> {
        if index < 1 {
            return None;
        }
        self.trsfs.get((index - 1) as usize)
    }

    /// `Parameter(P)` (`cxx:680-698`).
    pub fn parameter(&self, p: i32) -> Option<f64> {
        let path = self.path.as_ref()?;
        if p == 1 {
            return Some(path.first_parameter());
        }
        if p == self.nb_sections {
            return Some(path.last_parameter());
        }
        let u1 = path.first_parameter();
        let u2 = path.last_parameter();
        Some(
            ((self.nb_sections - p) as f64 * u1 + (p - 1) as f64 * u2)
                / (self.nb_sections - 1) as f64,
        )
    }

    /// `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `NbSections()`.
    pub fn nb_sections(&self) -> i32 {
        self.nb_sections
    }

    /// `Path()` (`GeomFill_SweepSectionGenerator.hxx`).
    pub fn path(&self) -> Option<&GeomBSplineCurve> {
        self.path.as_ref()
    }

    /// `FirstSect()`.
    pub fn first_section(&self) -> Option<&GeomBSplineCurve> {
        self.first_sect.as_ref()
    }

    /// `LastSect()`.
    pub fn last_section(&self) -> Option<&GeomBSplineCurve> {
        self.last_sect.as_ref()
    }
}

/// `if (Path->IsKind(STANDARD_TYPE(Geom_BSplineCurve))) down_cast(Copy()) else
/// GeomConvert::CurveToBSplineCurve(Path)` (`cxx:106-113`).
fn path_to_bspline(path: &dyn Curve) -> Result<GeomBSplineCurve, String> {
    if is_bspline(path) {
        return bspline_copy(path).map_err(|e| format!("GeomFill_SweepSectionGenerator: {e:?}"));
    }
    curve_to_bspline_curve(path, ParameterisationType::TgtThetaOver2)
        .map_err(|e| format!("GeomFill_SweepSectionGenerator: {e:?}"))
}

/// `if (Sect->IsKind(STANDARD_TYPE(Geom_BSplineCurve))) down_cast(Copy()) else
/// GeomConvert::CurveToBSplineCurve(Sect, Convert_QuasiAngular)`
/// (`cxx:142-152`).
fn section_to_bspline(sect: &dyn Curve) -> Result<GeomBSplineCurve, String> {
    if is_bspline(sect) {
        return bspline_copy(sect).map_err(|e| format!("GeomFill_SweepSectionGenerator: {e:?}"));
    }
    curve_to_bspline_curve(sect, ParameterisationType::QuasiAngular)
        .map_err(|e| format!("GeomFill_SweepSectionGenerator: {e:?}"))
}
