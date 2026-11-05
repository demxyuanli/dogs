//! `GCPnts_QuasiUniformDeflection` instantiated on `Adaptor2d_Curve2d`.
//!
//! Source: `ModelingData/TKGeomBase/GCPnts/GCPnts_QuasiUniformDeflection.cxx`,
//! a template over `TheCurve` that OCCT instantiates twice
//! (`GCPnts_QuasiUniformDeflection::Initialize(const Adaptor2d_Curve2d&, ...)`
//! and the `Adaptor3d_Curve` overload, `cxx:562-570`). The `Adaptor3d_Curve`
//! half is [`occt_geom::gcpnts::QuasiUniformDeflection`]; this is the
//! `Adaptor2d_Curve2d` half, needed by `BRepTopAdaptor_FClass2d`'s
//! re-discretization loop (`BRepTopAdaptor_FClass2d.cxx:349-436`).
//!
//! The algorithm is identical branch for branch: every geometric quantity it
//! compares (chord length, deflection, `Resolution`) is a distance, and the 2D
//! distances are the same numbers the 3D instantiation would compute on the
//! `z = 0` embedding, so the point counts match the template.
//!
//! Translation notes:
//! - `TheCurve::GetType()` is read through the port's `Curve2d` queries:
//!   [`Curve2d::is_line`] (`GeomAbs_Line`), [`Curve2d::gp_circ2d`]
//!   (`GeomAbs_Circle`), [`Curve2d::bspline2d`] (`GeomAbs_BSplineCurve`) and
//!   [`Curve2d::bezier_nb_poles`] (`GeomAbs_BezierCurve`), which is how
//!   `Geom2dAdaptor_Curve::GetType` classifies them.
//! - `BSplCLib::Hunt` (`BSplCLib.cxx:74-107`) is [`hunt`], same dichotomy as
//!   the 3D half.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::{epsilon, CONFUSION};

use crate::curve::Curve2d;

/// `MyMaxQuasiFleshe = 2000` (`cxx:24`).
const MY_MAX_QUASI_FLECHE: i32 = 2000;

/// `GCPnts_DeflectionType` (`GCPnts_DeflectionType.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DeflectionType {
    Linear,
    Circular,
    Curved,
    DefComposite,
}

/// `BSplCLib::Hunt(theArray, theX, theXPos)` (`BSplCLib.cxx:74-107`) on a
/// 1-based array.
fn hunt(array: &[f64], x: f64) -> i32 {
    let lower = 1i32;
    let upper = array.len() as i32;
    if array.is_empty() {
        return 0;
    }
    let at = |i: i32| array[(i - 1) as usize];
    if at(lower) > x {
        return lower - 1;
    }
    if at(upper) < x {
        return upper + 1;
    }
    let mut x_pos = lower;
    if array.len() <= 1 {
        return x_pos;
    }
    let mut hi = upper;
    while hi - x_pos != 1 {
        let mid = (hi + x_pos) / 2;
        if at(mid) < x {
            x_pos = mid;
        } else {
            hi = mid;
        }
    }
    x_pos
}

/// `GetDefType` (`cxx:322-348`) for the 2D adaptor.
fn get_def_type(c: &dyn Curve2d) -> DeflectionType {
    // `NbIntervals(GeomAbs_C1) > 1` — `C1` is 2 in the crate's encoding.
    if c.nb_intervals(2) > 1 {
        return DeflectionType::DefComposite;
    }
    if c.is_line() {
        return DeflectionType::Linear;
    }
    if c.gp_circ2d().is_some() {
        return DeflectionType::Circular;
    }
    if let Some(nb_poles) = c.bspline2d().map(|b| b.nb_poles()) {
        return if nb_poles == 2 {
            DeflectionType::Linear
        } else {
            DeflectionType::Curved
        };
    }
    // A trimmed / reparametrised wrapper reports `GeomAbs_BSplineCurve` through
    // `is_bspline2d` but exposes no `Geom2d_BSplineCurve`; its pole count is
    // not reachable through `Curve2d`, so the `NbPoles() == 2` shortcut is
    // skipped and the curve takes the generic arm.
    if c.is_bspline2d() {
        return DeflectionType::Curved;
    }
    if c.is_bezier2d() {
        return if c.bezier_nb_poles() == Some(2) {
            DeflectionType::Linear
        } else {
            DeflectionType::Curved
        };
    }
    DeflectionType::Curved
}

/// `PerformLinear` (`cxx:276-291`).
fn perform_linear(
    c: &dyn Curve2d,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt2d>,
    u1: f64,
    u2: f64,
) -> bool {
    params.push(u1);
    points.push(c.d0(u1));
    params.push(u2);
    points.push(c.d0(u2));
    true
}

/// `PerformCircular` (`cxx:296-316`).
fn perform_circular(
    c: &dyn Curve2d,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt2d>,
    deflection: f64,
    u1: f64,
    u2: f64,
) -> bool {
    let radius = c.gp_circ2d().map(|c| c.radius()).unwrap_or(1.0);
    let mut angle = (1.0 - (deflection / radius)).max(0.0);
    angle = 2.0 * angle.acos();
    let mut nb_points = ((u2 - u1) / angle) as i32;
    nb_points += 2;
    angle = (u2 - u1) / (nb_points - 1) as f64;
    let mut u = u1;
    for _ in 1..=nb_points {
        params.push(u);
        points.push(c.d0(u));
        u += angle;
    }
    true
}

/// `QuasiFleche` with derivatives (`cxx:150-184`).
#[allow(clippy::too_many_arguments)]
fn quasi_fleche_d1(
    c: &dyn Curve2d,
    deflection2: f64,
    u_deb: f64,
    p_deb: &GpPnt2d,
    v_deb: &GpVec2d,
    u_fin: f64,
    p_fin: &GpPnt2d,
    v_fin: &GpVec2d,
    nb_min: i32,
    eps: f64,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt2d>,
    nb_calls: &mut i32,
) {
    *nb_calls += 1;
    if *nb_calls >= MY_MAX_QUASI_FLECHE {
        return;
    }
    let pts_length = points.len() as i32;
    if *nb_calls > 100 && pts_length < 2 {
        return;
    }

    let mut u_delta = u_fin - u_deb;
    let (p_delta, v_delta);
    if nb_min > 2 {
        u_delta /= (nb_min - 1) as f64;
        let (p, v) = c.d1(u_deb + u_delta);
        p_delta = p;
        v_delta = v;
    } else {
        p_delta = *p_fin;
        v_delta = *v_fin;
    }

    // Square length of the chord.
    let norme = p_deb.coord.subtracted(&p_delta.coord).square_modulus();
    let mut fleche = 0.0;
    let mut is_fleche_ok = false;
    if norme > eps && norme > 16.0 * deflection2 {
        // Evaluation of the deflection by interpolation.
        let n1 = v_deb.coord.square_modulus();
        let n2 = v_delta.coord.square_modulus();
        if n1 > eps && n2 > eps {
            let normediff = match (v_deb.coord.normalized(), v_delta.coord.normalized()) {
                (Ok(a), Ok(b)) => a.subtracted(&b).square_modulus(),
                _ => 0.0,
            };
            if normediff > eps {
                fleche = normediff * norme / 64.0;
                is_fleche_ok = true;
            }
        }
    }

    let pmid = GpPnt2d::from_xy(p_deb.coord.added(&p_delta.coord).multiply(0.5));
    let p_verif = c.d0(u_deb + u_delta * 0.5);
    let fleche_mid_mid = pmid.square_distance(&p_verif);

    if is_fleche_ok {
        if fleche_mid_mid > norme / 4.0 + deflection2 {
            fleche = fleche_mid_mid;
        }
    } else {
        fleche = fleche_mid_mid;
    }

    if fleche < deflection2 {
        params.push(u_deb + u_delta);
        points.push(p_delta);
    } else {
        quasi_fleche_d1(
            c,
            deflection2,
            u_deb,
            p_deb,
            v_deb,
            u_deb + u_delta * 0.5,
            &p_verif,
            &c.d1(u_deb + u_delta * 0.5).1,
            2,
            eps,
            params,
            points,
            nb_calls,
        );
        quasi_fleche_d1(
            c,
            deflection2,
            u_deb + u_delta * 0.5,
            &p_verif,
            &c.d1(u_deb + u_delta * 0.5).1,
            u_deb + u_delta,
            &p_delta,
            &v_delta,
            2,
            eps,
            params,
            points,
            nb_calls,
        );
    }

    if nb_min > 2 {
        let remaining = nb_min - (points.len() as i32 - pts_length);
        quasi_fleche_d1(
            c,
            deflection2,
            u_deb + u_delta,
            &p_delta,
            &v_delta,
            u_fin,
            p_fin,
            v_fin,
            remaining,
            eps,
            params,
            points,
            nb_calls,
        );
    }
    *nb_calls -= 1;
}

/// `QuasiFleche` without derivatives (`cxx:189-273`).
#[allow(clippy::too_many_arguments)]
fn quasi_fleche(
    c: &dyn Curve2d,
    deflection2: f64,
    u_deb: f64,
    p_deb: &GpPnt2d,
    u_fin: f64,
    p_fin: &GpPnt2d,
    nb_min: i32,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt2d>,
    nb_calls: &mut i32,
) {
    *nb_calls += 1;
    if *nb_calls >= MY_MAX_QUASI_FLECHE {
        return;
    }
    let pts_length = points.len() as i32;
    if *nb_calls > 100 && pts_length < 2 {
        return;
    }

    let mut u_delta = u_fin - u_deb;
    let p_delta;
    if nb_min > 2 {
        u_delta /= (nb_min - 1) as f64;
        p_delta = c.d0(u_deb + u_delta);
    } else {
        p_delta = *p_fin;
    }

    let pmid = GpPnt2d::from_xy(p_deb.coord.added(&p_delta.coord).multiply(0.5));
    let p_verif = c.d0(u_deb + u_delta * 0.5);
    let fleche = pmid.square_distance(&p_verif);
    if fleche < deflection2 {
        params.push(u_deb + u_delta);
        points.push(p_delta);
    } else {
        quasi_fleche(
            c,
            deflection2,
            u_deb,
            p_deb,
            u_deb + u_delta * 0.5,
            &p_verif,
            2,
            params,
            points,
            nb_calls,
        );
        quasi_fleche(
            c,
            deflection2,
            u_deb + u_delta * 0.5,
            &p_verif,
            u_deb + u_delta,
            &p_delta,
            2,
            params,
            points,
            nb_calls,
        );
    }

    if nb_min > 2 {
        let remaining = nb_min - (points.len() as i32 - pts_length);
        quasi_fleche(
            c,
            deflection2,
            u_deb + u_delta,
            &p_delta,
            u_fin,
            p_fin,
            remaining,
            params,
            points,
            nb_calls,
        );
    }
    *nb_calls -= 1;
}

/// `PerformCurve` (`cxx:353-410`).
#[allow(clippy::too_many_arguments)]
fn perform_curve(
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt2d>,
    c: &dyn Curve2d,
    deflection: f64,
    u1: f64,
    u2: f64,
    epsilon_: f64,
    continuity: u8,
) -> bool {
    let nbmin = 2i32;
    let mut nb_call_qf = 0i32;

    if continuity <= 1 {
        let pdeb = c.d0(u1);
        params.push(u1);
        points.push(pdeb);

        let pfin = c.d0(u2);
        quasi_fleche(
            c,
            deflection * deflection,
            u1,
            &pdeb,
            u2,
            &pfin,
            nbmin,
            params,
            points,
            &mut nb_call_qf,
        );
    } else {
        let (pdeb, ddeb) = c.d1(u1);
        params.push(u1);
        points.push(pdeb);

        let decreased = u2 - epsilon(u2) * 10.0;
        let (pfin, dfin) = c.d1(decreased);
        quasi_fleche_d1(
            c,
            deflection * deflection,
            u1,
            &pdeb,
            &ddeb,
            u2,
            &pfin,
            &dfin,
            nbmin,
            epsilon_ * epsilon_,
            params,
            points,
            &mut nb_call_qf,
        );
    }
    true
}

/// `PerformComposite` (`cxx:415-459`).
#[allow(clippy::too_many_arguments)]
fn perform_composite(
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt2d>,
    c: &dyn Curve2d,
    deflection: f64,
    u1: f64,
    u2: f64,
    epsilon_: f64,
    continuity: u8,
) -> bool {
    // Consistency with Intervals.
    let nb_intervals = c.nb_intervals(3);
    let ti = c.parameter_intervals(3);
    let p_index = hunt(&ti, u1);

    // Iterate by continuous segments.
    let mut ua = u1;
    let mut an_index = p_index;
    loop {
        let ub = if an_index + 1 <= ti.len() as i32 {
            u2.min(ti[an_index as usize])
        } else {
            u2
        };
        if !perform_curve(params, points, c, deflection, ua, ub, epsilon_, continuity) {
            return false;
        }
        an_index += 1;
        if an_index > nb_intervals || u2 < ti[(an_index - 1) as usize] {
            return true;
        }
        // Remove the last point to avoid duplication.
        params.pop();
        points.pop();
        ua = ub;
    }
}

/// `GCPnts_QuasiUniformDeflection` over a 2D adaptor.
pub struct QuasiUniformDeflection2d {
    done: bool,
    deflection: f64,
    cont: u8,
    params: Vec<f64>,
    points: Vec<GpPnt2d>,
}

impl Default for QuasiUniformDeflection2d {
    fn default() -> Self {
        Self::new()
    }
}

impl QuasiUniformDeflection2d {
    /// `GCPnts_QuasiUniformDeflection()` (`cxx:472-477`).
    pub fn new() -> Self {
        Self {
            done: false,
            deflection: 0.0,
            cont: 1,
            params: Vec::new(),
            points: Vec::new(),
        }
    }

    /// `Initialize(theC, theDeflection, theU1, theU2, theContinuity)`
    /// (`cxx:574-623`).
    pub fn initialize(&mut self, c: &dyn Curve2d, deflection: f64, u1: f64, u2: f64, continuity: u8) {
        self.cont = if continuity > 1 { 2 } else { 0 };
        self.deflection = deflection;
        self.done = false;
        self.params.clear();
        self.points.clear();

        let epsilon_ = c.resolution(CONFUSION).min(1.0e50);
        let def_type = get_def_type(c);
        let lower = u1.min(u2);
        let upper = u1.max(u2);

        if (def_type == DeflectionType::Curved || def_type == DeflectionType::DefComposite)
            && (c.is_bspline2d() || c.is_bezier2d())
        {
            let max_par = c.first_parameter().abs().max(c.last_parameter().abs());
            if epsilon_ < epsilon(max_par) {
                return;
            }
        }

        match def_type {
            DeflectionType::Linear => {
                self.done = perform_linear(c, &mut self.params, &mut self.points, lower, upper);
            }
            DeflectionType::Circular => {
                self.done = perform_circular(
                    c,
                    &mut self.params,
                    &mut self.points,
                    deflection,
                    lower,
                    upper,
                );
            }
            DeflectionType::Curved => {
                self.done = perform_curve(
                    &mut self.params,
                    &mut self.points,
                    c,
                    deflection,
                    lower,
                    upper,
                    epsilon_,
                    self.cont,
                );
            }
            DeflectionType::DefComposite => {
                self.done = perform_composite(
                    &mut self.params,
                    &mut self.points,
                    c,
                    deflection,
                    lower,
                    upper,
                    epsilon_,
                    self.cont,
                );
            }
        }
    }

    /// `Initialize(theC, theDeflection, theContinuity)` (`cxx:532-538`).
    pub fn initialize_whole(&mut self, c: &dyn Curve2d, deflection: f64, continuity: u8) {
        self.initialize(c, deflection, c.first_parameter(), c.last_parameter(), continuity);
    }

    /// `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbPoints()`.
    pub fn nb_points(&self) -> i32 {
        self.points.len() as i32
    }

    /// `Parameter(Index)` (`GCPnts_QuasiUniformDeflection.hxx`), 1-based.
    pub fn parameter(&self, index: i32) -> Option<f64> {
        if !self.done || index < 1 {
            return None;
        }
        self.params.get((index - 1) as usize).copied()
    }

    /// `Value(Index)` (`cxx:463-467`), 1-based.
    pub fn value(&self, index: i32) -> Option<GpPnt2d> {
        if !self.done || index < 1 {
            return None;
        }
        self.points.get((index - 1) as usize).copied()
    }

    /// `Deflection()`.
    pub fn deflection(&self) -> f64 {
        self.deflection
    }
}
