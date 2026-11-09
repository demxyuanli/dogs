//! `GCPnts_QuasiUniformDeflection` — the recursive deflection sampler the
//! sweep-section generator uses to size its section count.
//!
//! Source: `ModelingData/TKGeomBase/GCPnts/GCPnts_QuasiUniformDeflection.cxx:44-623`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - The `template <class TheCurve>` over `Adaptor3d_Curve` /
//!   `Adaptor2d_Curve2d` is the private [`SampleCurve`] trait: it is
//!   implemented for every `Curve` and for the 2D [`Curve2dRange`] wrapper,
//!   which is `Geom2dAdaptor_Curve(C, U1, U2)` (`Geom2dAdaptor_Curve.cxx:
//!   272-291`). Both 2D overloads (`cxx:496-507`, `:521-528`, `:542-548`,
//!   `:562-568`) therefore run the same template body.
//! - `Value(theC, u)` is [`Curve::d0`], `D1(theC, u, P, V)` is
//!   [`Curve::d1`]; `theC.Resolution(Precision::Confusion())` is
//!   [`Curve::resolution`]. For the 2D adaptor, `Value` / `D1`
//!   (`cxx:31-52`) lift the `gp_Pnt2d` / `gp_Vec2d` into `z = 0`.
//! - `GetType()` is resolved from the adaptor queries the port's `Curve` trait
//!   already exposes (`gp_line` / `gp_circ` / `bspline_poles` / `bezier_poles`),
//!   which is how `GeomAdaptor_Curve::GetType` (`GeomAdaptor_Curve.cxx:960-984`)
//!   classifies them. The 2D twin reads `is_line` / `gp_circ2d` /
//!   `bspline2d().nb_poles()` / `bezier_nb_poles()`, which is
//!   `Geom2dAdaptor_Curve::GetType` (`Geom2dAdaptor_Curve.cxx:288-330`).
//! - `GeomAbs_Shape` continuity is the crate's `u8` encoding (see
//!   `curve.rs:174-181`: `0 = C0`, `1 = G1`, `2 = C1`, …); `myCont` keeps
//!   OCCT's `> GeomAbs_G1` test (`cxx:580`).
//! - `NCollection_Sequence` becomes `Vec`, with OCCT's 1-based `Value(i)`
//!   access kept as [`QuasiUniformDeflection::value`] / `parameter` taking a
//!   1-based index.
//! - `BSplCLib::Hunt` (`BSplCLib.cxx:74-107`) is [`hunt`], same dichotomy.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{epsilon, CONFUSION};
use occt_geom2d::Curve2d;

use crate::curve::Curve;

/// `MyMaxQuasiFleshe = 2000` (`cxx:24`).
const MY_MAX_QUASI_FLECHE: i32 = 2000;

/// The `GCPnts_DeflectionType` cases the sampler switches on
/// (`GCPnts_DeflectionType.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DeflectionType {
    /// `GCPnts_Linear`.
    Linear,
    /// `GCPnts_Circular`.
    Circular,
    /// `GCPnts_Curved`.
    Curved,
    /// `GCPnts_DefComposite`.
    DefComposite,
}

/// The adaptor interface the template body of
/// `GCPnts_QuasiUniformDeflection` reads (`cxx:26-52`, `:322-348`, `:587`).
///
/// OCCT instantiates `TheCurve` with `Adaptor3d_Curve` and
/// `Adaptor2d_Curve2d`; the port uses this trait for the pair.
trait SampleCurve {
    /// `Value(theC, u)` (`cxx:26-36`).
    fn d0(&self, u: f64) -> GpPnt;
    /// `D1(theC, u, P, V)` (`cxx:39-52`).
    fn d1(&self, u: f64) -> (GpPnt, GpVec);
    fn first_parameter(&self) -> f64;
    fn last_parameter(&self) -> f64;
    /// `theC.Resolution(Precision::Confusion())` (`cxx:587`).
    fn resolution(&self, r3d: f64) -> f64;
    /// `theC.Circle().Radius()` (`cxx:304`).
    fn circle_radius(&self) -> Option<f64>;
    /// `GetDefType(theC)` (`cxx:322-348`).
    fn def_type(&self) -> DeflectionType;
    /// `theC.NbIntervals(S)` (`cxx:326`, `:424`).
    fn nb_intervals(&self, s: u8) -> i32;
    /// `theC.Intervals(aTI, S)` (`cxx:427`).
    fn intervals(&self, s: u8) -> Vec<f64>;
    /// `theC.GetType()` is `GeomAbs_BSplineCurve` or `GeomAbs_BezierCurve`
    /// (`cxx:596`).
    fn knot_based(&self) -> bool;
}

/// `QuasiFleche`'s `GetDefType` over a 3D adaptor (`cxx:322-348`).
fn deflection_type_of_curve<T: Curve + ?Sized>(c: &T) -> DeflectionType {
    // `NbIntervals(GeomAbs_C1) > 1` - `C1` is 2 in the crate's encoding.
    if c.nb_intervals(2) > 1 {
        return DeflectionType::DefComposite;
    }
    if c.gp_line().is_some() {
        return DeflectionType::Linear;
    }
    if c.gp_circ().is_some() {
        return DeflectionType::Circular;
    }
    if let Some(poles) = c.bspline_poles() {
        return if poles.len() == 2 {
            DeflectionType::Linear
        } else {
            DeflectionType::Curved
        };
    }
    if let Some(poles) = c.bezier_poles() {
        return if poles.len() == 2 {
            DeflectionType::Linear
        } else {
            DeflectionType::Curved
        };
    }
    DeflectionType::Curved
}

impl<T: Curve + ?Sized> SampleCurve for T {
    fn d0(&self, u: f64) -> GpPnt {
        Curve::d0(self, u)
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        Curve::d1(self, u)
    }
    fn first_parameter(&self) -> f64 {
        Curve::first_parameter(self)
    }
    fn last_parameter(&self) -> f64 {
        Curve::last_parameter(self)
    }
    fn resolution(&self, r3d: f64) -> f64 {
        Curve::resolution(self, r3d)
    }
    fn circle_radius(&self) -> Option<f64> {
        Curve::circle_radius(self)
    }
    fn def_type(&self) -> DeflectionType {
        deflection_type_of_curve(self)
    }
    fn nb_intervals(&self, s: u8) -> i32 {
        Curve::nb_intervals(self, s)
    }
    fn intervals(&self, s: u8) -> Vec<f64> {
        Curve::parameter_intervals(self, s)
    }
    fn knot_based(&self) -> bool {
        Curve::bspline_poles(self).is_some() || Curve::bezier_poles(self).is_some()
    }
}

/// `Geom2dAdaptor_Curve(C, UFirst, ULast)` (`Geom2dAdaptor_Curve.cxx:272-291`):
/// a 2D curve restricted to `[UFirst, ULast]`, lifted into `z = 0`.
///
/// `load` keeps the requested range and, for a `Geom2d_TrimmedCurve`, switches
/// the *type* to the basis curve (`:288-291`); [`Curve2d::trimmed_basis`] is
/// that unwrap.
pub struct Curve2dRange<'a> {
    curve: &'a dyn Curve2d,
    basis: &'a dyn Curve2d,
    first: f64,
    last: f64,
}

impl<'a> Curve2dRange<'a> {
    /// `Geom2dAdaptor_Curve::Load(C, UFirst, ULast)` (`cxx:272-291`).
    pub fn new(curve: &'a dyn Curve2d, first: f64, last: f64) -> Self {
        let basis = curve.trimmed_basis().unwrap_or(curve);
        Self {
            curve,
            basis,
            first,
            last,
        }
    }
}

impl SampleCurve for Curve2dRange<'_> {
    fn d0(&self, u: f64) -> GpPnt {
        // `Value(const Adaptor2d_Curve2d&, u)` (`cxx:31-36`).
        let p = self.curve.d0(u);
        GpPnt::new(p.x(), p.y(), 0.0)
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        // `D1(const Adaptor2d_Curve2d&, u, P, V)` (`cxx:44-52`).
        let (p, v) = self.curve.d1(u);
        (
            GpPnt::new(p.x(), p.y(), 0.0),
            GpVec::new(v.x(), v.y(), 0.0),
        )
    }
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn resolution(&self, ruv: f64) -> f64 {
        self.basis.resolution(ruv)
    }
    fn circle_radius(&self) -> Option<f64> {
        self.basis.gp_circ2d().map(|c| c.radius())
    }
    fn def_type(&self) -> DeflectionType {
        // `GetDefType` over the 2D adaptor (`cxx:322-348`).
        if self.basis.nb_intervals(2) > 1 {
            return DeflectionType::DefComposite;
        }
        if self.basis.is_line() {
            return DeflectionType::Linear;
        }
        if self.basis.gp_circ2d().is_some() {
            return DeflectionType::Circular;
        }
        if let Some(poles) = self.basis.bspline2d().map(|b| b.nb_poles()) {
            return if poles == 2 {
                DeflectionType::Linear
            } else {
                DeflectionType::Curved
            };
        }
        if let Some(poles) = self.basis.bezier_nb_poles() {
            return if poles == 2 {
                DeflectionType::Linear
            } else {
                DeflectionType::Curved
            };
        }
        DeflectionType::Curved
    }
    fn nb_intervals(&self, s: u8) -> i32 {
        self.basis.nb_intervals(s)
    }
    fn intervals(&self, s: u8) -> Vec<f64> {
        self.basis.parameter_intervals(s)
    }
    fn knot_based(&self) -> bool {
        self.basis.bspline2d().is_some() || self.basis.bezier_nb_poles().is_some()
    }
}

/// `BSplCLib::Hunt(theArray, theX, theXPos)` (`BSplCLib.cxx:74-107`) on a
/// 1-based array: `theXPos` such that `theArray(theXPos) <= theX <
/// theArray(theXPos + 1)`, `Lower - 1` below the range and `Upper + 1` above it.
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

/// `PerformLinear` (`cxx:276-291`).
fn perform_linear<C: SampleCurve + ?Sized>(
    c: &C,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
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
fn perform_circular<C: SampleCurve + ?Sized>(
    c: &C,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
    deflection: f64,
    u1: f64,
    u2: f64,
) -> bool {
    let radius = c.circle_radius().unwrap_or(1.0);
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

/// `GetDefType` (`cxx:322-348`).
fn get_def_type<C: SampleCurve + ?Sized>(c: &C) -> DeflectionType {
    c.def_type()
}

/// `QuasiFleche` with the `D1` overload (`cxx:56-185`).
#[allow(clippy::too_many_arguments)]
fn quasi_fleche_d1<C: SampleCurve + ?Sized>(
    c: &C,
    deflection2: f64,
    u_deb: f64,
    p_deb: &GpPnt,
    v_deb: &GpVec,
    u_fin: f64,
    p_fin: &GpPnt,
    v_fin: &GpVec,
    nb_min: i32,
    eps: f64,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
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
        p_delta = p_fin.clone();
        v_delta = v_fin.clone();
    }

    // Square length of the chord.
    let norme = GpVec::from_pnts(p_deb, &p_delta).square_magnitude();
    let mut fleche = 0.0;
    let mut is_fleche_ok = false;
    if norme > eps && norme > 16.0 * deflection2 {
        // Evaluation of the deflection by interpolation.
        let n1 = v_deb.square_magnitude();
        let n2 = v_delta.square_magnitude();
        if n1 > eps && n2 > eps {
            let normediff = v_deb
                .normalized()
                .coord
                .subtracted(&v_delta.normalized().coord)
                .square_modulus();
            if normediff > eps {
                fleche = normediff * norme / 64.0;
                is_fleche_ok = true;
            }
        }
    }

    let pmid = GpPnt::from_xyz(&p_deb.coord.added(&p_delta.coord).multiply(0.5));
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
            u_deb + u_delta,
            &p_delta,
            &v_delta,
            3,
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
fn quasi_fleche<C: SampleCurve + ?Sized>(
    c: &C,
    deflection2: f64,
    u_deb: f64,
    p_deb: &GpPnt,
    u_fin: f64,
    p_fin: &GpPnt,
    nb_min: i32,
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
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
        p_delta = p_fin.clone();
    }

    let pmid = GpPnt::from_xyz(&p_deb.coord.added(&p_delta.coord).multiply(0.5));
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
fn perform_curve<C: SampleCurve + ?Sized>(
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
    c: &C,
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
        points.push(pdeb.clone());

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
        points.push(pdeb.clone());

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
fn perform_composite<C: SampleCurve + ?Sized>(
    params: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
    c: &C,
    deflection: f64,
    u1: f64,
    u2: f64,
    epsilon_: f64,
    continuity: u8,
) -> bool {
    // Consistency with Intervals.
    let nb_intervals = c.nb_intervals(3);
    let ti = c.intervals(3);
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

/// `GCPnts_QuasiUniformDeflection` (`GCPnts_QuasiUniformDeflection.hxx`).
pub struct QuasiUniformDeflection {
    done: bool,
    deflection: f64,
    cont: u8,
    params: Vec<f64>,
    points: Vec<GpPnt>,
}

impl Default for QuasiUniformDeflection {
    fn default() -> Self {
        Self::new()
    }
}

impl QuasiUniformDeflection {
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
    /// (`cxx:552-558` -> `initialize`, `cxx:574-623`) for a 3D adaptor.
    pub fn initialize(&mut self, c: &dyn Curve, deflection: f64, u1: f64, u2: f64, continuity: u8) {
        self.initialize_curve(c, deflection, u1, u2, continuity);
    }

    /// `Initialize(const Adaptor2d_Curve2d&, theDeflection, theU1, theU2,
    /// theContinuity)` (`cxx:562-568` -> `initialize`, `cxx:574-623`).
    ///
    /// `u1` / `u2` are the `UFirst` / `ULast` of the
    /// `Geom2dAdaptor_Curve(C, UFirst, ULast)` the reference passes in.
    pub fn initialize_2d(
        &mut self,
        c: &dyn Curve2d,
        u1: f64,
        u2: f64,
        deflection: f64,
        continuity: u8,
    ) {
        self.initialize_curve(&Curve2dRange::new(c, u1, u2), deflection, u1, u2, continuity);
    }

    /// `initialize(theC, theDeflection, theU1, theU2, theContinuity)`
    /// (`cxx:574-623`).
    fn initialize_curve<C: SampleCurve + ?Sized>(
        &mut self,
        c: &C,
        deflection: f64,
        u1: f64,
        u2: f64,
        continuity: u8,
    ) {
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
            && c.knot_based()
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
    pub fn initialize_whole(&mut self, c: &dyn Curve, deflection: f64, continuity: u8) {
        self.initialize(c, deflection, c.first_parameter(), c.last_parameter(), continuity);
    }

    /// `Initialize(const Adaptor2d_Curve2d&, theDeflection, theContinuity)`
    /// (`cxx:542-548`).
    pub fn initialize_2d_whole(&mut self, c: &dyn Curve2d, deflection: f64, continuity: u8) {
        self.initialize_2d(
            c,
            c.first_parameter(),
            c.last_parameter(),
            deflection,
            continuity,
        );
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
    pub fn value(&self, index: i32) -> Option<GpPnt> {
        if !self.done || index < 1 {
            return None;
        }
        self.points.get((index - 1) as usize).cloned()
    }

    /// `Deflection()`.
    pub fn deflection(&self) -> f64 {
        self.deflection
    }
}
