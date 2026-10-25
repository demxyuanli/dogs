//! `ChFi3d_Boite` — the min/max UV box of the fillet to intersect.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:242-284`
//! (OCCT 8.0.0):
//!
//! - `ChFi3d_Boite(p1, p2, mu, Mu, mv, Mv)` (`cxx:242-253`),
//! - `ChFi3d_Boite(p1, p2, p3, p4, Du, Dv, mu, Mu, mv, Mv)` (`cxx:259-284`).
//!
//! Translation note: OCCT answers six out-parameters by reference; the port
//! answers the same values in the two plain structs below.

use occt_core::gp::GpPnt2d;

/// The `(mu, Mu, mv, Mv)` answered by `ChFi3d_Boite(p1, p2, ...)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boite {
    /// `mu`.
    pub u_min: f64,
    /// `Mu`.
    pub u_max: f64,
    /// `mv`.
    pub v_min: f64,
    /// `Mv`.
    pub v_max: f64,
}

/// The `(Du, Dv, mu, Mu, mv, Mv)` answered by `ChFi3d_Boite(p1, p2, p3, p4,
/// ...)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Boite4 {
    /// `Du`.
    pub du: f64,
    /// `Dv`.
    pub dv: f64,
    /// `mu`.
    pub u_min: f64,
    /// `Mu`.
    pub u_max: f64,
    /// `mv`.
    pub v_min: f64,
    /// `Mv`.
    pub v_max: f64,
}

/// `ChFi3d_Boite(p1, p2, mu, Mu, mv, Mv)` (`cxx:242-253`).
pub fn ch_fi3d_boite(p1: &GpPnt2d, p2: &GpPnt2d) -> Boite {
    Boite {
        u_min: p1.x().min(p2.x()),
        u_max: p1.x().max(p2.x()),
        v_min: p1.y().min(p2.y()),
        v_max: p1.y().max(p2.y()),
    }
}

/// `ChFi3d_Boite(p1, p2, p3, p4, Du, Dv, mu, Mu, mv, Mv)` (`cxx:259-284`).
pub fn ch_fi3d_boite4(p1: &GpPnt2d, p2: &GpPnt2d, p3: &GpPnt2d, p4: &GpPnt2d) -> Boite4 {
    let u_min = p1.x().min(p2.x()).min(p3.x().min(p4.x()));
    let u_max = p1.x().max(p2.x()).max(p3.x().max(p4.x()));
    let v_min = p1.y().min(p2.y()).min(p3.y().min(p4.y()));
    let v_max = p1.y().max(p2.y()).max(p3.y().max(p4.y()));
    Boite4 {
        du: u_max - u_min,
        dv: v_max - v_min,
        u_min,
        u_max,
        v_min,
        v_max,
    }
}
